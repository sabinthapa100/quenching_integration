# glauber.py
# Optical Glauber model with AA and pA support (CMS-style)
# Reproduces Mathematica notebook logic by Sabin Thapa & M. Strickland
# Author: (you)
#
# Qualities:
# - PEP8, docstrings, type hints, single-responsibility classes
# - Same formulas/limits as your Mathematica notebook (0806.1116, 1304.0901)
# - Robust units, progress prints, sanity checks, and centrality tables
# - Reusable API for eloss + (later) pT broadening
#
# Dependencies: numpy, matplotlib
# (No SciPy required. Gauss–Legendre quadratures + bilinear interpolation included.)

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable, Tuple, Dict, List, Literal, Optional
import math
import numpy as np
import matplotlib.pyplot as plt
import os


# ----------------------------- Utilities ------------------------------------


def _leggauss(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return Gauss–Legendre nodes/weights on [-1, 1]."""
    return np.polynomial.legendre.leggauss(n)


def _gl_integrate_1d(f: Callable[[float], float],
                     a: float, b: float, n: int = 64) -> float:
    """
    1D Gauss–Legendre integration on [a,b].

    IMPORTANT: f is evaluated *pointwise* (scalar), so callers do not need
    to vectorize their integrand. This prevents accidental broadcasting when
    f itself performs an internal multi-dimensional integral (e.g. over x,y).
    """
    x, w = _leggauss(n)
    xm = 0.5 * (b - a)
    xc = 0.5 * (b + a)
    total = 0.0
    for xi, wi in zip(x, w):
        val = f(float(xm * xi + xc))  # ensure scalar
        total += wi * float(np.asarray(val))
    return float(xm * total)


def _gl_integrate_2d(f: Callable[[np.ndarray, np.ndarray], np.ndarray],
                     ax: float, bx: float, ay: float, by: float,
                     nx: int = 48, ny: int = 48) -> float:
    """
    2D Gauss–Legendre integration on a rectangle.

    REQUIREMENT: bounds (ax,bx,ay,by) must be *scalars*.
    We construct a weight matrix W of shape (nx, ny) and multiply elementwise
    by f(X, Y) of the same shape.
    """
    for name, v in (("ax", ax), ("bx", bx), ("ay", ay), ("by", by)):
        if np.ndim(v) != 0:
            raise ValueError(
                f"_gl_integrate_2d: bound {name} must be scalar, got shape {np.shape(v)}"
            )

    x, wx = _leggauss(nx)
    y, wy = _leggauss(ny)
    xm, xc = 0.5 * (bx - ax), 0.5 * (bx + ax)
    ym, yc = 0.5 * (by - ay), 0.5 * (by + ay)

    X = xm * x[:, None] + xc                    # (nx, ny)
    Y = ym * y[None, :] + yc                    # (nx, ny)
    W = (xm * ym) * (wx[:, None] * wy[None, :])  # (nx, ny)

    F = f(X, Y)
    return float(np.sum(W * F))


def _interp1_linear(x: np.ndarray, y: np.ndarray,
                    xq: np.ndarray | float) -> np.ndarray | float:
    """Fast 1D linear interpolation with clamping on edges."""
    return np.interp(xq, x, y, left=y[0], right=y[-1])


def _bilinear(z: np.ndarray, x: np.ndarray, y: np.ndarray,
              xq: np.ndarray, yq: np.ndarray) -> np.ndarray:
    """
    Bilinear interpolation on a rect grid.

    Assumes z.shape == (nx, ny) with x along axis 0 and y along axis 1:
    z[ix, iy] ≡ z(x[ix], y[iy]).
    """
    xq = np.asarray(xq)
    yq = np.asarray(yq)
    xi = np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2)
    yi = np.clip(np.searchsorted(y, yq) - 1, 0, len(y) - 2)
    x0, x1 = x[xi], x[xi + 1]
    y0, y1 = y[yi], y[yi + 1]
    z00 = z[xi, yi]
    z10 = z[xi + 1, yi]
    z01 = z[xi, yi + 1]
    z11 = z[xi + 1, yi + 1]
    tx = (xq - x0) / np.maximum(x1 - x0, 1e-12)
    ty = (yq - y0) / np.maximum(y1 - y0, 1e-12)
    return ((1 - tx) * (1 - ty) * z00 +
            tx * (1 - ty) * z10 +
            (1 - tx) * ty * z01 +
            tx * ty * z11)


# --------------------------- Physics constants -------------------------------

MB_TO_FM2 = 0.1            # 1 mb = 0.1 fm^2
FM2_TO_MB = 10.0           # 1 fm^2 = 10 mb
DEFAULT_N0 = 0.17          # fm^-3 central density (Pb/Au standard)
HBARC = 0.1973269804       # GeV*fm (not used here but kept for completeness)

# √s-dependent σ_nn (mb) and Woods–Saxon diffuseness d (fm)
_SIGMA_NN_BY_ROOTS_MB = {
    200.0: 42.0,
    2760.0: 62.0,
    5023.0: 67.6,
    8160.0: 71.0,
}

_DIFFUSENESS_BY_ROOTS = {
    200.0: 0.535,
    2760.0: 0.549,
    5023.0: 0.549,
    8160.0: 0.549,
}


# --------------------------- Data classes ------------------------------------


@dataclass(frozen=True)
class WoodsSaxon:
    """
    Woods–Saxon nuclear density.

    n(r) = n0 / (1 + exp((r - Rn)/d)),   Rn = 1.12 A^(1/3) - 0.86 A^(-1/3)
    T_A(ρ) = ∫ dz n(√(ρ^2 + z^2))  (units: fm^-2)
    """
    A: int
    n0: float = DEFAULT_N0
    d_fm: float = 0.549
    rmax_fm: float = 20.0
    dr_fm: float = 0.05
    zmax_fm: float = 20.0
    nz: int = 200  # GL nodes for z-integration

    def radius_rn(self) -> float:
        a13 = self.A ** (1 / 3)
        return 1.12 * a13 - 0.86 / a13

    def n_of_r(self, r: np.ndarray) -> np.ndarray:
        Rn = self.radius_rn()
        return self.n0 / (1.0 + np.exp((r - Rn) / self.d_fm))

    def tabulate_T_of_r(self) -> Tuple[np.ndarray, np.ndarray]:
        """Tabulate T_A(ρ) on ρ ∈ [0, rmax] using nz-point GL in z ∈ [-zmax, zmax]."""
        r_grid = np.arange(0.0, self.rmax_fm + 1e-12, self.dr_fm)
        z_nodes, z_w = _leggauss(self.nz)
        zm = self.zmax_fm  # map [-1,1] -> [-zmax,zmax] is z = zm * x
        Z = zm * z_nodes  # (nz,)
        T_vals = []
        for r in r_grid:
            rr = np.sqrt(r * r + Z * Z)
            T_vals.append(zm * np.sum(z_w * self.n_of_r(rr)))
        return r_grid, np.array(T_vals)


@dataclass
class ProtonProfile:
    """
    Generalized exponential (Weibull-like) proton transverse profile:

    T_p(ρ) = m / [2π r_p^2 Γ(2/m)] * exp[- (ρ / r_p)^m], normalized to 1.
    """
    m: float = 1.85
    r_p_fm: float = 0.975

    def norm_const(self) -> float:
        return self.m / (2.0 * math.pi * self.r_p_fm ** 2 * math.gamma(2.0 / self.m))

    def T_p(self, rho: np.ndarray) -> np.ndarray:
        return self.norm_const() * np.exp(-(rho / self.r_p_fm) ** self.m)


@dataclass
class SystemSpec:
    """
    Collision system and beam energy.

    system: 'AA' (A+A) or 'pA' (p+A)
    """
    system: Literal["AA", "pA"]
    roots_GeV: float
    A: int
    sigma_nn_mb: Optional[float] = None
    diffuseness_fm: Optional[float] = None

    def resolve(self) -> Tuple[float, float]:
        """Return (σ_nn in mb, diffuseness d in fm) with √s defaults when unspecified."""
        sigma = self.sigma_nn_mb if self.sigma_nn_mb is not None else \
            _SIGMA_NN_BY_ROOTS_MB.get(self.roots_GeV, 67.6)
        dval = self.diffuseness_fm if self.diffuseness_fm is not None else \
            _DIFFUSENESS_BY_ROOTS.get(self.roots_GeV, 0.549)
        return float(sigma), float(dval)


# ----------------------- Optical Glauber model (OOP) -------------------------


class OpticalGlauber:
    """
    Optical Glauber machinery for AA and pA:
      - builds T_A(r) (Woods–Saxon) once, reuses via bilinear interpolation
      - computes T_AA(b), T_pA(b)
      - computes N_part(b), N_coll(b)
      - maps centrality percentiles -> b-edges, and bin-averages
        (⟨b⟩, ⟨N_part⟩, ⟨N_coll⟩), with inelastic weights
    """

    def __init__(self,
                 spec: SystemSpec,
                 rmax_fm: float = 20.0,
                 dr_fm: float = 0.05,
                 zmax_fm: float = 20.0,
                 nz_z: int = 200,
                 xylim_fm: float = 15.0,
                 nx: int = 120,
                 ny: int = 120,
                 nx_pa_window: int = 120,
                 ny_pa_window: int = 120,
                 pa_x_half_width_fm: float = 4.0,
                 pa_y_half_width_fm: float = 10.0,
                 verbose: bool = True) -> None:
        """
        Parameters mirror your Mathematica notebook:
          - T_A(r): z integration ±20 fm (nz=200)
          - T_AA(b): x,y ∈ [-15,15] fm
          - T_pA(b): x ∈ [b-4,b+4], y ∈ [-10,10] fm
        """
        self.spec = spec
        self.sigma_nn_mb, self.d_fm = spec.resolve()
        self.verbose = verbose

        # Woods–Saxon & tabulated T_A(ρ)
        self.ws = WoodsSaxon(
            A=spec.A, d_fm=self.d_fm,
            rmax_fm=rmax_fm, dr_fm=dr_fm,
            zmax_fm=zmax_fm, nz=nz_z
        )
        if verbose:
            print(f"[Glauber] Building T_A(r) table: A={spec.A}, d={self.d_fm:.3f} fm, "
                  f"r ∈ [0,{rmax_fm}] fm, dr={dr_fm} fm, zmax={zmax_fm} fm, nz={nz_z}")
        self.r_grid, self.T_r = self.ws.tabulate_T_of_r()

        # 2D grid for T_A(x,y) via bilinear interpolation
        self.xylim_fm = float(xylim_fm)
        self.nx, self.ny = int(nx), int(ny)
        self.x_grid = np.linspace(-self.xylim_fm, self.xylim_fm, self.nx)
        self.y_grid = np.linspace(-self.xylim_fm, self.xylim_fm, self.ny)
        X, Y = np.meshgrid(self.x_grid, self.y_grid, indexing="ij")  # (nx, ny)
        R = np.sqrt(X ** 2 + Y ** 2)
        self.T_xy = _interp1_linear(self.r_grid, self.T_r, R)  # (nx, ny)

        if verbose:
            # nested trapezoid (np.trapezoid) for ∫T_A d^2x
            integ_check = np.trapezoid(np.trapezoid(self.T_xy, self.y_grid, axis=1),
                                       self.x_grid, axis=0)
            print(f"[Glauber] ∫ T_A(x,y) d^2x ≈ {integ_check:.3f} (target A = {self.spec.A})")

        # Proton profile (normalized)
        self.proton = ProtonProfile()

        # pA integration windows (per notebook)
        self.pa_x_hw = float(pa_x_half_width_fm)
        self.pa_y_hw = float(pa_y_half_width_fm)
        self.nx_pa = int(nx_pa_window)
        self.ny_pa = int(ny_pa_window)

        # b-grid for T(b), cross-sections, and centrality mapping
        self.bmax_fm = 20.0
        self.db_fm = 0.1
        self.b_grid = np.arange(0.0, self.bmax_fm + 1e-12, self.db_fm)

        # Tabulate T_AA(b), T_pA(b) with progress
        if verbose:
            print("[Glauber] Tabulating T_AA(b) on b-grid…")
        self.TAA_b = np.empty_like(self.b_grid, dtype=float)
        for i, b in enumerate(self.b_grid):
            if verbose and (i % max(1, len(self.b_grid) // 10) == 0):
                print(f"  • T_AA: {i:3d}/{len(self.b_grid) - 1}   b = {b:5.2f} fm")
            self.TAA_b[i] = self._TAA_of_b(b)

        if verbose:
            print("[Glauber] Tabulating T_pA(b) on b-grid…")
        self.TpA_b = np.empty_like(self.b_grid, dtype=float)
        for i, b in enumerate(self.b_grid):
            if verbose and (i % max(1, len(self.b_grid) // 10) == 0):
                print(f"  • T_pA: {i:3d}/{len(self.b_grid) - 1}   b = {b:5.2f} fm")
            self.TpA_b[i] = self._TpA_of_b(b)

        # Total cross sections (mb) & cumulative fractions
        if verbose:
            print("[Glauber] Computing total σ_tot and cumulative fractions…")
        self.sigma_AA_tot_mb = self._sigma_tot_mb(kind="AA")
        self.sigma_pA_tot_mb = self._sigma_tot_mb(kind="pA")
        self.cum_AA = self._cumulative_fraction(kind="AA")
        self.cum_pA = self._cumulative_fraction(kind="pA")

        if verbose:
            print(f"[Glauber] σ_tot^AA ≈ {self.sigma_AA_tot_mb:.2f} mb,  "
                  f"σ_tot^pA ≈ {self.sigma_pA_tot_mb:.2f} mb")

    # ----------------------- Core field evaluators ---------------------------

    def T_A(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Thickness function T_A(x,y) via bilinear interpolation on pretabulated grid."""
        return _bilinear(self.T_xy, self.x_grid, self.y_grid, x, y)

    def _TAA_of_b(self, b_fm: float) -> float:
        """
        T_AA(b) = ∫ d^2x T_A(x+b/2,y) T_A(x-b/2,y)  (fm^-2),
        integrated over x,y ∈ [-15,15] fm.
        """
        bx = b_fm / 2.0

        def f(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
            Ta = self.T_A(X + bx, Y)
            Tb = self.T_A(X - bx, Y)
            return Ta * Tb

        return _gl_integrate_2d(f, -self.xylim_fm, self.xylim_fm,
                                -self.xylim_fm, self.xylim_fm,
                                nx=self.nx, ny=self.ny)

    def _TpA_of_b(self, b_fm: float) -> float:
        """
        T_pA(b) = ∫ d^2x T_A(x+b/2,y) T_p(x-b/2,y)  (fm^-2),
        with window x ∈ [b−4,b+4], y ∈ [−10,10] fm (as in the notebook).
        """
        bx = b_fm / 2.0
        xa, xb = b_fm - self.pa_x_hw, b_fm + self.pa_x_hw
        ya, yb = -self.pa_y_hw, self.pa_y_hw

        def f(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
            Ta = self.T_A(X + bx, Y)
            rho = np.sqrt((X - bx) ** 2 + Y ** 2)
            Tp = self.proton.T_p(rho)
            return Ta * Tp

        return _gl_integrate_2d(f, xa, xb, ya, yb, nx=self.nx_pa, ny=self.ny_pa)

    # ------------------------ Participants & Collisions ----------------------

    def N_coll_of_b(self, b_fm: float, kind: Literal["AA", "pA"]) -> float:
        """N_coll(b) = σ_nn(fm^2) * T_kind(b).  σ_nn(fm^2) = σ_nn(mb) × 0.1"""
        sigma_fm2 = self.sigma_nn_mb * MB_TO_FM2
        T = _interp1_linear(self.b_grid,
                            self.TAA_b if kind == "AA" else self.TpA_b, b_fm)
        return float(sigma_fm2 * T)

    def N_part_of_b_AA(self, b_fm: float) -> float:
        """
        AA participants:
          n_part(x,y,b) = T_A(x+b/2,y) * [1 - (1 - σ_nn T_A(x-b/2,y)/A)^A]
                        + (x+b/2 ↔ x-b/2)
          N_part(b) = ∫ d^2x n_part
        """
        bx = b_fm / 2.0
        sigma_fm2 = self.sigma_nn_mb * MB_TO_FM2
        A = float(self.spec.A)

        def f(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
            Ta = self.T_A(X + bx, Y)
            Tb = self.T_A(X - bx, Y)
            one_minus_a = np.power(1.0 - sigma_fm2 * Tb / A, A)
            one_minus_b = np.power(1.0 - sigma_fm2 * Ta / A, A)
            return Ta * (1.0 - one_minus_a) + Tb * (1.0 - one_minus_b)

        return _gl_integrate_2d(f, -self.xylim_fm, self.xylim_fm,
                                -self.xylim_fm, self.xylim_fm,
                                nx=self.nx, ny=self.ny)

    def N_part_of_b_pA(self, b_fm: float) -> float:
        """
        pA participants:
          n_part^pA = T_A(x+b/2,y) * [σ_nn T_p(x-b/2,y)]
                    + T_p(x-b/2,y) * [1 - (1 - σ_nn T_A(x+b/2,y)/A)^A]
          N_part(b) = ∫ d^2x n_part^pA
        Window: x ∈ [b−5,b+5], y ∈ [−15,15] fm (per notebook).
        """
        bx = b_fm / 2.0
        sigma_fm2 = self.sigma_nn_mb * MB_TO_FM2
        A = float(self.spec.A)
        xa, xb = b_fm - 5.0, b_fm + 5.0
        ya, yb = -15.0, 15.0

        def f(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
            Ta = self.T_A(X + bx, Y)
            rho = np.sqrt((X - bx) ** 2 + Y ** 2)
            Tp = self.proton.T_p(rho)
            term1 = Ta * (sigma_fm2 * Tp)
            one_minus = np.power(1.0 - sigma_fm2 * Ta / A, A)
            term2 = Tp * (1.0 - one_minus)
            return term1 + term2

        return _gl_integrate_2d(f, xa, xb, ya, yb, nx=self.nx_pa, ny=self.ny_pa)

    # ------------------------ Total cross sections ---------------------------

    def _sigma_tot_mb(self, kind: Literal["AA", "pA"]) -> float:
        """
        Total inelastic cross section (mb):
          σ_tot = 2π ∫_0^∞ b [1 - exp(-σ_nn T(b))] db  × 10
        """
        sigma_fm2 = self.sigma_nn_mb * MB_TO_FM2
        T = self.TAA_b if kind == "AA" else self.TpA_b
        integrand = self.b_grid * (1.0 - np.exp(-T * sigma_fm2))
        val_fm2 = 2.0 * math.pi * np.trapezoid(integrand, self.b_grid)
        return float(val_fm2 * FM2_TO_MB)

    def _cumulative_fraction(self, kind: Literal["AA", "pA"]) -> np.ndarray:
        """
        CDF on the b-grid:
          sigFrac(bmax) = [2π ∫_0^{bmax} b (1 - e^{-σ_nn T(b)}) db × 10] / σ_tot
        """
        sigma_fm2 = self.sigma_nn_mb * MB_TO_FM2
        T = self.TAA_b if kind == "AA" else self.TpA_b
        integrand = self.b_grid * (1.0 - np.exp(-T * sigma_fm2))
        cum_fm2 = 2.0 * math.pi * np.cumsum(integrand) * self.db_fm
        sig_tot_mb = self.sigma_AA_tot_mb if kind == "AA" else self.sigma_pA_tot_mb
        return (cum_fm2 * FM2_TO_MB) / max(sig_tot_mb, 1e-12)

    # --------------------------- Centrality mapping --------------------------

    def b_from_percentile(self, c: float, kind: Literal["AA", "pA"]) -> float:
        """Return b_max (fm) for a cross-section percentile c ∈ [0,1] via inverse CDF."""
        c = float(np.clip(c, 0.0, 1.0))
        cum = self.cum_AA if kind == "AA" else self.cum_pA
        return float(_interp1_linear(cum, self.b_grid, c))

    def bin_sigma_mb(self, cmin: float, cmax: float, kind: Literal["AA", "pA"]) -> float:
        """σ_bin (mb) between [cmin, cmax]."""
        bmin = self.b_from_percentile(cmin, kind)
        bmax = self.b_from_percentile(cmax, kind)
        sigma_fm2 = self.sigma_nn_mb * MB_TO_FM2
        T = self.TAA_b if kind == "AA" else self.TpA_b

        def w(b: float) -> float:
            Tb = _interp1_linear(self.b_grid, T, b)
            return float(b * (1.0 - math.exp(-Tb * sigma_fm2)))

        val_fm2 = 2.0 * math.pi * _gl_integrate_1d(w, bmin, bmax, n=32)
        return float(val_fm2 * FM2_TO_MB)

    def avg_b(self, cmin: float, cmax: float, kind: Literal["AA", "pA"]) -> float:
        """⟨b⟩ (fm) in a centrality bin using inelastic weights."""
        bmin = self.b_from_percentile(cmin, kind)
        bmax = self.b_from_percentile(cmax, kind)
        sigma_fm2 = self.sigma_nn_mb * MB_TO_FM2
        T = self.TAA_b if kind == "AA" else self.TpA_b

        def num(bb: float) -> float:
            Tb = _interp1_linear(self.b_grid, T, bb)
            return float(bb * bb * (1.0 - np.exp(-Tb * sigma_fm2)))

        def den(bb: float) -> float:
            Tb = _interp1_linear(self.b_grid, T, bb)
            return float(bb * (1.0 - np.exp(-Tb * sigma_fm2)))

        num_fm3 = 2.0 * math.pi * _gl_integrate_1d(num, bmin, bmax, n=32)
        den_fm2 = 2.0 * math.pi * _gl_integrate_1d(den, bmin, bmax, n=32)
        return float(num_fm3 / max(den_fm2, 1e-12))

    def avg_Npart(self, cmin: float, cmax: float, kind: Literal["AA", "pA"]) -> float:
        """⟨N_part⟩ in a centrality bin with the same inelastic weight as σ_bin."""
        bmin = self.b_from_percentile(cmin, kind)
        bmax = self.b_from_percentile(cmax, kind)
        sigma_fm2 = self.sigma_nn_mb * MB_TO_FM2
        T = self.TAA_b if kind == "AA" else self.TpA_b

        def Np(bb: float) -> float:
            return (self.N_part_of_b_AA(bb) if kind == "AA"
                    else self.N_part_of_b_pA(bb))

        def num(bb: float) -> float:
            Tb = _interp1_linear(self.b_grid, T, bb)
            return float(bb * Np(bb) * (1.0 - np.exp(-Tb * sigma_fm2)))

        def den(bb: float) -> float:
            Tb = _interp1_linear(self.b_grid, T, bb)
            return float(bb * (1.0 - np.exp(-Tb * sigma_fm2)))

        num_val = 2.0 * math.pi * _gl_integrate_1d(num, bmin, bmax, n=24)
        den_val = 2.0 * math.pi * _gl_integrate_1d(den, bmin, bmax, n=24)
        return float(num_val / max(den_val, 1e-12))

    # ----------------------------- Public API --------------------------------

    def centrality_table(self,
                         edges: Iterable[float],
                         kind: Literal["AA", "pA"]) -> List[Dict[str, float]]:
        """
        Build centrality table for adjacent percentile edges:

          returns list of dicts with keys:
          {'cmin','cmax','bmin','bmax','b_avg','Npart_avg','Ncoll_avg','sigma_bin_mb'}
        """
        edges = [float(c) for c in edges]
        if not (min(edges) >= 0.0 and max(edges) <= 1.0):
            raise ValueError("centrality edges must lie within [0,1]")
        if any(edges[i] >= edges[i + 1] for i in range(len(edges) - 1)):
            raise ValueError("centrality edges must be strictly increasing")

        out: List[Dict[str, float]] = []
        for i in range(len(edges) - 1):
            cmin, cmax = edges[i], edges[i + 1]
            bmin = self.b_from_percentile(cmin, kind)
            bmax = self.b_from_percentile(cmax, kind)
            b_avg = self.avg_b(cmin, cmax, kind)
            Np_avg = self.avg_Npart(cmin, cmax, kind)

            sigma_fm2 = self.sigma_nn_mb * MB_TO_FM2
            T = self.TAA_b if kind == "AA" else self.TpA_b

            def Ncoll(bb: float) -> float:
                return self.N_coll_of_b(bb, kind=kind)

            def num(bb: float) -> float:
                Tb = _interp1_linear(self.b_grid, T, bb)
                return float(bb * Ncoll(bb) * (1.0 - np.exp(-Tb * sigma_fm2)))

            def den(bb: float) -> float:
                Tb = _interp1_linear(self.b_grid, T, bb)
                return float(bb * (1.0 - np.exp(-Tb * sigma_fm2)))

            num_val = 2.0 * math.pi * _gl_integrate_1d(num, bmin, bmax, n=24)
            den_val = 2.0 * math.pi * _gl_integrate_1d(den, bmin, bmax, n=24)
            Ncoll_avg = float(num_val / max(den_val, 1e-12))
            sig_bin_mb = self.bin_sigma_mb(cmin, cmax, kind)

            out.append(dict(cmin=cmin, cmax=cmax,
                            bmin=bmin, bmax=bmax, b_avg=b_avg,
                            Npart_avg=Np_avg, Ncoll_avg=Ncoll_avg,
                            sigma_bin_mb=sig_bin_mb))
        return out

    def min_bias(self, kind: Literal["AA", "pA"]) -> Dict[str, float]:
        """
        Minimum-bias (0–100%) averages using the same inelastic weight:
          returns dict with keys: cmin,cmax,bmin,bmax,b_avg,Npart_avg,Ncoll_avg,sigma_bin_mb, sigma_tot_mb
        """
        mb_row = self.centrality_table([0.0, 1.0], kind=kind)[0]
        # add explicit σ_tot from trapezoid integral (should match sigma_bin_mb very closely)
        mb_row["sigma_tot_mb"] = self.sigma_AA_tot_mb if kind == "AA" else self.sigma_pA_tot_mb
        return mb_row

    # ----------------------------- Plot helpers ------------------------------

    def plot_summary(self,
                     kind: Literal["AA", "pA"],
                     edges: Iterable[float] = (0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0),
                     show: bool = True,
                     savepath: Optional[str] = None) -> None:
        """
        Reproduce key notebook figures:
          - T_kind(b) vs b
          - N_part(b) and N_coll(b) vs b
          - b vs centrality (edges & ⟨b⟩)
          - ⟨N_part⟩ vs centrality
        """
        edges = list(edges)
        if kind == "AA":
            T_b = self.TAA_b
            title = "AA"
            npart_of_b = np.vectorize(self.N_part_of_b_AA)
        else:
            T_b = self.TpA_b
            title = "pA"
            npart_of_b = np.vectorize(self.N_part_of_b_pA)

        ncoll_of_b = np.vectorize(lambda bb: self.N_coll_of_b(bb, kind=kind))
        tab = self.centrality_table(edges, kind=kind)

        fig, axs = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)

        ax = axs[0, 0]
        ax.plot(self.b_grid, T_b, lw=2)
        ax.set_xlabel("b [fm]")
        ax.set_ylabel(r"$T_{" + title + r"}(b)$ [fm$^{-2}$]")
        ax.set_title(f"{title}: thickness vs impact parameter")

        ax = axs[0, 1]
        ax.plot(self.b_grid, npart_of_b(self.b_grid), lw=2, label=r"$N_\mathrm{part}(b)$")
        ax.plot(self.b_grid, ncoll_of_b(self.b_grid), lw=2, label=r"$N_\mathrm{coll}(b)$")
        ax.set_xlabel("b [fm]")
        ax.set_title(f"{title}: participants & collisions vs b")
        ax.legend(frameon=False)

        # b vs centrality edges and ⟨b⟩
        ax = axs[1, 0]
        c_mid = [0.5 * (t["cmin"] + t["cmax"]) for t in tab]
        b_avg = [t["b_avg"] for t in tab]
        b_lo = [t["bmin"] for t in tab]
        b_hi = [t["bmax"] for t in tab]
        ax.fill_between(c_mid, b_lo, b_hi, alpha=0.2, label="b edges")
        ax.plot(c_mid, b_avg, lw=2, label=r"$\langle b \rangle$")
        ax.set_xlabel("centrality fraction c")
        ax.set_ylabel("b [fm]")
        ax.set_title(f"{title}: b edges and ⟨b⟩ vs centrality")
        ax.legend(frameon=False)

        # ⟨N_part⟩ vs centrality
        ax = axs[1, 1]
        Np_avg = [t["Npart_avg"] for t in tab]
        ax.plot(c_mid, Np_avg, lw=2)
        ax.set_xlabel("centrality fraction c")
        ax.set_ylabel(r"$\langle N_\mathrm{part} \rangle$")
        ax.set_title(f"{title}: ⟨N_part⟩ vs centrality")

        supt = (f"Optical Glauber summary ({title})  •  "
                f"σ_nn={self.sigma_nn_mb:.1f} mb, A={self.spec.A}, d={self.d_fm:.3f} fm")
        fig.suptitle(supt, fontsize=12)

        if savepath:
            plt.savefig(savepath, dpi=200, bbox_inches="tight")
        if show:
            plt.show()

    # ------------------------------ Debug / I/O ------------------------------

    def sanity_report(self, kind: Literal["AA", "pA"]) -> None:
        """Prints normalization checks and spot values."""
        print("\n[Glauber] Sanity report")
        print("-----------------------")
        integ_TA = np.trapezoid(np.trapezoid(self.T_xy, self.y_grid, axis=1),
                                self.x_grid, axis=0)
        print(f"∫ T_A(x,y) d^2x ≈ {integ_TA:.3f} (target A = {self.spec.A})")
        print(f"σ_tot^AA = {self.sigma_AA_tot_mb:.2f} mb")
        print(f"σ_tot^pA = {self.sigma_pA_tot_mb:.2f} mb")
        for b in (0.0, 5.0, 10.0):
            if kind == "AA":
                npart = self.N_part_of_b_AA(b)
                ncoll = self.N_coll_of_b(b, kind="AA")
                T = _interp1_linear(self.b_grid, self.TAA_b, b)
            else:
                npart = self.N_part_of_b_pA(b)
                ncoll = self.N_coll_of_b(b, kind="pA")
                T = _interp1_linear(self.b_grid, self.TpA_b, b)
            print(f"b={b:4.1f} fm:  T(b)={T:.5f} fm^-2,  N_part={npart:.3f},  N_coll={ncoll:.3f}")
        print()

    def min_bias_report(self, kind: Literal["AA", "pA"]) -> None:
        """Print a compact minimum-bias (0–100%) summary."""
        mb = self.min_bias(kind)
        cstr = f"{mb['cmin']:.2f}-{mb['cmax']:.2f}"
        print("[Glauber] Minimum-bias (0–100%)")
        print("------------------------------")
        print(f" cmin-cmax: {cstr}")
        print(f" bmin-bmax: {mb['bmin']:.2f}–{mb['bmax']:.2f} fm;  <b> = {mb['b_avg']:.3f} fm")
        print(f" <Npart>  : {mb['Npart_avg']:.3f}")
        print(f" <Ncoll>  : {mb['Ncoll_avg']:.3f}")
        print(f" σ_bin    : {mb['sigma_bin_mb']:.2f} mb")
        print(f" σ_tot    : {mb['sigma_tot_mb']:.2f} mb\n")

    def debug_print(self, kind: Literal["AA", "pA"],
                    c_edges: Iterable[float] = (0, .1, .2, .4, .6, .8, 1.0)) -> None:
        """Print table mirroring the Mathematica grid output."""
        tab = self.centrality_table(c_edges, kind=kind)
        hdr = (f"{'cmin-cmax':>12} | {'bmin':>6} {'bmax':>6} {'<b>':>6} | "
               f"{'<Npart>':>9} {'<Ncoll>':>9} | {'σ_bin [mb]':>10}")
        print(hdr)
        print("-" * len(hdr))
        for t in tab:
            cstr = f"{t['cmin']:.2f}-{t['cmax']:.2f}"
            print(f"{cstr:>12} | {t['bmin']:6.2f} {t['bmax']:6.2f} {t['b_avg']:6.2f} | "
                  f"{t['Npart_avg']:9.3f} {t['Ncoll_avg']:9.3f} | {t['sigma_bin_mb']:10.2f}")

    def export_tsv(self,
                   outdir: str,
                   kind: Literal["AA", "pA"],
                   cdf_step: float = 1e-3) -> None:
        """
        Export TSVs analogous to your Mathematica outputs:
          - bvscData.tsv: {c, bmax(c)}
          - npartvsbData.tsv: {b, N_part(b)}  [pA: N_part^pA; AA: N_part^AA]
          - nbinvsbData.tsv: {b, N_coll(b)}   [i.e., σ_nn/10 * T_kind(b)]
        """
        os.makedirs(outdir, exist_ok=True)
        # b vs N_part and N_coll
        b = self.b_grid
        if kind == "AA":
            npart = np.array([self.N_part_of_b_AA(bb) for bb in b])
            T = self.TAA_b
        else:
            npart = np.array([self.N_part_of_b_pA(bb) for bb in b])
            T = self.TpA_b
        ncoll = (self.sigma_nn_mb * MB_TO_FM2) * T

        np.savetxt(os.path.join(outdir, "npartvsbData.tsv"),
            np.column_stack([b, npart]), fmt="%.6f", delimiter="\t")
        np.savetxt(os.path.join(outdir, "nbinvsbData.tsv"),
            np.column_stack([b, ncoll]), fmt="%.6f", delimiter="\t")

        # c vs bmax(c)
        cs = np.arange(0.0, 1.0 + 1e-12, cdf_step)
        bmax = np.array([self.b_from_percentile(float(c), kind) for c in cs])
        np.savetxt(os.path.join(outdir, "bvscData.tsv"),
            np.column_stack([cs, bmax]), fmt="%.6f", delimiter="\t")

# ------------------------------ CLI / main ----------------------------------


def _suggest_A_and_sigma(roots: float, system: str) -> Tuple[int, float]:
    """
    Sensible defaults:
      - AuAu@200 GeV: A=197, σ_nn from map
      - PbPb / pPb at LHC: A=208, σ_nn from map
    """
    if system in ("AA", "pA"):
        A = 197 if abs(roots - 200.0) < 1e-6 else 208
        sigma = _SIGMA_NN_BY_ROOTS_MB.get(roots, 67.6)
        return A, sigma
    raise ValueError("system must be 'AA' or 'pA'.")


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="Optical Glauber (AA / pA) — CMS style")
    ap.add_argument("--system", type=str, default="pA", choices=["AA", "pA"],
                    help="Collision system kind.")
    ap.add_argument("--roots", type=float, default=5023.0,
                    help="√s_NN in GeV. Presets: 200, 2760, 5023, 8160.")
    ap.add_argument("--A", type=int, default=None,
                    help="Mass number of nucleus (default: 197 at 200 GeV, else 208).")
    ap.add_argument("--sigma-nn", type=float, default=None,
                    help="σ_nn in mb (defaults by √s).")
    ap.add_argument("--d", type=float, default=None,
                    help="Woods–Saxon diffuseness d in fm (defaults by √s).")
    ap.add_argument("--edges", type=float, nargs="+",
                    default=[0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0],
                    help="Centrality edges (fractions in [0,1]).")
    ap.add_argument("--export", type=str, default=None,
                    help="Directory to export bvsc / npartvsb / nbinvsb TSVs.")
    ap.add_argument("--no-show", action="store_true",
                    help="Do not show figures (still computes).")
    ap.add_argument("--save", type=str, default=None,
                    help="Save summary figure to this path.")
    ap.add_argument("--verbose", action="store_true",
                    help="Verbose progress prints.")
    args = ap.parse_args()

    # Defaults for A and σ_nn if not provided
    if args.A is None or args.sigma_nn is None:
        A_def, sig_def = _suggest_A_and_sigma(args.roots, args.system)
        A = A_def if args.A is None else args.A
        sigma_mb = sig_def if args.sigma_nn is None else args.sigma_nn
    else:
        A, sigma_mb = args.A, args.sigma_nn

    spec = SystemSpec(system=args.system, roots_GeV=args.roots, A=A,
                      sigma_nn_mb=sigma_mb, diffuseness_fm=args.d)

    gl = OpticalGlauber(spec, verbose=args.verbose)
    gl.sanity_report(kind=args.system)

    # Minimum-bias quick summary
    gl.min_bias_report(kind=args.system)

    # Centrality table + plots
    edges = tuple(args.edges)
    gl.debug_print(kind=args.system, c_edges=edges)
    gl.plot_summary(kind=args.system, edges=edges, show=not args.no_show, savepath=args.save)

    # Optional TSV export matching Mathematica
    if args.export:
        gl.export_tsv(args.export, kind=args.system)


if __name__ == "__main__":
    main()
