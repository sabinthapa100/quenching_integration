# quenching_fast_parts.py (final, robust)
# -----------------------------------------------------------------------------
# Factorized R_pA(y,pT) into energy-loss and pT-broadening pieces
# using the same AP kernel and α_s policy as in quenching_fast, but with a
# *safe* σ_pp accessor that avoids any Torch GPU indexing pitfalls.
#
# R_loss(y,pT)  = ∫ d(δy) P̂_A(z) · [σ_pp(y+δy,pT) / σ_pp(y,pT)]
# R_broad(y,pT) = ⟨ σ_pp(y, |p⃗_T − Δp⃗_T^A|) / σ_pp(y,pT) ⟩_φ
# R_fact ≈ R_loss · R_broad
#
# Design:
#  • Uses quenching_fast.PhatA_t and ΔpT tapering for physics consistency.
#  • σ_pp is accessed through SafeSigmaPP (bilinear over Torch table or plain callable).
#  • Deterministic Gauss–Legendre in δy and φ, edge-stable near δy_max→0 and L→ℓ_p.
#  • Includes simple binning utilities and quick plotting helpers.
# -----------------------------------------------------------------------------
from __future__ import annotations
from dataclasses import replace
from typing import Callable, Optional, Tuple, Literal
import math, numpy as np

# ---- import core physics kernels from your fast module ----
from quenching_fast import (
    QuenchParams,
    TorchSigmaPPTable,
    PhatA_t,              # AP kernel (torch-based)
    _dpt_from_xL_t,       # tapered ΔpT (torch-based)
    y_max,                # kinematics
)

# Torch (only for AP kernel eval; σ_pp sampling is CPU-safe here)
try:
    import torch
    _HAS_TORCH = True
except Exception:
    _HAS_TORCH = False

# numerics / consts ------------------------------------------------------------
M_PROTON = 0.938
Z_FLOOR  = 1e-12
DY_EPS   = 1e-6
LOG2     = math.log(2.0)

# -----------------------------------------------------------------------------
# GL nodes (NumPy) and φ-nodes
# -----------------------------------------------------------------------------
def _gl_nodes_np(a: float, b: float, n: int):
    x, w = np.polynomial.legendre.leggauss(int(n))
    xm, xc = 0.5*(b-a), 0.5*(b+a)
    return xc + xm*x, xm*w

# -----------------------------------------------------------------------------
# σ_pp accessor – *robust* CPU-side bilinear interpolant
# -----------------------------------------------------------------------------
# --- replace your SafeSigmaPP class with this version ---

class SafeSigmaPP:
    def __init__(self, src):
        self.is_grid = False
        self.f = None

        # Torch-backed table
        if _HAS_TORCH and isinstance(src, TorchSigmaPPTable):
            self.y  = src.y.detach().cpu().numpy().astype(float)
            self.pt = src.pt.detach().cpu().numpy().astype(float)
            self.Z  = src.Z.detach().cpu().numpy().astype(float)
            self.is_grid = True
            return

        # Generic callable (SigmaPPTable or function)
        if callable(src):
            self.f = src
            return

        # Helpful message for common mistake
        if isinstance(src, (int, float)):
            raise TypeError(
                "Got a scalar (e.g. sigma_nn) where a pp double-differential "
                "cross-section is required. Pass TorchSigmaPPTable or a "
                "callable σ_pp(y,pT[,roots])."
            )

        raise TypeError("Unsupported σ_pp source; need TorchSigmaPPTable or callable")

    # internal: call a 2-arg or 3-arg callable safely
    def _call_sigma_callable(self, y: float, pT: float, roots_GeV: float) -> float:
        # Try (y,pT,roots), fallback to (y,pT)
        try:
            return float(self.f(y, pT, roots_GeV))
        except TypeError:
            return float(self.f(y, pT))

    def __call__(self, y: float, pT: float, roots_GeV: float) -> float:
        if self.is_grid:
            # CPU-side bilinear interpolation with interior clamping
            yv = float(np.clip(y,  self.y[0],  self.y[-1]  - 1e-12))
            pv = float(np.clip(pT, self.pt[0], self.pt[-1] - 1e-12))
            i = int(np.clip(np.searchsorted(self.y,  yv) - 1, 0, self.y.size  - 2))
            j = int(np.clip(np.searchsorted(self.pt, pv) - 1, 0, self.pt.size - 2))
            y0, y1 = self.y[i],  self.y[i+1]
            p0, p1 = self.pt[j], self.pt[j+1]
            ty = 0.0 if y1 == y0 else (yv - y0)/(y1 - y0)
            tp = 0.0 if p1 == p0 else (pv - p0)/(p1 - p0)
            z00 = self.Z[i, j];   z01 = self.Z[i, j+1]
            z10 = self.Z[i+1, j]; z11 = self.Z[i+1, j+1]
            z0 = (1.0 - tp)*z00 + tp*z01
            z1 = (1.0 - tp)*z10 + tp*z11
            return float((1.0 - ty)*z0 + ty*z1)
        else:
            return self._call_sigma_callable(y, pT, roots_GeV)

    def vec_y(self, y_vec, pT: float, roots_GeV: float):
        y_arr = np.asarray(y_vec, float)
        return np.array([self(y, pT, roots_GeV) for y in y_arr], float)

    def vec_p(self, y: float, p_vec, roots_GeV: float):
        p_arr = np.asarray(p_vec, float)
        return np.array([self(y, p, roots_GeV) for p in p_arr], float)


def make_safe_sigma_pp(src):
    return src if isinstance(src, SafeSigmaPP) else SafeSigmaPP(src)



def make_safe_sigma_pp(src):
    return src if isinstance(src, SafeSigmaPP) else SafeSigmaPP(src)

# -----------------------------------------------------------------------------
# small-x helper for coherence bound
# -----------------------------------------------------------------------------
def xA0_from_L(L_fm: float) -> float:
    return 1.0/(2.0*M_PROTON*max(L_fm, 1e-12))

def xA_for_side_A(LA_fm: float, mT: float, roots_GeV: float, y: float) -> float:
    return min(xA0_from_L(LA_fm), (mT/roots_GeV) * math.exp(-y))

# -----------------------------------------------------------------------------
# Adaptive Ny (matches quenching_fast logic)
# -----------------------------------------------------------------------------
def _Ny_from_dymax(dym: float) -> int:
    if dym < 0.02:  return 64
    if dym < 0.05:  return 48
    if dym < 0.10:  return 40
    return 32

# -----------------------------------------------------------------------------
# Energy-loss-only piece (robust; σ_pp via SafeSigmaPP, AP via torch)
# -----------------------------------------------------------------------------
def R_loss_pA(
    *, P, roots_GeV: float, qp: QuenchParams,
    y: float, pT: float, table_or_callable,
    Ny: Optional[int] = None, mapping: Literal["exp","linear"] = "exp",
    use_torch: bool = True,
) -> float:
    mT = float(P.mT(pT))
    dym = max(0.0, min(LOG2, y_max(roots_GeV, mT) - y))
    if dym <= DY_EPS:
        return 1.0

    # σ_pp accessor (CPU-safe)
    sigpp = make_safe_sigma_pp(table_or_callable)

    # Denominator σ_pp(y,pT)
    sig_den = sigpp(y, pT, roots_GeV)
    if not (sig_den > 0.0 and math.isfinite(sig_den)):
        return 1.0

    # δy nodes
    Ny_loc = _Ny_from_dymax(dym) if Ny is None else int(Ny)
    if mapping == "exp":
        u, wu = _gl_nodes_np(-30.0, math.log(max(dym, 1e-300)), Ny_loc)
        dy = np.exp(u)
        z  = np.expm1(dy).clip(min=Z_FLOOR)
        # AP kernel via torch per node (stable & consistent with quenching_fast)
        Ph = []
        if _HAS_TORCH:
            for zi in z:
                Ph.append(float(PhatA_t(torch.tensor([zi],dtype=torch.float64), mT,
                                         torch.tensor([xA_for_side_A(qp.LA_fm, mT, roots_GeV, y)],dtype=torch.float64),
                                         qp, pT=pT)[0].item()))
        else:
            # should not happen in your setup; keep graceful fallback
            for zi in z:
                Ph.append(0.0)
        Ph = np.asarray(Ph, float)
        # σ ratio
        sig_num = sigpp.vec_y(y + dy, pT, roots_GeV)
        ratio   = sig_num / max(sig_den, 1e-300)
        return float(np.sum(wu * np.exp(u) * Ph * ratio))
    else:
        dy, wy = _gl_nodes_np(0.0, dym, Ny_loc)
        z  = np.expm1(dy).clip(min=Z_FLOOR)
        Ph = []
        if _HAS_TORCH:
            for zi in z:
                Ph.append(float(PhatA_t(torch.tensor([zi],dtype=torch.float64), mT,
                                         torch.tensor([xA_for_side_A(qp.LA_fm, mT, roots_GeV, y)],dtype=torch.float64),
                                         qp, pT=pT)[0].item()))
        else:
            for zi in z:
                Ph.append(0.0)
        Ph = np.asarray(Ph, float)
        sig_num = sigpp.vec_y(y + dy, pT, roots_GeV)
        ratio   = sig_num / max(sig_den, 1e-300)
        return float(np.sum(wy * Ph * ratio))

# -----------------------------------------------------------------------------
# Broadening-only piece (robust; σ_pp via SafeSigmaPP, ΔpT via torch)
# -----------------------------------------------------------------------------
def R_broad_pA(
    *, P, roots_GeV: float, qp: QuenchParams,
    y: float, pT: float, table_or_callable,
    Nphi: int = 32, use_torch: bool = True,
) -> float:
    mT = float(P.mT(pT))
    sigpp = make_safe_sigma_pp(table_or_callable)
    sig_den = sigpp(y, pT, roots_GeV)
    if not (sig_den > 0.0 and math.isfinite(sig_den)):
        return 1.0

    # ΔpT^A with taper (torch)
    if _HAS_TORCH:
        xA = xA_for_side_A(qp.LA_fm, mT, roots_GeV, y)
        dpta = float(_dpt_from_xL_t(qp, torch.tensor([xA],dtype=torch.float64), qp.LA_fm,
                                    hard=qp.use_hard_cronin)[0].item())
    else:
        dpta = 0.0

    # φ-average on CPU
    phi, wphi = _gl_nodes_np(0.0, 2.0*math.pi, int(Nphi))
    pphi = np.sqrt((pT - dpta*np.cos(phi))**2 + (dpta*np.sin(phi))**2)
    sig  = sigpp.vec_p(y, pphi, roots_GeV)
    return float(np.sum((sig / max(sig_den, 1e-300)) * (wphi/(2.0*math.pi))))

# -----------------------------------------------------------------------------
# Combined (product) – returns tuple
# -----------------------------------------------------------------------------
def R_factored_pA(
    *, P, roots_GeV: float, qp: QuenchParams,
    y: float, pT: float, table_or_callable,
    Ny: Optional[int] = None, Nphi: int = 32,
    mapping: Literal["exp","linear"] = "exp",
    use_torch: bool = True,
) -> tuple[float, float, float]:
    Rl = R_loss_pA(P=P, roots_GeV=roots_GeV, qp=qp, y=y, pT=pT,
                   table_or_callable=table_or_callable, Ny=Ny,
                   mapping=mapping, use_torch=use_torch)
    Rb = R_broad_pA(P=P, roots_GeV=roots_GeV, qp=qp, y=y, pT=pT,
                    table_or_callable=table_or_callable, Nphi=Nphi,
                    use_torch=use_torch)
    return Rl, Rb, Rl*Rb

# -----------------------------------------------------------------------------
# Convenience scans vs pT and vs y
# -----------------------------------------------------------------------------
def parts_vs_pT(
    *, P, roots_GeV: float, qp: QuenchParams, table_or_callable,
    p_grid: np.ndarray, y: float,
    Ny_loss: Optional[int] = None, Nphi_broad: int = 32,
    mapping: Literal["exp","linear"] = "exp",
    use_torch: bool = True,
):
    Rl, Rb, Rf = [], [], []
    for p in p_grid:
        rL = R_loss_pA(P=P, roots_GeV=roots_GeV, qp=qp, y=float(y), pT=float(p),
                       table_or_callable=table_or_callable, Ny=Ny_loss,
                       mapping=mapping, use_torch=use_torch)
        rB = R_broad_pA(P=P, roots_GeV=roots_GeV, qp=qp, y=float(y), pT=float(p),
                        table_or_callable=table_or_callable, Nphi=Nphi_broad,
                        use_torch=use_torch)
        Rl.append(rL); Rb.append(rB); Rf.append(rL*rB)
    return np.array(Rl,float), np.array(Rb,float), np.array(Rf,float)


def parts_vs_y(
    *, P, roots_GeV: float, qp: QuenchParams, table_or_callable,
    y_grid: np.ndarray, pT: float,
    Ny_loss: Optional[int] = None, Nphi_broad: int = 32,
    mapping: Literal["exp","linear"] = "exp",
    use_torch: bool = True,
):
    Rl, Rb, Rf = [], [], []
    for y in y_grid:
        rL = R_loss_pA(P=P, roots_GeV=roots_GeV, qp=qp, y=float(y), pT=float(pT),
                       table_or_callable=table_or_callable, Ny=Ny_loss,
                       mapping=mapping, use_torch=use_torch)
        rB = R_broad_pA(P=P, roots_GeV=roots_GeV, qp=qp, y=float(y), pT=float(pT),
                        table_or_callable=table_or_callable, Nphi=Nphi_broad,
                        use_torch=use_torch)
        Rl.append(rL); Rb.append(rB); Rf.append(rL*rB)
    return np.array(Rl,float), np.array(Rb,float), np.array(Rf,float)

# -----------------------------------------------------------------------------
# Binning utilities (y,pT) – parts averaged over a window
# -----------------------------------------------------------------------------
def _weights_over_pT(
    *, P, roots_GeV: float, qp: QuenchParams, table_or_callable,
    p_nodes: np.ndarray, y_wref: float,
    wkind: Literal["flat","pp","pA"] = "pp",
    Nphi_weight: int = 12,
    use_torch: bool = True,
) -> np.ndarray:
    sigpp = make_safe_sigma_pp(table_or_callable)
    if wkind == "flat":
        return np.ones_like(p_nodes, float)
    if wkind == "pp":
        sig = sigpp.vec_p(float(y_wref), p_nodes, roots_GeV)
        return np.asarray(sig, float) * np.maximum(p_nodes, 1e-12)
    # pA weight ≈ σ_pp · R_broad (fast, δy-independent)
    out = []
    for pj in p_nodes:
        s0 = sigpp(float(y_wref), float(pj), roots_GeV)
        rb = R_broad_pA(P=P, roots_GeV=roots_GeV, qp=qp, y=float(y_wref), pT=float(pj),
                        table_or_callable=sigpp, Nphi=Nphi_weight, use_torch=use_torch)
        out.append(rb * s0 * max(float(pj), 1e-12))
    return np.array(out, float)


def rpa_parts_binned(
    *, P, roots_GeV: float, qp: QuenchParams, table_or_callable,
    y_range: Tuple[float,float], pt_range: Tuple[float,float],
    Ny_bin: int = 24, Npt_bin: int = 48,
    weight_kind: Literal["flat","pp","pA","auto"] = "auto",
    weight_ref_y: Optional[float] = 0.0,
    Ny_loss: Optional[int] = None, Nphi_broad: int = 24,
    mapping: Literal["exp","linear"] = "exp",
    use_torch: bool = True,
) -> dict[str, float]:
    yl, yr = map(float, y_range); pl, pr = map(float, pt_range)
    y_nodes, y_w = _gl_nodes_np(yl, yr, Ny_bin)
    p_nodes, p_w = _gl_nodes_np(pl, pr, Npt_bin)
    wkind = ("pA" if weight_kind == "auto" else weight_kind)
    y_wref = float(0.0 if weight_ref_y is None else weight_ref_y)
    Wp = _weights_over_pT(P=P, roots_GeV=roots_GeV, qp=qp, table_or_callable=table_or_callable,
                          p_nodes=np.asarray(p_nodes,float), y_wref=y_wref,
                          wkind=wkind if wkind in {"pp","pA"} else "pp",
                          Nphi_weight=12, use_torch=use_torch)
    accL = accB = accF = accW = 0.0
    for yi, wy in zip(y_nodes, y_w):
        for pj, wp, wj in zip(p_nodes, p_w, Wp):
            rL, rB, rF = R_factored_pA(P=P, roots_GeV=roots_GeV, qp=qp, y=float(yi), pT=float(pj),
                                       table_or_callable=table_or_callable,
                                       Ny=Ny_loss, Nphi=Nphi_broad, mapping=mapping, use_torch=use_torch)
            w = wy * wp * float(wj)
            accL += w * rL; accB += w * rB; accF += w * rF; accW += w
    if accW <= 0:
        return dict(R_loss=accL, R_broad=accB, R_fact=accF)
    return dict(R_loss=accL/accW, R_broad=accB/accW, R_fact=accF/accW)

# -----------------------------------------------------------------------------
# Centrality-binned (optical Glauber): one L_eff per bin, parts averaged
# -----------------------------------------------------------------------------
_MB_TO_FM2 = 0.1

def _optical_bin_weight_pA(glauber, c0_percent: float, c1_percent: float, n_sub: int = 1200) -> float:
    bmin = float(glauber.b_from_percentile(c0_percent/100.0, kind="pA"))
    bmax = float(glauber.b_from_percentile(c1_percent/100.0, kind="pA"))
    if bmax <= bmin:
        return 0.0
    b_sub = np.linspace(bmin, bmax, n_sub)
    TpA_sub = np.interp(b_sub, np.asarray(glauber.b_grid, float), np.asarray(glauber.TpA_b, float))
    sigma_fm2 = float(glauber.spec.sigma_nn_mb) * _MB_TO_FM2
    pinel = 1.0 - np.exp(-sigma_fm2 * np.maximum(TpA_sub, 0.0))
    integrand = 2.0 * math.pi * b_sub * pinel
    numer_fm2 = float(np.trapezoid(integrand, b_sub))
    sigma_tot_fm2 = float(glauber.sigma_pA_tot_mb) * _MB_TO_FM2
    return numer_fm2 / max(sigma_tot_fm2, 1e-30)


def rpa_parts_centrality_binned(
    *, P, roots_GeV: float, qp_base: QuenchParams, table_or_callable,
    glauber, cent_edges_or_bins,
    y_range: Tuple[float,float], pt_range: Tuple[float,float],
    Ny_bin: int = 18, Npt_bin: int = 36,
    weight_kind: Literal["flat","pp","pA","auto"] = "auto",
    weight_ref_y: Optional[float] = 0.0,
    Ny_loss: Optional[int] = None, Nphi_broad: int = 24,
    mapping: Literal["exp","linear"] = "exp",
    use_torch: bool = True,
) -> dict:
    arr = np.asarray(cent_edges_or_bins, dtype=float)
    if arr.ndim == 1:
        edges = np.unique(np.sort(arr))
        bins_list = [(int(edges[i]), int(edges[i+1])) for i in range(edges.size-1)]
    else:
        bins_list = [(int(a), int(b)) for (a,b) in cent_edges_or_bins]
    labels  = [f"{a}-{b}%" for (a,b) in bins_list]
    centers = 0.5*np.array([a+b for (a,b) in bins_list], float)

    w = np.array([_optical_bin_weight_pA(glauber, a, b) for (a, b) in bins_list], float)
    w = w / max(w.sum(), 1e-30)

    try:
        L_by = glauber.leff_bins_pA(bins_list, method="optical")
    except TypeError:
        L_by = glauber.leff_bins_pA(bins_list)
    Leff = np.array([float(L_by[tag]) for tag in labels], float)

    Rloss, Rbroad, Rfact = [], [], []
    for L in Leff:
        qp = replace(qp_base, LA_fm=float(L))
        parts = rpa_parts_binned(P=P, roots_GeV=roots_GeV, qp=qp, table_or_callable=table_or_callable,
                                 y_range=y_range, pt_range=pt_range,
                                 Ny_bin=Ny_bin, Npt_bin=Npt_bin,
                                 weight_kind=weight_kind, weight_ref_y=weight_ref_y,
                                 Ny_loss=Ny_loss, Nphi_broad=Nphi_broad,
                                 mapping=mapping, use_torch=use_torch)
        Rloss.append(parts['R_loss']); Rbroad.append(parts['R_broad']); Rfact.append(parts['R_fact'])

    Rloss = np.asarray(Rloss, float); Rbroad = np.asarray(Rbroad, float); Rfact = np.asarray(Rfact, float)
    Rmb_loss  = float(np.sum(w * Rloss))
    Rmb_broad = float(np.sum(w * Rbroad))
    Rmb_fact  = float(np.sum(w * Rfact))

    return dict(
        centers=centers.astype(float),
        labels=np.array(labels),
        Leff=Leff.astype(float),
        weights=w.astype(float),
        R_loss=Rloss, R_broad=Rbroad, R_fact=Rfact,
        R_minbias_loss=Rmb_loss,
        R_minbias_broad=Rmb_broad,
        R_minbias_fact=Rmb_fact,
    )

# -----------------------------------------------------------------------------
# Quick plotting helpers (optional)
# -----------------------------------------------------------------------------
def quick_plot_parts_vs_pT(ax, p_grid, Rl, Rb, Rf, label_prefix=""):
    ax.plot(p_grid, Rl, label=(label_prefix+" loss" if label_prefix else "loss"))
    ax.plot(p_grid, Rb, label=(label_prefix+" broad" if label_prefix else "broad"))
    ax.plot(p_grid, Rf, label=(label_prefix+" both" if label_prefix else "both"))
    ax.set_xlabel(r"$p_T$ [GeV]"); ax.set_ylabel(r"$R_{pA}$")


def quick_plot_parts_vs_y(ax, y_grid, Rl, Rb, Rf, label_prefix=""):
    ax.plot(y_grid, Rl, label=(label_prefix+" loss" if label_prefix else "loss"))
    ax.plot(y_grid, Rb, label=(label_prefix+" broad" if label_prefix else "broad"))
    ax.plot(y_grid, Rf, label=(label_prefix+" both" if label_prefix else "both"))
    ax.set_xlabel(r"$y$"); ax.set_ylabel(r"$R_{pA}$")
