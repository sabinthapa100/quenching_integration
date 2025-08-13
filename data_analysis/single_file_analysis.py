# single_file_analysis.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from data_handler import DataFile


Number = Union[int, float]
PathLike = Union[str, Path]


def _nearest_value(array: Sequence[Number], target: Number) -> Number:
    arr = np.asarray(array, dtype=float)
    idx = int(np.nanargmin(np.abs(arr - float(target))))
    return float(arr[idx])


@dataclass
class SingleFileAnalyzer:
    """Analyze and plot **one** TSV/CSV file (one centrality class).

    Parameters
    ----------
    data_file : DataFile
        Loaded DataFile object (from `data_handler.DataFile.from_path`).
    """

    data_file: DataFile

    @classmethod
    def from_path(cls, path: PathLike, value_floor: float = 1e-29) -> "SingleFileAnalyzer":
        return cls(DataFile.from_path(path, value_floor=value_floor))

    @property
    def df(self) -> pd.DataFrame:
        return self.data_file.data

    @property
    def meta(self):
        return self.data_file.meta

    # -------------------------
    # Summary helpers
    # -------------------------
    def summary(self, mode: str = "auto") -> dict:
        m, e, n = self.data_file.mean(mode=mode)
        return dict(mean=m, mean_error=e, n_points=n, meta=self.meta)

    def slice_vs_y(self, pt: Optional[Number] = None) -> pd.DataFrame:
        """Return a tidy slice at a given pT (nearest if not exact)."""
        df = self.df
        if pt is None:
            # choose mid pT
            pt = float(np.nanmedian(df["pt"])) if len(df) else np.nan
        pt_levels = np.sort(df["pt"].unique())
        pt_used = _nearest_value(pt_levels, pt)
        sub = df[df["pt"] == pt_used].sort_values("y")
        return sub.assign(pt_requested=pt, pt_used=pt_used).reset_index(drop=True)

    def slice_vs_pt(self, y: Optional[Number] = None) -> pd.DataFrame:
        """Return a tidy slice at a given rapidity y (nearest if not exact)."""
        df = self.df
        if y is None:
            y = float(np.nanmedian(df["y"])) if len(df) else np.nan
        y_levels = np.sort(df["y"].unique())
        y_used = _nearest_value(y_levels, y)
        sub = df[df["y"] == y_used].sort_values("pt")
        return sub.assign(y_requested=y, y_used=y_used).reset_index(drop=True)

    # -------------------------
    # Plots (always 1 chart per function)
    # -------------------------
    def plot_heatmap(self, log_values: bool = False, ax: Optional[plt.Axes] = None):
        """Heatmap of VALUES over (y, pT).

        Parameters
        ----------
        log_values : bool
            If True, plot log(value). Non-positive values are masked.
        ax : Optional[plt.Axes]
        """
        grid = self.data_file.to_grid()
        y_vals = grid.index.to_numpy(dtype=float)
        pt_vals = grid.columns.to_numpy(dtype=float)
        Z = grid.to_numpy(dtype=float)

        if log_values:
            with np.errstate(divide="ignore", invalid="ignore"):
                Z = np.where(Z > 0, np.log(Z), np.nan)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        # imshow expects [xmin, xmax, ymin, ymax]
        extent = [pt_vals.min(), pt_vals.max(), y_vals.min(), y_vals.max()]
        im = ax.imshow(
            Z,
            extent=extent,
            aspect="auto",
            origin="lower",
            interpolation="nearest",
        )
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("log(value)" if log_values else "value")

        ax.set_xlabel("pT [GeV]")
        ax.set_ylabel("y")
        title = f"{self.meta.get('particle')} {self.meta.get('kind')} — {self.meta.get('tag')}"
        ax.set_title(title)
        ax.grid(True, which="both", linestyle=":", alpha=0.3)
        return fig, ax

    def plot_y_curve(self, pt: Optional[Number] = None, with_errorband: bool = True, yscale: str = "linear", ax: Optional[plt.Axes] = None):
        """Plot value(y) at a chosen pT (nearest level used)."""
        s = self.slice_vs_y(pt=pt)
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        ax.plot(s["y"], s["value"], "-", label=f"pT={s['pt_used'].iloc[0]:g} GeV")
        if with_errorband and "error" in s:
            err = s["error"].to_numpy(dtype=float)
            if np.isfinite(err).any():
                v = s["value"].to_numpy(dtype=float)
                ylow = v - err
                yhigh = v + err
                ax.fill_between(s["y"], ylow, yhigh, alpha=0.2)

        ax.set_xlabel("y")
        ax.set_ylabel(self.meta.get("kind", "value"))
        ax.set_title(f"{self.meta.get('particle')} {self.meta.get('kind')} vs y — {self.meta.get('tag')}")
        ax.set_yscale(yscale)
        ax.grid(True, which="both", linestyle=":")
        ax.legend()
        return fig, ax, s

    def plot_pt_curve(self, y: Optional[Number] = None, with_errorband: bool = True, yscale: str = "linear", ax: Optional[plt.Axes] = None):
        """Plot value(pT) at a chosen y (nearest level used)."""
        s = self.slice_vs_pt(y=y)
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        ax.plot(s["pt"], s["value"], "-", label=f"y={s['y_used'].iloc[0]:g}")
        if with_errorband and "error" in s:
            err = s["error"].to_numpy(dtype=float)
            if np.isfinite(err).any():
                v = s["value"].to_numpy(dtype=float)
                ylow = v - err
                yhigh = v + err
                ax.fill_between(s["pt"], ylow, yhigh, alpha=0.2)

        ax.set_xlabel("pT [GeV]")
        ax.set_ylabel(self.meta.get("kind", "value"))
        ax.set_title(f"{self.meta.get('particle')} {self.meta.get('kind')} vs pT — {self.meta.get('tag')}")
        ax.set_yscale(yscale)
        ax.grid(True, which="both", linestyle=":")
        ax.legend()
        return fig, ax, s

    def plot_multi_y_curves(self, pts: Sequence[Number], with_errorband: bool = False, ax: Optional[plt.Axes] = None):
        """Plot value(y) for several pT choices on one axis (one chart)."""
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        for p in pts:
            s = self.slice_vs_y(pt=p)
            ax.plot(s["y"], s["value"], label=f"pT≈{s['pt_used'].iloc[0]:g}")
        ax.set_xlabel("y")
        ax.set_ylabel(self.meta.get("kind", "value"))
        ax.set_title(f"{self.meta.get('particle')} {self.meta.get('kind')} vs y — {self.meta.get('tag')}")
        ax.grid(True, which="both", linestyle=":")
        ax.legend(title="slices")
        return fig, ax

    def plot_multi_pt_curves(self, ys: Sequence[Number], with_errorband: bool = False, ax: Optional[plt.Axes] = None):
        """Plot value(pT) for several y choices on one axis (one chart)."""
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        for y in ys:
            s = self.slice_vs_pt(y=y)
            ax.plot(s["pt"], s["value"], label=f"y≈{s['y_used'].iloc[0]:g}")
        ax.set_xlabel("pT [GeV]")
        ax.set_ylabel(self.meta.get("kind", "value"))
        ax.set_title(f"{self.meta.get('particle')} {self.meta.get('kind')} vs pT — {self.meta.get('tag')}")
        ax.grid(True, which="both", linestyle=":")
        ax.legend(title="slices")
        return fig, ax