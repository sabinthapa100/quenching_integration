# data_plottting.py
from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt

from data_handler import DataFile, scan_output_dir

__all__ = ["rpA_vs_centrality", "rpA_vs_y", "rpA_vs_pt"]

def rpA_vs_centrality(output_dir: Union[str, 'Path'], particle: str = "JPsi", centrality_to_npart: Optional[Dict[Tuple[int, int], float]] = None, show_errorbars: bool = True, ax: Optional[plt.Axes] = None):
    files = scan_output_dir(output_dir, particle=particle, kind="RpA")
    if not files:
        raise FileNotFoundError(f"No RpA files found for particle={particle} in {output_dir}")

    data = []
    for f in files:
        mean, err, n = f.mean(mode="linear")
        if f.centrality is None:
            continue
        if centrality_to_npart:
            x = float(centrality_to_npart.get(f.centrality, np.nan))
        else:
            low, high = f.centrality
            x = 0.5 * (low + high)
        data.append(dict(centrality=f.centrality, x=x, mean=mean, err=err, n=n))

    data.sort(key=lambda d: d["x"])
    xs = [d["x"] for d in data]
    ys = [d["mean"] for d in data]
    es = [d["err"] if d["err"] is not None else 0.0 for d in data]

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    if show_errorbars and any(e > 0 for e in es):
        ax.errorbar(xs, ys, yerr=es, fmt="o-", capsize=3)
    else:
        ax.plot(xs, ys, "o-")

    ax.set_xlabel("<N_part>" if centrality_to_npart else "Centrality midpoint [%]")
    ax.set_ylabel("R_pA (〈over y, pT〉)")
    ax.set_title(f"{particle}: RpA vs centrality")
    ax.grid(True, which="both", linestyle=":")
    return fig, ax, data

def _choose_files_for_avg(files: List[DataFile]) -> List[DataFile]:
    mb = [f for f in files if f.centrality is None]
    return mb if mb else files

def rpA_vs_y(output_dir: Union[str, 'Path'], particle: str = "JPsi", pt_range: Optional[Tuple[float, float]] = None, average_over_all_centralities_if_no_minbias: bool = True, show_errorbars: bool = False, ax: Optional[plt.Axes] = None):
    files = scan_output_dir(output_dir, particle=particle, kind="RpA")
    if not files:
        raise FileNotFoundError(f"No RpA files found for particle={particle} in {output_dir}")

    chosen = _choose_files_for_avg(files) if average_over_all_centralities_if_no_minbias else files
    curves = [f.mean_vs_y(pt_range=pt_range, mode="linear") for f in chosen]
    ys_common = sorted(set.intersection(*[set(map(float, c["y"])) for c in curves]))
    vals = []
    errs = []
    for y in ys_common:
        v = []
        e2 = []
        for c in curves:
            row = c[c["y"] == y]
            if not row.empty:
                v.append(float(row["value"].iloc[0]))
                err = row["error"].iloc[0]
                if not np.isnan(err):
                    e2.append(float(err) ** 2)
        if v:
            vals.append(np.mean(v))
            errs.append(np.sqrt(np.sum(e2)) / max(1, len(v)) if e2 else 0.0)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    if show_errorbars and any(e > 0 for e in errs):
        ax.errorbar(ys_common, vals, yerr=errs, fmt="-")
    else:
        ax.plot(ys_common, vals, "-")

    label = "minbias" if any(f.centrality is None for f in chosen) else "〈over centrality〉"
    ax.set_xlabel("y")
    ax.set_ylabel("R_pA (〈over pT〉)")
    ax.set_title(f"{particle}: RpA vs y ({label})")
    ax.grid(True, which="both", linestyle=":")
    return fig, ax, {"y": ys_common, "mean": vals, "err": errs}

def rpA_vs_pt(output_dir: Union[str, 'Path'], particle: str = "JPsi", y_range: Optional[Tuple[float, float]] = None, average_over_all_centralities_if_no_minbias: bool = True, show_errorbars: bool = False, ax: Optional[plt.Axes] = None):
    files = scan_output_dir(output_dir, particle=particle, kind="RpA")
    if not files:
        raise FileNotFoundError(f"No RpA files found for particle={particle} in {output_dir}")

    chosen = _choose_files_for_avg(files) if average_over_all_centralities_if_no_minbias else files
    curves = [f.mean_vs_pt(y_range=y_range, mode="linear") for f in chosen]
    pts_common = sorted(set.intersection(*[set(map(float, c["pt"])) for c in curves]))
    vals = []
    errs = []
    for pt in pts_common:
        v = []
        e2 = []
        for c in curves:
            row = c[c["pt"] == pt]
            if not row.empty:
                v.append(float(row["value"].iloc[0]))
                err = row["error"].iloc[0]
                if not np.isnan(err):
                    e2.append(float(err) ** 2)
        if v:
            vals.append(np.mean(v))
            errs.append(np.sqrt(np.sum(e2)) / max(1, len(v)) if e2 else 0.0)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    if show_errorbars and any(e > 0 for e in errs):
        ax.errorbar(pts_common, vals, yerr=errs, fmt="-")
    else:
        ax.plot(pts_common, vals, "-")

    label = "minbias" if any(f.centrality is None for f in chosen) else "〈over centrality〉"
    ax.set_xlabel("pT [GeV]")
    ax.set_ylabel("R_pA (〈over y〉)")
    subtitle = f"{label}" + (f", y∈[{y_range[0]}, {y_range[1]}]" if y_range else "")
    ax.set_title(f"{particle}: RpA vs pT ({subtitle})")
    ax.grid(True, which="both", linestyle=":")
    return fig, ax, {"pt": pts_common, "mean": vals, "err": errs}