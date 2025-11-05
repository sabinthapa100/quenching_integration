#!/usr/bin/env python3
import argparse, os, glob, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

NUM_COLS_ALL = [
    "y","pt","muA","alphaA","muB","alphaB","dyAmax","dyBmax",
    "zA_med","PhA_med","zB_med","PhB_med",
    "zA_q25","PhA_q25","zA_q50","PhA_q50","zA_q75","PhA_q75",
    "zB_q25","PhB_q25","zB_q50","PhB_q50","zB_q75","PhB_q75",
]

def nice_label_from_path(path: str) -> str:
    base = os.path.basename(path)
    tag  = base.replace("diagnostics.csv","").rstrip("_")
    if tag.startswith("cent_"): tag = tag[len("cent_"):]
    if tag.startswith("minbias"): return "minbias"
    return tag.replace("_","").replace("-", "–") + "%"

def subplot_grid(n: int):
    cols = min(3, max(1, int(math.ceil(n**0.5))))
    rows = int(math.ceil(n/cols))
    return rows, cols

def load_numeric_csv(fp: str) -> pd.DataFrame:
    # tolerant to repeated headers & stray text
    df = pd.read_csv(fp, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    keep = [c for c in NUM_COLS_ALL if c in df.columns]
    df = df[keep].copy()
    for c in keep:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["muA","alphaA","y","pt"])
    return df

def shared_limits(files):
    all_muA, all_pt, all_y, all_a = [], [], [], []
    for fp in files:
        d = load_numeric_csv(fp)
        if not d.empty:
            all_muA.append(d["muA"].values)
            all_pt.append(d["pt"].values)
            all_y.append(d["y"].values)
            all_a.append(d["alphaA"].values)
    muA = np.concatenate(all_muA) if all_muA else np.array([0,1])
    pt  = np.concatenate(all_pt)  if all_pt  else np.array([0,1])
    y   = np.concatenate(all_y)   if all_y   else np.array([-5,5])
    a   = np.concatenate(all_a)   if all_a   else np.array([0.2,0.4])
    def pad(lo, hi, frac=0.05):
        span = hi - lo
        return lo - frac*span, hi + frac*span
    xpt  = pad(pt.min(),  pt.max())
    xy   = pad(y.min(),   y.max())
    xmuA = pad(max(muA.min(), 0.0), muA.max())
    xa   = pad(max(a.min(), 0.0), a.max())
    return xpt, xy, xmuA, xa

def place_left_colorbar_for_grid(fig, axes, mappable, label):
    # Place a single vertical colorbar hugging the LEFT of the entire grid.
    # Compute the bbox of all visible axes:
    bbs = [ax.get_position(fig) for ax in axes.ravel() if ax.get_visible()]
    left   = min(bb.x0 for bb in bbs)
    bottom = min(bb.y0 for bb in bbs)
    top    = max(bb.y1 for bb in bbs)
    height = top - bottom
    # Leave a thin margin to the left of the leftmost axes:
    cax_left = max(0.02, left - 0.035)
    cax = fig.add_axes([cax_left, bottom, 0.018, height])
    cb = fig.colorbar(mappable, cax=cax)
    cb.set_label(label)
    return cb

def pick_quench(df, side="A"):
    z_med, ph_med = f"z{side}_med", f"Ph{side}_med"
    z50, ph50     = f"z{side}_q50", f"Ph{side}_q50"
    if z_med in df and ph_med in df and not df[z_med].isna().all():
        return df[z_med].to_numpy(), df[ph_med].to_numpy()
    if z50 in df and ph50 in df and not df[z50].isna().all():
        return df[z50].to_numpy(), df[ph50].to_numpy()
    for q in ("q25","q75"):
        zq, phq = f"z{side}_{q}", f"Ph{side}_{q}"
        if zq in df and phq in df and not df[zq].isna().all():
            return df[zq].to_numpy(), df[phq].to_numpy()
    return np.array([]), np.array([])

# ---------------------------- PLOTS ---------------------------------

def plot_alphas_vs_muA(files, save_path):
    n = len(files)
    rows, cols = subplot_grid(n)
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows),
                             squeeze=False, constrained_layout=True)
    _, _, xmuA, xa = shared_limits(files)
    for ax, fp in zip(axes.ravel(), files):
        lab = nice_label_from_path(fp); df = load_numeric_csv(fp)
        if df.empty: ax.set_title(f"{lab} (no data)"); ax.axis("off"); continue
        ax.scatter(df["muA"], df["alphaA"], s=6, alpha=0.5)
        try:
            q = min(20, max(5, max(1, len(df)//200)))
            bins = pd.qcut(df["muA"], q=q, duplicates="drop")
            m = df.groupby(bins, observed=False)["alphaA"].median()
            x = [b.mid for b in m.index.categories]
            ax.plot(x, m.values, lw=2)
        except Exception:
            pass
        ax.set_xlim(*xmuA); ax.set_ylim(*xa)
        ax.set_xlabel(r"$\mu_A \equiv \Delta p_{T,A}$ [GeV]")
        ax.set_ylabel(r"$\alpha_s(\mu_A)$")
        ax.set_title(lab); ax.grid(True, ls=":")
    for ax in axes.ravel()[n:]: ax.axis("off")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path, dpi=200)

def plot_muA_vs_y(files, save_path):
    n = len(files)
    rows, cols = subplot_grid(n)
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows),
                             squeeze=False, constrained_layout=True)
    last_sc = None
    (_, xy, xmuA, _) = shared_limits(files)
    for ax, fp in zip(axes.ravel(), files):
        lab = nice_label_from_path(fp); df = load_numeric_csv(fp)
        if df.empty: ax.set_title(f"{lab} (no data)"); ax.axis("off"); continue
        last_sc = ax.scatter(df["y"], df["muA"], s=6, c=df["pt"], alpha=0.6)
        ax.set_xlim(*xy); ax.set_ylim(*xmuA)
        ax.set_xlabel(r"$y$"); ax.set_ylabel(r"$\mu_A$ [GeV]")
        ax.set_title(lab); ax.grid(True, ls=":")
    for ax in axes.ravel()[n:]: ax.axis("off")
    if last_sc is not None:
        place_left_colorbar_for_grid(fig, axes, last_sc, r"$p_T$ [GeV]")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path, dpi=200)

def plot_muA_vs_pt(files, save_path):
    n = len(files)
    rows, cols = subplot_grid(n)
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows),
                             squeeze=False, constrained_layout=True)
    last_sc = None
    (xpt, _, xmuA, _) = shared_limits(files)
    for ax, fp in zip(axes.ravel(), files):
        lab = nice_label_from_path(fp); df = load_numeric_csv(fp)
        if df.empty: ax.set_title(f"{lab} (no data)"); ax.axis("off"); continue
        last_sc = ax.scatter(df["pt"], df["muA"], s=6, c=df["y"], alpha=0.6)
        ax.set_xlim(*xpt); ax.set_ylim(*xmuA)
        ax.set_xlabel(r"$p_T$ [GeV]"); ax.set_ylabel(r"$\mu_A$ [GeV]")
        ax.set_title(lab); ax.grid(True, ls=":")
    for ax in axes.ravel()[n:]: ax.axis("off")
    if last_sc is not None:
        place_left_colorbar_for_grid(fig, axes, last_sc, r"$y$")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path, dpi=200)

def nearest_vals(values, targets):
    vals = np.array(sorted(np.unique(values)))
    chosen = []
    for t in targets:
        idx = np.argmin(np.abs(vals - t))
        chosen.append(vals[idx])
    return sorted(list(dict.fromkeys(chosen)))  # unique, sorted

def plot_quench_lines_vs_pt(files, save_path, side="A", ys=(-4,-2,0,2,4)):
    n = len(files)
    rows, cols = subplot_grid(n)
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows),
                             squeeze=False, constrained_layout=True)
    # global y-limits (in log10) for consistency
    logmins, logmaxs = [], []
    cache = []
    for fp in files:
        df = load_numeric_csv(fp); _, Ph = pick_quench(df, side=side)
        if df.empty or len(Ph)==0: cache.append(None); continue
        ys_sel = nearest_vals(df["y"].values, ys)
        curves = []
        for yy in ys_sel:
            sub = df.loc[np.isclose(df["y"], yy)]
            sub = sub.sort_values("pt")
            Z, Ph = pick_quench(sub, side=side)
            if len(Ph)==0: continue
            logP = np.log10(np.clip(Ph, 1e-30, None))
            curves.append((yy, sub["pt"].to_numpy(), logP))
            logmins.append(logP.min()); logmaxs.append(logP.max())
        cache.append(curves)
    vmin = (min(logmins) if logmins else -6.0); vmax = (max(logmaxs) if logmaxs else -1.0)

    for ax, curves, fp in zip(axes.ravel(), cache, files):
        lab = nice_label_from_path(fp)
        if not curves: ax.set_title(f"{lab} (no data)"); ax.axis("off"); continue
        for yy, xpt, logP in curves:
            ax.plot(xpt, logP, lw=1.6, label=f"y≈{yy:g}")
        ax.set_ylim(vmin, vmax)
        ax.set_xlabel(r"$p_T$ [GeV]")
        ax.set_ylabel(r"$\log_{10}\,P_{\hat{}}^{(%s)}$" % side)
        ax.set_title(lab); ax.grid(True, ls=":")
        ax.legend(fontsize=8, ncol=2, frameon=False)
    for ax in axes.ravel()[n:]: ax.axis("off")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path, dpi=200)

def plot_quench_lines_vs_y(files, save_path, side="A", pts=(0.5,1,2,5,10,20,40)):
    n = len(files)
    rows, cols = subplot_grid(n)
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows),
                             squeeze=False, constrained_layout=True)
    logmins, logmaxs = [], []
    cache = []
    for fp in files:
        df = load_numeric_csv(fp); _, Ph = pick_quench(df, side=side)
        if df.empty or len(Ph)==0: cache.append(None); continue
        pts_sel = nearest_vals(df["pt"].values, pts)
        curves = []
        for pp in pts_sel:
            sub = df.loc[np.isclose(df["pt"], pp)]
            sub = sub.sort_values("y")
            Z, Ph = pick_quench(sub, side=side)
            if len(Ph)==0: continue
            logP = np.log10(np.clip(Ph, 1e-30, None))
            curves.append((pp, sub["y"].to_numpy(), logP))
            logmins.append(logP.min()); logmaxs.append(logP.max())
        cache.append(curves)
    vmin = (min(logmins) if logmins else -6.0); vmax = (max(logmaxs) if logmaxs else -1.0)

    for ax, curves, fp in zip(axes.ravel(), cache, files):
        lab = nice_label_from_path(fp)
        if not curves: ax.set_title(f"{lab} (no data)"); ax.axis("off"); continue
        for pp, yy, logP in curves:
            ax.plot(yy, logP, lw=1.6, label=f"$p_T\\approx{pp:g}$")
        ax.set_ylim(vmin, vmax)
        ax.set_xlabel(r"$y$")
        ax.set_ylabel(r"$\log_{10}\,P_{\hat{}}^{(%s)}$" % side)
        ax.set_title(lab); ax.grid(True, ls=":")
        ax.legend(fontsize=8, ncol=2, frameon=False)
    for ax in axes.ravel()[n:]: ax.axis("off")
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    fig.savefig(save_path, dpi=200)

# --------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Plot ΔpT & running-coupling diagnostics for all centralities.")
    ap.add_argument("--output-dir", default="output", help="Directory with *diagnostics.csv files")
    ap.add_argument("--save-dir",   default=None,     help="Directory to save plots (default: <output-dir>/plots)")
    ap.add_argument("--ys",  type=float, nargs="*", default=[-4,-2,0,2,4],
                    help="Rapidity values to trace Ph(y,pT) lines (nearest available will be used).")
    ap.add_argument("--pts", type=float, nargs="*", default=[0.5,1,2,5,10,20,40],
                    help="pT values to trace Ph(y,pT) lines (nearest available will be used).")
    args = ap.parse_args()
    outdir  = args.output_dir
    savedir = args.save_dir or os.path.join(outdir, "plots")

    files = sorted(glob.glob(os.path.join(outdir, "*diagnostics.csv")))
    if not files:
        print("No diagnostics.csv files found in", outdir)
        return
    print(f"Found {len(files)} diagnostics files:")
    for f in files: print("  -", os.path.basename(f))
    os.makedirs(savedir, exist_ok=True)

    plot_alphas_vs_muA(files, os.path.join(savedir, "alphas_vs_muA_allbins.png"))
    plot_muA_vs_y    (files, os.path.join(savedir, "muA_vs_y_allbins.png"))
    plot_muA_vs_pt   (files, os.path.join(savedir, "muA_vs_pt_allbins.png"))

    # Quenching weights: scatter (already done previously) not requested now.
    # New: line plots for easier comparison:
    plot_quench_lines_vs_pt(files, os.path.join(savedir, "PhA_lines_vs_pt_allbins.png"),
                            side="A", ys=args.ys)
    plot_quench_lines_vs_y (files, os.path.join(savedir, "PhA_lines_vs_y_allbins.png"),
                            side="A", pts=args.pts)

    # If you also want B-side, uncomment:
    # plot_quench_lines_vs_pt(files, os.path.join(savedir, "PhB_lines_vs_pt_allbins.png"), side="B", ys=args.ys)
    # plot_quench_lines_vs_y (files, os.path.join(savedir, "PhB_lines_vs_y_allbins.png"), side="B", pts=args.pts)

    print("Saved plots to", savedir)

if __name__ == "__main__":
    main()
