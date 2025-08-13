# data_handler.py
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

__all__ = [
    "parse_filename",
    "read_tsv",
    "DataFile",
    "scan_output_dir",
]

_f_minbias = re.compile(r"^minbias_(?P<particle>[A-Za-z]+)_(?P<kind>.+)\.tsv$")
_f_central = re.compile(r"^cent_(?P<low>\d+)-(?P<high>\d+)_(?P<particle>[A-Za-z]+)_(?P<kind>.+)\.tsv$")

def parse_filename(path: Union[str, Path]) -> Dict[str, Union[str, Tuple[int, int], None]]:
    name = Path(path).name
    m = _f_minbias.match(name)
    if m:
        d = m.groupdict()
        return dict(tag="minbias", centrality=None, particle=d["particle"], kind=d["kind"])
    m = _f_central.match(name)
    if m:
        d = m.groupdict()
        low, high = int(d["low"]), int(d["high"])
        return dict(tag=f"cent_{low}-{high}", centrality=(low, high), particle=d["particle"], kind=d["kind"])
    return dict(tag=None, centrality=None, particle=None, kind=None)

def read_tsv(path: Union[str, Path], value_floor: float = 1e-29) -> pd.DataFrame:
    path = Path(path)
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python", comment="#")
    if df.shape[1] == 4:
        df.columns = ["y", "pt", "value", "error"]
    elif df.shape[1] == 3:
        df.columns = ["y", "pt", "value"]
        df["error"] = np.nan
    else:
        raise ValueError(f"Unexpected number of columns ({df.shape[1]}) in {path}")
    for c in ["y", "pt", "value", "error"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    if value_floor is not None and value_floor > 0:
        df.loc[df["value"].abs() <= value_floor, "value"] = np.nan
        df.loc[df["error"].abs() <= value_floor, "error"] = np.nan
    df = df.sort_values(["y", "pt"], kind="mergesort").reset_index(drop=True)
    return df

@dataclass
class DataFile:
    path: Path
    meta: Dict[str, Union[str, Tuple[int, int], None]]
    data: pd.DataFrame

    @classmethod
    def from_path(cls, path: Union[str, Path], value_floor: float = 1e-29) -> "DataFile":
        path = Path(path)
        meta = parse_filename(path)
        df = read_tsv(path, value_floor=value_floor)
        return cls(path=path, meta=meta, data=df)

    @property
    def particle(self) -> Optional[str]:
        return self.meta.get("particle")

    @property
    def kind(self) -> Optional[str]:
        return self.meta.get("kind")

    @property
    def centrality(self) -> Optional[Tuple[int, int]]:
        return self.meta.get("centrality")

    @property
    def tag(self) -> Optional[str]:
        return self.meta.get("tag")

    def y_values(self) -> np.ndarray:
        return np.sort(self.data["y"].unique())

    def pt_values(self) -> np.ndarray:
        return np.sort(self.data["pt"].unique())

    def to_grid(self) -> pd.DataFrame:
        grid = self.data.pivot_table(index="y", columns="pt", values="value", aggfunc="mean")
        grid = grid.sort_index(axis=0).sort_index(axis=1)
        return grid

    def _subset(self, y_range: Optional[Tuple[float, float]], pt_range: Optional[Tuple[float, float]]) -> pd.DataFrame:
        df = self.data
        if y_range is not None:
            y0, y1 = y_range
            df = df[(df["y"] >= min(y0, y1)) & (df["y"] <= max(y0, y1))]
        if pt_range is not None:
            p0, p1 = pt_range
            df = df[(df["pt"] >= min(p0, p1)) & (df["pt"] <= max(p0, p1))]
        return df

    def mean(self, y_range: Optional[Tuple[float, float]] = None, pt_range: Optional[Tuple[float, float]] = None, mode: str = "auto") -> Tuple[float, Optional[float], int]:
        df = self._subset(y_range, pt_range).copy()
        df = df.dropna(subset=["value"])
        n = len(df)
        if n == 0:
            return (np.nan, np.nan, 0)
        if mode == "auto":
            k = (self.kind or "").lower()
            mode = "log" if "cross-section" in k else "linear"
        if mode == "linear":
            m = float(df["value"].mean())
        elif mode == "log":
            vals = df["value"].to_numpy()
            vals = vals[vals > 0]
            if len(vals) == 0:
                return (np.nan, np.nan, 0)
            m = float(np.exp(np.log(vals).mean()))
        else:
            raise ValueError("mode must be 'auto', 'linear', or 'log'")
        if "error" in df and df["error"].notna().any():
            err = df["error"].fillna(0.0).to_numpy()
            mean_err = float(np.sqrt(np.sum(err ** 2)) / max(1, n))
        else:
            mean_err = None
        return (m, mean_err, n)

    def mean_vs_y(self, pt_range: Optional[Tuple[float, float]] = None, mode: str = "auto") -> pd.DataFrame:
        rows = []
        for y in self.y_values():
            m, e, n = self.mean(y_range=(y, y), pt_range=pt_range, mode=mode)
            rows.append(dict(y=float(y), value=m, error=e, n=n))
        return pd.DataFrame(rows).dropna(subset=["value"])

    def mean_vs_pt(self, y_range: Optional[Tuple[float, float]] = None, mode: str = "auto") -> pd.DataFrame:
        rows = []
        for pt in self.pt_values():
            m, e, n = self.mean(y_range=y_range, pt_range=(pt, pt), mode=mode)
            rows.append(dict(pt=float(pt), value=m, error=e, n=n))
        return pd.DataFrame(rows).dropna(subset=["value"])

def scan_output_dir(root: Union[str, Path], particle: Optional[str] = None, kind: Optional[str] = None) -> List[DataFile]:
    root = Path(root)
    files = [p for p in root.glob("*.tsv")]
    out: List[DataFile] = []
    for p in files:
        meta = parse_filename(p)
        if particle and meta.get("particle") != particle:
            continue
        if kind and meta.get("kind") != kind:
            continue
        out.append(DataFile.from_path(p))

    def cent_mid(df: DataFile) -> float:
        c = df.centrality
        if c is None:
            return 1e9
        low, high = c
        return 0.5 * (low + high)

    out.sort(key=lambda d: (cent_mid(d), d.path.name))
    return out