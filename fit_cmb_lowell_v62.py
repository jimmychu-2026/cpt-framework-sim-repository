import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np


@dataclass
class LowellPoint:
    ell: int
    cl_tt: float
    cl_tt_err: float
    cl_lcdm: float
    cl_residual: float
    sigma_cv: float
    weight: float
    mask_flag: int


@dataclass
class FitResult:
    label: str
    k: int
    f_lss: float
    ell_ir: float
    p: float
    chi2: float
    aic: float
    bic: float
    n_points: int


@dataclass
class ScanConfig:
    label: str
    k: int
    f_grid: np.ndarray
    ell_ir_grid: np.ndarray
    p_grid: np.ndarray
    fixed_p: Optional[float] = None
    fixed_ell_ir: Optional[float] = None


def load_lowell_csv(path: str | Path) -> List[LowellPoint]:
    points: List[LowellPoint] = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            points.append(
                LowellPoint(
                    ell=int(row["ell"]),
                    cl_tt=float(row["Cl_TT"]),
                    cl_tt_err=float(row["Cl_TT_err"]),
                    cl_lcdm=float(row["Cl_LCDM"]),
                    cl_residual=float(row["Cl_residual"]),
                    sigma_cv=float(row.get("sigma_cv", row["Cl_TT_err"])),
                    weight=float(row.get("weight", 1.0)),
                    mask_flag=int(row.get("mask_flag", 0)),
                )
            )
    return points


def suppression_kernel(ell: np.ndarray, f_lss: float, ell_ir: float, p: float) -> np.ndarray:
    return (f_lss ** 2) / (1.0 + (ell / ell_ir) ** p)


def model_cl(lcdm: np.ndarray, ell: np.ndarray, f_lss: float, ell_ir: float, p: float) -> np.ndarray:
    return lcdm * (1.0 - suppression_kernel(ell, f_lss, ell_ir, p))


def chi2(points: List[LowellPoint], f_lss: float, ell_ir: float, p: float) -> float:
    ell = np.array([pt.ell for pt in points if pt.mask_flag == 0], dtype=float)
    obs = np.array([pt.cl_tt for pt in points if pt.mask_flag == 0], dtype=float)
    err = np.array([pt.cl_tt_err for pt in points if pt.mask_flag == 0], dtype=float)
    lcdm = np.array([pt.cl_lcdm for pt in points if pt.mask_flag == 0], dtype=float)
    w = np.array([pt.weight for pt in points if pt.mask_flag == 0], dtype=float)
    mod = model_cl(lcdm, ell, f_lss, ell_ir, p)
    return float(np.sum(w * ((obs - mod) / err) ** 2))


def grid_search(points: List[LowellPoint], config: ScanConfig) -> FitResult:
    best = None
    best_params: Optional[Tuple[float, float, float]] = None

    for f_lss in config.f_grid:
        for ell_ir in config.ell_ir_grid:
            if config.fixed_ell_ir is not None and not math.isclose(ell_ir, config.fixed_ell_ir, rel_tol=0.0, abs_tol=1e-12):
                continue
            for p in config.p_grid:
                if config.fixed_p is not None and not math.isclose(p, config.fixed_p, rel_tol=0.0, abs_tol=1e-12):
                    continue
                c2 = chi2(points, f_lss, ell_ir, p)
                if best is None or c2 < best:
                    best = c2
                    best_params = (f_lss, ell_ir, p)

    assert best is not None and best_params is not None
    n = sum(1 for pt in points if pt.mask_flag == 0)
    aic = best + 2 * config.k
    bic = best + config.k * math.log(n)
    return FitResult(config.label, config.k, best_params[0], best_params[1], best_params[2], float(best), float(aic), float(bic), n)


def baseline_chi2(points: List[LowellPoint]) -> float:
    ell = np.array([pt.ell for pt in points if pt.mask_flag == 0], dtype=float)
    obs = np.array([pt.cl_tt for pt in points if pt.mask_flag == 0], dtype=float)
    err = np.array([pt.cl_tt_err for pt in points if pt.mask_flag == 0], dtype=float)
    lcdm = np.array([pt.cl_lcdm for pt in points if pt.mask_flag == 0], dtype=float)
    w = np.array([pt.weight for pt in points if pt.mask_flag == 0], dtype=float)
    return float(np.sum(w * ((obs - lcdm) / err) ** 2))


def summarize(points: List[LowellPoint], fits: List[FitResult]) -> str:
    base = baseline_chi2(points)
    lines = [f"n_points={sum(1 for pt in points if pt.mask_flag == 0)}", f"baseline_chi2={base:.6f}", ""]
    for fit in fits:
        lines.extend([
            f"[{fit.label}]",
            f"best_f_lss={fit.f_lss:.6f}",
            f"best_ell_ir={fit.ell_ir:.6f}",
            f"best_p={fit.p:.6f}",
            f"best_chi2={fit.chi2:.6f}",
            f"delta_chi2={fit.chi2 - base:.6f}",
            f"AIC={fit.aic:.6f}",
            f"BIC={fit.bic:.6f}",
            "",
        ])
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Fit mock low-ell TT data with V6.2 suppression kernel.")
    parser.add_argument("csv_path", help="Path to low-ell CSV file")
    args = parser.parse_args()

    pts = load_lowell_csv(args.csv_path)

    configs = [
        ScanConfig(
            label="3-parameter free scan",
            k=3,
            f_grid=np.linspace(0.0, 1.0, 400),
            ell_ir_grid=np.linspace(2.0, 10.0, 200),
            p_grid=np.linspace(1.0, 10.0, 181),
        ),
        ScanConfig(
            label="2-parameter reduced scan (fixed p = 5)",
            k=2,
            f_grid=np.linspace(0.0, 1.0, 400),
            ell_ir_grid=np.linspace(2.0, 10.0, 200),
            p_grid=np.array([5.0]),
            fixed_p=5.0,
        ),
        ScanConfig(
            label="1-parameter constrained template (fixed p = 5, ell_IR = 4)",
            k=1,
            f_grid=np.linspace(0.0, 1.0, 400),
            ell_ir_grid=np.array([4.0]),
            p_grid=np.array([5.0]),
            fixed_p=5.0,
            fixed_ell_ir=4.0,
        ),
    ]

    fits = [grid_search(pts, cfg) for cfg in configs]
    print(summarize(pts, fits))
