import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

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
    f_lss: float
    ell_ir: float
    p: float
    chi2: float
    aic: float
    bic: float
    n_points: int


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


def grid_search(points: List[LowellPoint], f_grid=None, ell_ir_grid=None, p_grid=None) -> FitResult:
    if f_grid is None:
        f_grid = np.linspace(0.0, 1.0, 400)
    if ell_ir_grid is None:
        ell_ir_grid = np.linspace(2.0, 10.0, 200)
    if p_grid is None:
        p_grid = np.linspace(1.0, 10.0, 181)

    best = None
    best_params = None
    for f_lss in f_grid:
        for ell_ir in ell_ir_grid:
            for p in p_grid:
                c2 = chi2(points, f_lss, ell_ir, p)
                if best is None or c2 < best:
                    best = c2
                    best_params = (f_lss, ell_ir, p)

    assert best is not None and best_params is not None
    n = sum(1 for pt in points if pt.mask_flag == 0)
    k = 3
    aic = best + 2 * k
    bic = best + k * math.log(n)
    return FitResult(best_params[0], best_params[1], best_params[2], float(best), float(aic), float(bic), n)


def summarize(points: List[LowellPoint], fit: FitResult) -> str:
    baseline = chi2(points, 0.0, fit.ell_ir, fit.p)
    return (
        f"n_points={fit.n_points}\n"
        f"best_f_lss={fit.f_lss:.6f}\n"
        f"best_ell_ir={fit.ell_ir:.6f}\n"
        f"best_p={fit.p:.6f}\n"
        f"baseline_chi2={baseline:.6f}\n"
        f"best_chi2={fit.chi2:.6f}\n"
        f"delta_chi2={fit.chi2 - baseline:.6f}\n"
        f"AIC={fit.aic:.6f}\n"
        f"BIC={fit.bic:.6f}\n"
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Fit mock low-ell TT data with V6.2 suppression kernel.")
    parser.add_argument("csv_path", help="Path to low-ell CSV file")
    args = parser.parse_args()

    pts = load_lowell_csv(args.csv_path)
    fit = grid_search(pts)
    print(summarize(pts, fit))
