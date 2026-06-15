"""
fit_cmb_lowell_v62.py
---------------------
Preliminary phenomenological prototype for a minimal viable CMB low-ell fit
within the V6.2 CPT-Curvature framework.

This script implements a simple IR-suppression model for the low-multipole
TT power spectrum and performs a coarse grid scan over three free parameters:
  f_LSS   - overall suppression amplitude linked to the LSS projection factor
  ell_IR  - IR transition multipole (pivot scale)
  p       - transition sharpness exponent

The suppression kernel is:
    S_ell = f_LSS**2 / (1 + (ell / ell_IR)**p)
    C_ell_model = C_ell_base * (1 - S_ell)

This is NOT a full Planck likelihood analysis and does not use CAMB or CLASS.
It is intended as a first-pass consistency check for the V6.2 hypothesis.

Outputs
-------
- Console: best-fit parameters, chi-squared values, delta chi-squared,
  plus simple AIC/BIC model-comparison summaries
- Figures saved to the output directory:
    lowell_spectrum.png         -- low-ell TT spectrum comparison
    suppression_kernel.png      -- S_ell kernel for best-fit parameters
    chisq_heatmap.png           -- chi-squared heatmap (f_LSS vs ell_IR)
    angular_correlation.png     -- C(theta) computed via Legendre polynomials
"""

import argparse
import math
import os

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval


# ---------------------------------------------------------------------------
# Baseline spectrum utilities
# ---------------------------------------------------------------------------

def generate_toy_baseline(ell_min: int, ell_max: int) -> np.ndarray:
    """Return a simple smooth decaying toy TT spectrum for ell in [ell_min, ell_max].

    The shape loosely mimics the Sachs-Wolfe plateau at low ell:
        C_ell_base ~ A * ell^(-0.1) * exp(-ell / 300)
    scaled so that C_2 ~ 1000 (arbitrary units consistent with mu K^2 order).
    """
    ells = np.arange(ell_min, ell_max + 1, dtype=float)
    A = 1000.0 * 2.0**0.1 * np.exp(2.0 / 300.0)  # normalise so C_2 ~ 1000
    cl_base = A * ells**(-0.1) * np.exp(-ells / 300.0)
    return cl_base


def load_observed_data(csv_path: str, ell_min: int, ell_max: int):
    """Load observed low-ell TT data from a CSV file.

    Expected columns: ell, cl_obs, [sigma]
    Returns (ells, cl_obs, sigma) arrays filtered to [ell_min, ell_max].
    """
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    mask = (data["ell"] >= ell_min) & (data["ell"] <= ell_max)
    ells = data["ell"][mask].astype(int)
    cl_obs = data["cl_obs"][mask]
    if "sigma" in data.dtype.names:
        sigma = data["sigma"][mask]
    else:
        # Cosmic-variance-inspired approximation
        sigma = np.sqrt(2.0 / (2.0 * ells + 1.0)) * cl_obs
    return ells, cl_obs, sigma


def cosmic_variance_sigma(ells: np.ndarray, cl: np.ndarray) -> np.ndarray:
    """Estimate uncertainty using the cosmic variance approximation.

    sigma_ell = sqrt(2 / (2*ell + 1)) * C_ell
    """
    return np.sqrt(2.0 / (2.0 * ells + 1.0)) * cl


# ---------------------------------------------------------------------------
# V6.2 suppression model
# ---------------------------------------------------------------------------

def suppression_kernel(ells: np.ndarray, f_lss: float, ell_ir: float, p: float) -> np.ndarray:
    """Compute the V6.2 IR suppression kernel.

    S_ell = f_LSS**2 / (1 + (ell / ell_IR)**p)
    """
    return f_lss**2 / (1.0 + (ells / ell_ir) ** p)


def apply_suppression(cl_base: np.ndarray, s_ell: np.ndarray) -> np.ndarray:
    """Apply the suppression kernel to the baseline spectrum.

    C_ell_model = C_ell_base * (1 - S_ell)
    Clamps output to zero to avoid unphysical negative values.
    """
    return np.maximum(cl_base * (1.0 - s_ell), 0.0)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def chi_squared(cl_model: np.ndarray, cl_obs: np.ndarray, sigma: np.ndarray) -> float:
    """Compute the simple chi-squared between model and observed spectrum."""
    return float(np.sum(((cl_model - cl_obs) / sigma) ** 2))


def information_criteria(chi2: float, n_points: int, n_params: int) -> tuple[float, float]:
    """Return (AIC, BIC) for a model.

    AIC = chi2 + 2k
    BIC = chi2 + k ln(N)
    where k is the number of fitted parameters and N is the number of data points.
    """
    aic = float(chi2 + 2.0 * n_params)
    bic = float(chi2 + n_params * math.log(n_points))
    return aic, bic


# ---------------------------------------------------------------------------
# Grid scan
# ---------------------------------------------------------------------------

def grid_scan(
    ells: np.ndarray,
    cl_base: np.ndarray,
    cl_obs: np.ndarray,
    sigma: np.ndarray,
    f_lss_vals: np.ndarray,
    ell_ir_vals: np.ndarray,
    p_vals: np.ndarray,
):
    """Perform a coarse grid scan over (f_LSS, ell_IR, p).

    Returns a dict with shape (n_f, n_ir, n_p) chi-squared array and
    the index and values of the best-fit point.
    """
    n_f = len(f_lss_vals)
    n_ir = len(ell_ir_vals)
    n_p = len(p_vals)
    chi2_grid = np.full((n_f, n_ir, n_p), np.inf)

    for i, f_lss in enumerate(f_lss_vals):
        for j, ell_ir in enumerate(ell_ir_vals):
            for k, p in enumerate(p_vals):
                s_ell = suppression_kernel(ells, f_lss, ell_ir, p)
                cl_model = apply_suppression(cl_base, s_ell)
                chi2_grid[i, j, k] = chi_squared(cl_model, cl_obs, sigma)

    # Best-fit location
    idx = np.unravel_index(np.argmin(chi2_grid), chi2_grid.shape)
    best = {
        "f_lss": float(f_lss_vals[idx[0]]),
        "ell_ir": float(ell_ir_vals[idx[1]]),
        "p": float(p_vals[idx[2]]),
        "chi2": float(chi2_grid[idx]),
        "idx": idx,
    }
    return chi2_grid, best


# ---------------------------------------------------------------------------
# Angular correlation function
# ---------------------------------------------------------------------------

def angular_correlation(cl_values: np.ndarray, ells: np.ndarray, theta_deg: np.ndarray) -> np.ndarray:
    """Compute the angular two-point correlation function C(theta).

    C(theta) = (1/4*pi) * sum_ell (2*ell+1) * C_ell * P_ell(cos(theta))

    Uses numpy's Legendre polynomial evaluation.
    """
    cos_theta = np.cos(np.radians(theta_deg))
    result = np.zeros_like(cos_theta)
    for ell, cl in zip(ells, cl_values):
        coeffs = np.zeros(int(ell) + 1)
        coeffs[-1] = (2.0 * ell + 1.0) * cl
        result += legval(cos_theta, coeffs)
    return result / (4.0 * np.pi)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_spectrum(ells, cl_base, cl_model, cl_obs, sigma, outdir: str):
    """Plot the low-ell TT power spectrum comparison."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(ells, cl_obs, yerr=sigma, fmt="o", color="black", label="Observed", zorder=5)
    ax.plot(ells, cl_base, "--", color="royalblue", label=r"Baseline $\Lambda$CDM")
    ax.plot(ells, cl_model, "-", color="tomato", label="V6.2 best-fit model")
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"$C_\ell$  [arbitrary units]")
    ax.set_title("Low-$\\ell$ TT Power Spectrum (V6.2 suppression)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    path = os.path.join(outdir, "lowell_spectrum.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_suppression_kernel(ells, f_lss, ell_ir, p, outdir: str):
    """Plot the suppression kernel S_ell for the best-fit parameters."""
    s_ell = suppression_kernel(ells.astype(float), f_lss, ell_ir, p)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ells, s_ell, color="darkorange", lw=2)
    ax.axvline(ell_ir, color="gray", linestyle=":", label=rf"$\ell_{{IR}}={ell_ir:.1f}$")
    ax.set_xlabel(r"Multipole $\ell$")
    ax.set_ylabel(r"$S_\ell$")
    ax.set_title(
        rf"V6.2 Suppression Kernel  ($f_{{LSS}}={f_lss:.3f}$, $\ell_{{IR}}={ell_ir:.1f}$, $p={p:.1f}$)"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)
    path = os.path.join(outdir, "suppression_kernel.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_chisq_heatmap(chi2_grid, f_lss_vals, ell_ir_vals, best, outdir: str):
    """Plot a chi-squared heatmap for the best p value."""
    k_best = best["idx"][2]
    slice_2d = chi2_grid[:, :, k_best]
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.pcolormesh(
        ell_ir_vals, f_lss_vals, slice_2d,
        shading="auto", cmap="viridis_r",
        vmax=np.nanpercentile(slice_2d, 80),
    )
    fig.colorbar(im, ax=ax, label=r"$\chi^2$")
    ax.scatter(
        [best["ell_ir"]], [best["f_lss"]],
        color="red", marker="*", s=200, zorder=5, label="Best fit"
    )
    ax.set_xlabel(r"$\ell_{IR}$")
    ax.set_ylabel(r"$f_{LSS}$")
    ax.set_title(rf"$\chi^2$ Heatmap  ($p = {best['p']:.1f}$)")
    ax.legend()
    path = os.path.join(outdir, "chisq_heatmap.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_angular_correlation(ells, cl_base, cl_model, outdir: str):
    """Plot the angular two-point correlation function C(theta)."""
    theta = np.linspace(0.0, 180.0, 360)
    c_base = angular_correlation(cl_base, ells, theta)
    c_model = angular_correlation(cl_model, ells, theta)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(theta, c_base, "--", color="royalblue", label=r"Baseline $\Lambda$CDM")
    ax.plot(theta, c_model, "-", color="tomato", label="V6.2 best-fit model")
    ax.axhline(0, color="black", lw=0.8, linestyle=":")
    ax.set_xlabel(r"$\theta$ [degrees]")
    ax.set_ylabel(r"$C(\theta)$  [arbitrary units]")
    ax.set_title("Angular Two-Point Correlation Function (low-$\\ell$ only)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    path = os.path.join(outdir, "angular_correlation.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Minimal phenomenological CMB low-ell fit for the V6.2 framework. "
            "Performs a grid scan over suppression parameters and outputs plots."
        )
    )
    parser.add_argument(
        "--csv", default=None,
        help="Path to CSV file with columns: ell, cl_obs[, sigma]. "
             "If omitted, a toy internal baseline is used as both baseline and mock observed data.",
    )
    parser.add_argument("--ell-min", type=int, default=2, help="Minimum multipole (default: 2)")
    parser.add_argument("--ell-max", type=int, default=30, help="Maximum multipole (default: 30)")

    # Scan ranges
    parser.add_argument("--f-lss-min", type=float, default=0.01, help="Min f_LSS (default: 0.01)")
    parser.add_argument("--f-lss-max", type=float, default=0.80, help="Max f_LSS (default: 0.80)")
    parser.add_argument("--f-lss-n", type=int, default=20, help="Number of f_LSS grid points (default: 20)")

    parser.add_argument("--ell-ir-min", type=float, default=2.0, help="Min ell_IR (default: 2.0)")
    parser.add_argument("--ell-ir-max", type=float, default=20.0, help="Max ell_IR (default: 20.0)")
    parser.add_argument("--ell-ir-n", type=int, default=19, help="Number of ell_IR grid points (default: 19)")

    parser.add_argument("--p-min", type=float, default=2.0, help="Min p exponent (default: 2.0)")
    parser.add_argument("--p-max", type=float, default=10.0, help="Max p exponent (default: 10.0)")
    parser.add_argument("--p-n", type=int, default=9, help="Number of p grid points (default: 9)")

    parser.add_argument(
        "--outdir", default=".",
        help="Directory for output plot files (default: current directory)",
    )
    return parser


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = build_parser().parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ell_min = args.ell_min
    ell_max = args.ell_max
    ells = np.arange(ell_min, ell_max + 1, dtype=float)

    # --- Baseline spectrum ---
    cl_base = generate_toy_baseline(ell_min, ell_max)

    # --- Observed data ---
    if args.csv is not None:
        ells_obs_int, cl_obs, sigma = load_observed_data(args.csv, ell_min, ell_max)
        ells = ells_obs_int.astype(float)
        base_all = {int(e): c for e, c in zip(np.arange(ell_min, ell_max + 1), generate_toy_baseline(ell_min, ell_max))}
        cl_base = np.array([base_all[int(e)] for e in ells])
        cl_obs_arr = cl_obs
        sigma_arr = sigma
        print(f"Loaded {len(ells)} data points from '{args.csv}'.")
    else:
        cl_obs_arr = cl_base.copy()
        sigma_arr = cosmic_variance_sigma(ells, cl_obs_arr)
        print("No CSV provided — using toy baseline as mock observed data.")
        print("Toy mode reminder: this is a smoke test only; chi-squared and AIC/BIC are not physically informative.")

    chi2_base = chi_squared(cl_base, cl_obs_arr, sigma_arr)
    n_points = len(ells)
    n_params_baseline = 0
    n_params_v62 = 3
    aic_base, bic_base = information_criteria(chi2_base, n_points, n_params_baseline)

    f_lss_vals = np.linspace(args.f_lss_min, args.f_lss_max, args.f_lss_n)
    ell_ir_vals = np.linspace(args.ell_ir_min, args.ell_ir_max, args.ell_ir_n)
    p_vals = np.linspace(args.p_min, args.p_max, args.p_n)

    print(
        f"\nGrid scan: f_LSS x ell_IR x p = "
        f"{args.f_lss_n} x {args.ell_ir_n} x {args.p_n} = "
        f"{args.f_lss_n * args.ell_ir_n * args.p_n} evaluations ..."
    )
    chi2_grid, best = grid_scan(ells, cl_base, cl_obs_arr, sigma_arr, f_lss_vals, ell_ir_vals, p_vals)

    s_best = suppression_kernel(ells, best["f_lss"], best["ell_ir"], best["p"])
    cl_best = apply_suppression(cl_base, s_best)

    aic_best, bic_best = information_criteria(best["chi2"], n_points, n_params_v62)

    print("\n" + "=" * 50)
    print("RESULTS")
    print("=" * 50)
    print(f"  ell range          : {ell_min} - {ell_max}")
    print(f"  Data points (N)    : {n_points}")
    print(f"  Baseline chi2      : {chi2_base:.3f}")
    print(f"  Best-fit chi2      : {best['chi2']:.3f}")
    print(f"  Delta chi2         : {best['chi2'] - chi2_base:.3f}")
    print(f"  Baseline AIC       : {aic_base:.3f}")
    print(f"  Best-fit AIC       : {aic_best:.3f}")
    print(f"  Delta AIC          : {aic_best - aic_base:.3f}")
    print(f"  Baseline BIC       : {bic_base:.3f}")
    print(f"  Best-fit BIC       : {bic_best:.3f}")
    print(f"  Delta BIC          : {bic_best - bic_base:.3f}")
    print(f"  Best-fit f_LSS     : {best['f_lss']:.4f}")
    print(f"  Best-fit ell_IR    : {best['ell_ir']:.4f}")
    print(f"  Best-fit p         : {best['p']:.4f}")
    print("=" * 50)

    print("\nSaving plots ...")
    plot_spectrum(ells, cl_base, cl_best, cl_obs_arr, sigma_arr, args.outdir)
    plot_suppression_kernel(ells, best["f_lss"], best["ell_ir"], best["p"], args.outdir)
    plot_chisq_heatmap(chi2_grid, f_lss_vals, ell_ir_vals, best, args.outdir)
    plot_angular_correlation(ells, cl_base, cl_best, args.outdir)
    print("\nDone.")


if __name__ == "__main__":
    main()
