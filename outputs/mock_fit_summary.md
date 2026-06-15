# Mock Fit Summary for V6.2 Low-ℓ Prototype

## Overview

This note summarizes three mock-data tests of the phenomenological low-ℓ suppression model implemented in `fit_cmb_lowell_v62.py`.

The goal is to compare:

1. A fully free **3-parameter** suppression model
2. A reduced **2-parameter** model with fixed `p = 5`
3. A constrained **1-parameter template** with fixed `p = 5` and fixed `ell_IR = 4`

The mock input data are provided in:

- `lowell_mock.csv`

These tests are **not** based on the full Planck low-ℓ likelihood. They are preliminary consistency checks using a simplified toy baseline and a simple chi-squared comparison.

---

## Suppression Model

The phenomenological V6.2 suppression kernel is:

```text
S_ell = f_LSS^2 / (1 + (ell / ell_IR)^p)
C_ell_model = C_ell_base * (1 - S_ell)
```

where:

- `f_LSS` is the overall low-ℓ suppression amplitude
- `ell_IR` is the IR transition scale
- `p` controls transition sharpness

---

## Test Results

| Model version | Effective parameters k | Best-fit f_LSS | Best-fit ell_IR | Best-fit p | Baseline chi2 | Best-fit chi2 | Delta chi2 | Delta AIC | Delta BIC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3-parameter free scan | 3 | 0.5921 | 4.0000 | 5.0000 | 4.201 | 0.353 | -3.847 | +2.153 | +6.255 |
| 2-parameter reduced scan (`p = 5`) | 2 | 0.5921 | 4.0000 | 5.0000 | 4.201 | 0.353 | -3.847 | +0.153 | +2.887 |
| 1-parameter constrained template (`p = 5`, `ell_IR = 4`) | 1 | 0.5921 | 4.0000 | 5.0000 | 4.201 | 0.353 | -3.847 | -1.847 | -0.480 |

---

## Interpretation

### 1. Three-parameter free model
The fully free model improves the fit substantially at the chi-squared level:

- `Delta chi2 = -3.847`

However, once the extra parameter cost is included, both AIC and BIC worsen:

- `Delta AIC = +2.153`
- `Delta BIC = +6.255`

This means that the free 3-parameter phenomenological version is **too flexible** relative to the improvement it buys.

### 2. Two-parameter reduced model
Fixing `p = 5` preserves the same best-fit location and the same fit improvement, while reducing the information-criterion penalty.

Even so:

- `Delta AIC = +0.153`
- `Delta BIC = +2.887`

This version is **close to competitive**, but still not clearly preferred over the baseline.

### 3. One-parameter constrained template
Fixing both:

- `p = 5`
- `ell_IR = 4`

and allowing only `f_LSS` to vary yields the same chi-squared improvement while reducing the complexity penalty enough to reverse the model-selection verdict:

- `Delta AIC = -1.847`
- `Delta BIC = -0.480`

This is the most promising result of the mock-data exercise.

---

## Main Conclusion

The mock tests suggest the following:

1. The V6.2 low-ℓ suppression kernel can reproduce an anomaly-like low-ℓ pattern.
2. The best-fit region is stable across the free, reduced, and constrained scans.
3. The strongest support appears when the suppression shape is treated as a **theory-motivated fixed template** and only the amplitude `f_LSS` is left free.

A concise statement is:

> The V6.2 low-ℓ suppression ansatz appears viable in a constrained predictive form, but not yet in an unconstrained multi-parameter phenomenological form.

---

## Important Limitations

These results should be interpreted cautiously:

- The data are mock, not a full observational likelihood.
- The baseline spectrum is a toy analytic spectrum, not CAMB/CLASS output.
- The chi-squared is simplified and does not include a full covariance matrix.
- The analysis uses only low-ℓ TT information.
- The apparent preference in the 1-parameter case is therefore **suggestive**, not conclusive.

---

## Recommended Next Steps

1. Treat the 1-parameter template (`p = 5`, `ell_IR = 4`) as the primary demonstrator.
2. Repeat the same test on a more realistic low-ℓ TT dataset.
3. Replace the toy baseline with a proper ΛCDM reference spectrum.
4. Extend the summary statistics to include angular-correlation suppression measures such as large-angle deficits.
5. Only after those steps consider a more complete likelihood-based analysis.
