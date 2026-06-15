# Mock Fit Summary for V6.2 Low-ℓ Prototype

## Overview

This note summarizes the mock-data tests of the phenomenological low-ℓ suppression model implemented in `fit_cmb_lowell_v62.py`.

The goal is to compare:

1. A fully free **3-parameter** suppression model
2. A reduced **2-parameter** model with fixed `p = 5`
3. A constrained **1-parameter template** with fixed `p = 5` and fixed `ell_IR = 4`

The mock input data are provided in:

- `lowell_mock.csv`
- `lowell_tt_mock_extended.csv`
- `lowell_tt_realistic.csv`

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

### Earlier mock baseline results

| Model version | Effective parameters k | Best-fit f_LSS | Best-fit ell_IR | Best-fit p | Baseline chi2 | Best-fit chi2 | Delta chi2 | Delta AIC | Delta BIC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3-parameter free scan | 3 | 0.5921 | 4.0000 | 5.0000 | 4.201 | 0.353 | -3.847 | +2.153 | +6.255 |
| 2-parameter reduced scan (`p = 5`) | 2 | 0.5921 | 4.0000 | 5.0000 | 4.201 | 0.353 | -3.847 | +0.153 | +2.887 |
| 1-parameter constrained template (`p = 5`, `ell_IR = 4`) | 1 | 0.5921 | 4.0000 | 5.0000 | 4.201 | 0.353 | -3.847 | -1.847 | -0.480 |

### Console result from the broader mock dataset

| Model version | Effective parameters k | Best-fit f_LSS | Best-fit ell_IR | Best-fit p | Baseline chi2 | Best-fit chi2 | Delta chi2 | Delta AIC | Delta BIC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3-parameter free scan | 3 | 0.273183 | 10.000000 | 1.150000 | 5.595672 | 0.288643 | -5.307028 | +6.289 | +10.391 |
| 2-parameter reduced scan (`p = 5`) | 2 | 0.248120 | 10.000000 | 5.000000 | 5.595672 | 1.634128 | -3.961543 | +5.634 | +8.369 |
| 1-parameter constrained template (`p = 5`, `ell_IR = 4`) | 1 | 0.258145 | 4.000000 | 5.000000 | 5.595672 | 4.356493 | -1.239179 | +6.356 | +7.724 |

### Shape implication from the latest fit

The latest console output indicates the best-fit region lies at:

- relatively **large `ell_IR`**
- relatively **small `p`** in the fully free scan

This points to a **broader, flatter suppression kernel** than the earlier toy example.

---

## Interpretation

### 1. Three-parameter free model
The fully free model improves the fit at the chi-squared level:

- `Delta chi2 = -5.307028`

However, once the extra parameter cost is included, both AIC and BIC worsen:

- `Delta AIC = +6.289`
- `Delta BIC = +10.391`

This means the free 3-parameter phenomenological version is still **too flexible** relative to the improvement it buys.

### 2. Two-parameter reduced model
Fixing `p = 5` changes the best-fit values and preserves a fit improvement:

- `Delta chi2 = -3.961543`

Even so, the information-criterion penalty remains substantial:

- `Delta AIC = +5.634`
- `Delta BIC = +8.369`

This version remains **disfavored** relative to the baseline for the current mock dataset.

### 3. One-parameter constrained template
Fixing both:

- `p = 5`
- `ell_IR = 4`

and allowing only `f_LSS` to vary gives a smaller improvement:

- `Delta chi2 = -1.239179`

The corresponding information-criterion penalties are still positive:

- `Delta AIC = +6.356`
- `Delta BIC = +7.724`

So, for this newer mock dataset, the constrained template is **not yet preferred** by AIC/BIC.

### 4. Comparison across datasets
The earlier toy mock favored a tightly constrained template, whereas the newer broader dataset favors a wider IR suppression shape in the free scan but does not yet provide information-criterion support for the template forms.

This indicates that the preferred shape is **dataset-dependent** and that the template family should be scanned over a broader range before drawing a stronger conclusion.

---

## Main Conclusion

The mock tests suggest the following:

1. The V6.2 low-ℓ suppression kernel can reproduce an anomaly-like low-ℓ pattern.
2. The best-fit region depends on the dataset, but the latest broader mock fit favors a wider and shallower suppression shape.
3. The earlier toy exercise supported a fixed-template form, whereas the latest broader dataset does not yet deliver information-criterion support for either the reduced or constrained forms.
4. The current evidence remains suggestive only; it does **not** establish preference over baseline in the latest mock dataset.

A concise statement is:

> The V6.2 low-ℓ suppression ansatz remains viable as a broad phenomenological template, but the latest mock fit prefers a wider and shallower suppression shape and does not yet deliver information-criterion support over the baseline.

---

## Important Limitations

These results should be interpreted cautiously:

- The data are mock, not a full observational likelihood.
- The baseline spectrum is a toy analytic spectrum, not CAMB/CLASS output.
- The chi-squared is simplified and does not include a full covariance matrix.
- The analysis uses only low-ℓ TT information.
- The apparent preference in the latest fit is therefore **suggestive**, not conclusive.

---

## Recommended Next Steps

1. Re-scan `p` and `ell_IR` over a broader range centered on the latest best-fit region.
2. Test whether a flatter kernel family can improve AIC/BIC.
3. Repeat the same test on a more realistic low-ℓ TT dataset.
4. Replace the toy baseline with a proper ΛCDM reference spectrum.
5. Extend the summary statistics to include angular-correlation suppression measures such as large-angle deficits.
6. Only after those steps consider a more complete likelihood-based analysis.
