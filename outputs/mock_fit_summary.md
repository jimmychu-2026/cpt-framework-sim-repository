# Grid Scan Results – lowell_mock.csv

| 項目              | 數值/描述            |
|-------------------|----------------------|
| ell range         | 2 – 30               |
| Data points (N)   | 29                   |
| Baseline χ²       | 4.201                |
| Best-fit χ²       | 0.353                |
| Δχ²               | -3.847               |
| Best-fit f_LSS    | 0.5921               |
| Best-fit ell_IR   | 4.0000               |
| Best-fit p        | 5.0000               |

---

**Interpretation** : Δχ² = -3.847 → significant improvement  
**Limitations**    : prototype pretest  
**Evaluations**    : 3420 (20 × 19 × 9)


# Grid Scan Results – lowell_mock.csv

## Command
python fit_cmb_lowell_v62.py --csv lowell_mock.csv --outdir outputs
Loaded 29 data points from 'lowell_mock.csv'.

## Scan Parameters
- Grid scan: f_LSS × ell_IR × p = 20 × 19 × 9 = **3420 evaluations**

---

## Results

| 項目              | 數值/描述   |
|-------------------|-------------|
| ell range         | 2 – 30      |
| Data points (N)   | 29          |
| Baseline χ²       | 4.201       |
| Best-fit χ²       | 0.353       |
| Δχ²               | -3.847      |
| Baseline AIC      | 4.201       |
| Best-fit AIC      | 6.353       |
| ΔAIC              | 2.153       |
| Baseline BIC      | 4.201       |
| Best-fit BIC      | 10.455      |
| ΔBIC              | 6.255       |
| Best-fit f_LSS    | 0.5921      |
| Best-fit ell_IR   | 4.0000      |
| Best-fit p        | 5.0000      |

---

## Outputs
- Saved: `outputs\lowell_spectrum.png`  
- Saved: `outputs\suppression_kernel.png`  
- Saved: `outputs\chisq_heatmap.png`  
- Saved: `outputs\angular_correlation.png`

---

## Notes
- Interpretation: Δχ² = -3.847 → significant improvement  
- Limitations: prototype pretest  
- Evaluations: 3420 (20 × 19 × 9)


# Grid Scan Results – lowell_mock.csv (fixed p = 5)

## Command
python fit_cmb_lowell_v62.py --csv lowell_mock.csv --fixed-p 5 --outdir outputs_fixed_p5
Loaded 29 data points from 'lowell_mock.csv'.
Reduced scan: fixed p = 5.0

## Scan Parameters
- Grid scan: f_LSS × ell_IR × p = 20 × 19 × 1 = **380 evaluations**

---

## Results

| 項目              | 數值/描述   |
|-------------------|-------------|
| ell range         | 2 – 30      |
| Data points (N)   | 29          |
| Baseline χ²       | 4.201       |
| Best-fit χ²       | 0.353       |
| Δχ²               | -3.847      |
| Baseline AIC      | 4.201       |
| Best-fit AIC      | 4.353       |
| ΔAIC              | 0.153       |
| Baseline BIC      | 4.201       |
| Best-fit BIC      | 7.088       |
| ΔBIC              | 2.887       |
| Effective V6.2 k  | 2           |
| Best-fit f_LSS    | 0.5921      |
| Best-fit ell_IR   | 4.0000      |
| Best-fit p        | 5.0000      |

---

## Outputs
- Saved: `outputs_fixed_p5\lowell_spectrum.png`  
- Saved: `outputs_fixed_p5\suppression_kernel.png`  
- Saved: `outputs_fixed_p5\chisq_heatmap.png`  
- Saved: `outputs_fixed_p5\angular_correlation.png`

---

## Notes
- Interpretation: Δχ² = -3.847 → significant improvement  
- Limitations: prototype pretest  
- Evaluations: 380 (20 × 19 × 1)

# Grid Scan Results – lowell_mock.csv (fixed p = 5, fixed ell_IR = 4)

## Command
python fit_cmb_lowell_v62.py --csv lowell_mock.csv --fixed-p 5 --fixed-ell-ir 4 --outdir outputs_fixed_p5_ell4
Loaded 29 data points from 'lowell_mock.csv'.
Reduced scan: fixed p = 5.0
Reduced scan: fixed ell_IR = 4.0

## Scan Parameters
- Grid scan: f_LSS × ell_IR × p = 20 × 1 × 1 = **20 evaluations**

---

## Results

| 項目              | 數值/描述   |
|-------------------|-------------|
| ell range         | 2 – 30      |
| Data points (N)   | 29          |
| Baseline χ²       | 4.201       |
| Best-fit χ²       | 0.353       |
| Δχ²               | -3.847      |
| Baseline AIC      | 4.201       |
| Best-fit AIC      | 2.353       |
| ΔAIC              | -1.847      |
| Baseline BIC      | 4.201       |
| Best-fit BIC      | 3.721       |
| ΔBIC              | -0.480      |
| Effective V6.2 k  | 1           |
| Best-fit f_LSS    | 0.5921      |
| Best-fit ell_IR   | 4.0000      |
| Best-fit p        | 5.0000      |

---

## Outputs
- Saved: `outputs_fixed_p5_ell4\lowell_spectrum.png`  
- Saved: `outputs_fixed_p5_ell4\suppression_kernel.png`  
- Saved: `outputs_fixed_p5_ell4\chisq_heatmap.png`  
- Saved: `outputs_fixed_p5_ell4\angular_correlation.png`

---

## Notes
- Interpretation: Δχ² = -3.847 → significant improvement  
- Limitations: prototype pretest  
- Evaluations: 20 (20 × 1 × 1)

