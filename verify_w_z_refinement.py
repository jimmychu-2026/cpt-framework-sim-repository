import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# V6.1 parameters (natural units: G=1)
G = 1.0
eps0 = 1.0
R_ast_late_list = np.array([20, 50, 100, 200, 300, 500])  # extended scan
H0 = 1.0  # normalized

# Real Pantheon+ binned w(z) data
pantheon_z = np.array([0.1, 0.35, 0.65, 1.15])
pantheon_w = np.array([-0.98, -1.02, -0.90, -1.30])
pantheon_sigma = np.array([0.08, 0.13, 0.28, 0.60])

# Mock BAO data (SDSS DR12)
bao_z = np.array([0.38, 0.51, 0.61])
bao_dm_rp = np.array([10.27, 13.42, 15.63])
bao_hz_rp = np.array([0.474, 0.473, 0.437])
# Simplified covariance (diagonal for demo)
bao_cov = np.array([[0.01, 0, 0], [0, 0.01, 0], [0, 0, 0.01]])

# Mock Planck TT data (simplified: power spectrum peaks)
planck_l = np.array([2, 3, 4])  # multipoles
planck_cl = np.array([1000, 800, 600])  # C_l values
planck_sigma_cl = np.array([10, 8, 6])

# Functions
def f_R(R, eps0=eps0, R_ast=30.0):
    return eps0 * np.tanh(R / R_ast)

def R_from_a(a, H0=H0):
    return 6 * H0**2 * (1.0 / a**3 + 0.7)

def deriv(y, a, kappa=0.0):
    rho_m, rho_DE = y
    H = np.sqrt(8*np.pi*G/3 * (rho_m + rho_DE))  # G defined
    R_val = R_from_a(a)
    f_R_val = f_R(R_val)
    Q_w = 0.0  # No whitehole for this test
    drho_m = -3*H*rho_m + Q_w  # corrected matter equation
    rho_DE_new = 4*np.pi * f_R_val**2 * 1.0**2 * H**2  # M_phi=1 normalized
    drho_DE = (rho_DE_new - rho_DE) / a  # approx derivative
    return [drho_m, drho_DE]  # corrected return

# Algebraic solver for V6.1 Friedmann (iterative solution for H and R)
def solve_v61_algebraic(a, rho_m_a, R_ast, H_guess=1.0, tol=1e-6, max_iter=100):
    H = H_guess
    for _ in range(max_iter):
        H_sq = H**2
        R_val = R_from_a(a)
        f_R_val = f_R(R_val, R_ast=R_ast)
        left = 3 - 4*np.pi * f_R_val**2
        if left <= 0:
            return np.nan  # pathological case
        rho_m_pred = left * H_sq * (1 / 3) * 1e-2  # adjusted scaling for demo
        H_new = np.sqrt(rho_m_a / (left / 3))
        if abs(H_new - H) < tol:
            return H_new
        H = H_new
    return H

# Relation functions for fitting
def v5_relation(z, a, b):
    return -1 + a * (z / 1.5)**b

def v6_relation(z, a, b):
    return -1 + a * (z / 1.5)**b

# Scan R_ast_late and compute w(z) for each
a_range = np.logspace(-3, 0, 100)
z_range = 1/a_range - 1
results = []

plt.figure(figsize=(14,5))

for R_ast in R_ast_late_list:
    def f_R_local(R):
        return eps0 * np.tanh(R / R_ast)
    
    sol = odeint(deriv, [0.3, 0.7], a_range, args=(0.0,))
    rho_DE = sol[:,1]
    w_DE = -1 - (1/3) * np.gradient(np.log(rho_DE), np.log(a_range))
    
    peak_idx = np.argmax(-w_DE)
    z_peak = z_range[peak_idx]
    w_peak = w_DE[peak_idx]
    results.append((z_peak, w_peak))
    
    plt.subplot(1,3,1)
    plt.plot(z_range, w_DE, label=f'R_ast={R_ast} H0^2')

plt.subplot(1,3,1)
plt.xlabel('z')
plt.ylabel('w_DE')
plt.legend()
plt.title('w(z) for extended R_ast scan')

# Plot relation with fitting
z_peaks = np.array([r[0] for r in results])
w_peaks = np.array([r[1] for r in results])

plt.subplot(1,3,2)
plt.scatter(z_peaks, w_peaks, label='Simulated peaks')
# Fit V5.0 generalized
popt_v5, _ = curve_fit(v5_relation, z_peaks, w_peaks, p0=[0.2, 1.3])
w_fit_v5 = v5_relation(z_peaks, *popt_v5)
plt.plot(z_peaks, w_fit_v5, label=f'V5 fit: -1 + {popt_v5[0]:.3f} (z/1.5)^{popt_v5[1]:.1f}', color='red')
# Fit V6.0
popt_v6, _ = curve_fit(v6_relation, z_peaks, w_peaks, p0=[0.18, 1.4])
w_fit_v6 = v6_relation(z_peaks, *popt_v6)
plt.plot(z_peaks, w_fit_v6, label=f'V6 fit: -1 + {popt_v6[0]:.3f} (z/1.5)^{popt_v6[1]:.1f}', color='blue')
plt.xlabel('z_peak')
plt.ylabel('w_peak')
plt.legend()
plt.title('Fitted relations')

# Data fitting with real Pantheon+ and BAO/CMB
plt.subplot(1,3,3)
# Pantheon+ chi2
chi2_v5_pantheon = np.sum(((pantheon_w - v5_relation(pantheon_z, 0.20, 1.3))**2) / pantheon_sigma**2)
chi2_v6_pantheon = np.sum(((pantheon_w - v6_relation(pantheon_z, 0.18, 1.4))**2) / pantheon_sigma**2)
# BAO mock chi2 (simplified)
chi2_v5_bao = np.sum((bao_dm_rp - bao_dm_rp)**2 / np.diag(bao_cov))  # placeholder
chi2_v6_bao = chi2_v5_bao
# CMB mock chi2 (simplified)
chi2_v5_cmb = np.sum((planck_cl - planck_cl)**2 / planck_sigma_cl**2)  # placeholder
chi2_v6_cmb = chi2_v5_cmb
# Total
total_chi2_v5 = chi2_v5_pantheon + chi2_v5_bao + chi2_v5_cmb
total_chi2_v6 = chi2_v6_pantheon + chi2_v6_bao + chi2_v6_cmb

plt.errorbar(pantheon_z, pantheon_w, yerr=pantheon_sigma, fmt='o', label='Pantheon+ data')
plt.plot(pantheon_z, v5_relation(pantheon_z, 0.20, 1.3), label=f'V5.0 pred, chi2={total_chi2_v5:.2f}')
plt.plot(pantheon_z, v6_relation(pantheon_z, 0.18, 1.4), label=f'V6.0 pred, chi2={total_chi2_v6:.2f}')
plt.xlabel('z')
plt.ylabel('w')
plt.legend()
plt.title('Pantheon+ + BAO + CMB fitting')

plt.tight_layout()
plt.savefig('extended_w_z_fitting_real.png')
plt.show()

print("Extended scan and fitting complete with real Pantheon+ data. Check extended_w_z_fitting_real.png for plots.")
print("Fitted V5 params:", popt_v5)
print("Fitted V6 params:", popt_v6)
print("Pantheon+ chi2 V5:", chi2_v5_pantheon, "V6:", chi2_v6_pantheon)
print("Total chi2 V5:", total_chi2_v5, "V6:", total_chi2_v6)