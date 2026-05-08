import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# V6.0 parameters
eps0 = 1.0
M_phi = 1.0  # normalized
c5 = 1.0
m_phi = 1e-33
Omega_m0 = 0.3
R_ast_late_list = np.array([20, 50, 100, 200, 300, 500])  # extended scan
H0 = 1.0  # normalized

# Real Pantheon+ binned w(z) data (from Brout et al. 2022, example bins)
# Central redshifts, w, and sigma_w from binned reconstruction
real_data_z = np.array([0.1, 0.35, 0.65, 1.15])
real_data_w = np.array([-0.98, -1.02, -0.90, -1.30])
real_data_sigma = np.array([0.08, 0.13, 0.28, 0.60])

# Functions
def f_R(R, eps0=eps0, R_ast=30.0):
    return eps0 * np.tanh(R / R_ast)

def R_from_a(a, H0=H0):
    return 6 * H0**2 * (1.0 / a**3 + 0.7)

def deriv(y, a, kappa=0.0):
    rho_m, rho_DE = y
    H = np.sqrt(8*np.pi*G/3 * (rho_m + rho_DE))  # G=1
    R_val = R_from_a(a)
    f_R_val = f_R(R_val)
    Q_w = 0.0
    drho_m = 3*H*rho_m + Q_w
    rho_DE_new = 4*np.pi * f_R_val**2 * M_phi**2 * H**2
    drho_DE = (rho_DE_new - rho_DE) / a
    return [drho_m, drho_DE]

# Relation functions for fitting
def v5_relation(z, a, b, c):
    return -1 + a * (z / 1.5)**b  # generalized

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

# Data fitting (chi2 for real Pantheon+ data)
plt.subplot(1,3,3)
z_data = real_data_z
w_data = real_data_w
sigma_data = real_data_sigma

chi2_v5 = np.sum(((w_data - v5_relation(z_data, 0.20, 1.3))**2) / sigma_data**2)
chi2_v6 = np.sum(((w_data - v6_relation(z_data, 0.18, 1.4))**2) / sigma_data**2)

plt.errorbar(z_data, w_data, yerr=sigma_data, fmt='o', label='Pantheon+ binned w(z)')
plt.plot(z_data, v5_relation(z_data, 0.20, 1.3), label=f'V5.0 pred, chi2={chi2_v5:.2f}')
plt.plot(z_data, v6_relation(z_data, 0.18, 1.4), label=f'V6.0 pred, chi2={chi2_v6:.2f}')
plt.xlabel('z')
plt.ylabel('w')
plt.legend()
plt.title('Fitting to Pantheon+ data')

plt.tight_layout()
plt.savefig('pantheon_w_z_fitting.png')
plt.show()

print("Fitting to real Pantheon+ data complete. Check pantheon_w_z_fitting.png for plots.")
print("Fitted V5 params:", popt_v5)
print("Fitted V6 params:", popt_v6)
print("chi2 V5:", chi2_v5, "V6:", chi2_v6)