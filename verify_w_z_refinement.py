import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# V6.0 parameters
eps0 = 1.0
M_phi = 1.0  # normalized
c5 = 1.0
m_phi = 1e-33
Omega_m0 = 0.3
R_ast_late_list = np.array([20, 50, 100, 200]) * 1.0  # in H0^2 units, scan
H0 = 1.0  # normalized

# Functions
def f_R(R, eps0=eps0, R_ast=30.0):
    return eps0 * np.tanh(R / R_ast)

def R_from_a(a, H0=H0):
    # Approx R = 6 H0^2 (1/a^3 + 0.7) for late universe
    return 6 * H0**2 * (1.0 / a**3 + 0.7)

def deriv(y, a, kappa=0.0):  # kappa for compatibility, but V6.0 uses algebraic rho_DE
    rho_m, rho_DE = y
    H = np.sqrt(8*np.pi*G/3 * (rho_m + rho_DE))  # G=1
    R_val = R_from_a(a)
    f_R_val = f_R(R_val)
    # For V6.0: rho_DE is algebraic, but simulate evolution approx
    Q_w = 0.0  # No whitehole for this test
    drho_m = 3*H*rho_m + Q_w
    # rho_DE evolves as per V6.0: track H^2
    rho_DE_new = 4*np.pi * f_R_val**2 * M_phi**2 * H**2
    drho_DE = (rho_DE_new - rho_DE) / a  # approx derivative
    return [drho_m, drho_DE]

# Scan R_ast_late and compute w(z) for each
a_range = np.logspace(-3, 0, 100)  # a from 0.001 to 1
z_range = 1/a_range - 1
results = []

plt.figure(figsize=(12,5))

for R_ast in R_ast_late_list:
    # Redefine f_R with this R_ast
    def f_R_local(R):
        return eps0 * np.tanh(R / R_ast)
    
    sol = odeint(deriv, [0.3, 0.7], a_range, args=(0.0,))
    rho_m = sol[:,0]
    rho_DE = sol[:,1]
    
    # Compute w(z)
    H_list = np.sqrt(8*np.pi/3 * (rho_m + rho_DE))
    w_eff = (rho_m + rho_DE) / (3 * rho_m + 4 * rho_DE) - 1  # approx
    w_DE = -1 - (1/3) * np.gradient(np.log(rho_DE), np.log(a_range))
    
    # Find peak
    peak_idx = np.argmax(-w_DE)  # w_DE most negative
    z_peak = z_range[peak_idx]
    w_peak = w_DE[peak_idx]
    results.append((z_peak, w_peak))
    
    plt.subplot(1,2,1)
    plt.plot(z_range, w_DE, label=f'R_ast={R_ast} H0^2, z_peak={z_peak:.1f}, w_peak={w_peak:.3f}')

# Plot V5.0 vs V6.0 relation
z_peaks = np.array([r[0] for r in results])
w_peaks = np.array([r[1] for r in results])

plt.subplot(1,2,1)
plt.xlabel('z')
plt.ylabel('w_DE')
plt.legend()
plt.title('w(z) for different R_ast_late')

plt.subplot(1,2,2)
plt.scatter(z_peaks, w_peaks, label='Simulated points')
# V5.0 relation: w_peak = -1 + 0.20 * (z_peak / 1.5)**1.3
z_fit = np.linspace(0.5, 3, 100)
w_v5 = -1 + 0.20 * (z_fit / 1.5)**1.3
plt.plot(z_fit, w_v5, label='V5.0: -1 + 0.20 (z/1.5)^1.3', color='red')
# V6.0 relation: w_peak = -1 + 0.18 * (z_peak / 1.5)**1.4
w_v6 = -1 + 0.18 * (z_fit / 1.5)**1.4
plt.plot(z_fit, w_v6, label='V6.0: -1 + 0.18 (z/1.5)^1.4', color='blue')
plt.xlabel('z_peak')
plt.ylabel('w_peak')
plt.legend()
plt.title('w_peak vs z_peak relation')

plt.tight_layout()
plt.savefig('w_z_refinement.png')
plt.show()

print("Verification complete. Check w_z_refinement.png for plots.")
print("Results:", results)