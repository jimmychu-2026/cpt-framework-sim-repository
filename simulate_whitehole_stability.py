import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# V6.1 parameters (natural units: G=1)
G = 1.0
eps0 = 1.0
R_ast_late = 30.0
M_phi = 1.0
c5 = 1.0
m_phi = 1e-33
Omega_m0 = 0.3
kappa_white_list = [0.0, 1e-10, 1e-9, 1e-8]

# Functions
def f_R(R, eps0=eps0, R_ast=R_ast_late):
    return eps0 * np.tanh(R / R_ast)

def R_from_a(a, H0=1.0):
    return 6 * H0**2 * (1.0 / a**3 + 0.7)

def Q_white(rho_DE, kappa, f_R_val, H):
    return -kappa * f_R_val * H * np.sqrt(np.abs(rho_DE))

def deriv(y, a, kappa):
    rho_m, rho_DE = y
    H = np.sqrt(8*np.pi*G/3 * (rho_m + rho_DE))  # G defined
    R_val = R_from_a(a)
    f_R_val = f_R(R_val)
    Q_w = Q_white(rho_DE, kappa, f_R_val, H)
    drho_m = -3*H*rho_m + Q_w  # corrected matter equation
    w_phi = -0.9  # defined before use
    drho_DE = 3*H*(1 + w_phi)*rho_DE - Q_w  # corrected
    return [drho_m, drho_DE]  # corrected return

# Stability check: estimate M_eff^2 from perturbation
def stability_check(kappa, a=1.0):
    H = np.sqrt(8*np.pi*G/3 * (0.3 + 0.7))  # G defined
    R_val = R_from_a(a)
    f_R_val = f_R(R_val)
    M_eff2 = m_phi**2 + (c5 * f_R_val / M_phi)**2 * kappa**2 * H**2
    c_s2 = 1 + kappa**2 * (f_R_val / M_phi)**2  # corrected identifier
    return M_eff2 > 0 and c_s2 > 0

# Simulation
a_range = np.linspace(0.001, 1.0, 100)
plt.figure(figsize=(12,5))

# Plot 1: rho_m evolution
plt.subplot(1,2,1)
for kappa in kappa_white_list:
    sol = odeint(deriv, [0.3, 0.7], a_range, args=(kappa,))
    plt.plot(a_range, sol[:,0], label=f'kappa={kappa}')
plt.title('Matter Density rho_m (with Whitehole Evaporation)')
plt.xlabel('a')
plt.ylabel('rho_m')
plt.legend()

# Plot 2: Stability summary
plt.subplot(1,2,2)
stable = [all(stability_check(kappa, a) for a in a_range) for kappa in kappa_white_list]
plt.bar(['0', '1e-10', '1e-9', '1e-8'], stable, color=['green' if s else 'red' for s in stable])
plt.title('Stability Check (M_eff^2 >0 & c_s^2 >0)')
plt.ylabel('Stable?')
plt.savefig('whitehole_stability.png')
plt.show()

print("Simulation complete. Check whitehole_stability.png for plots.")