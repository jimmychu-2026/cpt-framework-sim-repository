import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# V5.0 parameters
eps0 = 1.0
R_ast_late = 30.0
M_phi = 1.0
c5 = 1.0
m_phi = 1e-33
Omega_m0 = 0.3

# Functions
def f_R(R, eps0=eps0, R_ast=R_ast_late):
    return eps0 * np.tanh(R / R_ast)

def R_from_a(a, H0=1.0):
    return 6 * H0**2 * (1.0 / a**3 + 0.7)

def Q_white(rho_DE, kappa, f_R_val, H):
    return -kappa * f_R_val * H * np.sqrt(np.abs(rho_DE))

def deriv(y, a, kappa):
    rho_m, rho_DE = y
    H = np.sqrt(8*np.pi*G/3 * (rho_m + rho_DE))  # G=1
    R_val = R_from_a(a)
    f_R_val = f_R(R_val)
    Q_w = Q_white(rho_DE, kappa, f_R_val, H)
    drho_m = 3*H*rho_m + Q_w
    drho_DE = 3*H*(1 + w_phi)*rho_DE - Q_w
    w_phi = -0.9
    return [drho_m, rho_DE]

# Stability metrics
def stability_metrics(kappa, a_range):
    M_eff2_list = []
    c_s2_list = []
    for a in a_range:
        H = np.sqrt(8*np.pi*G/3 * (0.3 + 0.7))  # approx
        R_val = R_from_a(a)
        f_R_val = f_R(R_val)
        M_eff2 = m_phi**2 + (c5 * f_R_val / M_phi)**2 * kappa**2 * H**2
        c_s2 = 1 + kappa**2 * (f_R_val / M_phi)**2
        M_eff2_list.append(M_eff2)
        c_s2_list.append(c_s2)
    return np.array(M_eff2_list), np.array(c_s2_list)

# Comparison
kappa_no = 0.0  # No inverse Hawking
kappa_yes = 1e-10  # With inverse Hawking (small kappa)
a_range = np.linspace(0.001, 1.0, 100)

plt.figure(figsize=(14,5))

# Plot 1: rho_m comparison
plt.subplot(1,3,1)
for kappa, label in [(kappa_no, 'No Inverse Hawking'), (kappa_yes, 'With Inverse Hawking')]:
    sol = odeint(deriv, [0.3, 0.7], a_range, args=(kappa,))
    plt.plot(a_range, sol[:,0], label=label)
plt.title('rho_m Evolution')
plt.xlabel('a')
plt.ylabel('rho_m')
plt.legend()

# Plot 2: M_eff^2 comparison
plt.subplot(1,3,2)
M_no, _ = stability_metrics(kappa_no, a_range)
M_yes, _ = stability_metrics(kappa_yes, a_range)
plt.plot(a_range, M_no, label='No Inverse Hawking', color='blue')
plt.plot(a_range, M_yes, label='With Inverse Hawking', color='red')
plt.title('Effective Mass Squared M_eff^2')
plt.xlabel('a')
plt.ylabel('M_eff^2')
plt.legend()

# Plot 3: c_s^2 comparison
plt.subplot(1,3,3)
_, c_no = stability_metrics(kappa_no, a_range)
_, c_yes = stability_metrics(kappa_yes, a_range)
plt.plot(a_range, c_no, label='No Inverse Hawking', color='blue')
plt.plot(a_range, c_yes, label='With Inverse Hawking', color='red')
plt.title('Sound Speed Squared c_s^2')
plt.xlabel('a')
plt.ylabel('c_s^2')
plt.legend()

plt.tight_layout()
plt.savefig('compare_stability.png')
plt.show()

print("Comparison complete. Check compare_stability.png for plots.")