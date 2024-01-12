import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

def mag_and_angle_to_cmplx(mag, angle):
    return mag * np.exp(1j * angle)

# Parameters
fault_start = 1
fault_end = 1.05

fn = 60
H_gen = 3.5
X_gen = 0.2
X_ibb = 0.1
X_line = 0.65

# Values are initialized from loadflow
E_fd_gen = 1.075
E_fd_ibb = 1.033
P_m_gen = 1998/2200

# init states of variables
omega_gen_init = 0 # init state
delta_gen_init = np.deg2rad(45.9) # init state
delta_ibb_init = np.deg2rad(-5.0) # init state

# Number of steps for output
t_final = 5
steps = 1000
P_e = []

# Define the ODE system
def ODE_system(state, t, T_m_gen, P_e_gen):
    global H_gen
    omega, delta = state
    
    d_omega_dt = 1 / (2 * H_gen) * (T_m_gen(omega) - P_e_gen(delta, t))
    d_delta_dt = omega

    P_e.append(P_e_gen(delta, t))

    return [d_omega_dt, d_delta_dt]

# External functions for dependencies
def T_m_gen(omega):
    # Assuming a simple linear function for demonstration purposes
    return P_m_gen / (1 + omega)

def P_e_gen(delta, z):
    if fault_start <= t < fault_end: # sc_on = True
        y_adm = np.array([[(-1j / X_gen - 1j / X_line) + 1000000, 1j / X_line],[1j / X_line, -1j / X_line - 1j / X_ibb]])
    else: # sc_on = False
        y_adm = np.array([[-1j / X_gen - 1j / X_line, 1j / X_line],[1j / X_line, -1j / X_line - 1j / X_ibb]])

    # Calculate the inverse of the admittance matrix (Y^-1)
    y_inv = np.linalg.inv(y_adm)

    # Calculate current injections of the generator and the infinite busbar
    i_inj_gen = mag_and_angle_to_cmplx(E_fd_gen, delta) / (1j * X_gen)
    i_inj_ibb = mag_and_angle_to_cmplx(E_fd_ibb, delta_ibb_init) / (1j * X_ibb)

    # Calculate voltages at the bus by multiplying the inverse of the admittance matrix with the current injections
    v_bb_gen = y_inv[0, 0] * i_inj_gen + y_inv[0, 1] * i_inj_ibb
    v_bb_ibb = y_inv[1, 0] * i_inj_gen + y_inv[1, 1] * i_inj_ibb

    E_gen_cmplx = mag_and_angle_to_cmplx(E_fd_gen, delta)
    
    return (v_bb_gen * np.conj((E_gen_cmplx - v_bb_gen) / (1j * X_gen))).real

# Time vector
t = np.linspace(0, t_final, steps)

# Initial conditions
initial_conditions = [omega_gen_init, delta_gen_init]

# Integrate the ODE system
solution = odeint(ODE_system, initial_conditions, t, args=(T_m_gen,P_e_gen)).T

# Output
omega = solution[:, 0]
delta = solution[:, 1]

# print('The ODE system solution is:\n')
# print('t, omega, delta\n')
# for t, om, delt in zip(t, omega, delta):
#     print(f'{t}, {om}, {delt}')

plt.plot(t, omega, label='omega')
plt.legend()
plt.title('Solution of Blackbox.ai - Omega with odeint')
plt.show()

plt.plot(t, delta, label='delta')
plt.legend()
plt.title('Solution of Blackbox.ai - Delta with odeint')
plt.show()