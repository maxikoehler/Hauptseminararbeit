import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.integrate import odeint, RK45

# helping function for calculation with complex numbers
def mag_and_angle_to_cmplx(mag, angle):
    return mag * np.exp(1j * angle)

# variable setting
fn = 60
H_gen = 3.5
X_gen = 0.2
X_ibb = 0.1
X_line = 0.65

# Values are initialized from loadflow
E_fd_gen = 1.075
E_fd_ibb = 1.033
P_m_gen = 1998/2200
# P_m_gen = 0.6

omega_gen_init = 0
delta_gen_init = np.deg2rad(45.9)
delta_ibb_init = np.deg2rad(-5.0)

# simulation setups
t_start = 0
t_end = 200
t_step = 0.005

fault_start = 1
fault_end = t_end

def P_e(E_1, E_2, X, delta):
    P_e_gen = (E_1 * E_2) / X * np.sin(delta)
    return P_e_gen

def P_m(omega):
    global P_m_gen
    P_t = P_m_gen / (1 + omega)
    return P_m_gen

def ODE_system(state, t):
    global H_gen
    global E_fd_gen
    global E_fd_ibb
    global X_gen
    global X_line
    global fault_start
    global fault_end

    if fault_start <= t < fault_end:
        X = 0
        E_ibb = 0
    else:
        X = X_line
        E_ibb = E_fd_ibb

    omega, delta = state
    
    d_omega_dt = 1 / (2 * H_gen) * (P_m(omega) - P_e(E_fd_gen, E_ibb, (X_gen + X), (delta)))
    d_delta_dt = omega

    return [d_omega_dt, d_delta_dt]

def rk_system(t, state):
    return ODE_system(state, t)

def P_r_deg(delta, E1, E2, X):
    # determing the P_e curve under input in degrees
    return (E1*E2)/(X)*np.sin(delta) - P_m_gen

def P_t_deg(x):
    # determing the P_t curve under input in degrees
    return P_m_gen*np.ones(np.size(x))

def stability_eac_clearing(delta_0, delta_act, delta_max, E_p, E_k_pre, E_k_post, X_pre, X_post):
    # Compare the acceleration area until the given delta and compare it to the braking area left until the dynamic stability point is passed
    area_acc = sp.integrate.quad(P_r_deg, delta_0, delta_act, args=(E_p, E_k_post, X_post))
    area_dec = sp.integrate.quad(P_r_deg, delta_act, delta_max, args=(E_p, E_k_pre, X_pre))

    if abs(area_acc[0]) < abs(area_dec[0]): # True: stable, False: NOT stable 
        return True
    else:
        return False

def stability_eac_non_clearing(delta_0, delta_act, delta_max, E_p, E_k_post, X_post):
    # Compare the acceleration area until the given delta and compare it to the braking area left until the dynamic stability point is passed
    area_acc = sp.integrate.quad(P_r_deg, delta_0, delta_act, args=(E_p, E_k_post, X_post))
    area_dec = sp.integrate.quad(P_r_deg, delta_act, delta_max, args=(E_p, E_k_post, X_post))

    if abs(area_acc[0]) < abs(area_dec[0]): # True: stable, False: NOT stable 
        return True
    else:
        return False

def determine_cct(mode, t_sim, delta, delta_0, E_p, E_k_pre, E_k_post, X_pre, X_post):
    # t_sim and delta are result arrays
    # delta_0 is the initial angle delta of the stable system pre-fault

    # Save current time and delta at time point i; iterate through i to test any given time until stability can't be remained; delta_cc and t_cc is the angle and time at the last stable point
    if mode == 'c':
        delta_max = np.pi - delta_0

        i = 0
        t_cc, delta_cc = -1, -1
        while stability_eac_clearing(delta_0, delta[i], delta_max, E_p, E_k_pre, E_k_post, X_pre, X_post) and i < np.size(t_sim)-1:
            t_cc = t_sim[i]
            delta_cc = delta[i]
            i = i + 1

        if t_cc < 0:
            return False, -1, -1
        else:
            return True, t_cc, delta_cc
            
    elif mode == 'nc':
        delta_max = np.pi - np.arcsin(P_m_gen/P_e(E_p, E_k_post, X_post, np.pi/2))

        if P_m_gen > P_e(E_p, E_k_post, X_post, np.pi/2):
            return False

        for i in range(np.size(t_sim)):
            if stability_eac_non_clearing(delta_0, delta[i], delta_max, E_p, E_k_post, X_post) == False:
                return False
        return True

def plot_P_delta(mode, E_p, E_k_pre, E_k_post, X_pre, X_post):
    # plotting of P - delta curve and the relavant areas. Consideration of pre- and post-fault function as partially sin-function in [0, 2*pi], as well as constant power function of the turbine
    # marking the critical power angle delta_c

    # calculation of P_e_pre, P_e_post, and P_t
    x_deg = np.linspace(0, 180) # linear vector for plotting in deg
    x_rad = np.linspace(0, np.pi) # linear vector for calculation in rad
    P_e_pre = P_e(E_p, E_k_pre, X_pre, x_rad)
    P_e_post = P_e(E_p, E_k_post, X_post, x_rad)
    P_t = P_t_deg(x_rad)

    # setting up the plot and plotting P_e and P_t
    fig, ax = plt.subplots()
    ax.plot(x_deg, P_e_pre, '-', linewidth=2, label='$P_\mathrm{e}$ pre-fault')
    if P_e(E_p, E_k_post, X_post, np.pi/2)!= 0:
        ax.plot(x_deg, P_e_post, '-', linewidth=2, label='$P_\mathrm{e}$ post-fault')
    ax.plot(x_deg, P_t, '-', linewidth=2, label='$P_\mathrm{T}$ of the turbine')

    # Determing the boundary angles
    delta_0 = np.arcsin(P_m_gen/P_e(E_p, E_k_pre, X_pre, np.pi/2))

    ax.set_ylim(bottom=0)
    if mode == 'c':
        delta_max = np.pi - delta_0
        stability, t_cc, delta_c = determine_cct(mode, t_sim, delta, delta_0, E_p, E_k_pre, E_k_post, X_pre, X_post)

        delta_0_deg = int(np.rad2deg(delta_0))
        delta_max_deg = int(np.rad2deg(delta_max))
        delta_c_deg = int(np.rad2deg(delta_c))

        if stability:
            # Make the shaded region for area_acc
            ix1 = np.linspace(delta_0_deg, delta_c_deg)
            iy1 = P_e(E_p, E_k_post, X_post, np.deg2rad(ix1))
            ax.fill_between(ix1, iy1, P_m_gen, facecolor='0.9', edgecolor='0.5')

            # Make the shaded region for area_dec, https://matplotlib.org/stable/gallery/lines_bars_and_markers/fill_between_demo.html
            ix2 = np.linspace(delta_c_deg, delta_max_deg) # -> does this have to be in rad or in deg?
            iy2 = P_e(E_p, E_k_pre, X_pre, np.deg2rad(ix2))
            ax.fill_between(ix2, iy2, P_m_gen, facecolor='0.9', edgecolor='0.5')

            ax.set_xticks([0, 180, delta_0_deg, delta_c_deg, delta_max_deg], labels=['0', '180', '$\delta_\mathrm{0}$', '$\delta_\mathrm{c}$', '$\delta_\mathrm{max}$'])
        else:
            ix1 = np.linspace(delta_0_deg, delta_max_deg)
            iy1 = np.zeros(np.size(ix1))
            ax.fill_between(ix1, iy1, P_m_gen, facecolor='red', edgecolor='red', alpha=0.5)
    elif mode == 'nc':
        if P_m_gen > P_e(E_p, E_k_post, X_post, np.pi/2):
            delta_max = np.pi
        else:
            delta_max = np.pi - np.arcsin(P_m_gen/P_e(E_p, E_k_post, X_post, np.pi/2))
        
        stability = determine_cct(mode, t_sim, delta, delta_0, E_p, E_k_pre, E_k_post, X_pre, X_post)

        delta_0_deg = int(np.rad2deg(delta_0))
        delta_max_deg = int(np.rad2deg(delta_max))
        
        if stability:
            # Make the shaded region for area_acc
            ix1 = np.linspace(delta_0_deg, delta_max_deg) # -> does this have to be in rad or in deg?
            iy1 = P_e(E_p, E_k_post, X_post, np.deg2rad(ix1))
            ax.fill_between(ix1, iy1, P_m_gen, facecolor='0.9', edgecolor='0.5')

            ax.set_xticks([0, 180, delta_0_deg, delta_max_deg], labels=['0', '180', '$\delta_\mathrm{0}$', '$\delta_\mathrm{max}$'])
        else:
            # Make the shaded region for area_acc
            ix1 = np.linspace(delta_0_deg, delta_max_deg) # -> does this have to be in rad or in deg?
            iy1 = P_e(E_p, E_k_post, X_post, np.deg2rad(ix1))
            ax.fill_between(ix1, iy1, P_m_gen, facecolor='red', edgecolor='red', alpha=.5)

            ax.set_xticks([0, 180, delta_0_deg, delta_max_deg], labels=['0', '180', '$\delta_\mathrm{0}$', '$\delta_\mathrm{max}$'])

        
    # adding area descriptions
    # ax.text(0.5*(delta_0_deg + delta_c_deg), P_m_gen-0.1, r"$A_\mathrm{acc}$", horizontalalignment='center', fontsize=15)
    # ax.text(0.47*(delta_max_deg + delta_c_deg), P_m_gen+0.05, r"$A_\mathrm{dec}$", horizontalalignment='center', fontsize=15)
    # 0.3*(P_e(E_p, E_k_pre, X_pre, np.pi/2) - P_m_gen)

    # axis text and optical manipulation of the plot
    fig.text(0.9, 0.05, '$\delta$ in $^\circle$')
    fig.text(0.1, 0.9, '$P$ in $pu$')

    plt.legend()
    plt.show()
    return



if __name__ == "__main__":
    # setup simulation inputs 
    t_sim = np.arange(t_start, t_end, t_step)
    initial_conditions = [omega_gen_init, delta_gen_init]

    # solve ODE with python solver
    solution = odeint(ODE_system, initial_conditions, t_sim)

    # save solution
    omega = solution[:, 0]
    delta = solution[:, 1]

    determine_cct('c', t_sim, delta, delta_gen_init, E_fd_gen, E_fd_ibb, 0, (X_gen + X_line), X_gen)