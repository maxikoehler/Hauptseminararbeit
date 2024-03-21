############################
# Base module with the smib model.
# Input of the interested variables delta_0, E_bus, E_gen, P_m, ..., fault_start, fault_end, ...
# Export of stability, t_cc, delta_cc, t_sim, delta(t_sim), omega(t_sim), 
# Condsideration of TDS in just stable and just unstable regime
############################

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib as mpl
import numpy as np
import scipy as sp
from scipy.integrate import odeint

# redefining plot save parameters
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Charter"],
    "font.size": 10
})

# uncomment for updating savefig options for latex export
# mpl.use("pgf")

# helping function for calculation with complex numbers
def mag_and_angle_to_cmplx(mag, angle):
    return mag * np.exp(1j * angle)

def algebraic(delta_gen, fault_on):
    global E_fd_gen
    global E_fd_ibb
    global delta_ibb_init
    global X_gen, X_line, X_ibb, X_trans

    # If the SC is on, the admittance matrix is different.
    # The SC on busbar 0 is expressed in the admittance matrix as a very large admittance (1000000) i.e. a very small impedance.
    if fault_on:
        y_adm = np.array([X_fault,
                          [1j / X_line, -1j / X_line - 1j / X_ibb]])
    else:
        y_adm = np.array([[-1j / X_gen - 1j / X_line, 1j / X_line],
                          [1j / X_line, -1j / X_line - 1j / X_ibb]])

    # Calculate the inverse of the admittance matrix (Y^-1)
    y_inv = np.linalg.inv(y_adm)

    # Calculate current injections of the generator and the infinite busbar
    i_inj_gen = mag_and_angle_to_cmplx(E_fd_gen, delta_gen) / (1j * X_gen)
    i_inj_ibb = mag_and_angle_to_cmplx(E_fd_ibb, delta_ibb_init) / (1j * X_ibb)

    # Calculate voltages at the bus by multiplying the inverse of the admittance matrix with the current injections
    v_bb_gen = y_inv[0, 0] * i_inj_gen + y_inv[0, 1] * i_inj_ibb
    v_bb_ibb = y_inv[1, 0] * i_inj_gen + y_inv[1, 1] * i_inj_ibb

    return v_bb_gen

def P_e(delta, fault_on):
    # function for determing P_e WITHOUT algebraic help
    global X_gen
    global X_line
    global X_fault
    global E_fd_gen
    global E_fd_ibb

    if fault_on:
        X = 1
        E_ibb = 0
    else:
        X = X_gen + X_line + X_ibb
        E_ibb = E_fd_ibb
        
    P_e_gen = E_fd_gen * E_ibb / X * np.sin(delta)
    return P_e_gen

def P_e_alg(delta, fault_on):
    # function for determing P_e WITH algebraic help
    global E_fd_gen
    global X_gen
    
    v_bb_gen = algebraic(delta, fault_on)

    E_gen_complex = mag_and_angle_to_cmplx(E_fd_gen, delta)
    P_e_gen = (v_bb_gen * np.conj((E_gen_complex - v_bb_gen) / (1j * X_gen))).real
    return P_e_gen

def P_m(omega):
    # returning the torque of the generator, depending on the rotor speed
    global P_m_gen
    global omega_gen_init
    P_t = P_m_gen / (1 + (omega_gen_init + omega))
    return P_t

def get_max_delta(gen_parameters, sim_parameters, alg):
    init(gen_parameters, sim_parameters)

    area_acc = sp.integrate.quad(P_r_deg, delta_gen_init, delta_0_fault, args=(0, True, alg))
    area_dec = [0, 0]
    max_delta = delta_0_fault
    while abs(area_dec[0]) <= abs(area_acc[0]):
        area_dec = sp.integrate.quad(P_r_deg, delta_0_fault, max_delta, args=(0, True, alg))
        max_delta = max_delta + 0.01

    return max_delta

def get_delta_0(gen_parameters, sim_parameters, alg):
    init(gen_parameters, sim_parameters)
    x_rad = np.linspace(0,np.pi/2,360)
    delta = -1
    for x in x_rad:
        if abs(P_r_deg(x, 0, False, alg)) <= 0.01:
            delta = x

    return delta

# function for using odeint as ode-solver 
def ODE_system(state, t, fault_start, fault_end, alg):

    omega, delta = state

    global H_gen
    global E_fd_gen
    global E_fd_ibb
    global X_gen
    global X_line
    global fn

    if fault_start <= t < fault_end:
        fault_on = True
        # P_e_conv = P_e(E_fd_gen, 0, X_gen, delta)
    else:
        fault_on = False
        # P_e_conv = P_e(E_fd_gen, E_fd_ibb, X_gen + X_line, delta)

    # including time dependent solving of algebraic equations
    if alg:
        P_e_gen = P_e_alg(delta, fault_on)
    else:
        P_e_gen = P_e(delta, fault_on)
    
    d_omega_dt = 1 / (2 * H_gen) * (P_m(omega) - P_e_gen)
    d_delta_dt = omega * 2 * np.pi * fn

    return [d_omega_dt, d_delta_dt]

# functions for determing the critical clearing time
def P_r_deg(delta, omega, fault_on, alg):
    # determing the P_e curve under input in degrees
    if alg:
        P_r = P_e_alg(delta, fault_on) - P_m(omega)
    else:
        P_r = P_e(delta, fault_on) - P_m(omega)
    return P_r

def P_t_deg(x):
    # determing the P_t curve under input in degrees
    global P_m_gen

    return P_m_gen*np.ones(np.size(x))

def stability_eac(delta_0, delta_act, omega_act, delta_max, alg):
    # global delta_new, omega_new

    # Compare the acceleration area until the given delta and compare it to the braking area left until the dynamic stability point is passed
    area_acc = sp.integrate.quad(P_r_deg, delta_0, delta_act, args=(omega_act, True, alg))
    area_dec = sp.integrate.quad(P_r_deg, delta_act, delta_max, args=(omega_act, (not clearing), alg))

    if abs(area_acc[0]) < abs(area_dec[0]): # True: stable, False: NOT stable 
        return True
    else:
        return False

def determine_cct(t_sim, delta, omega, delta_0, alg):
    # t_sim and delta are result arrays
    # delta_0 is the initial angle delta of the stable system pre-fault

    # Save current time and delta at time point i; iterate through i to test any given time until stability can't be remained; delta_cc and t_cc is the angle and time at the last stable point
    global delta_max_fault
    if clearing:
        delta_max = np.pi - delta_0
    else:
        delta_max = delta_max_fault

    i = 0
    t_cc, delta_cc, omega_cc = -1, -1, -1

    while stability_eac(delta_0, delta[i], omega[i], delta_max, alg) and i < np.size(t_sim)-1 and delta[i] < delta_max_fault:
        t_cc = t_sim[i]
        delta_cc = delta[i]
        omega_cc = omega[i]
        i = i + 1

    if t_cc < 0:
        return False, -1, -1, -1
    else:
        if clearing:
            return True, t_cc, delta_cc, omega_cc
        else:
            return True, t_cc, delta_0_fault, omega_gen_init

# execution functions for simulation
def do_sim(gen_parameters, sim_parameters, alg):
    init(gen_parameters, sim_parameters)

    # setup simulation inputs 
    t_sim = np.arange(sim_parameters["t_start"], sim_parameters["t_end"], sim_parameters["t_step"])
    initial_conditions = [gen_parameters["omega_gen_init"], gen_parameters["delta_gen_init"]]

    delta_0 = gen_parameters["delta_gen_init"]
    # delta_max = np.pi - delta_0

    for i in range(1,4,1):
        if i == 1: # first TDS with no fault-clearing
            # solve ODE with python solver
            solution = odeint(ODE_system, initial_conditions, t_sim, args=(sim_parameters["fault_start"], sim_parameters["fault_end"], alg))
            stability, t_cc, delta_cc, omega_cc = determine_cct(t_sim, solution[:, 1], solution[:, 0], delta_0, alg)
        elif i == 2: # second TDS with fault clearing just right
            # solve ODE with python solver
            fault_end = t_cc - 5 * sim_parameters["t_step"]
            solution_stable = odeint(ODE_system, initial_conditions, t_sim, args=(sim_parameters["fault_start"], fault_end, alg))
        elif i == 3: # second TDS with fault clearing just NOT right
            fault_end = t_cc + 2 * sim_parameters["t_step"]
            solution_unstable = odeint(ODE_system, initial_conditions, t_sim, args=(sim_parameters["fault_start"], fault_end, alg))

    return stability, t_cc, delta_cc, t_sim, solution_stable, solution_unstable

def do_sim_simple(gen_parameters, sim_parameters, alg):
    init(gen_parameters, sim_parameters)

    # setup simulation inputs 
    t_sim = np.arange(t_start, t_end, t_step)
    initial_conditions = [omega_gen_init, delta_gen_init]

    delta_0 = delta_gen_init

    solution = odeint(ODE_system, initial_conditions, t_sim, args=(fault_start, fault_end, alg))
    stability, t_cc, delta_cc, omega_cc = determine_cct(t_sim, solution[:, 1], solution[:, 0], delta_0, alg)

    return stability, t_cc, delta_cc, t_sim, solution

def init(gen_parameters, sim_parameters):
    global fn, H_gen, X_gen, X_ibb, X_line, X_trans, X_fault, E_fd_gen, E_fd_ibb, P_m_gen, omega_gen_init, delta_gen_init, delta_ibb_init, t_start, t_end, t_step, fault_start, fault_end, clearing

    fn = gen_parameters["fn"]
    H_gen = gen_parameters["H_gen"]
    X_gen = gen_parameters["X_gen"]
    X_ibb = gen_parameters["X_ibb"]
    X_line = gen_parameters["X_line"]
    X_fault = gen_parameters["X_fault"]
    X_trans = gen_parameters["X_trans"]

    E_fd_gen = gen_parameters["E_fd_gen"]
    E_fd_ibb = gen_parameters["E_fd_ibb"]
    P_m_gen = gen_parameters["P_m_gen"]

    omega_gen_init = gen_parameters["omega_gen_init"]
    delta_gen_init = gen_parameters["delta_gen_init"]
    delta_ibb_init = gen_parameters["delta_ibb_init"]

    t_start = sim_parameters["t_start"]
    t_end = sim_parameters["t_end"]
    t_step = sim_parameters["t_step"]

    fault_start = sim_parameters["fault_start"]
    fault_end = sim_parameters["fault_end"]
    clearing = sim_parameters["clearing"]

    # assessment of delta_0 and delta_max in fault case
    x_rad = np.linspace(0,np.pi, 360)
    i = 0
    global delta_0_fault, delta_max_fault
    delta_0_fault = np.pi
    delta_max_fault = np.pi
    while i < np.size(x_rad)/2:
        if abs(P_r_deg(x_rad[i], 0, True, True)) < 0.01:
            delta_0_fault = x_rad[i]
            delta_max_fault = np.pi - delta_0_fault
        i = i + 1

    return

if __name__ == "__main__":
    # setup simulation inputs
    gen_parameters = {
        "fn":       50,
        "H_gen":    3.3,
        "X_gen":    0.2,
        "X_trans":  0.1,
        "X_ibb":    0.1,
        "X_line":   0.65,
        "X_fault":  0.0001,

        "E_fd_gen": 1.14,
        "E_fd_ibb": 1.0,
        "P_m_gen":  0.9,

        "omega_gen_init": 0,
        "delta_gen_init": np.deg2rad(48.59),
        "delta_ibb_init": np.deg2rad(0)
    }

    sim_parameters = {
        "t_start":      -1,
        "t_end":        5,
        "t_step":       0.001,
        
        "fault_start":  0,
        "fault_end":    5,
        "clearing":     True
    }

    gen_parameters["X_fault"] = [(-1j / gen_parameters["X_gen"] - 1j / gen_parameters["X_line"]) + 1000000, 1j / gen_parameters["X_line"]]

    # Execution of simulation
    alg = True
    stability, t_cc, delta_cc, t_sim, solution_stable, solution_unstable = do_sim(gen_parameters, sim_parameters, alg)

    # Evaluation of results
    print('t_cc:\t\t' + str(round(t_cc, 3)) + ' s')
    print('delta_cc:\t' + str(round(np.rad2deg(delta_cc), 1)) + ' deg')

    delta_stable = solution_stable[:,1]
    omega_stable = solution_stable[:,0]

    ##############################
    # Plot stable result
    ##############################
    fig, axs = plt.subplots(2, 1, figsize=(6,8), sharex=True)

    # determine the boundary angles
    delta_0 = delta_gen_init # delta_gen_init
    delta_max = np.pi - delta_0

    # calculation of P_e_pre, P_e_post, and P_t
    x_deg = np.linspace(0, 180) # linear vector for plotting in deg
    x_rad = np.linspace(0, np.pi) # linear vector for calculation in rad
    P_e_pre = P_e_alg(x_rad, False)
    P_e_post = P_e_alg(x_rad, True)
    P_t = P_t_deg(x_rad)

    plt.subplots_adjust(hspace=.0)

    ##############################
    # ax1
    ##############################
    axs[0].plot(x_deg, P_e_pre, '-', linewidth=2, label='$P_\mathrm{e}$ pre-fault')
    # axs[0].plot(x_deg, P_e_post, '-', linewidth=2, label='$P_\mathrm{e}$ post-fault')
    axs[0].plot(x_deg, P_t, '-', linewidth=2, label='$P_\mathrm{T}$ of the turbine')
    axs[0].set_ylim(bottom=0)
    delta_0_deg = np.rad2deg(delta_0)
    delta_max_deg = np.rad2deg(delta_max)
    delta_c_deg = np.rad2deg(delta_cc)
    # axs[0].set_xticks([0, 180, delta_0_deg, delta_c_deg, delta_max_deg], labels=['0', '180', '$\delta_\mathrm{0}$', '$\delta_\mathrm{c}$', '$\delta_\mathrm{max}$'])

    ix1 = np.linspace(delta_0_deg, delta_c_deg)
    iy1 = P_e_alg(np.deg2rad(ix1), True)
    axs[0].fill_between(ix1, iy1, P_m_gen, facecolor='0.9', edgecolor='0.5')

    # Make the shaded region for area_dec, https://matplotlib.org/stable/gallery/lines_bars_and_markers/fill_between_demo.html
    ix2 = np.linspace(delta_c_deg, delta_max_deg) # -> does this have to be in rad or in deg?
    iy2 = P_e_alg(np.deg2rad(ix2), False)
    axs[0].fill_between(ix2, iy2, P_m_gen, facecolor='0.9', edgecolor='0.5')
    axs[0].grid()
    axs[0].legend()
    axs[0].set_ylabel('power in p.u.')

    ##############################
    # ax2
    ##############################
    axs[1].plot(np.rad2deg(delta_stable), t_sim, label='delta')
    fig.gca().invert_yaxis()
    # axs[1].axhline(y=fault_end, linestyle='--', label='clearing of fault')
    axs[1].grid()
    axs[1].set_ylabel('time in s')
    axs[1].legend()
    plt.ylim(top=-.5)
    plt.xlim(left=0, right=180)
    plt.xlabel('power angle $\delta$ in deg')

    plt.suptitle('Stable scenario')
    plt.show()