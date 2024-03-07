############################
# simulation of fault 1
# t_cc, delta_cc and some plots
#
# scerario: full line fault (P_e = 0), clearing mode (around t_cc)
############################

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib as mpl
import numpy as np
import scipy as sp
from scipy.integrate import odeint

import smib_model as sim

# redefining plot save parameters
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Charter"],
    "font.size": 12
})

# uncomment for updating savefig options for latex export
# mpl.use("pgf")

def init(gen_parameters, sim_parameters):
    global fn, H_gen, X_gen, X_ibb, X_line, X_fault, E_fd_gen, E_fd_ibb, P_m_gen, omega_gen_init, delta_gen_init, delta_ibb_init, t_start, t_end, t_step, fault_start, fault_end, clearing

    fn = gen_parameters["fn"]
    H_gen = gen_parameters["H_gen"]
    X_gen = gen_parameters["X_gen"]
    X_ibb = gen_parameters["X_ibb"]
    X_line = gen_parameters["X_line"]
    X_fault = gen_parameters["X_fault"]

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

    return

if __name__ == "__main__":
    # setup simulation inputs
    gen_parameters = {
        "fn":       60,
        "H_gen":    3.5,
        "X_gen":    0.2,
        "X_ibb":    0.1,
        "X_line":   0.65,
        "X_fault":  0.0001,

        "E_fd_gen": 1.075,
        "E_fd_ibb": 1.033,
        # "P_m_gen":  0.3,
        "P_m_gen":  1998/2200,

        "omega_gen_init": 0,
        "delta_gen_init": np.deg2rad(50.9),
        "delta_ibb_init": np.deg2rad(0)
    }

    sim_parameters = {
        "t_start":      -1,
        "t_end":        1,
        "t_step":       0.001,
        
        "fault_start":  0,
        "fault_end":    1,
        "clearing":     True
    }

    gen_parameters["X_fault"] = [(-1j / gen_parameters["X_gen"] - 1j / gen_parameters["X_line"]) + 1000000, 1j / gen_parameters["X_line"]]

    init(gen_parameters, sim_parameters)

    # Execution of simulation
    alg = True
    stability, t_cc, delta_cc, t_sim, solution_stable, solution_unstable = sim.do_sim(gen_parameters, sim_parameters, alg)

    # Evaluation of results
    print('t_cc:\t\t' + str(round(t_cc, 3)) + ' s')
    print('delta_cc:\t' + str(round(np.rad2deg(delta_cc), 1)) + ' deg')

    delta_stable = solution_stable[:,1]
    delta_unstable = solution_unstable[:,1]

    ##############################
    # Plot unstable result
    ##############################
    fig, axs = plt.subplots(2, 1, figsize=(6,8), sharex=True)

    # determine the boundary angles
    delta_0 = delta_gen_init # delta_gen_init
    delta_max = np.pi - delta_0

    # calculation of P_e_pre, P_e_post, and P_t
    x_deg = np.linspace(0, 180) # linear vector for plotting in deg
    x_rad = np.linspace(0, np.pi) # linear vector for calculation in rad
    P_e_pre = sim.P_e_alg(x_rad, False)
    P_e_post = sim.P_e_alg(x_rad, True)
    P_t = sim.P_t_deg(x_rad)

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
    iy1 = sim.P_e_alg(np.deg2rad(ix1), True)
    axs[0].fill_between(ix1, iy1, P_m_gen, facecolor='0.9', edgecolor='0.5')

    # Make the shaded region for area_dec, https://matplotlib.org/stable/gallery/lines_bars_and_markers/fill_between_demo.html
    ix2 = np.linspace(delta_c_deg, delta_max_deg) # -> does this have to be in rad or in deg?
    iy2 = sim.P_e_alg(np.deg2rad(ix2), False)
    axs[0].fill_between(ix2, iy2, P_m_gen, facecolor='0.9', edgecolor='0.5')
    axs[0].grid()
    axs[0].legend()
    axs[0].set_ylabel('power in pu')

    ##############################
    # ax2
    ##############################
    axs[1].plot(np.rad2deg(delta_stable), t_sim, label='delta')
    fig.gca().invert_yaxis()
    axs[1].axhline(y=t_cc, linestyle='--', label='clearing of fault')
    axs[1].grid()
    axs[1].set_ylabel('time in s')
    axs[1].legend()
    plt.ylim(top=-.1)
    plt.xlim(left=0, right=180)
    plt.xlabel('power angle $\delta$ in deg')

    plt.suptitle('Stable scenario - fault 1')
    # plt.savefig('plots/fault1_stable.pgf')
    plt.show()
    # plt.close()

    ##############################
    # Plot UNstable result
    ##############################
    fig, axs = plt.subplots(2, 1, figsize=(6,8), sharex=True)

    # determine the boundary angles
    delta_0 = delta_gen_init # delta_gen_init
    delta_max = np.pi - delta_0

    # calculation of P_e_pre, P_e_post, and P_t
    x_deg = np.linspace(0, 180) # linear vector for plotting in deg
    x_rad = np.linspace(0, np.pi) # linear vector for calculation in rad
    P_e_pre = sim.P_e_alg(x_rad, False)
    P_e_post = sim.P_e_alg(x_rad, True)
    P_t = sim.P_t_deg(x_rad)

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
    iy1 = sim.P_e_alg(np.deg2rad(ix1), True)
    axs[0].fill_between(ix1, iy1, P_m_gen, facecolor='0.9', edgecolor='0.5')

    # Make the shaded region for area_dec, https://matplotlib.org/stable/gallery/lines_bars_and_markers/fill_between_demo.html
    ix2 = np.linspace(delta_c_deg, delta_max_deg) # -> does this have to be in rad or in deg?
    iy2 = sim.P_e_alg(np.deg2rad(ix2), False)
    axs[0].fill_between(ix2, iy2, P_m_gen, facecolor='0.9', edgecolor='0.5')
    axs[0].grid()
    axs[0].legend()
    axs[0].set_ylabel('power in pu')

    ##############################
    # ax2
    ##############################
    axs[1].plot(np.rad2deg(delta_unstable), t_sim, label='delta')
    fig.gca().invert_yaxis()
    axs[1].axhline(y=t_cc, linestyle='--', label='clearing of fault')
    axs[1].grid()
    axs[1].set_ylabel('time in s')
    axs[1].legend()
    plt.ylim(top=-.1)
    plt.xlim(left=0, right=180)
    plt.xlabel('power angle $\delta$ in deg')

    plt.suptitle('Unstable scenario - fault 1')
    # plt.savefig('plots/fault1_unstable.pgf')
    plt.show()
    # plt.close()