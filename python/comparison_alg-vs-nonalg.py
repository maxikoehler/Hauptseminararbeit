############################
# Comparison of calculation with P_e with algebraic equations vs. P_e = E*E/X *sin delta
#
# OR: comparison of solving the determine_cct with or without the TDS -> wrong speed/delta correlations and thus wrong Momentum (-> difference of T_m and P_e). Leading to false t_cc...
############################

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib as mpl
import numpy as np
import scipy as sp
from scipy.integrate import odeint

import smib_model as sim

# redefining plot displaying parameters
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Charter"],
    "font.size": 12
})

# uncomment for updating savefig options for latex export
# mpl.use("pgf")

def init(gen_parameters, sim_parameters):
    global fn, H_gen, X_gen, X_ibb, X_line, X_fault, E_fd_gen, E_fd_ibb, P_m_gen, omega_gen_init, delta_gen_init, delta_ibb_init, t_start, t_end, t_step, fault_start, fault_end

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
        "P_m_gen":  1998/2200,

        "omega_gen_init": 0,
        "delta_gen_init": np.deg2rad(50.9),
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

    init(gen_parameters, sim_parameters)

    # Execution of simulation
    alg = True
    stability_alg, t_cc_alg, delta_cc_alg, t_sim, solution_alg_stable, solution_alg_unstable = sim.do_sim(gen_parameters, sim_parameters, alg)

    print('t_cc:\t\t' + str(round(t_cc_alg, 3)) + ' s')
    print('delta_cc:\t' + str(round(np.rad2deg(delta_cc_alg), 1)) + ' deg')

    alg = False
    gen_parameters["delta_gen_init"] = sim.get_delta_0(gen_parameters, sim_parameters, alg)

    stability, t_cc, delta_cc, t_sim, solution_stable, solution_unstable = sim.do_sim(gen_parameters, sim_parameters, alg)

    print('t_cc:\t\t' + str(round(t_cc, 3)) + ' s')
    print('delta_cc:\t' + str(round(np.rad2deg(delta_cc), 1)) + ' deg')

    ##############################
    # TDS Plot of the different solutions
    ##############################
    delta_alg = solution_alg_stable[:,1]
    delta_non_alg = solution_stable[:,1]

    plt.plot(t_sim, np.rad2deg(delta_alg), label="algebraic")
    plt.plot(t_sim, np.rad2deg(delta_non_alg), label="non-algebraic")

    plt.legend()
    plt.grid()
    plt.show()
    # plt.savefig("plots/comparison_alg-vs-nonalg.pgf")
    # plt.close()

    ##############################
    # Different solutions: Plot of the different P_e's
    ##############################
    P_e = np.zeros(int((t_end-t_start)/t_step))
    for i in range(np.size(t_sim)-1):
        if fault_start <= t_sim[i] < t_cc:
            fault = True
        else:
            fault = False
        P_e[i] = sim.P_e(delta_non_alg[i], fault)

    P_e_alg = np.zeros(int((t_end-t_start)/t_step))
    for i in range(np.size(t_sim)-1):
        if fault_start <= t_sim[i] < t_cc_alg:
            fault = True
        else:
            fault = False
        P_e_alg[i] = sim.P_e_alg(delta_non_alg[i], fault)

    plt.plot(t_sim, P_e, label="power non-algebraic")
    plt.plot(t_sim, P_e_alg, label="power algebraic")
    plt.legend()
    plt.grid()
    plt.show()
    # plt.savefig("plots/comparison_alg-vs-nonalg.pgf")
    # plt.close()