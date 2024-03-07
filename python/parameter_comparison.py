############################
# Comparison of different machine parameters and operational points.
############################

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib as mpl
import numpy as np
import scipy as sp
from scipy.integrate import odeint
# import tikzplotlib as tp

import smib_model as sim

# redefining plot save parameterss
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Charter"],
    "font.size": 12
})

# # update savefig options for latex export
# mpl.use("pgf")

# comparing different H_gen in t_cc, delta_cc and their TDS when a fault is ON


# comparing different degree of utilization (P_e / P_e_max OR delta_P; direct correlated: delta_0)


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
        "X_fault":  0.001,

        "E_fd_gen": 1.075,
        "E_fd_ibb": 1.033,
        "P_m_gen":  1998/2200,

        "omega_gen_init": 0,
        "delta_gen_init": np.deg2rad(50.9),
        "delta_ibb_init": np.deg2rad(0)
    }

    sim_parameters = {
        "t_start":      0,
        "t_end":        5,
        "t_step":       0.001,
        
        "fault_start":  0,
        "fault_end":    5,
        "clearing":     True
    }

    gen_parameters["X_fault"] = [(-1j / gen_parameters["X_gen"] - 1j / gen_parameters["X_line"]) + 1000000, 1j / gen_parameters["X_line"]]

    # init(gen_parameters, sim_parameters)

    # Execution of simulations and savon in t_ccs array
    alg = True

    sim.init(gen_parameters, sim_parameters)
    P_e_max = sim.P_e_alg(np.pi/2, False)

    H = np.arange(0.1, 12.5, 0.5) # [2.5, 3, 3.5, 4, 4.5, 5]
    delta_P = [0.3, 0.5, 0.7, 0.9, 1.1]
    t_ccs = np.zeros((np.size(delta_P), np.size(H)))
    i = 0
    for delta_P_new in delta_P:
        gen_parameters["P_m_gen"] = delta_P_new
        gen_parameters["delta_gen_init"] = sim.get_delta_0(gen_parameters, sim_parameters, alg)
        z = 0
        for H_new in H:
            gen_parameters["H_gen"] = H_new

            stability, t_cc, delta_cc, t_sim, solution = sim.do_sim_simple(gen_parameters, sim_parameters, alg)
            print(t_cc)

            if stability:
                t_ccs[i, z] = t_cc
            else:
                t_ccs[i, z] = -1
            z = z + 1

        plt.plot(H, t_ccs[i,:], label=("$\Delta P =$ " + str(round(delta_P_new/P_e_max*100, 1)) + " $\%$"))
        i = i + 1

    plt.legend()
    plt.xlabel("$H_\mathrm{gen}$ in $\mathrm{s}$")
    plt.ylabel("$CCT$ in $\mathrm{s}$")
    # plt.title("Influence from $H_\mathrm{gen}$ and $\Delta P$ on CCT")
    plt.grid()
    # figure = plt.gcf() # get current figure
    # figure.set_size_inches(8, 6)
    # plt.savefig('plots/parameter_comparison.pgf', dpi=300)
    plt.show()
    # plt.close()