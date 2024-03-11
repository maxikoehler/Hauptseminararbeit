import numpy as np
import scipy as sp

################################
# some pre definitions
################################
gen_parameters = {
    "fn":       60,
    "H_gen":    3.5,
    "X_gen":    0.2,
    "X_ibb":    0.1,
    "X_line":   0.65,
    "X_fault":  0.0001,

    "E_gen": 1.075,
    "E_ibb": 1.033,
    # "P_m_gen":  0.3,
    "P_m_gen":  1998/2200,
    "P_e_max": 1.2,

    "omega_gen_init": 0,
    "delta_gen_init": np.deg2rad(50.9),
    "delta_ibb_init": np.deg2rad(0)
}

def P_e(delta):
    global P_e_max
    return P_e_max * np.sin(delta)

def init(gen_parameters):
    global fn, H_gen, X_gen, X_ibb, X_line, X_fault, E_gen, E_ibb, P_m_gen, omega_gen_init, delta_gen_init, delta_ibb_init, omega_0, delta_0, delta_max, P_e_max

    fn = gen_parameters["fn"]
    H_gen = gen_parameters["H_gen"]
    X_gen = gen_parameters["X_gen"]
    X_ibb = gen_parameters["X_ibb"]
    X_line = gen_parameters["X_line"]
    X_fault = gen_parameters["X_fault"]

    E_gen = gen_parameters["E_gen"]
    E_ibb = gen_parameters["E_ibb"]
    P_m_gen = gen_parameters["P_m_gen"]
    P_e_max = gen_parameters["P_e_max"]

    omega_gen_init = gen_parameters["omega_gen_init"]
    delta_gen_init = gen_parameters["delta_gen_init"]
    delta_ibb_init = gen_parameters["delta_ibb_init"]

    delta_0 = delta_gen_init - delta_ibb_init
    omega_0 = 2 * np.pi * fn
    delta_max = np.pi - delta_0

    return

################################
# calculation
################################
def calc():
    global H_gen, omega_0
    delta_cc = np.arccos(np.sin(delta_0) * (np.pi - 2 * delta_0) - np.cos(delta_0))
    t_cc = np.sqrt((4 * H_gen * (delta_cc - delta_0))/(omega_0 * abs(P_m_gen - P_e(0))))
    return delta_cc, t_cc

################################
# printing the resuts
################################
init(gen_parameters)

delta_cc, t_cc = calc()
print('t_cc:\t\t' + str(round(t_cc, 3)) + ' s')
print('delta_cc:\t' + str(round(np.rad2deg(delta_cc), 2)) + ' deg')