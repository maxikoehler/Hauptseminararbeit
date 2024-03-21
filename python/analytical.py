import numpy as np
import scipy as sp

################################
# some pre definitions
################################
gen_parameters1 = {
    "fn":       50,
    "H_gen":    3.3,
    "X_gen":    0.2,
    "X_ibb":    0.1,
    "X_line":   0.65,
    "X_fault":  0.0001,

    "E_gen": 1.14,
    "E_ibb": 1,
    "P_m_gen":  0.9,
    "P_e_max": 1.2,

    "omega_gen_init": 0,
    "delta_gen_init": np.deg2rad(48.6),
    "delta_ibb_init": np.deg2rad(0)
}

gen_parameters2 = {
    "fn":       50,
    "H_gen":    3.3,
    "X_gen":    0.2,
    "X_ibb":    0.1,
    "X_line":   0.65,
    "X_fault":  0.0001,

    "E_gen": 1.14,
    "E_ibb": 1,
    "P_m_gen":  0.9,
    "P_e_max": 1.2,

    "omega_gen_init": 0,
    "delta_gen_init": np.deg2rad(48.6),
    "delta_ibb_init": np.deg2rad(0)
}

def P_e(P_e_max_f, delta):
    return P_e_max_f * np.sin(delta)

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
def calc(P_e_max_f):
    global H_gen, omega_0, P_m_gen, P_e_max
    delta_cc = np.arccos((P_m_gen/(P_e_max - P_e_max_f)) * (np.pi - 2 * delta_0) - (P_e_max_f/(P_e_max - P_e_max_f)) * np.cos(delta_0) - (P_e_max/(P_e_max - P_e_max_f)) * np.cos(delta_0))
    t_cc = np.sqrt((4 * H_gen * (delta_cc - delta_0))/(omega_0 * abs(P_m_gen - P_e(P_e_max_f, delta_cc))))
    return delta_cc, t_cc

################################
# printing the resuts
################################

init(gen_parameters1)
delta_cc, t_cc = calc(0)
print('fault 1:')
print('t_cc:\t\t' + str(round(t_cc, 3)) + ' s')
print('delta_cc:\t' + str(round(np.rad2deg(delta_cc), 2)) + ' deg\n')

init(gen_parameters2)
delta_cc, t_cc = calc(0.808)
print('fault 2:')
print('t_cc:\t\t' + str(round(t_cc, 3)) + ' s')
print('delta_cc:\t' + str(round(np.rad2deg(delta_cc), 2)) + ' deg')