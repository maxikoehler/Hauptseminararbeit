import numpy as np
import scipy as sp

################################
# some pre definitions
################################

H_gen = 3.5
omega_0 = 2 * np.pi * 60
P_m = 1998/2200
P_e_max = 2200
X_line = 0.65
X_gen = 0.2
X_trans = 0.1
E_gen = 1.075
E_ibb = 1.033
delta_0 = np.deg2rad(50.9)
delta_max = np.pi - delta_0

def P_e(delta):
    global P_e_max
    return P_e_max * np.sin(delta)

################################
# calculation
################################
def calc():
    global H_gen, omega_0
    delta_cc = np.arccos(np.sin(delta_0) * (np.pi - 2 * delta_0) - np.cos(delta_0))
    t_cc = np.sqrt((4 * H_gen * (delta_cc - delta_0))/(omega_0 * abs(P_m - P_e(0))))
    return delta_cc, t_cc

################################
# printing the resuts
################################

delta_cc, t_cc = calc()
print('t_cc:\t\t' + str(round(t_cc, 3)) + ' s')
print('delta_cc:\t' + str(round(np.rad2deg(delta_cc), 2)) + ' deg')