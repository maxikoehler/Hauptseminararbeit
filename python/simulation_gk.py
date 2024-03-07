import matplotlib.pyplot as plt
import numpy as np
import scipy as sp


def mag_and_angle_to_cmplx(mag, angle):
    return mag * np.exp(1j * angle)


fn = 60

H_gen = 3.5
X_gen = 0.2
X_ibb = 0.1
X_line = 0.65

# Values are initialized from loadflow
E_fd_gen = 1.075
E_fd_ibb = 1.033
P_m_gen = 1998/2200

omega_gen_init = 0
delta_gen_init = np.deg2rad(45.9)
delta_ibb_init = np.deg2rad(-5.0)

v_bb_gen_init = mag_and_angle_to_cmplx(1.0, np.deg2rad(36.172))


def differential(omega, v_bb_gen, delta):
    # Calculate the electrical power extracted from the generator at its busbar.
    E_gen_cmplx = mag_and_angle_to_cmplx(E_fd_gen, delta)
    P_e_gen = (v_bb_gen * np.conj((E_gen_cmplx - v_bb_gen) / (1j * X_gen))).real

    # transform the constant mechanical energy into torque
    T_m_gen = P_m_gen / (1 + omega)

    # Differential equations of a generator according to Machowski
    domega_dt = 1 / (2 * H_gen) * (T_m_gen - P_e_gen)
    ddelta_dt = omega * 2 * np.pi * fn

    return domega_dt, ddelta_dt, P_e_gen, T_m_gen


def algebraic(delta_gen, sc_on):
    # If the SC is on, the admittance matrix is different.
    # The SC on busbar 0 is expressed in the admittance matrix as a very large admittance (1000000) i.e. a very small impedance.
    if sc_on:
        y_adm = np.array([[(-1j / X_gen - 1j / X_line) + 1000000, 1j / X_line],
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


def do_sim(fault_start, fault_end, sim_end):
    # Initialize the variables
    omega_gen = omega_gen_init
    delta_gen = delta_gen_init
    v_bb_gen = v_bb_gen_init

    # Define time. Here, the time step is 0.005 s and the simulation is 5 s long
    t = np.arange(0, sim_end, 0.001)
    res_omega = []
    res_delta = []
    res_P_e = []
    res_P_m = []

    for timestep in t:

        # Those lines cause a short circuit at t = 1 s until t = 1.05 s
        if fault_start <= timestep < fault_end:
            sc_on = True
        else:
            sc_on = False

        # Calculate the initial guess for the next step by executing the differential equations at the current step
        domega_dt_guess, ddelta_dt_guess, P_e_gen, P_m_gen = differential(omega_gen, v_bb_gen, delta_gen)
        omega_guess = omega_gen + domega_dt_guess * (t[1] - t[0])
        delta_guess = delta_gen + ddelta_dt_guess * (t[1] - t[0])

        v_bb_gen = algebraic(delta_guess, sc_on)

        # Calculate the differential equations with the initial guess
        domega_dt_guess2, ddelta_dt_guess2, P_e_gen2, P_m_gen2 = differential(omega_guess, v_bb_gen, delta_guess)

        domega_dt = (domega_dt_guess + domega_dt_guess2) / 2
        ddelta_dt = (ddelta_dt_guess + ddelta_dt_guess2) / 2

        omega_gen = omega_gen + domega_dt * (t[1] - t[0])
        delta_gen = delta_gen + ddelta_dt * (t[1] - t[0])

        v_bb_gen = algebraic(delta_gen, sc_on)


        # Save the results, so they can be plotted later
        res_omega.append(omega_gen)
        res_delta.append(delta_gen)
        res_P_e.append(P_e_gen)
        res_P_m.append(P_m_gen)

    # Convert the results to a numpy array
    res_omega = np.vstack(res_omega)
    res_delta = np.vstack(res_delta)
    P_e = np.vstack(res_P_e)
    P_m = np.vstack(res_P_m)
    return t, res_omega, res_delta, P_e, P_m

def stability_eac(delta_0, delta_act, delta_max):
    # Compare the acceleration area until the given delta and compare it to the braking area left until the dynamic stability point is passed
    area_acc = sp.integrate.quad(P_t_deg, delta_0, delta_act)
    area_dec = sp.integrate.quad(P_r_deg, delta_act, delta_max)

    if area_acc < area_dec: # True: stable, False: NOT stable 
        return True
    else:
        return False

def determine_max_stab(delta_0, delta, delta_max):
    # Save current time and delta at time point i; iterate through i to test any given time until stability can't be remained; delta_cc and t_cc is the angle and time at the last stable point
    i = 0

    t_cc = -1
    delta_cc = -1

    while stability_eac(delta_0, delta[i], delta_max):
        t_cc = t_sim[i]
        delta_cc = delta[i]
        i = i + 1
    
    return t_cc, delta_cc

if __name__ == '__main__':

    # Here the simulation is executed and the timesteps and corresponding results are returned.
    # In this example, the results are omega, delta, e_q_t, e_d_t, e_q_st, e_d_st of the generator and the IBB
    fault_start = 1
    fault_end = 1.05
    sim_end = 5

    t_sim, res, res_delta, P_e, P_m = do_sim(fault_start, fault_end, sim_end)

    # load the results from powerfactory for comparison
    delta_omega_pf = np.loadtxt('simulation_example_gk/pictures/powerfactory_data.csv', skiprows=1, delimiter=',')

    # Plot the results
    plt.plot(t_sim, np.rad2deg(res_delta[:, 0].real), label='delta_omega_gen_python')
    # plt.plot(delta_omega_pf[:, 0], delta_omega_pf[:, 1] - 1, label='delta_omega_gen_powerfactory')
    plt.legend()
    plt.grid()
    plt.title('Reaction of a generator to a short circuit')
    # plt.savefig('simulation_example_gk/pictures/short_circuit_improved.png')
    plt.show()

    fig, axs = plt.subplots(2)
    axs[0].plot(t_sim, P_e, label='electrical power generator')
    axs[0].plot(t_sim, P_m, label='mechanical power generator')
    axs[1].plot(t_sim, P_m-P_e, label='power difference')
    plt.grid()
    plt.show()