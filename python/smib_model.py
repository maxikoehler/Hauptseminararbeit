import matplotlib.pyplot as plt
import numpy as np

def mag_and_angle_to_cmplx(mag, angle):
    return mag * np.exp(1j * angle)

# Define the parameters of the system and starting values
fn = 60
H_gen = 3.5
X_gen = 0.2
X_ibb = 0.1
X_line = 0.65

E_fd_gen = 1.075
E_fd_ibb = 1.033
P_m_gen = 1998/2200

# init states of variables
omega_gen_init = 0 # init state
delta_gen_init = np.deg2rad(45.9) # init state
delta_ibb_init = np.deg2rad(-5.0) # init state
v_bb_gen_init = mag_and_angle_to_cmplx(1.0, np.deg2rad(36.172))

def init(t, fault, system_parameters, init_state):
    # -> wie global setzen?

    t_start, t_stop, timestep = t
    fault_type, fault_start = fault # fault_start, fault_end also neede for CCT-determination?
    fn, H_gen, X_gen, X_line, X_ibb, E_fd_gen, E_fd_ibb, P_m_gen = system_parameters
    omega_gen_init, delta_gen_init, delta_ibb_init, v_bb_gen_init = init_state
    return

def dynamic_system(init_state, t, fault_on):
    # global H_gen
    omega, delta = init_state
    
    d_omega_dt = 1 / (2 * H_gen) * (T_m_gen(omega) - P_e_gen(delta, fault_on))
    d_delta_dt = omega

    return [d_omega_dt, d_delta_dt]

def T_m_gen(omega):
    # Assuming a simple linear function for demonstration purposes
    return P_m_gen / (1 + omega)

def P_e_gen(delta, fault_on):
    if fault_on:
        y_adm = np.array([[(-1j / X_gen - 1j / X_line) + 1000000, 1j / X_line],[1j / X_line, -1j / X_line - 1j / X_ibb]])
    else:
        y_adm = np.array([[-1j / X_gen - 1j / X_line, 1j / X_line],[1j / X_line, -1j / X_line - 1j / X_ibb]])

    y_inv = np.linalg.inv(y_adm)

    i_inj_gen = mag_and_angle_to_cmplx(E_fd_gen, delta) / (1j * X_gen)
    i_inj_ibb = mag_and_angle_to_cmplx(E_fd_ibb, delta_ibb_init) / (1j * X_ibb)

    # Calculate voltages at the bus by multiplying the inverse of the admittance matrix with the current injections
    v_bb_gen = y_inv[0, 0] * i_inj_gen + y_inv[0, 1] * i_inj_ibb
    v_bb_ibb = y_inv[1, 0] * i_inj_gen + y_inv[1, 1] * i_inj_ibb

    E_gen_cmplx = mag_and_angle_to_cmplx(E_fd_gen, delta)
    
    P_e = (v_bb_gen * np.conj((E_gen_cmplx - v_bb_gen) / (1j * X_gen))).real
    return P_e

def stability_eac():
    # under construction!!! Aufruf zu jedem Zeitpunkt...

    # Berechne die aktuelle Beschleunigungsfläche -> integrate(P_t - P_e_fault)_delta-0^delta-akt
    # Berechne aktuell übrige Beschleunigungsfläche -> integrate(P_e_normal - P_t)_delta-akt^delta-max

    # vergleichen -> aktuell noch stabiler Betrieb oder nicht

    area_acc = 1
    area_dec = 1

    if area_acc < area_dec: # True: stable, False: NOT stable 
        return True
    else:
        return False

def solver_euler(t, fault_start):
    # Initialize the variables
    omega_gen = omega_gen_init
    delta_gen = delta_gen_init
    v_bb_gen = v_bb_gen_init

    # Define time. Here, the time step is 0.005 s and the simulation is 5 s long
    t = np.arange(0, 5, 0.005)
    x_result = []

    for timestep in t:

        # Those lines cause a short circuit at t = 1 s until t = 1.05 s
        if 1 <= timestep < 1.05:
            sc_on = True
        else:
            sc_on = False

        # Calculate the differences to the next step by executing the differential equations at the current step
        domega_dt, ddelta_dt = differential(omega_gen, v_bb_gen, delta_gen)
        omega_gen = omega_gen + domega_dt * (t[1] - t[0])
        delta_gen = delta_gen + ddelta_dt * (t[1] - t[0])

        v_bb_gen = algebraic(delta_gen, sc_on)

        # Save the results, so they can be plotted later
        x_result.append(omega_gen)

    # Convert the results to a numpy arrays
    res = np.vstack(x_result)

    # return variables
    t_sim = 1 # time vector
    omega = 1 # TDS: omega
    delta = 1 # TDS: delta
    P_e = 1 # TDS: P_e, terminal voltage at generator
    P_m = [] # TDS: P_m, Mechanical power of generator; -> const?; Moment(Torque) ist nicht const, da abhängig von omega
    t_cc = 1 # max. time under fault for stable clearing
    delta_cc = 1 # max. angle delta when clearing of fault results in stable operation and no cut-off or runaway

    return t_sim, omega, delta, P_e, P_m, t_cc, delta_cc

if __name__ == "__main__":
    print("Hello World!")