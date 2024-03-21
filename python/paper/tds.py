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

def P_e_alg(delta, fault_on):
    # function for determing P_e WITH algebraic help
    global E_fd_gen
    global X_gen
    
    v_bb_gen = algebraic(delta, fault_on)

    E_gen_complex = mag_and_angle_to_cmplx(E_fd_gen, delta)
    P_e_gen = (v_bb_gen * np.conj((E_gen_complex - v_bb_gen) / (1j * X_gen))).real
    return P_e_gen

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
        X = X_gen + X_line
        E_ibb = E_fd_ibb
        
    P_e_gen = E_fd_gen * E_ibb / X * np.sin(delta)
    return P_e_gen

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