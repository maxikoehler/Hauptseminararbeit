def stability_eac(delta_0, delta_act, omega_act, delta_max, alg):
    area_acc = sp.integrate.quad(P_r_deg, delta_0, delta_act, args=(omega_act, True, alg))
    area_dec = sp.integrate.quad(P_r_deg, delta_act, delta_max, args=(omega_act, (not clearing), alg))

    if abs(area_acc[0]) < abs(area_dec[0]):
        return True
    else:
        return False

def determine_cct(t_sim, delta, omega, delta_0, alg):
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