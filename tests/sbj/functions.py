from constants import * 
import numpy as np
from pymdo.core.variable import Variable


# Design variables #
def create_design_vars():

    design_vars = []
    design_vars_dict = {}
    design_vars_dict[z_cr_str] = [13e3, 18e3]
    design_vars_dict[M_cr_str] = [1.6, 2]
    design_vars_dict[S_w_str] = [100, 200]
    design_vars_dict[fwLE_w_str] = [np.deg2rad(45), np.deg2rad(70)]
    design_vars_dict[fwTE_w_str] = [np.deg2rad(-5), np.deg2rad(15)]
    design_vars_dict[tr_w_str] = [0.05, 0.5]
    design_vars_dict[tcr_w_str] = [0.04, 0.06]
    design_vars_dict[fwLE_v_str] = [np.deg2rad(45), np.deg2rad(70)]
    design_vars_dict[fwTE_v_str] = [np.deg2rad(-5), np.deg2rad(15)]
    design_vars_dict[tr_v_str] = [0.05, 0.5]
    design_vars_dict[tcr_v_str] = [0.06, 0.08]
    design_vars_dict[W_fuel_str] = [15e3, 30e3]

    for name in design_vars_dict:

        design_vars.append(Variable(name, 1, design_vars_dict[name][0], design_vars_dict[name][1]))

        value = (design_vars_dict[name][0] + design_vars_dict[name][1]) / 2.0

        design_vars_dict[name] = np.array([value])


    return design_vars_dict, design_vars

# Flight conditions #
def compute_temperature(z):
    if z <= 11e3:
        return T_ref+dTdz*z
    elif z <= 25e3:
        return 216.65


def compute_pressure(z, T):
    if z <= z11:
        return p_ref * (T / T_ref)**(-g/R/dTdz)
    else:
        return p11 * np.exp(-(g/R/T)*(z-z11))


def compute_density(z, T):
    if z <= z11:
        return rho_ref * (T/T_ref)**(-g/R/dTdz - 1)
    else:
        p = compute_pressure(z, T)
        return p/R/T


def compute_velocity(Mach, T):
    return Mach * np.sqrt(R * gamma * T)


def compute_cruising_conditions(M_cr, z_cr):
    T_cr = compute_temperature(z_cr)
    p_cr = compute_pressure(z_cr, T_cr)
    rho_cr = compute_density(z_cr,  T_cr)
    V_cr = compute_velocity(M_cr, T_cr)
    return T_cr, p_cr, rho_cr, V_cr


def ExecuteFlightConditions(input_values):
    M_cr = input_values[M_cr_str]
    z_cr = input_values[z_cr_str]
    T_cr, p_cr, rho_cr, V_cr = compute_cruising_conditions(
        M_cr[0],  z_cr[0])
    return {T_cr_str: np.array([T_cr]), p_cr_str: np.array([p_cr]), rho_cr_str: np.array([rho_cr]), V_cr_str: np.array([V_cr])}


# Wing geometry #
def compute_wing_chord_length(S_w, fwLE_w, fwTE_w, tr_w):
    return np.sqrt(S_w*(np.tan(fwLE_w)-np.tan(fwTE_w))/(1-tr_w**2))


def compute_wing_span(S_w, c_w, tr_w):
    return 2*S_w/(c_w*(1+tr_w))


def compute_wing_aspect_ratio(S_w, b_w):
    return b_w**2 / S_w


def compute_wing_sweep_angle_25(fwTE_w, AR_w, tr_w):
    # plus or minus?
    val = np.tan(fwTE_w) + 0.75 * 4/AR_w * (1-tr_w)/(1+tr_w)
    return np.arctan(val)


def compute_wing_exposed_surface(S_w, c_w, fwLE_w, fwTE_w):
    return S_w - D_f/2 * (2 * c_w - D_f/2*np.tan(fwLE_w) + D_f/2*np.tan(fwTE_w))


def compute_wing_wetted_surface(tcr_w, Sexp_w):
    return 2*(1 + 0.2 * tcr_w) * Sexp_w


def compute_wing_frontal_area():
    return 1


def compute_wing_geometry(S_w, fwLE_w, fwTE_w, tr_w, tcr_w):
    c_w = compute_wing_chord_length(S_w, fwLE_w, fwTE_w, tr_w)
    b_w = compute_wing_span(S_w, c_w, tr_w)
    AR_w = compute_wing_aspect_ratio(S_w, b_w)
    fw25_w = compute_wing_sweep_angle_25(fwTE_w, AR_w, tr_w)
    Sexp_w = compute_wing_exposed_surface(S_w, c_w, fwLE_w, fwTE_w)
    Swet_w = compute_wing_wetted_surface(tcr_w, Sexp_w)
    return c_w, b_w, AR_w, fw25_w, Sexp_w, Swet_w


def ExecuteWingGeometry(input_values):
    S_w = input_values[S_w_str]
    fwLE_w = input_values[fwLE_w_str]
    fwTE_w = input_values[fwTE_w_str]
    tr_w = input_values[tr_w_str]
    tcr_w = input_values[tcr_w_str]
    c_w, b_w, AR_w, fw25_w, Sexp_w, Swet_w = compute_wing_geometry(
        S_w[0], fwLE_w[0], fwTE_w[0], tr_w[0], tcr_w[0])
    return {c_w_str: np.array([c_w]),
            b_w_str: np.array([b_w]),
            AR_w_str: np.array([AR_w]),
            fw25_w_str: np.array([fw25_w]),
            Sexp_w_str: np.array([Sexp_w]),
            Swet_w_str: np.array([Swet_w])}


# Vertical Stablizer Geometry #
def compute_vertstab_surface(S_w):
    return 0.1 * S_w


def compute_vertstab_chord_length(S_v, fwLE_v, fwTE_v, tr_v):
    return np.sqrt(2 * S_v * (np.tan(fwLE_v) - np.tan(fwTE_v)) / (1-tr_v**2))


def compute_vertstab_span(S_v, c_v, tr_v, fwLE_v, fwTE_v):
    diftanv = np.tan(fwLE_v) - np.tan(fwTE_v)
    bspanv = c_v*(1.0-tr_v)/diftanv
    return bspanv
    # return 2 * S_v / c_v / (1 - tr_v**2)


def compute_vertstab_aspect_ratio(S_v, b_v):
    return 2.0 * b_v**2 / S_v


def compute_vertstab_wetted_surface(S_v, tcr_v):
    return 2*(1+0.2*tcr_v)*S_v


def compute_vertstab_sweep_angle_25(tr_v, fwTE_v, AR_v):
    e1 = 1.0
    e2 = 0.25
    part2 = (1.0 - tr_v)/(1.0 + tr_v)
    fv25 = np.arctan(np.tan(fwTE_v) + 4.0 / AR_v * part2 * (e1-e2))
    return fv25


def compute_vertstab_frontal_area():
    return 1


def compute_vertstaB_geometry(S_w, fwLE_v, fwTE_v, tcr_v, tr_v):
    S_v = compute_vertstab_surface(S_w)
    c_v = compute_vertstab_chord_length(S_v, fwLE_v, fwTE_v, tr_v)
    b_v = compute_vertstab_span(S_v, c_v, tr_v, fwLE_v, fwTE_v)
    Swet_v = compute_vertstab_wetted_surface(S_v, tcr_v)
    AR_v = compute_vertstab_aspect_ratio(S_v, b_v)
    fv25_v = compute_vertstab_sweep_angle_25(tr_v, fwTE_v, AR_v)
    return S_v, c_v, b_v, Swet_v, fv25_v


def ExecuteVerticalStabilizerGeometry(input_values):
    S_w = input_values[S_w_str]
    fwLE_v = input_values[fwLE_v_str]
    fwTE_v = input_values[fwTE_v_str]
    tcr_v = input_values[tcr_v_str]
    tr_v = input_values[tr_v_str]
    S_v, c_v, b_v, Swet_v, fw25_v = compute_vertstaB_geometry(
        S_w[0], fwLE_v[0], fwTE_v[0], tcr_v[0], tr_v[0])
    return {S_v_str: np.array([S_v]),
            c_v_str: np.array([c_v]),
            b_v_str: np.array([b_v]),
            Swet_v_str: np.array([Swet_v]),
            fw25_v_str: np.array([fw25_v])}

# Fuselage geometry #


def compute_length_1(Mach):
    return D_f/2/np.tan(0.3 * np.arcsin(1/Mach))


def compute_length_2():
    return 10


def compute_length_3(fuel_volume):
    return 4 * fuel_volume / np.pi / D_f**2


def compute_length_4():
    return D_f/2/np.tan(np.deg2rad(12))


def compute_fuselage_surface(l1, l2, l3, l4):
    surface = np.pi * D_f * \
        (l2 + l3 + 0.5 * (np.sqrt(D_f**2 / 4 + l1**2) + np.sqrt(D_f**2 / 4 + l4**2)))
    return surface


def compute_fuselage_length(l1, l2, l3, l4):
    return l1+l2+l3+l4


def compute_fuselage_frontal_area():
    return 1


def compute_fuselage_geometry(Mach, fuel_weight):
    l1 = compute_length_1(Mach)
    l2 = compute_length_2()
    l3 = compute_length_3(fuel_weight/rho_fuel)
    l4 = compute_length_4()
    return compute_fuselage_length(l1, l2, l3, l4), compute_fuselage_surface(l1, l2, l3, l4)


def ExecuteFuselageGeometry(input_values):
    W_fuel = input_values[W_fuel_str]
    Mach = input_values[M_cr_str]
    L_f, S_f = compute_fuselage_geometry(Mach[0], W_fuel[0])
    return {L_f_str: np.array([L_f]), S_f_str: np.array([S_f])}

# Weight #


def compute_takeoff_weight(zfw, fuel_weight):
    return zfw + fuel_weight


def compute_zero_fuel_weight(W_empty):
    return W_pc + W_op + W_empty


def compute_empty_weight(tow, W_w, W_f, W_v, W_eng):
    W_gear = 0.04 * tow
    W_fe = 0.08 * tow
    W_prop = 1.6 * W_eng
    return W_w + W_f + W_v + W_gear + N_eng * W_prop + W_fe


def compute_wing_weight(S_w, b_w, Sexp_w, tcr_w, fw25_w, tr_w, tow, zfw):
    # exposed or ref?
    val = 20.61 * Sexp_w
    term1 = N_ult * b_w**3 * np.sqrt(tow*zfw) * (1 + 2 * tr_w)
    term2 = tcr_w*np.cos(fw25_w)**2 * Sexp_w*(1+tr_w)
    val += 5.39 * 1e-6 * term1/term2
    return val


def compute_vert_stab_weight(S_v, b_v, Sexp_w, tcr_v, fw25_v, tow):
    # coeffs???
    #12.8 or 17.91
    term1 = N_ult*b_v**3 * (8.0 + 0.09012*tow/Sexp_w)
    term2 = (tcr_v*np.cos(fw25_v)**2)
    val = 17.91*S_v + 3.361*1e-4*term1/term2
    return val


def compute_single_engine_weight(P0):
    return 2.2 * 1e-2 * P0**0.985


def compute_fuselage_weight(S_f, L_f, p_cr, W_w, W_eng, zfw):
    T2000 = compute_temperature(2000)
    p_cab = compute_pressure(2000, T2000)
    W_prop = 1.6 * W_eng
    I_p = 1.03*1e-4 * (-p_cr + p_cab) * D_f  # pressure index
    I_b = 1.28*1e-4 * (zfw-W_w-N_wing_mounted_eng*W_prop) * \
        N_ult * L_f / D_f**2  # buckling index
    I_f = I_p if I_p > I_b else (I_p**2 + I_b**2)/(2*I_b)
    val = (5.1314 + 0.498 * I_f) * S_f
    return val


def compute_weight_iter(tow, W_fuel, S_w, b_w, Sexp_w, tcr_w, fw25_w, tr_w, S_v, b_v, tcr_v, fw25_v, P0, S_f, L_f, p_cr):
    zfw = tow - W_fuel
    W_w = compute_wing_weight(S_w, b_w, Sexp_w, tcr_w, fw25_w, tr_w, tow, zfw)
    W_v = compute_vert_stab_weight(S_v, b_v, Sexp_w, tcr_v, fw25_v, tow)
    W_eng = compute_single_engine_weight(P0)
    W_f = compute_fuselage_weight(S_f, L_f, p_cr, W_w, W_eng, zfw)
    W_empty = compute_empty_weight(tow, W_w, W_f, W_v, W_eng)
    zfw = compute_zero_fuel_weight(W_empty)
    tow = compute_takeoff_weight(zfw, W_fuel)
    return tow


def compute_weight(W_fuel, S_w, b_w, Sexp_w, tcr_w, fw25_w, tr_w, S_v, b_v, tcr_v, fw25_v, P0, S_f, L_f, p_cr):

    tow = W_fuel + 5000
    max_it = 500
    tol = 1e-8
    status = False
    for it in range(max_it):
        tow_new = compute_weight_iter(tow, W_fuel, S_w, b_w, Sexp_w, tcr_w, fw25_w, tr_w,
                                      S_v, b_v, tcr_v, fw25_v, P0, S_f, L_f, p_cr)
        if (np.abs(tow_new-tow)/tow_new <= tol):
            status = True
            break
        tow = tow_new
    if (status == False):
        print("TOW NOT CONVERGED")
    zfw = tow - W_fuel
    return tow, zfw


def ExecuteWeight(input_values):
    S_w = input_values[S_w_str]
    tcr_w = input_values[tcr_w_str]
    fw25_w = input_values[fw25_w_str]
    Sexp_w = input_values[Sexp_w_str]
    tr_w = input_values[tr_w_str]
    b_w = input_values[b_w_str]

    S_v = input_values[S_v_str]
    tcr_v = input_values[tcr_v_str]
    fw25_v = input_values[fw25_v_str]
    b_v = input_values[b_v_str]

    S_f = input_values[S_f_str]
    L_f = input_values[L_f_str]

    p_cr = input_values[p_cr_str]
    P0 = input_values[P0_str]
    W_fuel = input_values[W_fuel_str]
    tow, zfw = compute_weight(W_fuel[0], S_w[0], b_w[0], Sexp_w[0], tcr_w[0], fw25_w[0], tr_w[0],
                              S_v[0], b_v[0], tcr_v[0], fw25_v[0], P0[0], S_f[0], L_f[0], p_cr[0])

    return {tow_str: np.array([tow]), zfw_str: np.array([zfw])}


# Engine #
def compute_single_engine_power(CD, V_cr, rho_cr, S_w):
    return 0.5 * rho_cr * V_cr**2 * CD * S_w / N_eng


def compute_P_P0_ratio(p_cr, M_cr):
    return 0.6 * p_cr / p_ref * (1 + (gamma-1)/2 * M_cr**2)**(gamma/(gamma-1))


def compute_specific_fuel_consumption(M_cr, T_cr):
    return 2.83 * 1e-5 * (0.9 + 0.3 * M_cr) * np.sqrt(T_cr/T_ref)


def compute_engine_geometry(P0):
    D_eng = 0.00062 * P0**0.67
    L_eng = 0.025 * P0**0.46
    alin = 4.0*D_eng
    L_nac = alin + L_eng
    D_nac = 1.1 * D_eng
    frontal_area_eng = np.pi*D_nac*D_nac/4.0 - 0.540 * D_eng * D_eng
    return D_eng, L_eng, D_nac, L_nac, frontal_area_eng


def compute_engine_chars(p_cr, rho_cr, T_cr, M_cr,  S_w,  V_cr):
    PP0 = compute_P_P0_ratio(p_cr, M_cr)
    sfc = compute_specific_fuel_consumption(M_cr, T_cr)
    return PP0, sfc


def ExecuteEngine(input_values):
    p_cr = input_values[p_cr_str]
    rho_cr = input_values[rho_cr_str]
    T_cr = input_values[T_cr_str]
    M_cr = input_values[M_cr_str]
    S_w = input_values[S_w_str]
    V_cr = input_values[V_cr_str]
    CD = input_values[CD_str]
    P = compute_single_engine_power(CD[0], V_cr[0], rho_cr[0], S_w[0])
    P0 = P / compute_P_P0_ratio(p_cr[0], M_cr[0])
    sfc = compute_specific_fuel_consumption(M_cr[0], T_cr[0])
    D_eng, L_eng, D_nac, L_nac, fa_eng = compute_engine_geometry(P0)
    return {P_str: np.array([P]),
            P0_str: np.array([P0]),
            sfc_str: np.array([sfc]),
            L_nac_str: np.array([L_nac]),
            D_nac_str: np.array([D_nac]),
            fa_eng_str: np.array([fa_eng])}

# Aerodynamics #


def intterm(dfus,anloc,fle,fte,cr,x,knmax,kn,ktype):
    if ktype == 1: 
        an = dfus * 0.5 + (anloc - dfus * 0.5)*(kn-1) / float(knmax-1)
    elif ktype == 2: 
        an = anloc * (kn-1) / float(knmax-1)
    xle = an * np.tan(fle)
    xte = cr + an * np.tan(fte)
    ee = (x-xle)*(xte-x)/(xte-xle)
    return an, ee

def frontal_area(span, fle, fte, cr, threl, ktype):
    """
    b_w or b_v: span 
    fwLE_w or fwLE_v: fle
    fwTE_w or fwTE_v: fte
    c_w ot c_v: cr
    tcr_w or tcr_v: threl 
    ktype: 1 wing, 2 vtail 
    """
    armax = -1
    knmax = 100 
    kxmax = 100 
    x1 = 0
    x4 = cr
    if ktype == 1:
        x1 = D_f/2 * np.tan(fle)
        x4 = cr + span/2 * np.tan(fte)
    for kx in range(1, kxmax):
        x =  x1 + (x4-x1)*float(kx-1)/float(kxmax-1)
        anbb = min(span , x/np.tan(fle))
        if ktype == 1: 
            anbb = min(span*0.5, x/np.tan(fle))
        xtebb = cr + anbb*np.tan(fte)
        sf = 0.0
        for kn in range(1,knmax-1):
          aa, eea = intterm(D_f,anbb,fle,fte,cr,x,knmax,kn  ,ktype)
          bb, eeb = intterm(D_f,anbb,fle,fte,cr,x,knmax,kn+1,ktype)
          res = (bb-aa)*(eeb+eea) / 2.0
          sf = sf + 4.0 * threl * res 
        area = max(armax,sf)
        armax = area

    return area


def compute_drag_coeff(CD_0, CL, K):
    return CD_0 + K * CL**2


def compute_kappa(M_cr, AR_w, fwLE_w):
    #oswald = 4.61*(1 - 0.045 * AR_w**0.68)*np.cos(fwLE_w)**0.15 - 3.1
    if M_cr <= 1:
        return 1/np.pi/AR_w/1.0
    if M_cr > 1:
        return AR_w*(M_cr**2 - 1)/(4*AR_w*np.sqrt(M_cr**2 - 1)-2)*np.cos(fwLE_w)


def compute_drag_coeff_0(M_cr, T_cr, V_cr, rho_cr, L_f, L_nac, D_nac, fa_eng, c_w, c_v, tr_w,
                         tr_v, S_w, S_v, Swet_w, Swet_v, S_f, fwLE_w, tcr_w, tcr_v, b_w, b_v,
                         fwTE_w, fwLE_v, fwTE_v):
    CD_gear = 0.02
    CD_visc = compute_viscous_drag_coeff(M_cr, T_cr, V_cr, rho_cr, L_f, L_nac, D_nac, c_w,
                                         c_v, tr_w, tr_v, S_w,
                                         S_v, Swet_w, Swet_v, S_f)
    if M_cr < 1:
        return CD_gear + CD_visc
    if M_cr >= 1:
        frontal_area_w = 2.0 * frontal_area(b_w, fwLE_w, fwTE_w, c_w, tcr_w, 1)
        frontal_area_v = frontal_area(b_v, fwLE_v, fwTE_v, c_v, tcr_v, 2) 
        A_max = frontal_area_w + frontal_area_v + frontal_area_f + N_eng * fa_eng
        # print(frontal_area_w)
        # print(frontal_area_v)
        # print(frontal_area_f)
        # print(fa_eng)
        # print(f"AMAX: {A_max}")
        CD_wave = compute_wave_drag_coeff(M_cr, S_w, L_f, fwLE_w, A_max)
        return CD_visc + CD_wave

def compute_wave_drag_coeff(M_cr, S_w, L_f, fwLE_w, A_max):
    searshaack = 9.0 * np.pi / 2.0 * (A_max / L_f)**2
    term2 = 1.0 - 0.3860 * (M_cr - 1.20)**0.57 * \
        (1.0 - np.pi*(fwLE_w*180.0/np.pi)**0.77 / 100)
    cdw = E_wd * term2 * searshaack / S_w
    cdw = 0.75 * cdw
    return cdw


def compute_mean_chord(c, tr):
    return 2*c*(1+tr+tr**2)/(3+3*tr)


def compute_viscous_drag_coeff(M_cr, T_cr, V_cr, rho_cr, L_f, L_nac, D_nac, c_w, c_v, tr_w, tr_v,
                               S_w, S_v, Swet_w, Swet_v, S_f):

    Tw = T_cr * (1.0 + 0.178 * M_cr**2)
    air_visc = 17.15 * 1e-6 * (Tw/273.0)**1.5 * (273+110.4)/(Tw+110.4)
    r1 = 1.0 + 0.1151 * M_cr**2
    r2 = 0.5 * (1.0 + r1) * r1**1.5
    crat = r2**0.2 / r1

    Re_w = rho_cr * V_cr * compute_mean_chord(c_w, tr_w) / air_visc
    Re_v = rho_cr * V_cr * compute_mean_chord(c_v, tr_v) / air_visc
    Re_f = rho_cr * V_cr * L_f / air_visc
    Re_nac = rho_cr * V_cr * L_nac / air_visc
    Re_arr = [Re_w, Re_v, Re_f, Re_nac]
    S_wet = [Swet_w, Swet_v, S_f, N_eng * np.pi*L_nac*D_nac]
    S_ref = S_w

    cd_visc = 0.0
    for i in range(4):
        cd_visc += (0.074 * Re_arr[i]**(-0.2)) * crat * S_wet[i] / S_ref
    return cd_visc


def compute_lift_coeff(tow, rho_cr, V_cr, S_w):
    return 0.95 * g * tow / (0.5 * rho_cr * V_cr**2 * S_w)


def compute_lift_coeff_from_AoA(aoa, AR_w, fwLE_w):
    return 0.5 * np.pi * AR_w * np.cos(aoa) * np.sin(aoa)*(np.cos(aoa) + np.sin(aoa)*np.cos(aoa)/np.cos(fwLE_w) - 0.5*np.sin(aoa)/np.cos(fwLE_w))

def ExecuteAerodynaics(input_values):
    S_w = input_values[S_w_str]
    tcr_w = input_values[tcr_w_str]
    fwLE_w = input_values[fwLE_w_str]
    fwTE_w = input_values[fwTE_w_str]
    Swet_w = input_values[Swet_w_str]
    tr_w = input_values[tr_w_str]
    b_w = input_values[b_w_str]
    c_w = input_values[c_w_str]
    AR_w = input_values[AR_w_str]

    S_v = input_values[S_v_str]
    tcr_v = input_values[tcr_v_str]
    b_v = input_values[b_v_str]
    Swet_v = input_values[Swet_v_str]
    c_v = input_values[c_v_str]
    tr_v = input_values[tr_v_str]
    fwLE_v = input_values[fwLE_v_str]
    fwTE_v = input_values[fwTE_v_str]

    S_f = input_values[S_f_str]
    L_f = input_values[L_f_str]

    T_cr = input_values[T_cr_str]
    rho_cr = input_values[rho_cr_str]
    M_cr = input_values[M_cr_str]
    V_cr = input_values[V_cr_str]

    L_nac = input_values[L_nac_str]
    D_nac = input_values[D_nac_str]
    fa_eng = input_values[fa_eng_str]

    tow = input_values[tow_str]

    CD_0 = compute_drag_coeff_0(M_cr[0], T_cr[0], V_cr[0], rho_cr[0], L_f[0], L_nac[0], D_nac[0], fa_eng[0], c_w[0], c_v[0], tr_w[0], tr_v[0], S_w[0],
                                    S_v[0], Swet_w[0], Swet_v[0], S_f[0], fwLE_w[0], tcr_w[0], tcr_v[0], b_w[0], b_v[0], 
                                    fwTE_w[0], fwLE_v[0], fwTE_v[0])
    K = compute_kappa(M_cr[0], AR_w[0], fwLE_w[0])
    CL = compute_lift_coeff(tow[0], rho_cr[0], V_cr[0], S_w[0])
    CD = compute_drag_coeff(CD_0, CL, K)
    
    return {CL_str: np.array([CL]), CD_str: np.array([CD])}

# Performance #


def compute_range(V_cr, sfc, CL, CD, tow, W_fuel):
    W_ec = tow - 0.95 * W_fuel
    W_bc = 0.95 * tow
    Range = V_cr/g/sfc * CL/CD * np.log(W_bc/W_ec)
    return Range / 1e7


def compute_approach_velocity(S_w, tow, AR_w, fwLE_w):
    mlw = 0.8 * tow  # max_landing_weight
    CL_max = compute_lift_coeff_from_AoA(np.deg2rad(14), AR_w, fwLE_w)
    V_sep = np.sqrt(g*mlw/(0.5*rho_ref*S_w*CL_max))
    return 1.3 * V_sep


def ExecutePerformance(input_values):
    V_cr = input_values[V_cr_str]
    sfc = input_values[sfc_str]
    CL = input_values[CL_str]
    CD = input_values[CD_str]
    tow = input_values[tow_str]
    W_fuel = input_values[W_fuel_str]
    Range = compute_range(V_cr[0], sfc[0], CL[0], CD[0], tow[0], W_fuel[0])
    return {Range_str: np.array([Range])}
