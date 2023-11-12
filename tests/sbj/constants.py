""" CONSTANT VALUES """
g = 9.8
R = 287.04
gamma = 1.4
T_ref = 288.15
dTdz = -0.0065
p_ref = 101325
z11 = 11e3
T11 = T_ref + dTdz * z11
p11 = p_ref * (T11/T_ref) ** (-g/dTdz/R)
rho_ref = p_ref / T_ref / R
D_f = 2.15
frontal_area_f = (D_f/2)**2 * 3.14
rho_fuel = 807.5
N_eng = 2
N_wing_mounted_eng = 2
N_pass = 10 
N_crew = 3 
person_weight = 100
W_op = 210 
W_pc = 1500
N_ult = 4.5 
E_wd = 2.0
c1 = 0.00062
k1 = 0.67
c2 = 0.025
k2 = 0.46

"""
DESIGN VARIABLE NAMES 
"""
# Cruising altitude and Mach number
z_cr_str = "zcr"
M_cr_str = "Mcr"

# Wing geometry design variables
# w subscript denotes wing
S_w_str = "Sw"
fwLE_w_str = "ΛLEw"
fwTE_w_str = "ΛTEw"
tr_w_str = "λw"
tcr_w_str = "(t_c)w"

# Vertical stabilizer geometry design variables
# v subscript denotes vertical stabilizer
fwLE_v_str = "ΛLEv"
fwTE_v_str = "ΛTEv"
tr_v_str = "λv"
tcr_v_str = "(t_c)v"

# Fuel weight 
W_fuel_str = "Wfuel"

"""
COUPLING VARIABLE NAMES
"""
# Environemnt
p_cr_str = "pcr"
rho_cr_str = "rhocr"
T_cr_str = "Tcr"
V_cr_str = "Vcr"

# Engine
P_str = "T"
P0_str = "T0"
sfc_str = "SFC" 
L_nac_str = "Lnac"
D_nac_str = "Dnac"
fa_eng_str = "faeng"

# Geometry 
c_w_str = "crw"
b_w_str = "bw"
AR_w_str = "ARw"
fw25_w_str = "Λ25w"
Sexp_w_str = "Sexpw"
Swet_w_str = "Swetw"

S_v_str = "Sv"
c_v_str = "crv"
b_v_str = "bv"
fw25_v_str = "Λ25v"
Swet_v_str = "Swetv"

# f subscript denotes fuselage
S_f_str = "Sf"
L_f_str = "Lf"

# Weight
tow_str = "TOW"
zfw_str = "ZFW"

# Aerodynamics 
CL_str = "CL"
CD_str = "CD"

# Performance 
Range_str = "range"