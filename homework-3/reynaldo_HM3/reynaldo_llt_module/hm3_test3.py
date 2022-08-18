'''
Homework 3, item 2.3
Author: Reynaldo S. Lima
AP-266
'''

# IMPORTS
import numpy as np
from llt_module import llt_module as llt
from llt_module_b import llt_module_b as llt_b
from scipy.optimize import root


# GAMMA FUNCTION
def gamafunc(Gamma):
    X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
    alpha0 = np.zeros(4)
    chords = np.array([1.0, 1.0, 1.0, 1.0])
    cl0 = np.zeros(4)
    cla = np.array([6.283, 6.283, 6.283, 6.283])
    alpha = 5.498 * np.pi / 180
    Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
    rho = 1.0
    [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, Gamma, alpha0, chords, cl0, cla, Vinf, rho)
    return res


# INITIAL INPUT
Gama = np.array([1.0, 2.0, 2.0, 1.0])

# SOLUTION FOR GAMA
sol = root(gamafunc, Gama)
Gama = sol.x
X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
alpha0 = np.zeros(4)
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
alpha = 5.498 * np.pi / 180
Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
rho = 1.0
[res, Sref, CL, CD, L, D] = llt.tapenade_main(X, Gama.copy(), alpha0, chords, cl0, cla, Vinf, rho)


# ADJOINT METHOD

# FIRST, THE RESIDUAL FUNCTION

def adjEqResidual(resb):
    x_f = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
    gama_f = sol.x
    alpha0_f = np.zeros(4)
    chords_f = np.array([1.0, 1.0, 1.0, 1.0])
    cl0_f = np.zeros(4)
    cla_f = np.array([6.283, 6.283, 6.283, 6.283])
    alpha_f = 5.498 * np.pi / 180
    Vinf_f = np.array([np.cos(alpha_f), 0.0, np.sin(alpha_f)])
    rho_f = 1.0
    Srefb_f = 0.0
    CLb_f = 1.0  # Thus, we are calculating d()/dCL
    CDb_f = 0.0
    Lb_f = 0.0
    Db_f = 0.0
    res_f = np.zeros(4)
    Xb_f = np.zeros((3, 5), order='F')
    Gamab_f = np.zeros(4)
    alpha0b_f = np.zeros(4)
    chordsb_f = np.zeros(4)
    Sref_f = 0.0
    CL_f = 0.0
    CD_f = 0.0
    L_f = 0.0
    D_f = 0.0

    llt_b.tapenade_main_b(x_f, Xb_f, gama_f.copy(), Gamab_f, alpha0_f, alpha0b_f, chords_f, chordsb_f, cl0_f, cla_f, Vinf_f,
                          rho_f, res_f, resb.copy(), Sref_f, Srefb_f, CL_f, CLb_f, CD_f, CDb_f, L_f, Lb_f, D_f, Db_f)
    return Gamab_f

def adjEqResidual2(resb):
    x_f = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
    gama_f = sol.x
    alpha0_f = np.zeros(4)
    chords_f = np.array([1.0, 1.0, 1.0, 1.0])
    cl0_f = np.zeros(4)
    cla_f = np.array([6.283, 6.283, 6.283, 6.283])
    alpha_f = 5.498 * np.pi / 180
    Vinf_f = np.array([np.cos(alpha_f), 0.0, np.sin(alpha_f)])
    rho_f = 1.0
    Srefb_f = 1.0
    CLb_f = 0.0  # Thus, we are calculating d()/dCL
    CDb_f = 0.0
    Lb_f = 0.0
    Db_f = 0.0
    res_f = np.zeros(4)
    Xb_f = np.zeros((3, 5), order='F')
    Gamab_f = np.zeros(4)
    alpha0b_f = np.zeros(4)
    chordsb_f = np.zeros(4)
    Sref_f = 0.0
    CL_f = 0.0
    CD_f = 0.0
    L_f = 0.0
    D_f = 0.0

    llt_b.tapenade_main_b(x_f, Xb_f, gama_f.copy(), Gamab_f, alpha0_f, alpha0b_f, chords_f, chordsb_f, cl0_f, cla_f, Vinf_f,
                          rho_f, res_f, resb.copy(), Sref_f, Srefb_f, CL_f, CLb_f, CD_f, CDb_f, L_f, Lb_f, D_f, Db_f)
    return Gamab_f

def adjEqResidual3(resb):
    x_f = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
    gama_f = sol.x
    alpha0_f = np.zeros(4)
    chords_f = np.array([1.0, 1.0, 1.0, 1.0])
    cl0_f = np.zeros(4)
    cla_f = np.array([6.283, 6.283, 6.283, 6.283])
    alpha_f = 5.498 * np.pi / 180
    Vinf_f = np.array([np.cos(alpha_f), 0.0, np.sin(alpha_f)])
    rho_f = 1.0
    Srefb_f = 0.0
    CLb_f = 0.0
    CDb_f = 1.0 # Thus, we are calculating d()/dCD
    Lb_f = 0.0
    Db_f = 0.0
    res_f = np.zeros(4)
    Xb_f = np.zeros((3, 5), order='F')
    Gamab_f = np.zeros(4)
    alpha0b_f = np.zeros(4)
    chordsb_f = np.zeros(4)
    Sref_f = 0.0
    CL_f = 0.0
    CD_f = 0.0
    L_f = 0.0
    D_f = 0.0

    llt_b.tapenade_main_b(x_f, Xb_f, gama_f.copy(), Gamab_f, alpha0_f, alpha0b_f, chords_f, chordsb_f, cl0_f, cla_f, Vinf_f,
                          rho_f, res_f, resb.copy(), Sref_f, Srefb_f, CL_f, CLb_f, CD_f, CDb_f, L_f, Lb_f, D_f, Db_f)
    return Gamab_f

psi0 = 10 ** 6 * np.ones(4)
step2 = root(adjEqResidual, psi0)
step3 = np.zeros(4) #SREF DOES NOT NEED THIS STEP
step4 = root(adjEqResidual3, psi0)
print('-------------------------')
print(step2)
print('-------------------------')
print(' ')
print('-------------------------')
print(step3)
print('-------------------------')
print(' ')
print('-------------------------')
print(step4)
print('-------------------------')
print(' ')
print('-------------------------')

# IF TRUE, WE HAVE PSI
psi = step2.x
psi2 = step3
psi3 = step4.x
resb = psi.copy()
resb2 = psi2.copy()
resb3 = psi3.copy()

# START OF FIRST DERIVATIVE
X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
Gama = sol.x
alpha0 = np.zeros(4)
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
alpha = 5.498 * np.pi / 180
Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
rho = 1.0
Srefb = 0.0
CLb = 1.0  # Thus, we are calculating d()/dCL
CDb = 0.0
Lb = 0.0
Db = 0.0
res = np.zeros(4)
Xb = np.zeros((3, 5), order='F')
Gamab = np.zeros(4)
alpha0b = np.zeros(4)
chordsb = np.zeros(4)
Sref = 0.0
CL = 0.0
CD = 0.0
L = 0.0
D = 0.0
llt_b.tapenade_main_b(X, Xb, Gama.copy(), Gamab, alpha0, alpha0b, chords, chordsb, cl0, cla, Vinf,
                      rho, res, resb, Sref, Srefb, CL, CLb, CD, CDb, L, Lb, D, Db)
dCLdX = Xb
dCLdalpha = alpha0b
dCLdchords = chordsb
print('1st results: ')
print('dCLdX = ...')
print(dCLdX)
print('dCLda0 = ', dCLdalpha)
print('dCLdc = ', dCLdchords)
# END OF FIRST DERIVATIVE

# START OF 2ND DERIVATIVE
X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
Gama = sol.x
alpha0 = np.zeros(4)
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
alpha = 5.498 * np.pi / 180
Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
rho = 1.0
Srefb = 0.0
CLb = 1.0  # Thus, we are calculating d()/dCL
CDb = 0.0
Lb = 0.0
Db = 0.0
res = np.zeros(4)
Xb = np.zeros((3, 5), order='F')
Gamab = np.zeros(4)
alpha0b = np.zeros(4)
chordsb = np.zeros(4)
Sref = 0.0
CL = 0.0
CD = 0.0
L = 0.0
D = 0.0
llt_b.tapenade_main_b(X, Xb, Gama.copy(), Gamab, alpha0, alpha0b, chords, chordsb, cl0, cla, Vinf,
                      rho, res, resb2, Sref, Srefb, CL, CLb, CD, CDb, L, Lb, D, Db)
dCLdX = Xb
dCLdalpha = alpha0b
dCLdchords = chordsb
print(' ')
print('-------------------------')
print(' ')
print('2nd results: ')
print('dSrefdX = ...')
print(dCLdX)
print('dSrefda0 = ', dCLdalpha)
print('dSrefdc = ', dCLdchords)
print('-------------------------')
print(' ')
# END OF 2ND DERIVATIVE

# START OF 3RD DERIVATIVE
X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
Gama = sol.x
alpha0 = np.zeros(4)
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
alpha = 5.498 * np.pi / 180
Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
rho = 1.0
Srefb = 0.0
CLb = 1.0  # Thus, we are calculating d()/dCL
CDb = 0.0
Lb = 0.0
Db = 0.0
res = np.zeros(4)
Xb = np.zeros((3, 5), order='F')
Gamab = np.zeros(4)
alpha0b = np.zeros(4)
chordsb = np.zeros(4)
Sref = 0.0
CL = 0.0
CD = 0.0
L = 0.0
D = 0.0
llt_b.tapenade_main_b(X, Xb, Gama.copy(), Gamab, alpha0, alpha0b, chords, chordsb, cl0, cla, Vinf,
                      rho, res, resb3, Sref, Srefb, CL, CLb, CD, CDb, L, Lb, D, Db)
dCLdX = Xb
dCLdalpha = alpha0b
dCLdchords = chordsb
print('-------------------------')
print(' ')
print('3rd results: ')
print('dCDdX = ...')
print(dCLdX)
print('dCDda0 = ', dCLdalpha)
print('dCDdc = ', dCLdchords)
print('-------------------------')
print(' ')
# END OF 3RD DERIVATIVE