'''
Homework, item 5.4 & 5.5
Author: Reynaldo S. Lima
AP-266
'''

# IMPORTS
import numpy as np
from asa_module import asa_module as asa
from asa_module_b import asa_module_b as asab
from scipy.optimize import root

# Now, the add. parameters
Xfem = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
alpha0 = np.zeros(4)
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
Vinfm = 60
alpha = 5 * np.pi / 180
Vinf = np.array([Vinfm * np.cos(alpha), 0.0, Vinfm * np.sin(alpha)])
rhoAir = 1.225
R = np.array([0.1, 0.1, 0.1, 0.1])
t = np.array([0.005, 0.005, 0.005, 0.005])
E = 73.1E9
rhoMat = 2780.0
sigmaY = 324.0E6
pKS = 200.0
conIDs = np.array([5, 6], dtype='int32')
CD0 = 0.0270
fixedMass = 700.0
g = 9.8
Endurance = 4.0 * 60.0 * 60.0
TSFC = 0.5 / 3600.0
loadFactor = 3.0 * 1.5

# First, lets analyse the resfunc
def resfunc(stateVars):
    n = np.int((len(stateVars) - 2) / 3)
    Gama = stateVars[0:n] / 0.1
    d = stateVars[n:3 * n + 2] / 100.0
    # General constants for this test
    Xfem = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')

    alpha0 = np.zeros(4)
    chords = np.array([1.0, 1.0, 1.0, 1.0])
    cl0 = np.zeros(4)
    cla = np.array([6.283, 6.283, 6.283, 6.283])
    Vinfm = 60
    alpha = 5 * np.pi / 180
    Vinf = np.array([Vinfm * np.cos(alpha), 0.0, Vinfm * np.sin(alpha)])
    rhoAir = 1.225
    R = np.array([0.1, 0.1, 0.1, 0.1])
    t = np.array([0.005, 0.005, 0.005, 0.005])
    E = 73.1E9
    rhoMat = 2780.0
    sigmaY = 324.0E6
    pKS = 200.0
    conIDs = np.array([5, 6], dtype='int32')
    CD0 = 0.0270
    fixedMass = 700.0
    g = 9.8
    Endurance = 4.0 * 60.0 * 60.0
    TSFC = 0.5 / 3600.0
    loadFactor = 3.0 * 1.5

    resllt, resfem, liftExcess, margins, KSmargin, FB, Weight, Sref, CL = asa.asa_analysis(Gama, alpha0, chords,
                                                                                           Xfem, R, t, d,
                                                                                           cl0, cla, Vinf, rhoAir,
                                                                                           E, rhoMat, sigmaY, pKS,
                                                                                           conIDs,
                                                                                           CD0, fixedMass, g, Endurance,
                                                                                           TSFC, loadFactor)

    res = np.hstack([resllt, resfem])
    return res


Gama0 = np.array([80.0, 80.0, 80.0, 80.0])
d0 = 1E-3 * np.ones(10)
stateVars = np.concatenate((Gama0, d0))
sol = root(resfunc, stateVars)
stateVars = sol.x
n = np.int((len(stateVars) - 2) / 3)
Gama = stateVars[0:n] / 0.1
d = stateVars[n:3 * n + 2] / 100.0
print('----------------------- RES TEST ------------------------')
print('---------------------------------------------------------')
print(' ')
print('Gama = ', Gama)
print(' ')
print('d = ', d)
print('---------------------------------------------------------')


# Then, let's see the adjoint problem

def adjfunc(resb):
    n = np.int((len(resb) - 2) / 3)
    reslltb = resb[0:n]
    resfemb = resb[n:3 * n + 2]

    # Output values initialized as 0.0
    resllt = np.zeros(4)
    resfem = np.zeros(10)
    liftexcess = 0.0
    margins = np.zeros(8)
    ksmargin = 0.0
    fb = 0.0
    weight = 0.0
    sref = 0.0
    cl = 0.0

    # Input seeds initialized as 0.0
    tb = np.zeros(4)
    db = np.zeros(10)
    alpha0b = np.zeros(4)
    Gamab = np.zeros(4)

    # Seeds of the outputs
    liftexcessb = 0.0
    marginsb = np.zeros(8)
    ksmarginb = 1.0
    fbb = 0.0
    weightb = 0.0
    print(reslltb)
    asab.asa_analysis_b(Gama, Gamab, alpha0, alpha0b, chords, Xfem, R, t, tb, d, db, cl0, cla, Vinf,
                        rhoAir, E, rhoMat,
                        sigmaY, pKS, conIDs, CD0, fixedMass, g, Endurance, TSFC, loadFactor, resllt,
                        reslltb.copy(), resfem, resfemb.copy(), liftexcess, liftexcessb, margins, marginsb,
                        ksmargin, ksmarginb, fb, fbb, weight, weightb, sref, cl)

    stateVarsb = np.hstack([Gamab, db])
    return stateVarsb


resb0 = np.array([10 ** 6, 10 ** 6, 10 ** 6, 10 ** 6, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
sol = root(adjfunc, resb0)
resb = sol.x
psi = resb[0:n]
lambd = resb[n:3 * n + 2]
reslltb = psi
resfemb = lambd
print(' ')
print('----------------------- ADJ TEST ------------------------')
print('---------------------------------------------------------')
print(' ')
print('psi = ', psi)
print(' ')
print('lambda =', lambd)
print(' ')

# Output values initialized as 0.0
resllt = np.zeros(4)
resfem = np.zeros(10)
liftexcess = 0.0
margins = np.zeros(8)
ksmargin = 0.0
fb = 0.0
weight = 0.0
sref = 0.0
cl = 0.0

# Input seeds initialized as 0.0
tb = np.zeros(4)
db = np.zeros(10)
alpha0b = np.zeros(4)
Gamab = np.zeros(4)

# Seeds of the outputs
liftexcessb = 0.0
marginsb = np.zeros(8)
ksmarginb = 1.0
fbb = 0.0
weightb = 0.0

asab.asa_analysis_b(Gama, Gamab, alpha0, alpha0b, chords, Xfem, R, t, tb, d, db, cl0, cla, Vinf,
                        rhoAir, E, rhoMat,
                        sigmaY, pKS, conIDs, CD0, fixedMass, g, Endurance, TSFC, loadFactor, resllt,
                        reslltb.copy(), resfem, resfemb.copy(), liftexcess, liftexcessb, margins, marginsb,
                        ksmargin, ksmarginb, fb, fbb, weight, weightb, sref, cl)

print('dKSda0 =', alpha0b)
print(' ')
print('dKSdt = ', tb)
print('---------------------------------------------------------')