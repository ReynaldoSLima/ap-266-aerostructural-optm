'''
Homework, item 5.6
Author: Reynaldo S. Lima
AP-266
'''

# IMPORTS
import numpy as np
from asa_module import asa_module as asa
from asa_module_b import asa_module_b as asab
from scipy.optimize import root, minimize
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
from colour import Color
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)

# Specific to the problem parameters
nP = 20
AR = 6.0
sRef = 16.0
bRef = np.sqrt(AR * sRef)
cRef = sRef / bRef
chords = np.ones(nP) * cRef
Xfem = np.array([np.zeros(nP + 1), np.linspace(0.0, bRef, num=nP + 1), np.zeros(nP + 1)], order='F')
conIDs = np.array([nP + 1, nP + 2], dtype='int32')

# Additional parameters
cl0 = np.zeros(nP)
cla = np.ones(nP) * 6.283
Vinfm = 60
alpha = 5 * np.pi / 180
Vinf = np.array([Vinfm * np.cos(alpha), 0.0, Vinfm * np.sin(alpha)])
rhoAir = 1.225
R = np.ones(nP) * 0.05 * cRef
E = 73.1E9
rhoMat = 2780.0
sigmaY = 324.0E6
pKS = 200.0
CD0 = 0.0270
fixedMass = 700.0
g = 9.8
Endurance = 4.0 * 60.0 * 60.0
TSFC = 0.5 / 3600.0
loadFactor = 3.0 * 1.5

of = []
ah = []
th = []
L1 = []
KS = []
def callback(x):
    of.append(objFunc(x))
    ah.append(x[0:nP])
    th.append(x[nP:2 * nP] / 100.0)
    L1.append(conFunEq(x))
    KS.append(conFunIneq(x))

def solve_asa(desVars):
    nP = np.int(len(desVars) / 2)
    alpha0 = desVars[0:nP]
    t = desVars[nP:2 * nP] / 100.0

    # Defining our resfunc
    def resfunc(stateVars):
        Gama = stateVars[0:nP] / 0.1
        d = stateVars[nP:3 * nP + 2] / 100.0

        # Running asa_analysis
        resllt, resfem, liftExcess, margins, KSmargin, FB, Weight, Sref, CL = asa.asa_analysis(Gama, alpha0, chords,
                                                                                               Xfem, R, t, d,
                                                                                               cl0, cla, Vinf, rhoAir,
                                                                                               E, rhoMat, sigmaY, pKS,
                                                                                               conIDs,
                                                                                               CD0, fixedMass, g,
                                                                                               Endurance,
                                                                                               TSFC, loadFactor)
        res = np.hstack([resllt, resfem])
        return res

    # Finding which Gama and d satisfies the problem for our given alpha0/t
    Gama0 = np.ones(nP) * 8.0
    d0 = 1E-1 * np.ones(2 * (nP + 1))
    stateVars = np.concatenate((Gama0, d0))
    sol = root(resfunc, stateVars)
    stateVars = sol.x
    Gama = stateVars[0:nP] / 0.1
    d = stateVars[nP:3 * nP + 2] / 100.0
    resllt, resfem, liftExcess, margins, KSmargin, FB, Weight, Sref, CL = asa.asa_analysis(Gama, alpha0, chords,
                                                                                           Xfem, R, t, d,
                                                                                           cl0, cla, Vinf, rhoAir,
                                                                                           E, rhoMat, sigmaY, pKS,
                                                                                           conIDs,
                                                                                           CD0, fixedMass, g,
                                                                                           Endurance,
                                                                                           TSFC, loadFactor)
    return resllt, resfem, liftExcess, margins, KSmargin, FB, Weight, Sref, CL, Gama


def compute_grads(desVars, liftExcessb, KSmarginb, FBb, Weightb):
    nP = np.int(len(desVars) / 2)
    alpha0 = desVars[0:nP]
    t = desVars[nP:2 * nP] / 100.0

    a = liftExcessb
    b = KSmarginb
    c = FBb
    w = Weightb

    # Defining our resfunc
    def resfunc(stateVars):
        nP = np.int((len(stateVars) - 2) / 3)

        Gama = stateVars[0:nP] / 0.1
        d = stateVars[nP:3 * nP + 2] / 100.0

        resllt, resfem, liftExcess, margins, KSmargin, FB, Weight, Sref, CL = asa.asa_analysis(Gama, alpha0, chords,
                                                                                               Xfem, R, t, d,
                                                                                               cl0, cla, Vinf, rhoAir,
                                                                                               E, rhoMat, sigmaY, pKS,
                                                                                               conIDs,
                                                                                               CD0, fixedMass, g,
                                                                                               Endurance,
                                                                                               TSFC, loadFactor)

        res = np.hstack([resllt, resfem])
        return res

    # Finding which Gama and d satisfies the problem for our given alpha0/t
    Gama0 = np.ones(nP) * 8.0
    d0 = 1E-1 * np.ones(2 * (nP + 1))
    stateVars = np.concatenate((Gama0, d0))
    sol = root(resfunc, stateVars)
    stateVars = sol.x
    Gama = stateVars[0:nP] / 0.1
    d = stateVars[nP:3 * nP + 2] / 100.0

    # Solving the aerostructural adjoint:

    def adjfunc(resb):
        reslltb = resb[0:nP]
        resfemb = resb[nP:3 * nP + 2]

        # Output values initialized as 0.0
        resllt = np.zeros(nP)
        resfem = np.zeros(2 * (nP + 1))
        liftexcess = 0.0
        margins = np.zeros(2 * nP)
        ksmargin = 0.0
        fb = 0.0
        weight = 0.0
        sref = 0.0
        cl = 0.0

        # Input seeds initialized as 0.0
        tb = np.zeros(nP)
        db = np.zeros(2 * (nP + 1))
        alpha0b = np.zeros(nP)
        Gamab = np.zeros(nP)

        # Seeds of the outputs
        marginsb = np.zeros(2 * nP)

        liftExcessb = a
        KSmarginb = b
        FBb = c
        Weightb = w

        asab.asa_analysis_b(Gama.copy(), Gamab, alpha0.copy(), alpha0b, chords, Xfem, R, t.copy(), tb, d.copy(), db, cl0
                            , cla, Vinf, rhoAir, E, rhoMat,
                            sigmaY, pKS, conIDs, CD0, fixedMass, g, Endurance, TSFC, loadFactor, resllt,
                            reslltb.copy(), resfem, resfemb.copy(), liftexcess, liftExcessb, margins, marginsb,
                            ksmargin, KSmarginb, fb, FBb, weight, Weightb, sref, cl)
        stateVarsb = np.hstack([Gamab, db])
        return stateVarsb

    resb0 = np.concatenate((10 ** 6 * np.ones(nP), np.ones(2 * (nP + 1))))
    sol = root(adjfunc, resb0)
    resb = sol.x
    psi = resb[0:nP]
    lambd = resb[nP:3 * nP + 2]
    reslltb = psi
    resfemb = lambd
    # Output values initialized as 0.0
    resllt = np.zeros(nP)
    resfem = np.zeros(2 * (nP + 1))
    liftexcess = 0.0
    margins = np.zeros(2 * nP)
    ksmargin = 0.0
    fb = 0.0
    weight = 0.0
    sref = 0.0
    cl = 0.0

    # Input seeds initialized as 0.0
    tb = np.zeros(nP)
    db = np.zeros(2 * (nP + 1))
    alpha0b = np.zeros(nP)
    Gamab = np.zeros(nP)

    # Seeds of the outputs
    marginsb = np.zeros(2 * nP)

    liftExcessb = a
    KSmarginb = b
    FBb = c
    Weightb = w

    asab.asa_analysis_b(Gama, Gamab, alpha0, alpha0b, chords, Xfem, R, t, tb, d, db, cl0
                        , cla, Vinf, rhoAir, E, rhoMat,
                        sigmaY, pKS, conIDs, CD0, fixedMass, g, Endurance, TSFC, loadFactor, resllt,
                        reslltb.copy(), resfem, resfemb.copy(), liftexcess, liftExcessb, margins, marginsb,
                        ksmargin, KSmarginb, fb, FBb, weight, Weightb, sref, cl)
    dYdt = tb / 100.0
    dYda0 = alpha0b
    return np.hstack([dYda0, dYdt])


def objFunc(desVars):
    output = solve_asa(desVars)
    return output[5]


def objFuncGrad(desVars):
    dFBdVars = compute_grads(desVars, 0.0, 0.0, 1.0, 0.0)
    return dFBdVars


def conFunEq(desVars):
    output = solve_asa(desVars)
    return output[2]


def conFunEqGrad(desVars):
    dDLdVars = compute_grads(desVars, 1.0, 0.0, 0.0, 0.0)
    return dDLdVars


def conFunIneq(desVars):
    output = solve_asa(desVars)
    return output[4]


def conFunIneqGrad(desVars):
    dKSdVars = compute_grads(desVars, 0.0, 1.0, 0.0, 0.0)
    return dKSdVars


con1 = {'type': 'eq',
        'fun': conFunEq,
        'jac': conFunEqGrad}

con2 = {'type': 'ineq',
        'fun': conFunIneq,
        'jac': conFunIneqGrad}

# Let's stop for a second to admire our start point:

desVars0 = np.concatenate((np.ones(nP) * 0.0, np.ones(nP) * 0.5))
a00 = desVars0[0:nP]
t00 = desVars0[nP:2 * nP] / 100.0
np.set_printoptions(precision=4)
output0 = solve_asa(desVars0)
print('-------------------- OPTIMIZATION -----------------------')
print('-------------------- START POINT ------------------------')
print(' ')
print('---------------------- OBJ FUNC. ------------------------')
print('FB =', "%.4f" % output0[5])
print(' ')
print('---------------------- DES VARS. ------------------------')
print('a0 =', a00)
print(' ')
print('t0 =', t00)
print(' ')
print('---------------------- CON FUNC. ------------------------')
print('DL =', "%.4f" % output0[2])
print(' ')
print('KS =', "%.4f" % output0[4])
print(' ')
print('---------------------------------------------------------')
print('--------------------- END  POINT ------------------------')
bounds1 = [None, None] * nP
bounds2 = [0.1, None] * nP
bounds = np.concatenate([bounds1, bounds2])
bounds = bounds.reshape((2 * nP, 2))
# Run optimizer
result = minimize(objFunc, desVars0, jac=objFuncGrad,
                  constraints=[con1, con2],
                  bounds=bounds,
                  method='slsqp',
                  callback=callback)


# Final output
desVars = result.x
output = solve_asa(desVars)
a0 = desVars[0:nP]
t0 = desVars[nP:2 * nP] / 100.0
Gama = output[9]
print(' ')
print('---------------------- OBJ FUNC. ------------------------')
print('FB =', "%.4f" % output[5])
print(' ')
print('---------------------- DES VARS. ------------------------')
print('a0 =', a0)
print(' ')
print('t0 =', t0)
print(' ')
print('---------------------- CON FUNC. ------------------------')
print('DL =', "%.4f" % output[2])
print(' ')
print('KS =', "%.4f" % output[4])
print(' ')
print('---------------------------------------------------------')
print(' ')
# Weight fractions:
Weight = output[6]
ratioFB = output[5] / Weight
ratioFM = fixedMass *g / Weight
ratioSM = 1 - ratioFB - ratioFM
print('-------------------- WEIGHT RATIO -----------------------')
print('---------------------------------------------------------')
print(' ')
print('FB =', "%.2f" % (ratioFB * 100), '%')
print(' ')
print('FM =', "%.2f" % (ratioFM * 100), '%')
print(' ')
print('SM =', "%.2f" % (ratioSM * 100), '%')
print('---------------------------------------------------------')

margins = output[3]
yM = np.zeros(nP * 2)
prev = 1
for x in range(0, 2 * nP):
    if x % 2 == 0:
        yM[x] = Xfem[1, np.int(x / 2)]
        prev = np.int(x / 2 + 1)
    else:
        yM[x] = Xfem[1, prev]

# Ploting:
fig, ax = plt.subplots()
ax.plot(yM, margins, color='black', linewidth=2.0)

# Customizing the plot:

ax.tick_params(axis="y", labelsize=10)
ax.tick_params(axis="x", labelsize=10)
plt.xlabel('Wingspan', fontsize=20)
plt.ylabel('Margins', fontsize=20)
ax.set_facecolor('white')
ax.patch.set_facecolor('white')
plt.title('Margins, FB optimization', fontsize=20)
plt.show()


# Elliptical circulation
CL = output[8]
gama0 = 2 * Vinfm * sRef * CL / (np.pi * bRef)
y = np.linspace(0.0, bRef, num=nP)
Gamae = 2 * gama0 * np.sqrt(y / bRef - (y ** 2) / bRef ** 2)

# Ploting:
fig, ax = plt.subplots()
ax.plot(y, Gamae, label='Elliptical', color='mediumblue', ls='--', linewidth=1.5)
ax.plot(y, Gama, label='FB Minimization', color='black', linewidth=2.0)

# Customizing the plot:

ax.tick_params(axis="y", labelsize=10)
ax.tick_params(axis="x", labelsize=10)
plt.xlabel('Wingspan', fontsize=20)
plt.ylabel('Circulation', fontsize=20)
legend = ax.legend(loc='best')
ax.set_facecolor('white')
ax.patch.set_facecolor('white')
plt.title('Circulation distribution', fontsize=20)
plt.show()

# Ploting:
fig, ax = plt.subplots()
black = Color("black")
for i in range(-1, len(of)+1):
    if i == -1:
        plt.plot(output0[5], 0.0, marker='d', markersize=20.0,
                 color=[0.0, 0.0, 0.0], markeredgecolor='red')
    elif i == len(of):
        print('oi')
        plt.plot(output[5], 0.0, marker='o', markersize=20.0 - 10 * (len(of) - i) / len(of),
                 color=[0.0, 0.0, i * 1.0 / len(of)], markeredgecolor='red')
    else:
        plt.plot(of[i], 0.0, marker='o', markersize=20.0 - 10 * (len(of) - i) / len(of),
                 color=[0.0, 0.0, i * 1.0 / len(of)], ls='None')


ax.set_xlim(1450, 2010)
# Customizing the plot:
ax.get_yaxis().set_visible(False)
ax.tick_params(axis="x", labelsize=10)
plt.xlabel('FB value', fontsize=12)
ax.set_facecolor('white')
ax.patch.set_facecolor('white')
plt.title('FB value evolution, from black to blue', fontsize=20)
plt.show()

# Ploting:
fig, ax = plt.subplots(1, 2)
black = Color("black")
for i in range(-1, len(of)+1):
    if i == -1:
        ax[0].plot(output0[2], 0.0, marker='d', markersize=20.0,
                   color=[0.0, 0.0, 0.0], markeredgecolor='red')
        ax[1].plot(output0[4], 0.0, marker='d', markersize=20.0,
                   color=[0.0, 0.0, 0.0], markeredgecolor='red')
    elif i == len(of):
        ax[0].plot(output[2], 0.0, marker='o', markersize=20.0 - 10 * (len(of) - i) / len(of),
                   color=[0.0, 0.0, i * 1.0 / len(of)], markeredgecolor='red')
        ax[1].plot(output[4], 0.0, marker='o', markersize=20.0 - 10 * (len(of) - i) / len(of),
                   color=[0.0, 0.0, i * 1.0 / len(of)], markeredgecolor='red')
    else:
        ax[0].plot(L1[i], 0.0, marker='o', markersize=20.0 - 10 * (len(of) - i) / len(of),
                   color=[0.0, 0.0, i * 1.0 / len(of)], ls='None')
        ax[1].plot(KS[i], 0.0, marker='o', markersize=20.0 - 10 * (len(of) - i) / len(of),
                   color=[0.0, 0.0, i * 1.0 / len(of)], ls='None')


#ax.set_xlim(1890, 2010)
# Customizing the plot:

ax[0].get_yaxis().set_visible(False)
ax[0].tick_params(axis="x", labelsize=10)
ax[0].set_xlabel('$\Delta L$ value', fontsize=12)
ax[1].get_yaxis().set_visible(False)
ax[1].tick_params(axis="x", labelsize=10)
ax[1].set_xlabel('$KS$ value', fontsize=12)
fig.suptitle('Restrictions evolution, from black to blue', fontsize=20)
plt.show()

fig, ax = plt.subplots(1, 2)
black = Color("black")
for i in range(-1, len(of)+1):
    if i == -1:
        ax[0].plot(y, a00, marker='d', markersize=20.0,
                   color=[0.0, 0.0, 0.0], markeredgecolor='red')
        ax[1].plot(y, t00, marker='d', markersize=20.0,
                   color=[0.0, 0.0, 0.0], markeredgecolor='red')
    elif i == len(of):
        ax[0].plot(y, a0, marker='o', markersize=20.0 - 10 * (len(of) - i) / len(of),
                   color=[0.0, 0.0, i * 1.0 / len(of)], markeredgecolor='red')
        ax[1].plot(y, t0, marker='o', markersize=20.0 - 10 * (len(of) - i) / len(of),
                   color=[0.0, 0.0, i * 1.0 / len(of)], markeredgecolor='red')
    else:
        ax[0].plot(y, ah[i], marker='o', markersize=20.0 - 10 * (len(of) - i) / len(of),
                   color=[0.0, 0.0, i * 1.0 / len(of)], ls='None')
        ax[1].plot(y, th[i], marker='o', markersize=20.0 - 10 * (len(of) - i) / len(of),
                   color=[0.0, 0.0, i * 1.0 / len(of)], ls='None')



#ax.set_xlim(1890, 2010)
# Customizing the plot:

ax[0].tick_params(axis="y", labelsize=10)
ax[0].tick_params(axis="x", labelsize=10)
ax[1].set_ylabel('$t$ value', fontsize=12)
ax[0].set_xlabel('Wingspan', fontsize=12)

ax[1].tick_params(axis="y", labelsize=10)
ax[1].tick_params(axis="x", labelsize=10)
ax[0].set_ylabel('$a_0$ value', fontsize=12)
ax[1].set_xlabel('Wingspan', fontsize=12)
fig.suptitle('Statevars evolution, from black to blue', fontsize=20)
plt.show()