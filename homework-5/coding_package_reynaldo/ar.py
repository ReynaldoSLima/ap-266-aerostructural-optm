'''
Homework, item 5.8
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
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)

nP = 20

def fb_min(AR):
    # Specific to the problem parameters
    nP = 20
    sRef = 16.0
    bRef = np.sqrt(AR * sRef)
    print('FB, AR = ', AR, ', bRef =', bRef)
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

    print('-------------------- OPTIMIZATION -----------------------')
    print('------------------------- FB ----------------------------')
    desVars0 = np.concatenate((np.ones(nP) * 0.0, np.ones(nP) * 0.5))
    np.set_printoptions(precision=4)
    bounds1 = [None, None] * nP
    bounds2 = [0.1, None] * nP
    bounds = np.concatenate([bounds1, bounds2])
    bounds = bounds.reshape((2 * nP, 2))
    # Run optimizer
    result = minimize(objFunc, desVars0, jac=objFuncGrad,
                      constraints=[con1, con2],
                      bounds=bounds,
                      method='slsqp')


    # Final output
    desVars = result.x
    output = solve_asa(desVars)
    a0 = desVars[0:nP]
    t0 = desVars[nP:2 * nP] / 100.0
    Gama = output[9]
    print('AR =', AR)
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

    print('FB, AR = ', AR, ', CL =', output[8])
    return output

def we_min(AR):
    '''
    Homework, item 5.7
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
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)

    # Specific to the problem parameters
    nP = 20
    sRef = 16.0
    bRef = np.sqrt(AR * sRef)
    print('WE, AR = ', AR, ', bRef =', bRef)
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
                                                                                                   cl0, cla, Vinf,
                                                                                                   rhoAir,
                                                                                                   E, rhoMat, sigmaY,
                                                                                                   pKS,
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
                                                                                                   cl0, cla, Vinf,
                                                                                                   rhoAir,
                                                                                                   E, rhoMat, sigmaY,
                                                                                                   pKS,
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

            asab.asa_analysis_b(Gama.copy(), Gamab, alpha0.copy(), alpha0b, chords, Xfem, R, t.copy(), tb, d.copy(), db,
                                cl0
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
        return output[6]

    def objFuncGrad(desVars):
        dFBdVars = compute_grads(desVars, 0.0, 0.0, 0.0, 1.0)
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
    np.set_printoptions(precision=4)
    print('-------------------- OPTIMIZATION -----------------------')
    print('----------------------- WEIGHT --------------------------')
    bounds1 = [None, None] * nP
    bounds2 = [0.1, None] * nP
    bounds = np.concatenate([bounds1, bounds2])
    bounds = bounds.reshape((2 * nP, 2))
    # Run optimizer
    result = minimize(objFunc, desVars0, jac=objFuncGrad,
                      constraints=[con1, con2],
                      bounds=bounds,
                      method='slsqp')

    # Final output
    desVars = result.x
    output = solve_asa(desVars)
    a0 = desVars[0:nP]
    t0 = desVars[nP:2 * nP] / 100.0
    Gama = output[9]
    print('AR =', AR)
    print('---------------------- OBJ FUNC. ------------------------')
    print('W0 =', "%.4f" % output[6])
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
    print('WE, AR = ', AR, ', CL =', output[8])
    return output

fb1 = fb_min(6.0)
we1 = we_min(6.0)
fb2 = fb_min(10.0)
we2 = we_min(10.0)

# Weght ratio analysis:

fixedMass = 700.0
g = 9.8

FW = fixedMass * g
Weight = fb1[6]
FB = fb1[5]
ratioFB = FB/Weight
ratioFM = FW/Weight
ratioSM = 1 - ratioFB - ratioFM
print('-------------------- WEIGHT RATIO -----------------------')
print('---------------------- FB 6.0 ---------------------------')
print('W0 = ', Weight)
print('FB =', "%.2f" % (ratioFB * 100), '%')
print(' ')
print('FM =', "%.2f" % (ratioFM * 100), '%')
print(' ')
print('SM =', "%.2f" % (ratioSM * 100), '%')
print('---------------------------------------------------------')
Weight = fb2[6]
FB = fb2[5]
ratioFB = FB/Weight
ratioFM = FW/Weight
ratioSM = 1 - ratioFB - ratioFM
print('---------------------- FB 10.0 --------------------------')
print('W0 = ', Weight)
print('FB =', "%.2f" % (ratioFB * 100), '%')
print(' ')
print('FM =', "%.2f" % (ratioFM * 100), '%')
print(' ')
print('SM =', "%.2f" % (ratioSM * 100), '%')
print('---------------------------------------------------------')
Weight = we1[6]
FB = we1[5]
ratioFB = FB/Weight
ratioFM = FW/Weight
ratioSM = 1 - ratioFB - ratioFM
print('---------------------- WE 6.0 ---------------------------')
print('W0 = ', Weight)
print('FB =', "%.2f" % (ratioFB * 100), '%')
print(' ')
print('FM =', "%.2f" % (ratioFM * 100), '%')
print(' ')
print('SM =', "%.2f" % (ratioSM * 100), '%')
print('---------------------------------------------------------')
Weight = we2[6]
FB = we2[5]
ratioFB = FB/Weight
ratioFM = FW/Weight
ratioSM = 1 - ratioFB - ratioFM
print('---------------------- WE 10.0 --------------------------')
print('W0 = ', Weight)
print('FB =', "%.2f" % (ratioFB * 100), '%')
print(' ')
print('FM =', "%.2f" % (ratioFM * 100), '%')
print(' ')
print('SM =', "%.2f" % (ratioSM * 100), '%')
print('---------------------------------------------------------')


# Elliptical circulation
bRef = 9.797958971132712
bRef2 = 12.649110640673518
dx = bRef2 - bRef
dx = dx / 2
Vinfm = 60.0
sRef = 16.0
y = np.linspace(0.0, bRef, num=nP)
y2 = np.linspace(0.0, bRef2, num=nP)
CL = fb1[8]
gama0 = 2 * Vinfm * sRef * CL / (np.pi * bRef)
Gamae1 = 2 * gama0 * np.sqrt(y / bRef - (y ** 2) / bRef ** 2)

CL = fb2[8]
gama0 = 2 * Vinfm * sRef * CL / (np.pi * bRef2)
Gamae2 = 2 * gama0 * np.sqrt(y2 / bRef2 - (y2 ** 2) / bRef2 ** 2)

CL = we1[8]
gama0 = 2 * Vinfm * sRef * CL / (np.pi * bRef)
Gamae3 = 2 * gama0 * np.sqrt(y / bRef - (y ** 2) / bRef ** 2)

CL = we2[8]
gama0 = 2 * Vinfm * sRef * CL / (np.pi * bRef2)
Gamae4 = 2 * gama0 * np.sqrt(y2 / bRef2 - (y2 ** 2) / bRef2 ** 2)

# Ploting:
fig, ax = plt.subplots(2, 1)
ax[0].plot(y + dx, Gamae1, label='Elliptical, AR=6.0, FB min.', color='dimgray', ls='--', linewidth=1.0)
ax[1].plot(y2, Gamae2, label='Elliptical, AR=10.0, FB min.', color='dimgrey', ls='--', linewidth=1.0)
ax[0].plot(y + dx, Gamae3, label='Elliptical, AR=6.0, W0 min.', color='silver', ls='--', linewidth=1.0)
ax[1].plot(y2, Gamae4, label='Elliptical, AR=10.0, W0 min.', color='lightgrey', ls='--', linewidth=1.0)
ax[0].plot(y + dx, fb1[9], label='FB Minimization, AR = 6.0', color='darkblue', linewidth=2.0)
ax[1].plot(y2, fb2[9], label='FB Minimization, AR = 10.0', color='mediumblue', linewidth=2.0)
ax[0].plot(y + dx, we1[9], label='Weight Minimization, AR = 6.0', color='deepskyblue', linewidth=2.0)
ax[1].plot(y2, we2[9], label='Weight Minimization, AR = 10.0', color='dodgerblue', linewidth=2.0)

# Customizing the plot:

ax[0].grid()
ax[1].grid()
ax[0].tick_params(axis="y", labelsize=10)
ax[0].tick_params(axis="x", labelsize=10)
ax[1].tick_params(axis="y", labelsize=10)
ax[1].tick_params(axis="x", labelsize=10)
ax[0].set_xlim(0, bRef2)
ax[0].set_ylim(-2.0, 22.0)
ax[1].set_xlim(0, bRef2)
ax[1].set_ylim(-2.0, 22.0)
ax[1].set_xlabel('Wingspan', fontsize=20)
ax[0].set_ylabel('Circulation', fontsize=20)
ax[1].set_ylabel('Circulation', fontsize=20)
legend = ax[0].legend(loc='best')
legend2 = ax[1].legend(loc='best')
ax[0].set_title('Circulation distribution, AR 6.0', fontsize=20)
ax[1].set_title('Circulation distribution, AR 10.0', fontsize=20)
plt.show()