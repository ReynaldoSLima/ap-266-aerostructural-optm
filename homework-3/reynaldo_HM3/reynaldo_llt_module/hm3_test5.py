'''
Homework 3, item 2.5
Author: Reynaldo S. Lima
AP-266, references to class06_session.pdf
'''

## IMPORTS
import numpy as np
from llt_module import llt_module as llt
from llt_module_b import llt_module_b as llt_b
from scipy.optimize import root, minimize
import matplotlib.pyplot as plt
import matplotlib


## EXECUTION

def output(n_vortex, max):
    # Define Gama
    n = int(n_vortex / 2)
    nf = max
    G0 = np.concatenate([np.linspace(1.0, nf, n), np.linspace(nf, 1.0, n)])

    # Define alpha0
    a00 = np.zeros(n_vortex)

    # Define other variables
    X = np.array([np.zeros(n_vortex + 1), np.linspace(0.0, 8.0, num=n_vortex + 1), np.zeros(n_vortex + 1)], order='F')
    chords = np.ones(n_vortex)
    cl0 = np.zeros(n_vortex)
    cla = 2 * np.pi * np.ones(n_vortex)
    alpha = 5.498 * np.pi / 180
    Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
    rho = 1.0

    def objfun(a0):
        # Define function that computes residuals for given Gama, then get its CD coef.
        def resfunc(G):
            # Define other variables
            [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G.copy(), a0, chords, cl0, cla, Vinf, rho)

            return res

        # Solve the LLT physical problem
        sol = root(resfunc, G0)

        # Get solved Gama
        G = sol.x

        # Compute CD
        [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G.copy(), a0, chords, cl0, cla, Vinf, rho)
        return CD

    def objfungrad(a0):
        # Define function that computes residuals for given Gama, then get its CD coef.
        def resfunc(G):
            [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G.copy(), a0, chords, cl0, cla, Vinf, rho)

            return res

        # Solve the LLT physical problem
        sol = root(resfunc, G0)

        # Get solved Gama
        G = sol.x

        # Initialize output variables for reverse calls
        Xb = np.zeros((3, n_vortex + 1), order='F')
        Gb = np.zeros(n_vortex)
        a0b = np.zeros(n_vortex)
        chordsb = np.zeros(n_vortex)
        res = np.zeros(n_vortex)
        Sref = 0.0
        CL = 0.0
        CD = 0.0
        L = 0.0
        D = 0.0

        # Define function to compute the residual of the adjoint equation
        def adjfunc(resb):
            # Define seeds
            Srefb = 0.0
            CLb = 0.0
            CDb = 1.0  # Thus, we are calculating d()/dCD
            Lb = 0.0
            Db = 0.0
            # Initialize output variables for reverse calls
            Xb = np.zeros((3, n_vortex + 1), order='F')
            Gb = np.zeros(n_vortex)
            a0b = np.zeros(n_vortex)
            chordsb = np.zeros(n_vortex)
            res = np.zeros(n_vortex)
            Sref = 0.0
            CL = 0.0
            CD = 0.0
            L = 0.0
            D = 0.0
            # Run reverse mode
            llt_b.tapenade_main_b(X, Xb, G, Gb, a0, a0b, chords, chordsb, cl0, cla, Vinf,
                                  rho, res, resb.copy(), Sref, Srefb, CL, CLb, CD, CDb, L, Lb, D, Db)

            # Compute the residual of the adjoint equation
            adj_res = Gb

            return adj_res

        # Solve the adjoint problem
        psi0 = (10 ** 7) * np.ones(n_vortex)
        sol = root(adjfunc, psi0)
        psi = sol.x
        resb = psi

        # Initialize derivative seeds
        Srefb = 0.0
        CLb = 0.0
        CDb = 1.0  # Thus, we are calculating dCD/da0
        Lb = 0.0
        Db = 0.0

        # Run reverse mode
        llt_b.tapenade_main_b(X, Xb, G, Gb, a0, a0b, chords, chordsb, cl0, cla, Vinf,
                              rho, res, resb.copy(), Sref, Srefb, CL, CLb, CD, CDb, L, Lb, D, Db)
        # Get total derivative
        dCDda0 = a0b
        return dCDda0

    def confun(a0):
        # Define function that computes residuals for given Gama, then get its CL coef.
        def resfunc(G):
            [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G.copy(), a0, chords, cl0, cla, Vinf, rho)

            return res

        # Solve the LLT physical problem
        sol = root(resfunc, G0)

        # Get solved Gama
        G = sol.x


        # Compute CD
        [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G, a0, chords, cl0, cla, Vinf, rho)

        con = CL - 0.5

        return con

    def confungrad(a0):
        # Define function that computes residuals for given Gama, then get its CL coef.
        def resfunc(G):
            [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G.copy(), a0, chords, cl0, cla, Vinf, rho)

            return res

        # Solve the LLT physical problem
        sol = root(resfunc, G0)

        # Get solved Gama
        G = sol.x

        # Compute CD
        [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G, a0, chords, cl0, cla, Vinf, rho)

        # Initialize output variables for reverse calls
        Xb = np.zeros((3, n_vortex + 1), order='F')
        Gb = np.zeros(n_vortex)
        a0b = np.zeros(n_vortex)
        chordsb = np.zeros(n_vortex)
        res = np.zeros(n_vortex)
        Sref = 0.0
        CL = 0.0
        CD = 0.0
        L = 0.0
        D = 0.0

        # Define function to compute the residual of the adjoint equation
        def adjfunc(resb):
            # Define seeds
            Srefb = 0.0
            CLb = 1.0  # Thus, we are calculating d()/dCL
            CDb = 0.0
            Lb = 0.0
            Db = 0.0
            # Initialize output variables for reverse calls
            Xb = np.zeros((3, n_vortex + 1), order='F')
            Gb = np.zeros(n_vortex)
            a0b = np.zeros(n_vortex)
            chordsb = np.zeros(n_vortex)
            res = np.zeros(n_vortex)
            Sref = 0.0
            CL = 0.0
            CD = 0.0
            L = 0.0
            D = 0.0
            # Run reverse mode
            llt_b.tapenade_main_b(X, Xb, G, Gb, a0, a0b, chords, chordsb, cl0, cla, Vinf,
                                  rho, res, resb.copy(), Sref, Srefb, CL, CLb, CD, CDb, L, Lb, D, Db)

            # Compute the residual of the adjoint equation
            adj_res = Gb

            return adj_res

        # Solve the adjoint problem
        psi0 = (10 ** 7) * np.ones(n_vortex)
        sol = root(adjfunc, psi0)
        psi = sol.x
        resb = psi
        # Initialize derivative seeds
        Srefb = 0.0
        CLb = 1.0  # Thus, we are calculating dCL/da0
        CDb = 0.0
        Lb = 0.0
        Db = 0.0

        # Run reverse mode
        llt_b.tapenade_main_b(X, Xb, G, Gb, a0, a0b, chords, chordsb, cl0, cla, Vinf,
                              rho, res, resb.copy(), Sref, Srefb, CL, CLb, CD, CDb, L, Lb, D, Db)
        # Get total derivative
        dCLda0 = a0b

        congrad = dCLda0
        return congrad

    ## OPTIMIZATION

    # Create list of constraints

    con1 = {'type': 'eq',
            'fun': confun,
            'jac': confungrad}

    # Set lower bounds for areas
    bounds = [[None, None]] * n_vortex
    # Run optimizer
    result = minimize(objfun, a00, jac=objfungrad,
                      constraints=[con1],
                      bounds=bounds,
                      method='slsqp')

    # Compute function values at the optimum
    a0 = result.x
    out2 = result.nfev
    CD = objfun(a0)
    CL = confun(a0) + 0.5
    '''
    print(' ')
    print('--------------------------------------')
    print(' ')
    print('--------------------------------------')
    print(result)
    '''
    print(' ')
    print('--------------------------------------')
    print(' ')
    print('--------------------------------------')
    print('a0 = ', a0)
    print('CD = ', CD)
    print('CL = ', CL)

    def resfunc_aux(G):
        [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G.copy(), a0, chords, cl0, cla, Vinf, rho)
        return res

    # Solve the LLT physical problem

    sol = root(resfunc_aux, G0)

    # Get solved Gama
    out = sol.x


    return out, out2, CD


CD_ELIPTICAL = 8.0 * 0.25 / (3.1415 * 64.0)
out_4, nfev_4, CD_4 = output(4, 1.0)
out_8, nfev_8, CD_8 = output(8, 1.0)
out_16, nfev_16, CD_16 = output(16, 1.0)
