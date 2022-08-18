'''
Homework 3, item 2.5 and 2.6
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
out_40, nfev_40, CD_40 = output(40, 1.0)

# For printing
out_12, nfev_12, CD_12 = output(12, 1.0)
out_20, nfev_20, CD_20 = output(20, 1.0)
out_24, nfev_24, CD_24 = output(24, 1.0)
out_28, nfev_28, CD_28 = output(28, 1.0)
out_32, nfev_32, CD_32 = output(32, 1.0)
out_36, nfev_36, CD_36 = output(36, 1.0)
out_2, nfev_2, CD_2 = output(2, 1.0)
out_6, nfev_6, CD_6 = output(6, 1.0)
out_10, nfev_10, CD_10 = output(10, 1.0)
out_14, nfev_14, CD_14 = output(14, 1.0)
out_18, nfev_18, CD_18 = output(18, 1.0)
out_22, nfev_22, CD_22 = output(22, 1.0)
out_26, nfev_26, CD_26 = output(26, 1.0)
out_30, nfev_30, CD_30 = output(30, 1.0)
out_34, nfev_34, CD_34 = output(34, 1.0)
out_38, nfev_38, CD_38 = output(38, 1.0)

# Plots
y = np.linspace(0.0, 8.0, num=100)
y_4 = np.linspace(0.0, 8.0, num=4)
y_8 = np.linspace(0.0, 8.0, num=8)
y_16 = np.linspace(0.0, 8.0, num=16)
y_40 = np.linspace(0.0, 8.0, num=40)
gama0 = 2 * 8.0 * 1.0 * 0.5 / (3.1415 * 8.0)
gama = 2 * gama0 * np.sqrt(y / 8.0 - (y ** 2) / 64.0)

fig, ax = plt.subplots()
ax.plot(y, gama, label='Elliptical', ls='--', linewidth=1.5)
ax.plot(y_4, out_4, label='4 Vortices', linewidth=2.0)
ax.plot(y_8, out_8, label='8 Vortices', linewidth=2.0)
ax.plot(y_16, out_16, label='16 Vortices', linewidth=2.0)
ax.plot(y_40, out_40, label='40 Vortices', linewidth=2.0)
ax.tick_params(axis="y", labelsize=10)
ax.tick_params(axis="x", labelsize=10)
plt.xlabel('Span', fontsize=20)
plt.ylabel('Circulation', fontsize=20)
ax.tick_params(axis='x', width=5)
legend = ax.legend(loc='best')
ax.set_facecolor('white')
ax.patch.set_facecolor('white')
plt.xticks(np.arange(-1, 10, 2.0))
plt.title('Circulation distribution', fontsize=20)
plt.show()

x = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40])
z = np.array([CD_2, CD_4, CD_6, CD_8, CD_10, CD_12, CD_14, CD_16, CD_18, CD_20, CD_22, CD_24, CD_26, CD_28, CD_30, CD_32
                 , CD_34, CD_36, CD_38, CD_40])
#x = np.array([4, 8, 16, 40])
#z = np.array([CD_4, CD_8, CD_16, CD_40])
fig, ax = plt.subplots()
ax.plot(np.linspace(0.0, 40.0, num=100), CD_ELIPTICAL * np.ones(100), label='Elliptical', ls='--', linewidth=1.5)
ax.plot(x, z, marker='o', ls='None')
ax.tick_params(axis="y", labelsize=10)
ax.tick_params(axis="x", labelsize=10)
plt.xlabel('N_Vortex', fontsize=20)
plt.ylabel('CD', fontsize=20)
ax.tick_params(axis='x', width=5)
ax.set_facecolor('white')
ax.patch.set_facecolor('white')
plt.xticks(np.arange(0, 40, 2.0))
plt.title('CD deppending on Number of var.', fontsize=20)
plt.show()

x = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40])
z = np.array(
    [nfev_2, nfev_4, nfev_6, nfev_8, nfev_10, nfev_12, nfev_14, nfev_16, nfev_18, nfev_20, nfev_22, nfev_24, nfev_26
        , nfev_28, nfev_30, nfev_32, nfev_34, nfev_36, nfev_38, nfev_40])
fig, ax = plt.subplots()
ax.plot(x, z, marker='o', ls='None')
ax.tick_params(axis="y", labelsize=10)
ax.tick_params(axis="x", labelsize=10)
plt.xlabel('N_Vortex', fontsize=20)
plt.ylabel('NFEV', fontsize=20)
ax.tick_params(axis='x', width=5)
ax.set_facecolor('white')
ax.patch.set_facecolor('white')
plt.xticks(np.arange(0, 40, 2.0))
plt.title('Number of function evaluations', fontsize=20)
plt.show()
