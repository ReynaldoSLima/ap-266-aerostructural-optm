'''
Homework 3, item 2.3
Author: Reynaldo S. Lima
AP-266, references to class06_session.pdf
'''

## IMPORTS
import numpy as np
from llt_module import llt_module as llt
from llt_module_b import llt_module_b as llt_b
from scipy.optimize import root, minimize

## EXECUTION

# Define number of vortices
n_vortex = 4

# Define Gama
G0 = np.array([1.0, 2.0, 2.0, 1.0])

# Define alpha0
a0 = np.zeros(4)

# Define other variables
X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
alpha = 5.498 * np.pi / 180
Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
rho = 1.0


## OBJECTIVE FUNCTION
def objfun(a0):
    # Define function that computes residuals for given Gama, then get its CD coef.
    def resfunc(G):
        # Define other variables
        X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
        chords = np.array([1.0, 1.0, 1.0, 1.0])
        cl0 = np.zeros(4)
        cla = np.array([6.283, 6.283, 6.283, 6.283])
        alpha = 5.498 * np.pi / 180
        Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
        rho = 1.0
        [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G, a0, chords, cl0, cla, Vinf, rho)

        return res

    # Solve the LLT physical problem
    sol = root(resfunc, G0)

    # Get solved Gama
    G = sol.x

    # Define other variables
    X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
    chords = np.array([1.0, 1.0, 1.0, 1.0])
    cl0 = np.zeros(4)
    cla = np.array([6.283, 6.283, 6.283, 6.283])
    alpha = 5.498 * np.pi / 180
    Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
    rho = 1.0

    # Compute CD
    [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G, a0, chords, cl0, cla, Vinf, rho)
    return CD


def objfungrad(a0):
    # Define function that computes residuals for given Gama, then get its CD coef.
    def resfunc(G):
        # Define other variables
        X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
        chords = np.array([1.0, 1.0, 1.0, 1.0])
        cl0 = np.zeros(4)
        cla = np.array([6.283, 6.283, 6.283, 6.283])
        alpha = 5.498 * np.pi / 180
        Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
        rho = 1.0
        [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G, a0, chords, cl0, cla, Vinf, rho)

        return res

    # Solve the LLT physical problem
    sol = root(resfunc, G0)

    # Get solved Gama
    G = sol.x

    # Initialize output variables for reverse calls
    Xb = np.zeros((3, 5), order='F')
    Gb = np.zeros(4)
    a0b = np.zeros(4)
    chordsb = np.zeros(4)
    res = np.zeros(4)
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
        Xb = np.zeros((3, 5), order='F')
        Gb = np.zeros(4)
        a0b = np.zeros(4)
        chordsb = np.zeros(4)
        res = np.zeros(4)
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
    psi0 = (10 ** 6) * np.ones(4)
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
        # Define other variables
        X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
        chords = np.array([1.0, 1.0, 1.0, 1.0])
        cl0 = np.zeros(4)
        cla = np.array([6.283, 6.283, 6.283, 6.283])
        alpha = 5.498 * np.pi / 180
        Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
        rho = 1.0
        [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G, a0, chords, cl0, cla, Vinf, rho)

        return res

    # Solve the LLT physical problem
    sol = root(resfunc, G0)

    # Get solved Gama
    G = sol.x

    # Define other variables
    X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
    chords = np.array([1.0, 1.0, 1.0, 1.0])
    cl0 = np.zeros(4)
    cla = np.array([6.283, 6.283, 6.283, 6.283])
    alpha = 5.498 * np.pi / 180
    Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
    rho = 1.0

    # Compute CD
    [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G, a0, chords, cl0, cla, Vinf, rho)

    con = CL - 0.5

    return con


def confungrad(a0):
    # Define function that computes residuals for given Gama, then get its CL coef.
    def resfunc(G):
        # Define other variables
        X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
        chords = np.array([1.0, 1.0, 1.0, 1.0])
        cl0 = np.zeros(4)
        cla = np.array([6.283, 6.283, 6.283, 6.283])
        alpha = 5.498 * np.pi / 180
        Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
        rho = 1.0
        [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G, a0, chords, cl0, cla, Vinf, rho)

        return res

    # Solve the LLT physical problem
    sol = root(resfunc, G0)

    # Get solved Gama
    G = sol.x

    # Define other variables
    X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
    chords = np.array([1.0, 1.0, 1.0, 1.0])
    cl0 = np.zeros(4)
    cla = np.array([6.283, 6.283, 6.283, 6.283])
    alpha = 5.498 * np.pi / 180
    Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
    rho = 1.0

    # Compute CD
    [res, Sref, CL, CD, L, D] = llt.tapenade_main(X, G, a0, chords, cl0, cla, Vinf, rho)

    # Initialize output variables for reverse calls
    Xb = np.zeros((3, 5), order='F')
    Gb = np.zeros(4)
    a0b = np.zeros(4)
    chordsb = np.zeros(4)
    res = np.zeros(4)
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
        Xb = np.zeros((3, 5), order='F')
        Gb = np.zeros(4)
        a0b = np.zeros(4)
        chordsb = np.zeros(4)
        res = np.zeros(4)
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
    psi0 = (10 ** 6) * np.ones(4)
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
bounds = [[None, None]]*n_vortex

# Run optimizer
result = minimize(objfun, a0, jac=objfungrad,
                  constraints=[con1],
                  bounds=bounds,
                  method='slsqp')

# Compute function values at the optimum
a0 = result.x
CD = objfun(a0)
CL = confun(a0) + 0.5

print(result)
print(' ')
print('--------------------------------------')
print(' ')
print('--------------------------------------')
print('a0 = ', a0)
print('CD = ', CD)
print('CL = ', CL)
