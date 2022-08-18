""" Homework 2 : AP-266
    2.b)
    Author     : Reynaldo S. Lima """

import numpy as np
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

'''
Let's define functions for any n
'''


def rosenbrock(x):
    '''
    :type x: object
    :param x     : initial value
    :return      : Rosenbrock test function
    '''
    n = len(x)
    r = 0
    for i in range(0, n - 1):
        r = r + 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2
    return r


'''
------------------- Nelder-Mead -------------------------
'''
nm = np.zeros(21)
funm = np.zeros(21)
options = dict(maxiter=200000, fatol=1e-12, adaptive=True)
for i in range(1, 11):
    # Define initial guess
    x0 = np.zeros(2 * i)
    # Run optimization
    result = minimize(rosenbrock, x0, method='Nelder-Mead', options=options)
    nm[2 * i] = result.nfev
    funm[2 * i] = rosenbrock(result.x)

'''
------------------- Dif. Evol -------------------------
'''
ev = np.zeros(21)
fune = np.zeros(21)

bounds = np.array(
    [[-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4],
     [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4]])

for i in range(1, 11):
    # Define bounds
    # Run optimization
    result = differential_evolution(rosenbrock, bounds[0:2 * i, :], seed=1, maxiter=5000, polish=False, atol=1e-12)
    ev[2 * i] = result.nfev
    fune[2 * i] = rosenbrock(result.x)

'''
----------------------- CG -------------------------
'''
cg = np.zeros(21)
fcg = np.zeros(21)
options = dict(maxiter=200000)
for i in range(1, 11):
    # Define initial guess
    x0 = np.zeros(2 * i)
    # Run optimization
    result = minimize(rosenbrock, x0, method='CG', tol=1e-6, options=options)
    cg[2 * i] = result.nfev + result.njev
    fcg[2 * i] = rosenbrock(result.x)

'''
--------------------- BFGS -------------------------
'''
bfgs = np.zeros(21)
fbfgs = np.zeros(21)
options = dict(maxiter=200000)
for i in range(1, 11):
    # Define initial guess
    x0 = np.zeros(2 * i)
    # Run optimization
    result = minimize(rosenbrock, x0, method='BFGS', tol=1e-6, options=options)
    bfgs[2 * i] = result.nfev + result.njev
    fbfgs[2 * i] = rosenbrock(result.x)

np.savetxt("nelder-mead", (nm, funm))
np.savetxt("dif-evol", (ev, fune))
np.savetxt("grad-conj", (cg, fcg))
np.savetxt("bfgs", (bfgs, fbfgs))
