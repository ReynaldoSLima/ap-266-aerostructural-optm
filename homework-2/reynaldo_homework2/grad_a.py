""" Homework 2 : AP-266
    2.a)
    Author     : Reynaldo S. Lima """

# Import section
import numpy as np
from scipy.optimize import minimize

'''
Let's define functions for n = 2, 4 and 6
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
        r = r + 100*(x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2
    return r


CG = minimize(rosenbrock, np.zeros(2), method='CG')
BF = minimize(rosenbrock, np.zeros(2), method='BFGS')
print('For n = 2')
print(CG, BF)
CG = minimize(rosenbrock, np.zeros(4), method='CG')
BF = minimize(rosenbrock, np.zeros(4), method='BFGS')
print('For n = 4')
print(CG, BF)
CG = minimize(rosenbrock, np.zeros(6), method='CG')
BF = minimize(rosenbrock, np.zeros(6), method='BFGS')
print('For n = 6')
print(CG, BF)