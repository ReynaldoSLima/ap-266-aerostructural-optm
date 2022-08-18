'''
Homework 3, item 2.2
Author: Reynaldo S. Lima
AP-266
'''

# IMPORTS
import numpy as np
from llt_module import llt_module as llt
from llt_module_b import llt_module_b as llt_b
from scipy.optimize import root

# GAMMA FUNCTION
def resfunc(Gamma):
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

sol = root(resfunc, Gama)
print('Gama = ', sol.x)
print('------------------------------------------------------------------------------')
print(' ')
print('----------------------------------- outputs ----------------------------------')
Gama = sol.x
X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order = 'F')
alpha0 = np.zeros(4)
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
alpha = 5.498*np.pi/180
Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
rho = 1.0
[res, Sref, CL, CD, L, D] = llt.tapenade_main(X, Gama, alpha0, chords, cl0, cla, Vinf, rho)
print('res = ' , res)
print('Sref = ', Sref)
print('CL = ', CL)
print('CD = ', CD)
print('L = ', L)
print('D = ', D)