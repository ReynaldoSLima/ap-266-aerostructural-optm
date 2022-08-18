'''
Homework 3, item 2.1
Author: Reynaldo S. Lima
AP-266
'''

# IMPORTS
import numpy as np
from llt_module import llt_module as llt
from llt_module_b import llt_module_b as llt_b
from scipy.optimize import root

# LLT_MODULE EXECUTION
X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order = 'F')
Gama = np.array([1.0, 2.0, 2.0, 1.0])
alpha0 = np.zeros(4)
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
alpha = 5.498*np.pi/180
Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
rho = 1.0
print('---- LLT result ----')
[res_test, Sref, CL, CD, L, D] = llt.tapenade_main(X, Gama, alpha0, chords, cl0, cla, Vinf, rho)
print('res = ' , res_test)
print('Sref = ', Sref)
print('CL = ', CL)
print('CD = ', CD)
print('L = ', L)
print('D = ', D)

# LLT_MODULE_B EXECUTION
print('---- Reverse Diff. ----')
X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order = 'F')
Gama = np.array([1.0, 2.0, 2.0, 1.0])
alpha0 = np.zeros(4)
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
alpha = 5.498*np.pi/180
Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
rho = 1.0
resb = np.array([0.1, 0.2, 0.3, 0.5])
Srefb = -0.3
CLb = 1.3
CDb = -0.5
Lb = 3.0
Db = 0.02
res = np.zeros(4)
Xb = np.zeros((3,5), order = 'F')
Gamab = np.zeros(4)
alpha0b = np.zeros(4)
chordsb = np.zeros(4)
Sref = 0.0
CL = 0.0
CD = 0.0
L = 0.0
D = 0.0

output = llt_b.tapenade_main_b(X, Xb, Gama, Gamab, alpha0, alpha0b, chords, chordsb, cl0, cla, Vinf, rho, res, resb, Sref, Srefb, CL, CLb, CD, CDb, L, Lb, D, Db)
print('Xb = ...')
print(Xb)
print('Gamab = ', Gamab)
print('alpha0b = ', alpha0b)
print('chordsb = ', chordsb)

# LLT_MODULE_B EXECUTION: REQUESTED OUTPUT
print('---- Reverse Diff., requested output ----')
X = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order = 'F')
Gama = np.array([1.0, 2.0, 2.0, 1.0])
alpha0 = np.zeros(4)
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
alpha = 5.498*np.pi/180
Vinf = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
rho = 1.0
resb = np.array([5.0, -0.3, 0.0, 1.0])
Srefb = -0.3
CLb = 1.3
CDb = -0.5
Lb = 3.0
Db = 0.02
res = np.zeros(4)
Xb = np.zeros((3,5), order = 'F')
Gamab = np.zeros(4)
alpha0b = np.zeros(4)
chordsb = np.zeros(4)
Sref = 0.0
CL = 0.0
CD = 0.0
L = 0.0
D = 0.0

output = llt_b.tapenade_main_b(X, Xb, Gama, Gamab, alpha0, alpha0b, chords, chordsb, cl0, cla, Vinf, rho, res, resb, Sref, Srefb, CL, CLb, CD, CDb, L, Lb, D, Db)
print('Xb = ...')
print(Xb)
print('Gamab = ', Gamab)
print('alpha0b = ', alpha0b)
print('chordsb = ', chordsb)