'''
Homework, item 5.2
Author: Reynaldo S. Lima
AP-266
'''

# IMPORTS
import numpy as np
from asa_module import asa_module as asa

Xfem = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
Gama = np.array([80.0, 80.0, 80.0, 80.0])
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
d = 1E-3 * np.ones(10)
CD0 = 0.0270
fixedMass = 700.0
g = 9.8
Endurance = 4.0 * 60.0 * 60.0
TSFC = 0.5 / 3600.0
loadFactor = 3.0 * 1.5

resllt, resfem, liftExcess, margins, KSmargin, FB, Weight, Sref, CL = asa.asa_analysis(Gama, alpha0, chords,
                                                                                       Xfem, R, t, d,
                                                                                       cl0, cla, Vinf, rhoAir,
                                                                                       E, rhoMat, sigmaY, pKS, conIDs,
                                                                                       CD0, fixedMass, g, Endurance,
                                                                                       TSFC, loadFactor)

print('---------------------------------')
print(' ')
print('resllt = ', resllt)
print(' ')
print('resfem = ', resfem)
print(' ')
print('liftExcess = ', liftExcess)
print(' ')
print('margins = ', margins)
print(' ')
print('KSmargin = ', KSmargin)
print(' ')
print('FB = ', FB)
print(' ')
print('Weight = ', Weight)
print(' ')
print('Sref = ', Sref)
print(' ')
print('CL = ', CL)
print(' ')
print('---------------------------------')