'''
Homework, item 5.3
Author: Reynaldo S. Lima
AP-266
'''

# IMPORTS
import numpy as np
from asa_module import asa_module as asa
from asa_module_d import asa_module_d as asad
from asa_module_b import asa_module_b as asab

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

# FINITE DIFFERENCES TEST:

# SEEDS!!!!!
alpha0d = np.array([0.0, 0.1, -0.1, 0.2])
Gamad = np.array([-0.0, 0.5, 0.3, 0.2])
dd = 1E-4 * np.ones(10)
td = np.array([0.003, 0.002, -0.001, 0.001])
# H VALUE
h = 1E-7

resllt, resfem, liftExcess, margins, KSmargin, FB, Weight, Sref, CL = asa.asa_analysis(Gama, alpha0, chords,
                                                                                       Xfem, R, t, d,
                                                                                       cl0, cla, Vinf, rhoAir,
                                                                                       E, rhoMat, sigmaY, pKS, conIDs,
                                                                                       CD0, fixedMass, g, Endurance,
                                                                                       TSFC, loadFactor)

resllt1, resfem1, liftExcess1, margins1, KSmargin1, FB1, Weight1, Sref1, CL1 = asa.asa_analysis(Gama + h * Gamad, alpha0
                                                                                                + h * alpha0d, chords,
                                                                                                Xfem, R, t + h * td,
                                                                                                d + dd * h,
                                                                                                cl0, cla, Vinf, rhoAir,
                                                                                                E, rhoMat, sigmaY, pKS,
                                                                                                conIDs,
                                                                                                CD0, fixedMass, g,
                                                                                                Endurance,
                                                                                                TSFC, loadFactor)

FDllt = (resllt1 - resllt) / h
FDfem = (resfem1 - resfem) / h
FDlif = (liftExcess1 - liftExcess) / h
FDmar = (margins1 - margins) / h
FDKSm = (KSmargin1 - KSmargin) / h
FDFB = (FB1 - FB) / h
FDWe = (Weight1 - Weight) / h
# CALLING THE DIFF. CODE

resllt, reslltd, resfem, resfemd, liftexcess, liftexcessd, margins, marginsd, ksmargin, ksmargind, fb, fbd, weight, \
weightd, sref, cl = asad.asa_analysis_d(Gama, Gamad, alpha0, alpha0d, chords, Xfem, R, t, td, d, dd, cl0, cla, Vinf,
                                        rhoAir, E, rhoMat,
                                        sigmaY, pKS, conIDs, CD0, fixedMass, g, Endurance, TSFC, loadFactor)

out1 = FDllt - reslltd
out2 = FDfem - resfemd
out3 = FDlif - liftexcessd
out4 = FDmar - marginsd
out5 = FDKSm - ksmargind
out6 = FDFB - fbd
out7 = FDWe - weightd

print('------------------------ FD TEST ------------------------')
print('---------------------------------------------------------')
print('FD test results (0.0 expected)')
print('For reslltD we expect about 1.0')
print(' ')
print('reslltD = ', out1)
print(' ')
print('resfemD = ', out2)
print(' ')
print('liftExcessD = ', out3)
print(' ')
print('marginsD = ', out4)
print(' ')
print('KSmarginD = ', out5)
print(' ')
print('FBD = ', out6)
print(' ')
print('WeightD = ', out7)
print(' ')
print('---------------------------------------------------------')

# LET'S START THE SECOND TEST!!!
# DOT PRODUCT TEST

Xfem = np.array([[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 6.0, 8.0], [0.0, 0.0, 0.0, 0.0, 0.0]], order='F')
Gama = np.array([80.0, 80.0, 80.0, 80.0])
Gamab = np.zeros(4)
alpha0 = np.zeros(4)
alpha0b = np.zeros(4)
chords = np.array([1.0, 1.0, 1.0, 1.0])
cl0 = np.zeros(4)
cla = np.array([6.283, 6.283, 6.283, 6.283])
Vinfm = 60
alpha = 5 * np.pi / 180
Vinf = np.array([Vinfm * np.cos(alpha), 0.0, Vinfm * np.sin(alpha)])
rhoAir = 1.225
R = np.array([0.1, 0.1, 0.1, 0.1])
t = np.array([0.005, 0.005, 0.005, 0.005])
tb = np.zeros(4)
E = 73.1E9
rhoMat = 2780.0
sigmaY = 324.0E6
pKS = 200.0
conIDs = np.array([5, 6], dtype='int32')
d = 1E-3 * np.ones(10)
db = np.zeros(10)
CD0 = 0.0270
fixedMass = 700.0
g = 9.8
Endurance = 4.0 * 60.0 * 60.0
TSFC = 0.5 / 3600.0
loadFactor = 3.0 * 1.5

# CALLING THE REVERSE FUNCTION
resllt = np.zeros(4)
resfem = np.zeros(10)
liftexcess = 0.0
margins = np.zeros(8)
ksmargin = 0.0
fb = 0.0
weight = 0.0
sref = 0.0
cl = 0.0

# Seeds of the outputs (random)
reslltb = np.ones(4) * 0.3
resfemb = np.ones(10) * 0.2
liftexcessb = 50.5
marginsb = np.ones(8) * 0.1
ksmarginb = 102.0
fbb = 30.3
weightb = 10.0

output = asab.asa_analysis_b(Gama, Gamab, alpha0, alpha0b, chords, Xfem, R, t, tb, d, db, cl0, cla, Vinf,
                             rhoAir, E, rhoMat,
                             sigmaY, pKS, conIDs, CD0, fixedMass, g, Endurance, TSFC, loadFactor, resllt,
                             reslltb, resfem, resfemb, liftexcess, liftexcessb, margins, marginsb, ksmargin, ksmarginb,
                             fb,
                             fbb, weight, weightb, sref, cl)

# We can evaluate now the dot product
# Seeds of the outputs (random)
reslltb = np.ones(4) * 0.3
resfemb = np.ones(10) * 0.2
liftexcessb = 50.5
marginsb = np.ones(8) * 0.1
ksmarginb = 102.0
fbb = 30.3
weightb = 10.0

dotp = 1 + (- np.dot(reslltd, reslltb) - np.dot(resfemd, resfemb) - liftexcessd * liftexcessb - np.dot(marginsd, marginsb) - ksmargind * ksmarginb - fbd * fbb - weightd * weightb) / (np.dot(Gamad, Gamab) + np.dot(alpha0d, alpha0b) +
                                                                                                                                                                                       np.dot(td, tb) + np.dot(dd, db))

print('------------------------ DP TEST ------------------------')
print('---------------------------------------------------------')
print('Dot product results (0.0 expected)')
print('dotp = ', dotp)
print('---------------------------------------------------------')
