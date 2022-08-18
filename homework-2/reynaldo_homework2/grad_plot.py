""" Homework 2 : AP-266
    2.b)       : auxiliar plot file
    Author     : Reynaldo S. Lima """

import numpy as np
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
import matplotlib
plt.style.use('fivethirtyeight')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
nm = np.loadtxt('nelder-mead', np.dtype('d'))
ev = np.loadtxt('dif-evol', np.dtype('d'))
cg = np.loadtxt('grad-conj', np.dtype('d'))
bfgs = np.loadtxt('bfgs', np.dtype('d'))
"""
From the format files, each vector[:,0] gathers nfev (or nfev + njev) and vector[:,1] function value
"""

"""
First, assigning evi as a intermediate variable:
"""

ev1 = nm[0, :]
ev2 = ev[0, :]
ev3 = cg[0, :]
ev4 = bfgs[0, :]
x = range(0,21,2)
pr1 = ev1[::2]
pr2 = ev2[::2]
pr3 = ev3[::2]
pr4 = ev4[::2]

""""
Ploting (x, pri)
"""
fig, ax = plt.subplots()
ax.plot(x, pr1, label='Nelder-Mead', marker='+',ls='--', linewidth=2.0)
ax.plot(x, pr2, label='Evolutionary',marker=',',ls='-', linewidth=2.0)
ax.plot(x, pr3, label='Conjugate Gradient',marker='.',ls=':', linewidth=2.0)
ax.plot(x, pr4, label='BFGS',marker='>',ls='-.', linewidth=2.0)
ax.tick_params(axis="y", labelsize=20)
ax.tick_params(axis="x", labelsize=20)
plt.xlabel('Number of entries', fontsize=30)
plt.ylabel('Number of function evaluations', fontsize=30)
ax.tick_params(axis = 'x', width=5)
legend = ax.legend(loc='best')
ax.set_facecolor('white')
plt.yscale('log')
ax.patch.set_facecolor('white')
plt.xticks(np.arange(0, 20,2.0))
plt.title('Comparison of optimization algorithms - Log scale', fontsize=40)
plt.show()
