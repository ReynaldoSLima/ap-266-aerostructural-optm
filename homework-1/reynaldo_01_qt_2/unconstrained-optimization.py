# IMPORTS
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# EXECUTION

# Define objective function
def objfun(x):
    # Let's make the program for any value of n, then print for all values desired.
    n = len(x)
    f = 0
    for i in range(0, n - 1):
        f = f + 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2
    return f


# =========== NELDER-MEAD ============#
# nm counts the number of iterations
nm = np.zeros(21)
funm = np.zeros(21)

for i in range(1, 11):
    # Define initial guess
    x0 = np.zeros(2 * i)
    # Run optimization
    result = minimize(objfun, x0, method='Nelder-Mead', options={'maxiter': 1e10})
    if i < 4:
        print('NM. for n = ', 2 * i, result.x, result.message)
    nm[2 * i] = result.nfev
    funm[2 * i] = objfun(result.x)
    if i > 1:
        print(nm[2*(i-1)]/nm[2*i], (i-1)/i)

# =========== EVOL. METH. ===========#
# ev counts the number of iterations
ev = np.zeros(21)
fune = np.zeros(21)

bounds = np.array(
    [[-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4],
     [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4], [-4, 4]])

for i in range(1, 11):
    # Define bounds
    # Run optimization
    result = differential_evolution(objfun, bounds[0:2 * i, :], maxiter=int(1e10))
    # Print results
    if i < 4:
        print('Ev. for n = ', 2 * i, result.x, result.message)
    ev[2 * i] = result.nfev
    fune[2 * i] = objfun(result.x)
    if i > 3:
        print(ev[2*(i-1)] / ev[2*i], (i-1)/i, (i-1)**2/i**2)

# in order to know which order is Evolutionary:



fig, ax = plt.subplots()
ax.plot(nm[::2], label='Nelder-Mead')
ax.plot(ev[::2], label='Evolutionary')
legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
plt.title('Comparison of unconstrained optimization algorithms')
plt.show()

fig2, ax2 = plt.subplots()
ax2.plot(funm[::2], label='Nelder-Mead')
ax2.plot(fune[::2], label='Evolutionary')
legend2 = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
plt.title('Comparison of function values')
plt.show()
