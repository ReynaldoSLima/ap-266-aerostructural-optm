'''
Homework 4, item 1
Autor: Reynaldo S. Lima
AP-266
'''

# IMPORTS
import numpy as np
from scipy.optimize import root, minimize
# Imports de plot
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
# Função de acompanhamento da otimização:
X = []
Volume = []
q100 = []
def callback(x):
    fobj = objfun(x)
    fcon = confun1(x)
    Volume.append(fobj)
    X.append(x)
    q100.append(fcon)

# Função objetivo, volume
def objfun(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    L = 1.0

    # Primeiro vamos buscar a função dos resíduos, para que todo deslocamento obedeça o equilíbrio
    def resfunc(q):
        E = 7.0 * 10 ** 10
        K = 2 / 3 * E / L ** 3 * np.array(
            [[12 * (H1 ** 3 * T1 + H2 ** 3 * T2), 6 * L * (-H1 ** 3 * T1 + H2 ** 3 * T2), 6 * L * H2 ** 3 * T2],
             [6 * L * (-H1 ** 3 * T1 + H2 ** 3 * T2), 4 * L ** 2 * (H1 ** 3 * T1 + H2 ** 3 * T2),
              2 * L ** 2 * H2 ** 3 * T2],
             [6 * L * H2 ** 3 * T2, 2 * L ** 2 * H2 ** 3 * T2, 4 * L ** 2 * H2 ** 3 * T2]])
        F = np.array([10000.0, 0.0, -1833.3])
        Kd = K.dot(q)
        res = Kd - F
        return res

    # Chute inicial para os deslocamentos
    q0 = np.array([[0.0001], [0.0001], [0.0001]])
    # Resolver o problema estático:
    sol = root(resfunc, q0)

    # Escolhendo o vetor deslocamento
    q = sol.x
    # Podemos obter disto o nosso volume
    V = 4 * (T1 * H1 + T2 * H2) * L

    return V


# Gradiente da função objetivo
def objfungrad(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    L = 1.0

    # Primeiro vamos buscar a função dos resíduos, para que todo deslocamento obedeça o equilíbrio
    def resfunc(q):
        E = 7.0 * 10 ** 10
        K = 2 / 3 * E / L ** 3 * np.array(
            [[12 * (H1 ** 3 * T1 + H2 ** 3 * T2), 6 * L * (-H1 ** 3 * T1 + H2 ** 3 * T2), 6 * L * H2 ** 3 * T2],
             [6 * L * (-H1 ** 3 * T1 + H2 ** 3 * T2), 4 * L ** 2 * (H1 ** 3 * T1 + H2 ** 3 * T2),
              2 * L ** 2 * H2 ** 3 * T2],
             [6 * L * H2 ** 3 * T2, 2 * L ** 2 * H2 ** 3 * T2, 4 * L ** 2 * H2 ** 3 * T2]])
        F = np.array([10000.0, 0.0, -1833.3])
        Kd = K.dot(q)
        res = Kd - F
        return res

    # Chute inicial para os deslocamentos
    q0 = np.array([[0.0001], [0.0001], [0.0001]])
    # Resolver o problema estático:
    sol = root(resfunc, q0)

    # Escolhendo o vetor deslocamento
    q = sol.x
    # Podemos obter disto o nosso volume
    dVdA = 4 * L * np.array([T1, H1, T2, H2])

    return dVdA


# Primeira restrição, deslocamento q1
def confun1(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    L = 1.0

    # Primeiro vamos buscar a função dos resíduos, para que todo deslocamento obedeça o equilíbrio
    def resfunc(q):
        E = 7.0 * 10 ** 10
        K = 2 / 3 * E / L ** 3 * np.array(
            [[12 * (H1 ** 3 * T1 + H2 ** 3 * T2), 6 * L * (-H1 ** 3 * T1 + H2 ** 3 * T2), 6 * L * H2 ** 3 * T2],
             [6 * L * (-H1 ** 3 * T1 + H2 ** 3 * T2), 4 * L ** 2 * (H1 ** 3 * T1 + H2 ** 3 * T2),
              2 * L ** 2 * H2 ** 3 * T2],
             [6 * L * H2 ** 3 * T2, 2 * L ** 2 * H2 ** 3 * T2, 4 * L ** 2 * H2 ** 3 * T2]])
        F = np.array([10000.0, 0.0, -1833.3])
        Kd = K.dot(q)
        res = Kd - F
        return res

    # Chute inicial para os deslocamentos
    q0 = np.array([[0.0001], [0.0001], [0.0001]])
    # Resolver o problema estático:
    sol = root(resfunc, q0)

    # Escolhendo o vetor deslocamento
    q = sol.x
    q1 = q[0]
    return 0.001 - q1


# Gradiente da primeira restrição
def confun1grad(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    L = 1.0
    E = 7.0 * 10 ** 10
    K = 2 / 3 * E / L ** 3 * np.array(
        [[12 * (H1 ** 3 * T1 + H2 ** 3 * T2), 6 * L * (-H1 ** 3 * T1 + H2 ** 3 * T2), 6 * L * H2 ** 3 * T2],
         [6 * L * (-H1 ** 3 * T1 + H2 ** 3 * T2), 4 * L ** 2 * (H1 ** 3 * T1 + H2 ** 3 * T2),
          2 * L ** 2 * H2 ** 3 * T2],
         [6 * L * H2 ** 3 * T2, 2 * L ** 2 * H2 ** 3 * T2, 4 * L ** 2 * H2 ** 3 * T2]])

    # Primeiro vamos buscar a função dos resíduos, para que todo deslocamento obedeça o equilíbrio
    def resfunc(q):
        F = np.array([10000.0, 0.0, -1833.3])
        Kd = K.dot(q)
        res = Kd - F
        return res

    # Chute inicial para os deslocamentos
    q0 = np.array([[0.0001], [0.0001], [0.0001]])
    # Resolver o problema estático:
    sol = root(resfunc, q0)
    # Escolhendo o vetor deslocamento
    q = sol.x
    # Precisamos resolver K(dqDA) = - (dKdA)q
    dKdH1 = 2 * E / L ** 3 * np.array(
        [[12 * (H1 ** 2 * T1), 6 * L * (-H1 ** 2 * T1), 0.0],
         [6 * L * (-H1 ** 2 * T1), 4 * L ** 2 * (H1 ** 2 * T1), 0.0],
         [0.0, 0.0, 0.0]])
    dKdT1 = 2 / 3 * E / L ** 3 * np.array(
        [[12 * (H1 ** 3), 6 * L * (-H1 ** 3), 0.0],
         [6 * L * (-H1 ** 3), 4 * L ** 2 * (H1 ** 3), 0.0],
         [0.0, 0.0, 0.0]])
    dKdH2 = 2 * E / L ** 3 * np.array(
        [[12 * (H2 ** 2 * T2), 6 * L * (H2 ** 2 * T2), 6 * L * H2 ** 2 * T2],
         [6 * L * (H2 ** 2 * T2), 4 * L ** 2 * (H2 ** 2 * T2), 2 * L ** 2 * H2 ** 2 * T2],
         [6 * L * H2 ** 2 * T2, 2 * L ** 2 * H2 ** 2 * T2, 4 * L ** 2 * H2 ** 2 * T2]])
    dKdT2 = 2 / 3 * E / L ** 3 * np.array(
        [[12 * (H2 ** 3), 6 * L * (H2 ** 3), 6 * L * H2 ** 3],
         [6 * L * (H2 ** 3), 4 * L ** 2 * (H2 ** 3), 2 * L ** 2 * H2 ** 3],
         [6 * L * H2 ** 3, 2 * L ** 2 * H2 ** 3, 4 * L ** 2 * H2 ** 3]])

    dqdH1 = np.linalg.solve(K, -dKdH1.dot(q))
    dqdT1 = np.linalg.solve(K, -dKdT1.dot(q))
    dqdH2 = np.linalg.solve(K, -dKdH2.dot(q))
    dqdT2 = np.linalg.solve(K, -dKdT2.dot(q))

    dqdA = np.array([dqdH1, dqdT1, dqdH2, dqdT2])
    dqdA = np.transpose(dqdA)
    dq1dA = dqdA[0, :]
    return - dq1dA


# Segunda restrição, limite de H1
def confun2(HT):
    H1 = HT[0]

    return 1 - H1 / 0.15


# Gradiente da segunda restrição
def confun2grad(HT):
    dH1dA = np.array([1.0, 0.0, 0.0, 0.0])
    return - dH1dA / 0.15


# Terceira restrição, limite de H2
def confun3(HT):
    H2 = HT[2]

    return 1 - H2 / 0.15


# Gradiente da terceira restrição
def confun3grad(HT):
    dH2dA = np.array([0.0, 0.0, 1.0, 0.0])
    return - dH2dA / 0.15


# Criando lista das restrições, lembrando que no SLSQP do python, ineq é >= 0:

con1 = {'type': 'ineq',
        'fun': confun1,
        'jac': confun1grad}
con2 = {'type': 'ineq',
        'fun': confun2,
        'jac': confun2grad}
con3 = {'type': 'ineq',
        'fun': confun3,
        'jac': confun3grad}

bounds = [[0.0, None]] * 4

# Rodando otimização
A0 = np.array([0.1, 0.002, 0.1, 0.002])
result = minimize(objfun, A0, jac=objfungrad,
                  constraints=[con1, con2, con3],
                  bounds=bounds,
                  method='slsqp',
                  callback=callback)
print(' ')
print('Resultado da otimização')
print('-----------------------')
print(result)
print(' ')
print('-----------------------')

A = result.x
print('Vol. = ', objfun(A))
k = result.nit
print(' ')
print('-----------------------')
print(' ')

# Organizando os dados para gerar gráfico do progresso das variáveis
h1 = np.zeros((13,1))
t1 = np.zeros((13,1))
h2 = np.zeros((13,1))
t2 = np.zeros((13,1))
v = np.zeros((13, 1))
Q = np.zeros((13, 1))
v[0] = objfun(A0)
# os deslocamentos estarão em mm
Q[0] = (0.001 - confun1(A0))*1000
h1[0], t1[0], h2[0], t2[0] = A0[0], A0[1], A0[2], A0[3]
for i in range(k):
    h1[i+1], t1[i+1], h2[i+1], t2[i+1] = X[i]
    v[i+1] = Volume[i]
    Q[i+1] = (0.001 - q100[i])*1000

print(Q[12])

fig, ax = plt.subplots(2,2)
fig.suptitle("Variação do vetor área", fontsize=16)

ax[0,0].plot(h1, marker='o', color = 'crimson', ls='None')
ax[0,0].tick_params(axis="y", labelsize=10)
ax[0,0].tick_params(axis="x", labelsize=10)
ax[0,0].set_ylabel('$H_1$', fontsize=12)
ax[0,0].set_xlabel('$n_{iter.}$', fontsize=12)
ax[0,1].plot(t1, marker='o', color = 'crimson', ls='None')
ax[0,1].tick_params(axis="y", labelsize=10)
ax[0,1].tick_params(axis="x", labelsize=10)
ax[0,1].set_ylabel('$T_1$', fontsize=12)
ax[0,1].set_xlabel('$n_{iter.}$', fontsize=12)
ax[1,0].plot(h2, marker='o', color = 'crimson', ls='None')
ax[1,0].tick_params(axis="y", labelsize=10)
ax[1,0].tick_params(axis="x", labelsize=10)
ax[1,0].set_ylabel('$H_2$', fontsize=12)
ax[1,0].set_xlabel('$n_{iter.}$', fontsize=12)
ax[1,1].plot(t2, marker='o', color = 'crimson', ls='None')
ax[1,1].tick_params(axis="y", labelsize=10)
ax[1,1].tick_params(axis="x", labelsize=10)
ax[1,1].set_ylabel('$T_2$', fontsize=12)
ax[1,1].set_xlabel('$n_{iter.}$', fontsize=12)

ax[0,0].grid()
ax[1,0].grid()
ax[0,1].grid()
ax[1,1].grid()
fig.tight_layout()
fig.subplots_adjust(top=0.88)
plt.show()

fig, ax = plt.subplots(1)
fig.suptitle("Variação do volume", fontsize=16)

plt.plot(v, marker='o', color = 'crimson', ls='None')
plt.xlabel('$n_{iter.}$', fontsize=12)
plt.ylabel('$Vol.$', fontsize=12)
fig.tight_layout()
fig.subplots_adjust(top=0.88)

plt.grid()
plt.show()

fig, ax = plt.subplots(1)
fig.suptitle("Variação da restrição", fontsize=16)

plt.plot(Q, marker='o', color = 'crimson', ls='None')
plt.xlabel('$n_{iter.}$', fontsize=12)
plt.ylabel('$q_1.$', fontsize=12)
fig.tight_layout()
fig.subplots_adjust(top=0.88)

plt.grid()
plt.show()