'''
Homework 4, item 2
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
# Algumas constantes do problema:
A0 = np.array([0.1, 0.002, 0.1, 0.002])
H1, T1, H2, T2 = A0[0], A0[1], A0[2], A0[3]
L = 1.0
E = 7.0 * 10 ** 10
K = 2 / 3 * E / L ** 3 * np.array(
    [[12 * (H1 ** 3 * T1 + H2 ** 3 * T2), 6 * L * (-H1 ** 3 * T1 + H2 ** 3 * T2), 6 * L * H2 ** 3 * T2],
     [6 * L * (-H1 ** 3 * T1 + H2 ** 3 * T2), 4 * L ** 2 * (H1 ** 3 * T1 + H2 ** 3 * T2),
      2 * L ** 2 * H2 ** 3 * T2],
     [6 * L * H2 ** 3 * T2, 2 * L ** 2 * H2 ** 3 * T2, 4 * L ** 2 * H2 ** 3 * T2]])
F = np.array([10000.0, 0.0, -1833.3])
q0 = np.array([[0.001], [0.001], [0.001]])

# Função objetivo, volume
def objfun(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    L = 1.0
    # Podemos obter disto o nosso volume
    V = 4 * (T1 * H1 + T2 * H2) * L
    return V

# Gradiente da função objetivo
def objfungrad(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    L = 1.0
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
    q0 = np.array([[0.001], [0.001], [0.001]])
    # Resolver o problema estático:
    sol = root(resfunc, q0)

    # Escolhendo o vetor deslocamento
    q = sol.x
    q1 = q[0]
    q1n = q1 / 0.001
    return 1 - q1n

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
    q0 = np.array([[0.001], [0.001], [0.001]])
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
    dq1daN = dq1dA / 0.001
    return - dq1daN

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

# Definindo o OST como função da Janela
def OST(Janela):
    # Iniciando o vetor para não gastar tempo realocando memória

    obj = np.zeros((100, 4))
    obj[0, :] = A0
    v = np.zeros(100)
    v[0] = objfun(A0)
    p = 0
    q = np.zeros((100, 1))
    # Problema sequencial aproximado, montado num do, while(não converge)
    while True:
        Ap = obj[p, :]
        H1, T1, H2, T2 = Ap[0], Ap[1], Ap[2], Ap[3]
        Vp = objfun(Ap)
        v[p] = Vp
        g1 = confun1(Ap)
        q[p] = g1
        # Vamos identificar cada restrição por uma FLAG
        FLAG1 = True

        # O valor absoluto está forçando a sempre colocarmos todas as restrições (simplificando o problema)
        if abs(g1) > 1E-5:
            FLAG1 = False
        # Se todas as FLAGS são verdadeiras, o problema respeita todas as restrições
        if all([(v[p] - v[p - 1]) < (v[p - 1] * 0.01), FLAG1]):
            break
        if p >= 99:
            print('ERRO: NÃO CONVERGIU')
            break

        # Precisamos calcular o gradiente das restrições que não são 0
        dq1dA = np.zeros(4)
        dH1dA = np.zeros(4)
        dH2dA = np.zeros(4)
        dq1dA = confun1grad(Ap)

        # Precisamos do gradiente do volume
        dVdA = objfungrad(Ap)

        # Agora podemos resolver o problema aproximado
        def subobj(a):
            h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
            sobj = 4 * (t1 * h1 + t2 * h2)
            return sobj

        def subobjgrad(a):
            h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
            return 4 * L * np.array([t1, h1, t2, h2])

        def subcon1(a):
            dq1dH1 = dq1dA[0]
            dq1dT1 = dq1dA[1]
            dq1dH2 = dq1dA[2]
            dq1dT2 = dq1dA[3]
            h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
            scon1 = g1 + (h1 - H1) * H1 / h1 * dq1dH1 + (t1 - T1) * T1 / t1 * dq1dT1 + (h2 - H2) * H2 / h2 * dq1dH2 + \
                    (t2 - T2) * T2 / t2 * dq1dT2
            return scon1

        def subcon1grad(a):
            dq1dH1 = dq1dA[0]
            dq1dT1 = dq1dA[1]
            dq1dH2 = dq1dA[2]
            dq1dT2 = dq1dA[3]
            h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
            return (H1 / h1) ** 2 * dq1dH1 * np.array([1.0, 0.0, 0.0, 0.0]) \
                   + (T1 / t1) ** 2 * dq1dT1 * np.array([0.0, 1.0, 0.0, 0.0]) \
                   + (H2 / h2) ** 2 * dq1dH2 * np.array([0.0, 0.0, 1.0, 0.0]) \
                   + (T2 / t2) ** 2 * dq1dT2 * np.array([0.0, 0.0, 0.0, 1.0])


        # Criando lista das restrições, lembrando que no SLSQP do python, ineq é >= 0:
        con1 = {'type': 'ineq',
                'fun': subcon1,
                'jac': subcon1grad}

        max1 = np.amin([0.15, Ap[0] * (1 + Janela/100)])
        min1 = np.amax([0.0, Ap[0] * (1 - Janela/100)])
        max2 = (1 + Janela/100) * Ap[1]
        min2 = np.amax([0.0, (1 - Janela/100) * Ap[1]])
        max3 = np.amin([0.15, Ap[2] * (1 + Janela/100)])
        min3 = np.amax([0.0, Ap[2] * (1 - Janela/100)])
        max4 = (1 + Janela/100) * Ap[3]
        min4 = np.amax([0.0, (1 - Janela/100) * Ap[3]])
        bounds = [[min1, max1], [min2, max2], [min3, max3], [min4, max4]]
        result = minimize(subobj, Ap, constraints=[con1], jac = subobjgrad,
                          bounds=bounds,
                          method='slsqp')
        p = p + 1
        obj[p, :] = result.x

    print(' ')
    print('Resultado da otimização')
    print('-----------------------')
    print(obj[p])
    print(' ')
    print('-----------------------')

    A = obj[p]
    print('Vol. = ', objfun(A))
    print(' ')
    print('-----------------------')
    print('n loops = ', p)
    return q, obj, v

# Variando a janela para observar seu efeito nos resultados
q5, x5, v5 = OST(5)
q, x, v = OST(30)

# Os índices a seguir serviram apenas de testes intermediários
OST(65)
OST(100)

print('q1 = ', 1 - q[7])
#print(q[0:10])
#print(x[0:10, :])

fig, ax = plt.subplots(2,2)
fig.suptitle("Variação do vetor área", fontsize=16)

h1, t1, h2, t2 = x[0:7, 0], x[0:7, 1], x[0:7, 2], x[0:7, 3]
h15, t15, h25, t25 = x5[0:25, 0], x5[0:25, 1], x5[0:25, 2], x5[0:25, 3]
ax[0,0].plot(h1, marker='o', color = 'crimson', ls='None',label='30%')
ax[0,0].plot(h15, marker='o', color = 'yellow', ls='None', label='5%')
ax[0,0].tick_params(axis="y", labelsize=10)
ax[0,0].tick_params(axis="x", labelsize=10)
ax[0,0].set_ylabel('$H_1$', fontsize=12)
ax[0,0].set_xlabel('$n_{iter.}$', fontsize=12)
ax[0,0].legend(['30\%', '5\%'])
ax[0,1].plot(t1, marker='o', color = 'crimson', ls='None')
ax[0,1].plot(t15, marker='o', color = 'yellow', ls='None')
ax[0,1].tick_params(axis="y", labelsize=10)
ax[0,1].tick_params(axis="x", labelsize=10)
ax[0,1].set_ylabel('$T_1$', fontsize=12)
ax[0,1].set_xlabel('$n_{iter.}$', fontsize=12)
ax[0,1].legend(['30\%', '5\%'])
ax[1,0].plot(h2, marker='o', color = 'crimson', ls='None')
ax[1,0].plot(h25, marker='o', color = 'yellow', ls='None')
ax[1,0].tick_params(axis="y", labelsize=10)
ax[1,0].tick_params(axis="x", labelsize=10)
ax[1,0].set_ylabel('$H_2$', fontsize=12)
ax[1,0].set_xlabel('$n_{iter.}$', fontsize=12)
ax[1,0].legend(['30\%', '5\%'])
ax[1,1].plot(t2, marker='o', color = 'crimson', ls='None')
ax[1,1].plot(t25, marker='o', color = 'yellow', ls='None')
ax[1,1].tick_params(axis="y", labelsize=10)
ax[1,1].tick_params(axis="x", labelsize=10)
ax[1,1].set_ylabel('$T_2$', fontsize=12)
ax[1,1].set_xlabel('$n_{iter.}$', fontsize=12)
ax[1,1].legend(['30\%', '5\%'])
fig.tight_layout()
fig.subplots_adjust(top=0.88)
ax[0,0].grid()
ax[1,0].grid()
ax[0,1].grid()
ax[1,1].grid()
plt.show()

fig, ax = plt.subplots(1)
fig.suptitle("Variação do volume", fontsize=16)

plt.plot(v[0:7], marker='o', color = 'crimson', ls='None')
plt.plot(v5[0:25], marker='o', color = 'yellow', ls='None')
plt.xlabel('$n_{iter.}$', fontsize=12)
plt.ylabel('$Vol.$', fontsize=12)
fig.tight_layout()
fig.subplots_adjust(top=0.88)
plt.legend(['30\%', '5\%'])
plt.grid()
plt.show()

fig, ax = plt.subplots(1)
fig.suptitle("Variação da restrição", fontsize=16)

plt.plot(1-q[0:7], marker='o', color = 'crimson', ls='None')
plt.plot(1-q5[0:25], marker='o', color = 'yellow', ls='None')
plt.xlabel('$n_{iter.}$', fontsize=12)
plt.ylabel('$q_1.$', fontsize=12)
fig.tight_layout()
fig.subplots_adjust(top=0.88)
plt.legend(['30\%', '5\%'])
plt.grid()
plt.show()

