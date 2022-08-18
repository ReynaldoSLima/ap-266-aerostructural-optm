'''
Homework 4, item 3
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

# Valor inicial de área
A0 = np.array([0.1, 0.002, 0.1, 0.002])
H1, T1, H2, T2 = A0[0], A0[1], A0[2], A0[3]

# Função de acompanhamento da otimização:
X = []
V = []
g1n = []
g2n = []
g1c = []
g2c = []
def callback(x):
    fobj = objfun(x)
    V.append(fobj)
    X.append(x)
    g1n.append(conGNL1(x))
    g2n.append(conGNL2(x))
    g1c.append(conGCR1(x))
    g2c.append(conGCR2(x))

# Função objetivo, volume
def objfun(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    L = 1.0
    # Podemos obter disto o nosso volume
    V = 4 * (T1 * H1 + T2 * H2) * L
    return V

# gradiente da função objetivo
def objfungrad(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    L = 1.0
    # Podemos obter disto o nosso volume
    dVdA = 4 * L * np.array([T1, H1, T2, H2])

    return dVdA

# Função que resolve os displacements para o FEM
def disp(HT):
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
    q2 = q[1]
    q3 = q[2]
    return q1, q2, q3

# Gradiente dos displacements
def dispGrad(HT):
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
    dq2dA = dqdA[1, :]
    dq3dA = dqdA[2, :]
    return dq1dA, dq2dA, dq3dA

# Condição de tensão normal, elemento 1, elo 1
def conGNL1(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    W1 = 4 / 3 * H1 ** 2 * T1
    L = 1.0
    E = 7.0 * 10 ** 10
    # Matriz AUX refere-se aos termos constantes das matrizes de tensão 4x4 (para cada elemento)
    AUX = E / L ** 3 * np.array([[12, 6 * L, - 12, 6 * L],
                                 [6 * L, 4 * L ** 2, - 6 * L, 2 * L ** 2],
                                 [-12, -6 * L, 12, - 6 * L],
                                 [6 * L, 2 * L ** 2, - 6 * L, 4 * L ** 2]])
    # Termos constantes de força:
    F1 = np.array([5000.0, 10000.0 / 12, 5000.0, -10000.0 / 12])
    I1 = 2 / 3 * H1 ** 3 * T1
    K1 = I1 * AUX
    q1, q2, q3 = disp(HT)
    u1 = np.array([0, 0, q1, q2])
    S1 = K1.dot(u1) - F1
    M1 = S1[1]
    sigma1 = M1 / W1
    Y = - 9.0E7
    if sigma1 > 0:
        Y = 9.0E7
    r = sigma1 / Y
    return 1.0 - r

# Gradiente da restrição gNL1
def conGNL1grad(HT):
    Y = - 9.0E7
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    F1 = np.array([5000.0, 10000.0 / 12, 5000.0, -10000.0 / 12])
    W1 = 4 / 3 * H1 ** 2 * T1
    L = 1.0
    E = 7.0 * 10 ** 10
    # Matriz AUX refere-se aos termos constantes das matrizes de tensão 4x4 (para cada elemento)
    AUX = E / L ** 3 * np.array([[12, 6 * L, - 12, 6 * L],
                                 [6 * L, 4 * L ** 2, - 6 * L, 2 * L ** 2],
                                 [-12, -6 * L, 12, - 6 * L],
                                 [6 * L, 2 * L ** 2, - 6 * L, 4 * L ** 2]])
    # Termos constantes de força:
    I1 = 2 / 3 * H1 ** 3 * T1
    dI1dH1 = 2 * H1 ** 2 * T1
    dI1dT1 = 2 / 3 * H1 ** 3
    K1 = I1 * AUX
    dK1dH1 = dI1dH1 * AUX
    dK1dT1 = dI1dT1 * AUX
    q1, q2, q3 = disp(HT)
    dq1dA, dq2dA, dq3dA = dispGrad(HT)
    u1 = np.array([0, 0, q1, q2])
    du1dH1 = np.array([0, 0, dq1dA[0], dq2dA[0]])
    du1dT1 = np.array([0, 0, dq1dA[1], dq2dA[1]])
    du1dH2 = np.array([0, 0, dq1dA[2], dq2dA[2]])
    du1dT2 = np.array([0, 0, dq1dA[3], dq2dA[3]])
    dS1dH1 = dK1dH1.dot(u1) + K1.dot(du1dH1)
    dS1dT1 = dK1dT1.dot(u1) + K1.dot(du1dT1)
    dS1dH2 = K1.dot(du1dH2)
    dS1dT2 = K1.dot(du1dT2)
    dM1dH1 = dS1dH1[1]
    dM1dT1 = dS1dT1[1]
    dM1dH2 = dS1dH2[1]
    dM1dT2 = dS1dT2[1]
    S1 = K1.dot(u1) - F1
    M1 = S1[1]
    if M1 > 0:
        Y = 9.0E7
    dW1dH1 = 8 / 3 * H1 * T1
    dW1dT1 = 4 / 3 * H1 ** 2
    dg1dH1 = 1 / Y * (- dW1dH1 * M1 / W1 ** 2 + 1 / W1 * dM1dH1)
    dg1dT1 = 1 / Y * (- dW1dT1 * M1 / W1 ** 2 + 1 / W1 * dM1dT1)
    dg1dH2 = 1 / Y * (1 / W1 * dM1dH2)
    dg1dT2 = 1 / Y * (1 / W1 * dM1dT2)
    dg1dA = np.array([dg1dH1, dg1dT1, dg1dH2, dg1dT2])
    return - dg1dA

# Condição de tensão normal, elemento 2, elo 1
def conGNL2(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    W2 = 4 / 3 * H2 ** 2 * T2
    L = 1.0
    E = 7.0 * 10 ** 10
    # Matriz AUX refere-se aos termos constantes das matrizes de tensão 4x4 (para cada elemento)
    AUX = E / L ** 3 * np.array([[12, 6 * L, - 12, 6 * L],
                                 [6 * L, 4 * L ** 2, - 6 * L, 2 * L ** 2],
                                 [-12, -6 * L, 12, - 6 * L],
                                 [6 * L, 2 * L ** 2, - 6 * L, 4 * L ** 2]])
    # Termos constantes de força:
    F2 = np.array([5000.0, 10000.0 / 12, 5000.0, -10000.0 / 12])
    I2 = 2 / 3 * H2 ** 3 * T2
    K2 = I2 * AUX
    q1, q2, q3 = disp(HT)
    u2 = np.array([q1, q2, 0, q3])
    S2 = K2.dot(u2) - F2
    M2 = S2[1]
    sigma2 = M2 / W2
    Y = 9.0E7
    if sigma2 < 0:
        Y = - 9.0E7
    r = sigma2 / Y
    return 1.0 - r

# Gradiente da restrição gNL2
def conGNL2grad(HT):
    Y = 9.0E7
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    F2 = np.array([5000.0, 10000.0 / 12, 5000.0, -10000.0 / 12])
    W2 = 4 / 3 * H2 ** 2 * T2
    L = 1.0
    E = 7.0 * 10 ** 10
    # Matriz AUX refere-se aos termos constantes das matrizes de tensão 4x4 (para cada elemento)
    AUX = E / L ** 3 * np.array([[12, 6 * L, - 12, 6 * L],
                                 [6 * L, 4 * L ** 2, - 6 * L, 2 * L ** 2],
                                 [-12, -6 * L, 12, - 6 * L],
                                 [6 * L, 2 * L ** 2, - 6 * L, 4 * L ** 2]])
    # Termos constantes de força:
    I2 = 2 / 3 * H2 ** 3 * T2
    dI2dH2 = 2 * H2 ** 2 * T2
    dI2dT2 = 2 / 3 * H2 ** 3
    K2 = I2 * AUX
    dK2dH2 = dI2dH2 * AUX
    dK2dT2 = dI2dT2 * AUX
    q1, q2, q3 = disp(HT)
    dq1dA, dq2dA, dq3dA = dispGrad(HT)
    u2 = np.array([q1, q2, 0, q3])
    du2dH1 = np.array([dq1dA[0], dq2dA[0], 0, dq3dA[0]])
    du2dT1 = np.array([dq1dA[1], dq2dA[1], 0, dq3dA[1]])
    du2dH2 = np.array([dq1dA[2], dq2dA[2], 0, dq3dA[2]])
    du2dT2 = np.array([dq1dA[3], dq2dA[3], 0, dq3dA[3]])
    dS2dH1 = K2.dot(du2dH1)
    dS2dT1 = K2.dot(du2dT1)
    dS2dH2 = dK2dH2.dot(u2) + K2.dot(du2dH2)
    dS2dT2 = dK2dT2.dot(u2) + K2.dot(du2dT2)
    dM2dH1 = dS2dH1[1]
    dM2dT1 = dS2dT1[1]
    dM2dH2 = dS2dH2[1]
    dM2dT2 = dS2dT2[1]
    S2 = K2.dot(u2) - F2
    M2 = S2[1]
    if M2 < 0:
        Y = - 9.0E7
    dW2dH2 = 8 / 3 * H2 * T2
    dW2dT2 = 4 / 3 * H2 ** 2
    dg2dH1 = 1 / Y * (1 / W2 * dM2dH1)
    dg2dT1 = 1 / Y * (1 / W2 * dM2dT1)
    dg2dH2 = 1 / Y * (- dW2dH2 * M2 / (W2 ** 2) + 1 / W2 * dM2dH2)
    dg2dT2 = 1 / Y * (- dW2dT2 * M2 / (W2 ** 2) + 1 / W2 * dM2dT2)
    dg2dA = np.array([dg2dH1, dg2dT1, dg2dH2, dg2dT2])
    return -dg2dA

# Condição de não flambagem, elemento 1, elo 1
def conGCR1(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    W1 = 4 / 3 * H1 ** 2 * T1
    L = 1.0
    E = 7.0 * 10 ** 10
    v = 0.3
    # Matriz AUX refere-se aos termos constantes das matrizes de tensão 4x4 (para cada elemento)
    AUX = E / L ** 3 * np.array([[12, 6 * L, - 12, 6 * L],
                                 [6 * L, 4 * L ** 2, - 6 * L, 2 * L ** 2],
                                 [-12, -6 * L, 12, - 6 * L],
                                 [6 * L, 2 * L ** 2, - 6 * L, 4 * L ** 2]])
    # Termos constantes de força:
    F1 = np.array([5000.0, 10000.0 / 12, 5000.0, -10000.0 / 12])
    I1 = 2 / 3 * H1 ** 3 * T1
    K1 = I1 * AUX
    q1, q2, q3 = disp(HT)
    u1 = np.array([0, 0, q1, q2])
    S1 = K1.dot(u1) - F1
    M1 = S1[1]
    sigma1 = M1 / W1
    YCR = - E * 4 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T1 ** 2) / H1 ** 2
    if sigma1 > 0:
        YCR = E * 4 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T1 ** 2) / H1 ** 2
    r = sigma1 / YCR
    return 1.0 - r

# Gradiente da condição gCR1
def conGCR1grad(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    F1 = np.array([5000.0, 10000.0 / 12, 5000.0, -10000.0 / 12])
    W1 = 4 / 3 * H1 ** 2 * T1
    L = 1.0
    v = 0.3
    E = 7.0 * 10 ** 10
    # Matriz AUX refere-se aos termos constantes das matrizes de tensão 4x4 (para cada elemento)
    AUX = E / L ** 3 * np.array([[12, 6 * L, - 12, 6 * L],
                                 [6 * L, 4 * L ** 2, - 6 * L, 2 * L ** 2],
                                 [-12, -6 * L, 12, - 6 * L],
                                 [6 * L, 2 * L ** 2, - 6 * L, 4 * L ** 2]])
    # Termos constantes de força:
    I1 = 2 / 3 * H1 ** 3 * T1
    dI1dH1 = 2 * H1 ** 2 * T1
    dI1dT1 = 2 / 3 * H1 ** 3
    K1 = I1 * AUX
    dK1dH1 = dI1dH1 * AUX
    dK1dT1 = dI1dT1 * AUX
    q1, q2, q3 = disp(HT)
    dq1dA, dq2dA, dq3dA = dispGrad(HT)
    u1 = np.array([0, 0, q1, q2])
    du1dH1 = np.array([0, 0, dq1dA[0], dq2dA[0]])
    du1dT1 = np.array([0, 0, dq1dA[1], dq2dA[1]])
    du1dH2 = np.array([0, 0, dq1dA[2], dq2dA[2]])
    du1dT2 = np.array([0, 0, dq1dA[3], dq2dA[3]])
    dS1dH1 = dK1dH1.dot(u1) + K1.dot(du1dH1)
    dS1dT1 = dK1dT1.dot(u1) + K1.dot(du1dT1)
    dS1dH2 = K1.dot(du1dH2)
    dS1dT2 = K1.dot(du1dT2)
    dM1dH1 = dS1dH1[1]
    dM1dT1 = dS1dT1[1]
    dM1dH2 = dS1dH2[1]
    dM1dT2 = dS1dT2[1]
    S1 = K1.dot(u1) - F1
    M1 = S1[1]
    sigma1 = M1 / W1
    YCR = - E * 4 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T1 ** 2) / H1 ** 2
    dYCRdH1 = E * 8 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T1 ** 2) / H1 ** 3
    dYCRdT1 = - E * 8 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * T1 / H1 ** 2
    dW1dH1 = 8 / 3 * H1 * T1
    dW1dT1 = 4 / 3 * H1 ** 2
    if sigma1 > 0:
        YCR = E * 4 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T1 ** 2) / H1 ** 2
        dYCRdH1 = - E * 8 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T1 ** 2) / H1 ** 3
        dYCRdT1 = E * 8 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * T1 / H1 ** 2
    dg1dH1 = 1 / YCR * (- dW1dH1 * M1 / W1 ** 2 + 1 / W1 * dM1dH1) - 1 / YCR ** 2 * sigma1 * dYCRdH1
    dg1dT1 = 1 / YCR * (- dW1dT1 * M1 / W1 ** 2 + 1 / W1 * dM1dT1) - 1 / YCR ** 2 * sigma1 * dYCRdT1
    dg1dH2 = 1 / YCR * (1 / W1 * dM1dH2)
    dg1dT2 = 1 / YCR * (1 / W1 * dM1dT2)
    dgCR1dA = np.array([dg1dH1, dg1dT1, dg1dH2, dg1dT2])
    return - dgCR1dA

# Condição de não flambagem, elemento 2, elo 1
def conGCR2(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    W2 = 4 / 3 * H2 ** 2 * T2
    L = 1.0
    E = 7.0 * 10 ** 10
    v = 0.3
    # Matriz AUX refere-se aos termos constantes das matrizes de tensão 4x4 (para cada elemento)
    AUX = E / L ** 3 * np.array([[12, 6 * L, - 12, 6 * L],
                                 [6 * L, 4 * L ** 2, - 6 * L, 2 * L ** 2],
                                 [-12, -6 * L, 12, - 6 * L],
                                 [6 * L, 2 * L ** 2, - 6 * L, 4 * L ** 2]])
    # Termos constantes de força:
    F2 = np.array([5000.0, 10000.0 / 12, 5000.0, -10000.0 / 12])
    I2 = 2 / 3 * H2 ** 3 * T2
    K2 = I2 * AUX
    q1, q2, q3 = disp(HT)
    u2 = np.array([q1, q2, 0, q3])
    S2 = K2.dot(u2) - F2
    M2 = S2[1]
    sigma2 = - M2 / W2
    YCR = - E * 4 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T2 ** 2) / H2 ** 2
    if sigma2 > 0:
        YCR = E * 4 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T2 ** 2) / H2 ** 2
    r = sigma2 / YCR
    return 1.0 - r

# Gradiente da
def conGCR2grad(HT):
    H1, T1, H2, T2 = HT[0], HT[1], HT[2], HT[3]
    F2 = np.array([5000.0, 10000.0 / 12, 5000.0, -10000.0 / 12])
    W2 = 4 / 3 * H2 ** 2 * T2
    L = 1.0
    v = 0.3
    E = 7.0 * 10 ** 10
    # Matriz AUX refere-se aos termos constantes das matrizes de tensão 4x4 (para cada elemento)
    AUX = E / L ** 3 * np.array([[12, 6 * L, - 12, 6 * L],
                                 [6 * L, 4 * L ** 2, - 6 * L, 2 * L ** 2],
                                 [-12, -6 * L, 12, - 6 * L],
                                 [6 * L, 2 * L ** 2, - 6 * L, 4 * L ** 2]])
    # Termos constantes de força:
    I2 = 2 / 3 * H2 ** 3 * T2
    dI2dH2 = 2 * H2 ** 2 * T2
    dI2dT2 = 2 / 3 * H2 ** 3
    K2 = I2 * AUX
    dK2dH2 = dI2dH2 * AUX
    dK2dT2 = dI2dT2 * AUX
    q1, q2, q3 = disp(HT)
    dq1dA, dq2dA, dq3dA = dispGrad(HT)
    u2 = np.array([q1, q2, 0, q3])
    du2dH1 = np.array([dq1dA[0], dq2dA[0], 0, dq3dA[0]])
    du2dT1 = np.array([dq1dA[1], dq2dA[1], 0, dq3dA[1]])
    du2dH2 = np.array([dq1dA[2], dq2dA[2], 0, dq3dA[2]])
    du2dT2 = np.array([dq1dA[3], dq2dA[3], 0, dq3dA[3]])
    dS2dH1 = K2.dot(du2dH1)
    dS2dT1 = K2.dot(du2dT1)
    dS2dH2 = dK2dH2.dot(u2) + K2.dot(du2dH2)
    dS2dT2 = dK2dT2.dot(u2) + K2.dot(du2dT2)
    dM2dH1 = dS2dH1[1]
    dM2dT1 = dS2dT1[1]
    dM2dH2 = dS2dH2[1]
    dM2dT2 = dS2dT2[1]
    S2 = K2.dot(u2) - F2
    M2 = S2[1]
    sigma2 = - M2 / W2
    YCR = - E * 4 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T2 ** 2) / H2 ** 2
    dYCRdH2 = E * 8 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T2 ** 2) / H2 ** 3
    dYCRdT2 = - E * 8 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * T2 / H2 ** 2
    dW2dH2 = 8 / 3 * H2 * T2
    dW2dT2 = 4 / 3 * H2 ** 2
    if sigma2 > 0:
        YCR = E * 4 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T2 ** 2) / H2 ** 2
        dYCRdH2 = - E * 8 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * (T2 ** 2) / H2 ** 3
        dYCRdT2 = E * 8 * (np.pi ** 2) / 12 * 1 / (1 - v ** 2) * T2 / H2 ** 2
    dg2dH1 = - 1 / YCR * (1 / W2 * dM2dH1)
    dg2dT1 = - 1 / YCR * (1 / W2 * dM2dT1)
    dg2dH2 = - 1 / YCR * (- dW2dH2 * (M2) / (W2 ** 2) + 1 / W2 * dM2dH2) - 1 / YCR ** 2 * sigma2 * dYCRdH2
    dg2dT2 = - 1 / YCR * (- dW2dT2 * (M2) / (W2 ** 2) + 1 / W2 * dM2dT2) - 1 / YCR ** 2 * sigma2 * dYCRdT2
    dgCR2dA = np.array([dg2dH1, dg2dT1, dg2dH2, dg2dT2])
    return - dgCR2dA


def conH1(HT):
    H1 = HT[0]

    return 1 - H1 / 0.3


def conH1grad(HT):
    dH1dA = np.array([1.0, 0.0, 0.0, 0.0])
    return - dH1dA / 0.3


def conH2(HT):
    H2 = HT[2]

    return 1 - H2 / 0.3


def conH2grad(HT):
    dH2dA = np.array([0.0, 0.0, 1.0, 0.0])
    return - dH2dA / 0.3


# Criando lista das restrições, lembrando que no SLSQP do python, ineq é >= 0:
con1 = {'type': 'ineq',
        'fun': conGNL1,
        'jac': conGNL1grad}
con2 = {'type': 'ineq',
        'fun': conGNL2,
        'jac': conGNL2grad}
con3 = {'type': 'ineq',
        'fun': conGCR1,
        'jac': conGCR1grad}
con4 = {'type': 'ineq',
        'fun': conGCR2,
        'jac': conGCR2grad}
con5 = {'type': 'ineq',
        'fun': conH1,
        'jac': conH1grad}
con6 = {'type': 'ineq',
        'fun': conH2,
        'jac': conH2grad}

bounds = [[0.0003, None]] * 4

# Rodando otimização inicial
A0 = np.array([0.1, 0.002, 0.1, 0.002])
result = minimize(objfun, A0, jac=objfungrad,
                  constraints=[con1, con2, con3, con4, con5, con6],
                  bounds=bounds,
                  method='slsqp',
                  callback=callback)
print(' ')
print('Resultado da otimização')
print('-----------------------')
print(' ')
print(result)
print(' ')
print('-----------------------')
A = result.x
k = result.nit
print('Vol. = ', objfun(A))
print('qt = ', disp(A))
print('A = ', A)
print('gN1 = ', conGNL1(A))
print('gN2 = ', conGNL2(A))
print('gCR1 = ', conGCR1(A))
print('gCR2 = ', conGCR2(A))

print(' ')
print('-----------------------')
obj = np.zeros((50, 4))
obj[0, :] = A0
v = np.zeros(50)
v[0] = objfun(A0)
p = 0
L = 1.0

# Vamos guardar o histórico de restrições para o OSA
G1 = np.zeros((6, 1))
G2 = np.zeros((6, 1))
G3 = np.zeros((6, 1))
G4 = np.zeros((6 ,1))

# As restrições calculadas foram negativas
G1[0] = -(conGNL1(A0))
G2[0] = -(conGNL2(A0))
G3[0] = -(conGCR1(A0))
G4[0] = -(conGCR2(A0))
# Problema sequencial aproximado, montado num do, while(não converge)
# Para esse problema não utilizo os limitantes de H1 e H2 como restrições, mas sim nas janelas
while True:
    Ap = obj[p, :]
    H1, T1, H2, T2 = Ap[0], Ap[1], Ap[2], Ap[3]
    Vp = objfun(Ap)
    v[p] = Vp
    g1 = conGNL1(Ap)
    g2 = conGNL2(Ap)
    g3 = conGCR1(Ap)
    g4 = conGCR2(Ap)
    G1[p+1] = -g1
    G2[p+1] = -g2
    G3[p+1] = -g3
    G4[p+1] = -g4
    # Vamos identificar cada restrição por uma FLAG
    FLAG1 = True
    FLAG2 = True
    FLAG3 = True
    FLAG4 = True

    # O valor absoluto está forçando a sempre colocarmos todas as restrições (simplificando o problema)

    if abs(g1) > 1E-5:
        FLAG1 = False
    if g2 < - 1E-5:
        FLAG2 = False
    if abs(g3) > 1E-5:
        FLAG3 = False
    if abs(g4) > 1E-5:
        FLAG4 = False

    # Se todas as FLAGS são verdadeiras, o problema respeita todas as restrições
    if all([(v[p] - v[p - 1]) < (v[p - 1] * 0.01), FLAG1, FLAG2, FLAG3, FLAG4]):
        break
    if p >= 49:
        print('ERRO: NÃO CONVERGIU')
        break

    # Precisamos calcular o gradiente das restrições que não são 0
    dgN1dA = np.zeros(4)
    dgN2dA = np.zeros(4)
    dgC1dA = np.zeros(4)
    dgC2dA = np.zeros(4)
    if not FLAG1:
        dgN1dA = conGNL1grad(Ap)
    if not FLAG2:
        dgN2dA = conGNL2grad(Ap)
    if not FLAG3:
        dgC1dA = conGCR1grad(Ap)
    if not FLAG4:
        dgC2dA = conGCR2grad(Ap)

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
        dgdH1 = dgN1dA[0]
        dgdT1 = dgN1dA[1]
        dgdH2 = dgN1dA[2]
        dgdT2 = dgN1dA[3]
        h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
        scon1 = g1 + (h1 - H1) * H1 / h1 * dgdH1 + (t1 - T1) * T1 / t1 * dgdT1 + (h2 - H2) * H2 / h2 * dgdH2 + \
                (t2 - T2) * T2 / t2 * dgdT2
        return scon1


    def subcon1grad(a):
        dgdH1 = dgN1dA[0]
        dgdT1 = dgN1dA[1]
        dgdH2 = dgN1dA[2]
        dgdT2 = dgN1dA[3]
        h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
        return (H1 / h1) ** 2 * dgdH1 * np.array([1.0, 0.0, 0.0, 0.0]) \
               + (T1 / t1) ** 2 * dgdT1 * np.array([0.0, 1.0, 0.0, 0.0]) \
               + (H2 / h2) ** 2 * dgdH2 * np.array([0.0, 0.0, 1.0, 0.0]) \
               + (T2 / t2) ** 2 * dgdT2 * np.array([0.0, 0.0, 0.0, 1.0])


    def subcon2(a):
        dgdH1 = dgN2dA[0]
        dgdT1 = dgN2dA[1]
        dgdH2 = dgN2dA[2]
        dgdT2 = dgN2dA[3]
        h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
        scon2 = g2 + (h1 - H1) * H1 / h1 * dgdH1 + (t1 - T1) * T1 / t1 * dgdT1 + (h2 - H2) * H2 / h2 * dgdH2 + \
                (t2 - T2) * T2 / t2 * dgdT2

        return scon2


    def subcon2grad(a):
        dgdH1 = dgN2dA[0]
        dgdT1 = dgN2dA[1]
        dgdH2 = dgN2dA[2]
        dgdT2 = dgN2dA[3]
        h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
        return (H1 / h1) ** 2 * dgdH1 * np.array([1.0, 0.0, 0.0, 0.0]) \
               + (T1 / t1) ** 2 * dgdT1 * np.array([0.0, 1.0, 0.0, 0.0]) \
               + (H2 / h2) ** 2 * dgdH2 * np.array([0.0, 0.0, 1.0, 0.0]) \
               + (T2 / t2) ** 2 * dgdT2 * np.array([0.0, 0.0, 0.0, 1.0])


    def subcon3(a):
        dgdH1 = dgC1dA[0]
        dgdT1 = dgC1dA[1]
        dgdH2 = dgC1dA[2]
        dgdT2 = dgC1dA[3]
        h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
        scon3 = g3 + (h1 - H1) * H1 / h1 * dgdH1 + (t1 - T1) * T1 / t1 * dgdT1 + (h2 - H2) * H2 / h2 * dgdH2 + \
                (t2 - T2) * T2 / t2 * dgdT2
        return scon3


    def subcon3grad(a):
        dgdH1 = dgC1dA[0]
        dgdT1 = dgC1dA[1]
        dgdH2 = dgC1dA[2]
        dgdT2 = dgC1dA[3]
        h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
        return (H1 / h1) ** 2 * dgdH1 * np.array([1.0, 0.0, 0.0, 0.0]) \
               + (T1 / t1) ** 2 * dgdT1 * np.array([0.0, 1.0, 0.0, 0.0]) \
               + (H2 / h2) ** 2 * dgdH2 * np.array([0.0, 0.0, 1.0, 0.0]) \
               + (T2 / t2) ** 2 * dgdT2 * np.array([0.0, 0.0, 0.0, 1.0])


    def subcon4(a):
        dgdH1 = dgC2dA[0]
        dgdT1 = dgC2dA[1]
        dgdH2 = dgC2dA[2]
        dgdT2 = dgC2dA[3]
        h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
        scon4 = g4 + (h1 - H1) * H1 / h1 * dgdH1 + (t1 - T1) * T1 / t1 * dgdT1 + (h2 - H2) * H2 / h2 * dgdH2 + \
                (t2 - T2) * T2 / t2 * dgdT2

        return scon4


    def subcon4grad(a):
        dgdH1 = dgC2dA[0]
        dgdT1 = dgC2dA[1]
        dgdH2 = dgC2dA[2]
        dgdT2 = dgC2dA[3]
        h1, t1, h2, t2 = a[0], a[1], a[2], a[3]
        return (H1 / h1) ** 2 * dgdH1 * np.array([1.0, 0.0, 0.0, 0.0]) \
               + (T1 / t1) ** 2 * dgdT1 * np.array([0.0, 1.0, 0.0, 0.0]) \
               + (H2 / h2) ** 2 * dgdH2 * np.array([0.0, 0.0, 1.0, 0.0]) \
               + (T2 / t2) ** 2 * dgdT2 * np.array([0.0, 0.0, 0.0, 1.0])


    if not FLAG1:
        con1 = {'type': 'ineq',
                'fun': subcon1,
                'jac': subcon1grad}
    if not FLAG2:
        con2 = {'type': 'ineq',
                'fun': subcon2,
                'jac': subcon2grad}
    if not FLAG3:
        con3 = {'type': 'ineq',
                'fun': subcon3,
                'jac': subcon3grad}
    if not FLAG4:
        con4 = {'type': 'ineq',
                'fun': subcon4,
                'jac': subcon4grad}

    if not FLAG1:
        if not FLAG2:
            if not FLAG3:
                if not FLAG4:
                    const = [con1, con2, con3, con4]
                else:
                    const = [con1, con2, con3]
            else:
                if not FLAG4:
                    const = [con1, con2, con4]
                else:
                    const = [con1, con2]
        else:
            if not FLAG3:
                if not FLAG4:
                    const = [con1, con3, con4]
                else:
                    const = [con1, con3]
            else:
                if not FLAG4:
                    const = [con1, con4]
                else:
                    const = [con1]
    else:
        if not FLAG2:
            if not FLAG3:
                if not FLAG4:
                    const = [con2, con3, con4]
                else:
                    const = [con2, con3]
            else:
                if not FLAG4:
                    const = [con2, con4]
                else:
                    const = [con2]
        else:
            if not FLAG3:
                if not FLAG4:
                    const = [con3, con4]
                else:
                    const = [con3]
            else:
                if not FLAG4:
                    const = [con4]

    max1 = np.amin([0.3, Ap[0] * 1.3])
    min1 = np.amax([0.0, Ap[0] * 0.7])
    max2 = 1.3 * Ap[1]
    min2 = np.amax([0.0, 0.7 * Ap[1]])
    max3 = np.amin([0.3, Ap[2] * 1.3])
    min3 = np.amax([0.0, Ap[2] * 0.7])
    max4 = 1.3 * Ap[3]
    min4 = np.amax([0.0, 0.7 * Ap[3]])
    bounds = [[min1, max1], [min2, max2], [min3, max3], [min4, max4]]
    result = minimize(subobj, Ap, constraints=const, jac=subobjgrad,
                      bounds=bounds,
                      method='slsqp')
    p = p + 1
    obj[p, :] = result.x

print(' ')
print('Resultado da otimização')
print('-----------------------')
print(obj[p])
A = obj[p]
print(' ')
print('-----------------------')

print('Vol. = ', objfun(A))
print('qt = ', disp(A))
print('A = ', A)
print('gN1 = ', conGNL1(A))
print('gN2 = ', conGNL2(A))
print('gCR1 = ', conGCR1(A))
print('gCR2 = ', conGCR2(A))
print(' ')
print('-----------------------')

# Se passar do máximo de iterações, ele não convergiu
if p == 49:
    print('NÃO CONVERGIU:')
print('n loops = ', p)

# Vamos aos plots da análise
h_1 = np.zeros((10,1))
t_1 = np.zeros((10,1))
h_2 = np.zeros((10,1))
t_2 = np.zeros((10,1))
vol = np.zeros((10, 1))
Q1 = np.zeros((10, 1))
Q2 = np.zeros((10, 1))
Q3 = np.zeros((10, 1))
Q4 = np.zeros((10, 1))
vol[0] = objfun(A0)
# As restrições calculadas foram negativas
Q1[0] = -(conGNL1(A0))
Q2[0] = -(conGNL2(A0))
Q3[0] = -(conGCR1(A0))
Q4[0] = -(conGCR2(A0))
h_1[0], t_1[0], h_2[0], t_2[0] = A0[0], A0[1], A0[2], A0[3]
for i in range(k):
    h_1[i+1], t_1[i+1], h_2[i+1], t_2[i+1] = X[i]
    vol[i+1] = V[i]
    Q1[i+1] = -(g1n[i])
    Q2[i+1] = -(g2n[i])
    Q3[i+1] = -(g1c[i])
    Q4[i+1] = -(g2c[i])

fig, ax = plt.subplots(2,2)
fig.suptitle("Variação das restrições", fontsize=16)

ax[0,0].plot(Q1, marker='o', color = 'yellow', ls='None')
ax[0,0].plot(G1, marker='o', color = 'crimson', ls='None')
ax[0,0].tick_params(axis="y", labelsize=10)
ax[0,0].tick_params(axis="x", labelsize=10)
ax[0,0].set_ylabel('$g^1_N$', fontsize=12)
ax[0,0].set_xlabel('$n_{iter.}$', fontsize=12)
ax[0,0].legend(['SLSQP', 'OSA'])
ax[0,1].plot(Q2, marker='o', color = 'yellow', ls='None')
ax[0,1].plot(G2, marker='o', color = 'crimson', ls='None')
ax[0,1].tick_params(axis="y", labelsize=10)
ax[0,1].tick_params(axis="x", labelsize=10)
ax[0,1].set_ylabel('$g^2_{N}$', fontsize=12)
ax[0,1].set_xlabel('$n_{iter.}$', fontsize=12)
ax[0,1].legend(['SLSQP', 'OSA'])
ax[1,0].plot(Q3, marker='o', color = 'yellow', ls='None')
ax[1,0].plot(G3, marker='o', color = 'crimson', ls='None')
ax[1,0].tick_params(axis="y", labelsize=10)
ax[1,0].tick_params(axis="x", labelsize=10)
ax[1,0].set_ylabel('$g^1_{CR}$', fontsize=12)
ax[1,0].set_xlabel('$n_{iter.}$', fontsize=12)
ax[1,0].legend(['SLSQP', 'OSA'])
ax[1,1].plot(Q4, marker='o', color = 'yellow', ls='None')
ax[1,1].plot(G4, marker='o', color = 'crimson', ls='None')
ax[1,1].tick_params(axis="y", labelsize=10)
ax[1,1].tick_params(axis="x", labelsize=10)
ax[1,1].set_ylabel('$g^{2}_{CR}$', fontsize=12)
ax[1,1].set_xlabel('$n_{iter.}$', fontsize=12)
ax[1,1].legend(['SLSQP', 'OSA'])
ax[0,0].grid()
ax[1,0].grid()
ax[0,1].grid()
ax[1,1].grid()
fig.tight_layout()
fig.subplots_adjust(top=0.88)
plt.show()
fig, ax = plt.subplots(1)
fig.suptitle("Variação do volume", fontsize=16)

plt.plot(v[0:5], marker='o', color = 'crimson', ls='None')
plt.plot(vol, marker='o', color = 'yellow', ls='None')
plt.xlabel('$n_{iter.}$', fontsize=12)
plt.ylabel('$Vol.$', fontsize=12)
fig.tight_layout()
fig.subplots_adjust(top=0.88)
plt.legend(['OSA', 'SLSQP'])
plt.grid()
plt.show()
