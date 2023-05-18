import numpy as np
import math

from approximation import lennardJones, approx
from scipy.signal import argrelextrema
import scipy.constants as cst


def k(E, approx_y, x):
    m=cst.electron_mass
    print(222, approx_y)
    E0 = np.array([])
    kcar = np.array([], dtype=complex)
    b = np.array([])#position changement
    E0 = np.append(E0, approx_y[0])
    j = 1
    a=0
    for i in range(len(approx_y)-1):
        a+=1
        if approx_y[i] != approx_y[i+1]:
            E0 = np.append(E0, approx_y[i+1])
            print(111, E0, approx_y[i+1])
            kt = 2*m*(-E0[j-1]+E)/cst.hbar**2
            j+=1
            kcar= np.append(kcar, kt)
            b= np.append(b, x[i])
            a=0
    return kcar, b
def Mx(k, b):
    print(k, b)
    M = np.array([[math.exp(k*b), math.exp(-k*b)],[k*math.exp(k*b), -k*math.exp(k*b)]], dtype = complex)
    return M
def M(n, E, Eo, sigma):
    x_min, x_max = (0.4, 1.25)
    x = np.linspace(x_min, x_max, 1000)
    approx_y = approx(lennardJones, x_min, x_max, 1, Eo, sigma)
    K, c = k(E, approx_y, x)
    M = np.zeros((2,2))
    print(K, c, n)
    for i in range(n):
        Mi = Mx(K[n-i-1], c[n-i-1])
        if (i == 0):
            k1 = K[i]
            M = Mi
        else:
            M = np.linalg.inv(Mi)*M
    return M, k1, x[0]


def nbener(n, Eo, sigma):
    Etest = np.linspace(-Eo, 0, 1000)
    s = 0
    Mt = np.array([])
    for e in Etest:
        m, k1, a = M(n, e, Eo, sigma)
        Mt = np.append(Mt, m[0][0]-m[0][1]*(1/math.tan(k1*a)))
    e = Mt[argrelextrema(Mt, np.less)[0]]
    s = len(e)
    return s, e
def q3():
    n=1
    Eo = 0.45* 1.602176565e-19
    sigma = 0.45e-9
    return nbener(n, Eo, sigma)
def q4():
    n = 1
    sigma = np.arange(0.1, 1.05, 0.05)
    E0 = np.arange(1, 10.5, 0.5)
    soln = np.zeros((len(sigma), len(E0)))#matrice taille (i,j)
    sole = np.zeros((len(sigma), len(E0)))#matrice taille (i,j)
    a=b=0
    for i in sigma:
        for j in E0:
            s, e = nbener(n, i, j) #plus bas e
            soln[b] [a] = s
            sole[b] [a] = np.min(e)
            a+=1
        b+=1
    
    #reponse sous plot(brouillon), heat map une n et autre e
    
q3()
