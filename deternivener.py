import numpy as np
import math

from approximation import approx
from lennardJones import lennardJones
from scipy.signal import argrelextrema
import scipy.constants as cst

#Retourne le K**2 associé à chaque puits, et une matrice B contenant sa position (x)
def kVal(E, approx_y, x):
    m = cst.electron_mass
    
    #print(222, approx_y)
    
    E0 = np.array([])
    E0 = np.append(E0, approx_y[0])             #Potentiel minimum dans la région courante
    kSquared = np.array([], dtype = complex)    #Matrice des K carrés
    B = np.array([])                            #position de début / fin de puits (x)
    
    j = 1                                       #Nombre du puits (ou nombre du E0)
    
    for i in range(len(approx_y) - 1):

        if approx_y[i] != approx_y[i+1]:
            E0 = np.append(E0, approx_y[i+1])
            
            print(111, E0, approx_y[i+1])
            
            kCurrent = (2*m*(E - E0[j-1])) / (cst.hbar**2) #Calcul du coefficient K**2 par sa définition
            j += 1
            kSquared = np.append(kSquared, kCurrent)
            B = np.append(B, x[i])

    return kSquared, B

def Mx(k, b):
    print(k, b)
    
    M = np.array([[math.exp(k*b), math.exp(-k*b)],[k*math.exp(k*b), -k*math.exp(k*b)]], dtype = complex)
    return M

def M(n, E, Eo, sigma):
    x_min, x_max = (0.4, 1.25)
    x = np.linspace(x_min, x_max, 1000)
    approx_y = approx(lennardJones, x_min, x_max, 1, Eo, sigma)
    K, c = kVal(E, approx_y, x)
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
    n = 1
    Eo = 0.45*1.602176565**(-19)
    sigma = 0.45**(-9)
    
    return nbener(n, Eo, sigma)

def q4():
    n = 1
    sigma = np.arange(0.1, 1.05, 0.05)
    E0 = np.arange(1, 10.5, 0.5)
    soln = np.zeros((len(sigma), len(E0)))#matrice taille (i,j)
    sole = np.zeros((len(sigma), len(E0)))#matrice taille (i,j)
    a = b = 0

    for i in sigma:
        for j in E0:
            s, e = nbener(n, i, j) #plus bas e
            soln[b] [a] = s
            sole[b] [a] = np.min(e)
            a+=1
        b+=1
    
    #reponse sous plot(brouillon), heat map une n et autre e
    
q3()
