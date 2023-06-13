import numpy as np

from approximation import approx
from lennardJones import lennardJones
from scipy.signal import argrelextrema
import scipy.constants as cst
import cmath
import matplotlib.pyplot as plt

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
            kCurrent = (2*m*(E - E0[j-1])) / (cst.hbar**2) #Calcul du coefficient K**2 par sa définition
            j += 1
            kSquared = np.append(kSquared, kCurrent)
            B = np.append(B, x[i])
    return kSquared, B

def Mx(k_square, b):# calcul de une matrice Mi, avec un K donnée
    k = np.sqrt(k_square)
    M = np.array([[cmath.exp(k*b), cmath.exp(-k*b)],[k*cmath.exp(k*b), -k*cmath.exp(k*b)]], dtype = complex)
    return M

def M(n, E, Eo, sigma):#calcule de la matrice M finale
    x_min, x_max = (0.4, 1.25)
    x = np.linspace(x_min, x_max, 1000)
    approx_y = approx(lennardJones, x_min, x_max, n, Eo, sigma)
    K, c = kVal(E, approx_y, x)
    M = np.zeros((2,2))
    x=x*10**(-9)#passage en mètres
    approx_y = approx_y*cst.e#passage en joule
    for i in range(n):
        Mi = Mx(K[n-i-1], c[n-i-1])
        if (i == 0):
            k1 = K[i]
            M = Mi
        else:
            M = np.linalg.inv(Mi)*M

    return M, k1, x_min


def nbener(n, Eo, sigma): #calcul du nombre d'états liés et de la valeur d'énergie de ceux-ci
    Etest = np.linspace(-Eo*10**(-9), -1.12264629e-21, 100)
    s = 0
    Mt = np.array([], dtype= complex)
    i=0
    for e in Etest:
        print(i)
        m, k1, a = M(n, e, Eo, sigma)
        Mt = (np.append(Mt,- np.tan(k1*a)*m[0][0]-m[0][1])) #critere de continuité, np.tan(k1*a)
        i+=1
    e = Etest[argrelextrema(Mt, np.less)[0]]
    s = len(e)
    plt.plot(Etest, Mt)
    plt.show()
    
    return s, (e/cst.e)

def q3():
    n = 1
    Eo = 10
    sigma = 1
    
    return nbener(n, Eo, sigma)

def q4():
    n = 1
    sigma = np.arange(0.1, 1.05, 0.05)#en nm
    E0 = np.arange(1, 10.5, 0.5)#en ev
    sigma = sigma*10**(-9)#en m
    E0 = E0*cst.e#en j
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
    
print(q3())
