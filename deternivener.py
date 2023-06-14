import numpy as np

from approximation import approx
from lennardJones import lennardJones
from scipy.signal import argrelextrema
from heatMap import generate_heatmap
import scipy.constants as cst
import cmath
import matplotlib.pyplot as plt

#Retourne le K**2 associé à chaque puits, et une matrice B contenant sa position (x)
def kVal(E, approx_y, x, n):
    m = cst.electron_mass
    
    #print(222, approx_y)
    E0 = np.array([])
    approx_y[0]=approx_y[1]
    E0 = np.append(E0, approx_y[0])             #Potentiel minimum dans la région courante
    kSquared = np.array([], dtype = complex)    #Matrice des K carrés
    B = np.array([])                            #position de début / fin de puits (x)
    B = np.append(B, x[0])

    j = 1                                       #Nombre du puits (ou nombre du E0)
    
    for i in range(len(approx_y) - 1):
        if approx_y[i] != approx_y[i+1]:
            E0 = np.append(E0, approx_y[i+1])
            kCurrent = (2*m*(-E + E0[j-1]*cst.e) )/ (cst.hbar**2) #Calcul du coefficient K**2 par sa définition
            j += 1
            kSquared = np.append(kSquared, kCurrent)
            B = np.append(B, x[i])
    
    while(j < n+1): #Dans le cas ou les derniers puits sont négligeables (n > 50)
        B = np.append(B, B[len(B)-1])
        kSquared = np.append(kSquared, kCurrent)
        j += 1
    
    kCurrent = (2*m*(-E + E0[-1]*cst.e) )/ (cst.hbar**2)
    kSquared = np.append(kSquared, kCurrent)
    #print(B, kSquared, 111)
    return kSquared, B

def Mx(k_square, b):# calcul de une matrice Mi, avec un K donnée
    k = np.sqrt(k_square)
    b = b*10**(-9)
    M = np.array([[cmath.exp(k*b), cmath.exp(-k*b)],[k*cmath.exp(k*b), -k*cmath.exp(-k*b)]], dtype = complex)
    #print(M)
    return M

def M(n, E, Eo, sigma):#calcule de la matrice M finale
    x_min, x_max = (sigma-10**-2, sigma*5)
    x = np.linspace(x_min, x_max, 1000)
    approx_y = approx(lennardJones, x_min, x_max, n, Eo, sigma)
    K, b = kVal(E, approx_y, x, n)
    M = np.identity(2, dtype = complex)
    #print(x)
    #x=x*10**(-9)#passage en mètres
    #approx_y = approx_y*cst.e#passage en joule
    k1 = K[0]
    for i in range(n):
        #print(i, n,K, b)
        Mi1 = Mx(K[i], b[i+1])
        Mi2 = Mx(K[i+1], b[i+1])
        #print(b[n-i])
        
        M = np.matmul(M,np.linalg.inv(Mi1))
        M = np.matmul(M,Mi2)
    return M, k1, x_min


def nbEner(n, Eo, sigma): #calcul du nombre d'états liés et de la valeur d'énergie de ceux-ci
    Etest = np.linspace(-Eo*cst.e, 0, 100)
    s = 0
    Mt = np.array([], dtype= complex)
    i=0
    for e in Etest:
        #print(i)
        m, k1, a = M(n, e, Eo, sigma)
        #print( -np.tan(k1*a)**-1*m[0][1])
        Mt = (np.append(Mt,m[0][0] - np.tan(k1*a)**-1*m[0][1])) #critere de continuité, np.tan(k1*a)
        i+=1
    e = Etest[argrelextrema(Mt, np.less)[0]]
    nombreEnergies = len(e)
    Energies = -e/cst.e
    #plt.plot(Etest/cst.e, Mt)
    #plt.ylim(top = 1, bottom = -0.5)
    #plt.show()
    
    return nombreEnergies, Energies

def q3():
    n = 100
    Eo = 10
    sigma = 1
    
    return nbEner(n, Eo, sigma)

def q4():
    n = 1
    sigma = np.arange(0.1, 1.05, 0.05)#en nm
    E0 = np.arange(1, 10.5, 0.5)#en ev
     #sigma = sigma*10**(-9)#en m
    #E0 = E0*cst.e#en j
    soln = np.zeros((len(sigma), len(E0)))#matrice taille (i,j)
    sole = np.zeros((len(sigma), len(E0)))#matrice taille (i,j)
    a = 0
    b = 0

    for i in sigma:
        for j in E0:
            print(j)
            s, e = nbEner(n, i, j) #plus bas e
            soln[a] [b] = s
            sole[a] [b] = np.min(e)
            a+=1
        b+=1
        a = 0
    
    generate_heatmap(E0, sigma, sole)
    generate_heatmap(E0, sigma, soln)

    print(soln)
    #reponse sous plot(brouillon), heat map une n et autre e
    
#print(q3())
#q4()