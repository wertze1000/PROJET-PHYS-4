from scipy.optimize import fsolve
import time 
from deternivener import nbEner
import numpy as np

#PARAMETERS:
Eo = 7      #eV
sigma = 0.7 #nm

#Premier niveau d'énergie(pour Eo = 7, sigma = 0,7)
E100 = min(nbEner(100,Eo,sigma)[1]) #[eV]
print("Pour n=100", E100)

#FUNCTIONS:
def epsilon(n, E100): #RELATIVE ERROR TO FIRST ALLOWED ENERGY LEVEL
    EN = min(nbEner(n,Eo,sigma)[1])
    print("Pour n=",n,EN)
    return (EN - E100)/E100

def t(n): #elapsed time for execution of nbener (Schrödinger 1D) for a given n
    start = time.time()
    nbEner(n, Eo, sigma) #function deternivener for given n
    end = time.time()

    return end - start 

def T(n): #comparison of elapsed time to case n = 100
    return t(n)/t(100)

def p(n): #function describing convergence (to optimize)
    return 4*np.abs(epsilon(n)) + np.abs(T(n))

def optimization(p):
    nopt = fsolve(p)
    return nopt

#test:
n=1

for i in range(1,10,1):
    eps = epsilon(i,E100)
    print("Erreur", eps)