from scipy.optimize import fsolve
import time 
from deternivener import nbener

#PARAMETERS:
Eo = 7      #eV
sigma = 0.7 #nm

#FUNCTIONS:
def epsilon(E1, n): #RELATIVE ERROR TO FIRST ALLOWED ENERGY LEVEL
    return (E1(n) - E1(100))/E1(100)

def t(n): #elapsed time for execution of nbener (Schrödinger 1D) for a given n
    start = time.time()
    nbener(n, Eo, sigma) #function deternivener for given n
    end = time.time()

    return end - start 

def T(n): #comparison of elapsed time to case n = 100
    return t(n)/t(100)

def p(n): #function describing convergence (to optimize)
    return 4*epsilon(n) + T(n)

def optimization(p):
    nopt = fsolve(p)
    return nopt