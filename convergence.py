from scipy.optimize import fsolve
import time 
from deternivener import nbEner
import numpy as np
import matplotlib.pyplot as plt

#PARAMETERS:
Eo = 7      #eV
sigma = 0.7 #nm

#Premier niveau d'énergie(pour Eo = 7, sigma = 0,7)
start = time.time()
E100 = min(nbEner(100,Eo,sigma)[1]) #[eV]
end = time.time()

tE100 = end - start
print("Pour n=100", E100, "time", tE100, "s")

#FUNCTIONS:
def P(n, E100, tE100): 
    start = time.time()
    EN = min(nbEner(n,Eo,sigma)[1])
    end = time.time()
    tn = end - start

    print("Pour n=",n,EN, "time", tn, "s")
    epsilon = (EN - E100)/E100
    t = tn / tE100
    return 4*epsilon + t

nRange = range(1,100,1)
Plist = []

for i in nRange:
    Plist.append(P(i,E100, tE100))
    
plt.plot(nRange, Plist)
plt.ylabel("P(n)")
plt.xlabel("n")
plt.show()