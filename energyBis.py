from lennardJones import lennardJones
from approximation import approx
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
from scipy.signal import argrelextrema
#params:
Eo = 10     #eV
sigma = 1   #nm
n = 1 #wells
E = -5    #eV

#interval:
x_min, x_max = (sigma - 10**-2, sigma*5)

xInterval = np.linspace(x_min, x_max, 1000)

#approx:
yapprox = approx(lennardJones, x_min, x_max, n, Eo, sigma)

#functions:
def intervalBorders(y, x, n): #x is the interval in which y is defined (x values)
    #list storing start and end of interval, first interval should be [borders[0] = (a0; b0)]
    borders = []
    temp = []
    counter = 0

    for i in range(len(y) - 1):
        if(y[i] != y[i+1]):
            counter += 1
            temp.append(x[i])

    for j in range(len(temp) - 1):    
        borders.append((temp[j],temp[j+1]))
    
    if(counter < n+1): #The number of wells can be approximated for large values of n
        for k in range(n+1 - counter):
            borders.append((temp[len(temp) - 2], temp[len(temp) - 1]))
    
    return borders

def intervalV(y, x, interval, n): #V for interval n, in Joules
    borders = intervalBorders(y, x, n)

    if(interval >= n):
        V = 0
        return V

    if(interval < n):
        midWell = (((borders[interval][0] + borders[interval][1]) / 2)) #V[0] for interval 1
       
        for j in range(len(x)):
            if(midWell - x[j] < 10**(-4)):
                V = y[j]*cst.elementary_charge
                break  
    return V

def kValuesSquared(E, V): #JOULES
    kSqrd = (2*cst.electron_mass*(E-V))/(cst.hbar**2)
    #CAN BE NEGATIVE
    return kSqrd 

def matrixAddToProduct(Mi, M): #M
    MInverse = np.linalg.inv(M)
    return np.matmul(MInverse,Mi)
 
def matrixMcalculation(k, index, pivot, n, x, y): #I is the interval, x is xRange, y is approx
    #print("matrCalc",n)
    if(index > n - 1):
        index = index - 1
        b = intervalBorders(y, x, n)[index][pivot] *10**(-9) #to m
    else:
        b = intervalBorders(y, x, n)[index][pivot] *10**(-9) #to m
    
    #M contains the general solution to the wave function and its derivative (row 1, 2 respectively)
    #print("B for index", index, "=", b)
    Mi = np.array([[np.exp(k*b), np.exp(-k*b)],[k*np.exp(k*b), -k*np.exp(-k*b)]], dtype = complex)
    return Mi

def energyValues(n, E, x, y): #nb wells, energy, xInterval, approximation 
    E = E*cst.elementary_charge
    #calculate M final (accounting for all wells)
    k = []  #Array containig k values
    Mi = [] #Array containing the matrices for each interval
    Mtot = []

    for i in range(n + 1):
        #print("EVAL -> for interval =", i)
        V = intervalV(y, x, i, n)
        #print("     V =", V, "[J]","E", E,"[J]", V/cst.elementary_charge, E/cst.elementary_charge)
        kSqrd = kValuesSquared(E, V)
        k.append(np.array(kSqrd, dtype = complex)**0.5) #Using numpy to avoid NaN errors
        #print("     k^2 =", kSqrd)
    
    print("K VECTOR", k)
    #MI = 4 for 2 wells !!
    for j in range(n):
        if(j == 0):
            Mi.append(matrixMcalculation(k[j], j, 1, n, x, y)) #Initial matrix
        
        if(n == 1): #Specific case
            Mi.append(matrixMcalculation(k[j+1], j, 1, n, x, y))
            break

        if(j == n - 1): #End of well list case
            Mi.append(matrixMcalculation(k[j+1],j,1,n,x,y))

        else:
            Mi.append(matrixMcalculation(k[j+1], j+1, 0, n, x, y))

    #print("MI", Mi)
    Mtot = Mi[0]
    #for l in range(len(Mi)):
    #    if(len(Mi) != 1):
    #        Mtot = np.matmul(np.linalg.inv(Mi[l]), Minit)
    #        Minit = Mtot
    
    Mtot = np.matmul(Mtot, Mi[-1])
    
    #print(Mtot[0][0], Mtot[0][1])
    
    A = Mtot[0][0]
    B = Mtot[0][1]
    return A,B

def energyTester(Espace,x,y):
    
    V = intervalV(y,x,0,n)
    a = intervalBorders(y,x,n)[0][1]
    sol = np.empty((0))
    counter = 0

    for Etest in Espace:    
        A, B = energyValues(n,Etest, x, y)
        k = np.array(kValuesSquared(Etest, V), dtype = complex)**0.5
        ka = k*(a*10**(-9))
        print("iter nb", counter,"k", k, "ka", ka, "a",a*10**-9,"v",V)
        tempE = A*np.cos(ka) - B*np.sin(ka)
        print("ENERGY:",Etest, "tempE",tempE)
        sol = np.append(sol, tempE)
        counter += 1

    print(sol)
    states = Espace[argrelextrema(sol, np.less)[0]]
    print("Energies",states/cst.elementary_charge, "[eV]")
    plt.plot(Espace, np.abs(sol), color= "blue")
    #print("Efound", Efound)

#test
borders = intervalBorders(yapprox, xInterval, n)

#energyValues(n, E, xInterval , yapprox)  #The Eo parameter is in the approx function

energyTester(np.linspace(-Eo*cst.elementary_charge,0,1000),xInterval,yapprox)
print("borders",len(borders), borders[0])

#plt.axhline(y = E, color = "purple", label = "E")
#plt.plot(xInterval, yapprox, color = 'red', linestyle = "dotted")
#plt.plot(xInterval, lennardJones(xInterval, Eo, sigma), color = 'green')
#plt.legend()
#plt.ylabel("V(r) [eV]")
#plt.xlabel("r [nm]")
plt.show()