import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

#PARAMETERS
Eo = 0.45
sigma = 0.45

def lennardJones(r, Eo, sigma):
        lj = (4*Eo*((sigma/r)**12 - (sigma/r)**6))
        return lj 

def approx(arbitraryFunction, xmin, xmax, nWells, Eo, sigma):
    nWells+=1
    xSpace = np.linspace(xmin, xmax, 1000, retstep=True)    #Set of points whithin which the function should be welled, with a set precision
    step = xSpace[1]                                        #linspace gives the step thanks to the retstep argument
    functionMinimum = min(arbitraryFunction(xSpace[0], Eo, sigma))     #The minimum of the function within the desired space
    
    approx = []                                             #Empty list that will be filled with the calculations

    #Offsets the entire arbitrary function to avoid negative integrals (so the function is always positive)
    def f(x):                                               
         return arbitraryFunction(x, Eo, sigma) + np.abs(functionMinimum)
    
    #Dividing the total area calculated with an integral by the number of wells ( that should have an equal area )
    area = np.abs(integrate.quad(f, xmin, xmax))
    wellArea = area[0]/nWells
    
    #Setting the well interval
    tempMin = xmin
    i = tempMin

    #Tracking the current well and point (debug purposes)
    wellNb = 0 
    pointCount = 0
    
    #Giving the approximation an initial value at the leftmost point
    approx.append(f(xmin))

    #Well calculation
    while i < xmax:

        localArea = integrate.quad(f, tempMin, i)
        #print("i =", i, "tempMin", tempMin, "xmax", xmax, "localArea", localArea, "wellArea", wellArea, "wellNo", wellNb, "point", pointCount)
        
        if(np.abs(localArea[0]) >= np.abs(wellArea)):
            #print("localArea = ", localArea[0],">= wellArea =",wellArea, "for x range:(", xmin, ";", i, ")")
            j = tempMin
            localInterval = np.arange(j,i,step)
            localMin = min(f(localInterval))
            
            while j < i:
                    approx.append(localMin)
                    j = j+step
            tempMin = i
            wellNb += 1
            
            i = i+step

        #End of interval case
        elif(np.abs(i - xmax) <= 10**-3):
                
                j = tempMin
                
                while j < xmax:
                    approx.append(f(tempMin))
                    j = j+step
                
                break
        else:
                i = i + step
        pointCount += 1

    #reverse the initial offset
    return approx - np.abs(functionMinimum)

#EXECUTION
x_min, x_max = (0.4, 1.25)
x = np.linspace(x_min, x_max, 1000)
y = lennardJones(x, Eo, sigma)

#PLOT
plt.plot(x, y, color = "green",label = "Potentiel de Lennard-Jones")
plt.plot(x, approx(lennardJones, x_min, x_max, 1, Eo, sigma), color = "blue", linestyle = "dotted", label = "ordre 1")
plt.plot(x, approx(lennardJones, x_min, x_max, 10, Eo, sigma), color = "red", linestyle = "dashed", label = "ordre 10")

plt.axhline(y = 0, color = 'black')

plt.legend()
plt.ylabel("V(r) [eV]")
plt.xlabel("r [nm]")
plt.ylim(top = 1.5, bottom = -1)
plt.xlim(left = 0.1, right = 1.5)
plt.show()

plt.legend()
plt.ylabel("V(r) [eV]")
plt.xlabel("r [nm]")
plt.ylim(top = 1.5, bottom = -1)
plt.xlim(left = 0.1, right = 1.5)
plt.show()
