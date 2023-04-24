import numpy as np
import matplotlib.pyplot as plt

#Default values (to get a Lennard Jones potential close to the project's Figure 1)
Eo = 0.45
sigma = 0.45

def lennardJones(r):
        lj = (4*Eo*((sigma/r)**12 - (sigma/r)**6))
        return lj 
 
x_min, x_max = (0.4, 1.25)
x = np.linspace(x_min, x_max, 1000)
y = lennardJones(x)

def derivative(f,r):
        deltar = 0.001
        return (f(r+deltar) - f(r))/deltar

def approximation(f, xRange, nWells):
        df = derivative(f, xRange)
        intervals = np.array_split(xRange, nWells)
        #print(intervals)
        well = []
        for j in range(nWells):
                interval = intervals[j]
                wellLeftmost = interval[0]
                wellRighmost = interval[len(interval) - 1]
                wellMinVal = min(f(interval))
                print("minVal", wellMinVal)

                print("left", wellLeftmost,"right", wellRighmost, "min", wellMinVal)
                for i in interval:
                        if i == wellLeftmost:
                                well.append(f(wellLeftmost))
                        elif (i == wellRighmost):
                                well.append(f(wellRighmost))
                                #print("derivative was at", i, "=", derivative(lennardJones, i))
                        else :
                                well.append(wellMinVal)
        return well
approx1 = approximation(lennardJones, x, 1)
approx2 = approximation(lennardJones, x, 10)
#print(approx)

plt.plot(x, y, color = "green",label = "Potentiel de Lennard-Jones")
plt.plot(x, approx1, color = "blue", linestyle = "dotted", label = "ordre 1")
plt.plot(x, approx2, color = "red", linestyle = "dashed", label = "ordre 10")
#plt.plot(x, derivative(lennardJones,x))
plt.axhline(y = 0, color = 'black')

plt.legend()
plt.ylabel("V(r) [eV]")
plt.xlabel("r [nm]")
plt.ylim(top = 1.5, bottom = -1)
plt.xlim(left = 0.1, right = 1.5)
plt.show()

#poubelle
#or (np.abs(derivative(lennardJones, i)) <= 10**(-2)