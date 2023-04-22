import numpy as np
import matplotlib.pyplot as plt

#Default values (to get a Lennard Jones potential close to the project's Figure 1)
Eo = 0.45
sigma = 0.45

def lennardJones(r):
        lj = (4*Eo*((sigma/r)**12 - (sigma/r)**6))
        return lj 
 
""" x_min, x_max = (0.01, 2)
x = np.linspace(x_min, x_max, 1000)
y = lennardJones(x)


plt.plot(x, y)
plt.axhline(y = 0, color = 'black')
plt.ylabel("V(r) [eV]")
plt.xlabel("r [nm]")
plt.ylim(top = 1, bottom = -1)
plt.xlim(left = 0, right = 2)
plt.show()
 """