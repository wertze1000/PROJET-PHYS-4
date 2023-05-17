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
