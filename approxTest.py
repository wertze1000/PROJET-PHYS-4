import numpy as np
import matplotlib.pyplot as plt
from lennardJones import lennardJones
from approximation import approx
import scipy.constants as cst

#PARAMETERS
n = 1
Eo = 10 #eV
sigma = 0.15 #nm

#EXECUTION
x_min, x_max = (sigma - 10**-2, 1)
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

plt.show()