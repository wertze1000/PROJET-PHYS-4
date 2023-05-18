import numpy as np
import matplotlib.pyplot as plt
from lennardJones import lennardJones
from approximation import approx

#PARAMETERS
Eo = 0.45
sigma = 0.45

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