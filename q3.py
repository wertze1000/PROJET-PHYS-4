# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 19:39:02 2023

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
import math as m

 
largeur =0.4754754754755084*10**-9#pris dans fonction approx avec un print [m]
e0 = 5 #[eV]

r = np.sqrt(2*cst.electron_mass*e0*cst.e*largeur**2/cst.hbar**2)

x = np.linspace(0.00001, 10, 1000)
y = np.linspace(0.00001, 10, 1000)
for i in range(1000):
    y[i] = -x[i]/m.tan(x[i])

 
a =[]
b =[]
for i in range(1000):
    if x[i]< r:
        b = np.append(b,np.sqrt(r**2-x[i]**2))
        a = np.append(a,x[i])

xlin = x[:len(b)]
ylin = y[:len(b)]

idx = np.argwhere(np.diff(np.sign(ylin - b))).flatten()

plt.plot(x[idx], y[idx], 'ro')

K1 = largeur**-1*y[idx[0]]
k1 = largeur**-1*x[idx[0]]
K2 = largeur**-1*y[idx[2]]
k2 = largeur**-1*x[idx[2]]
E1 = K1**2/(2*cst.electron_mass)*cst.hbar**2
E2 = K2**2/(2*cst.electron_mass)*cst.hbar**2
e1 = E1/cst.e
e2 = E2/cst.e

plt.plot(x, y)
plt.plot(a, b, label = "$x^2 +y^2 = R^2$")
plt.axhline(y = 0, label="y = -x cotan(x)")
plt.ylabel("y = $Kb$")
plt.xlabel("x = kb")
plt.ylim(top = 8, bottom = 0)
plt.xlim(left = 0, right = 8)
plt.gca().set_aspect('equal')
plt.legend(loc="upper right")
plt.show()