# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 19:39:02 2023

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt
import math as m
#Default values (to get a Lennard Jones potential close to the project's Figure 1)
 

x = np.linspace(0.00001, 20, 1000)
y = np.linspace(0.0001, 20, 1000)
for i in range(1000):
    y[i] = -x[i]*m.tan(x[i])**(-1)
    
theta = np.linspace( 0 , 2 * np.pi , 1000 )
 
radius = 5
 
a = radius * np.cos( theta )
b = radius * np.sin( theta )
radius = 10
c = radius * np.cos( theta )
d = radius * np.sin( theta )
radius = 2.5
e = radius * np.cos( theta )
f = radius * np.sin( theta )
plt.plot(x, y)
plt.plot(a, b, label = "$x^2 +y^2 = R^2$ pour R = 5")
plt.plot(c, d, label = "$x^2 +y^2 = R^2$ pour R = 10")
plt.plot(e, f, label = "$x^2 +y^2 = R^2$ pour R = 2.5")
plt.axhline(y = 0, label="-x cotan(x)")
plt.ylabel("y = $K_2b$")
plt.xlabel("x = kb")
plt.ylim(top = 14, bottom = 0)
plt.xlim(left = 0, right = 16)
plt.gca().set_aspect('equal')
plt.legend(loc="upper right")
plt.show()