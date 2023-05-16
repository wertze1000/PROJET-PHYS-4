# -*- coding: utf-8 -*-
"""
Created on Mon May 15 17:21:00 2023

@author: User
"""
import numpy as np
import math

from lennardJones import lennardJones

import scipy.constants as cst

def piecewise_approximation(f, x_range, num_pieces):
    x_min, x_max = x_range
    x = np.linspace(x_min, x_max, 1000)
    y = f(x)
    approx_y = np.zeros_like(y)
    piece_width = (x_max - x_min) / num_pieces

    for i in range(num_pieces):
        x_piece = np.linspace(x_min + i*piece_width, x_min + (i+1)*piece_width, 1000//num_pieces)
        y_piece = f(x_piece)
        approx_y_piece = np.full_like(y_piece, np.mean(y_piece))
        indices = np.where((x >= x_min + i*piece_width) & (x <= x_min + (i+1)*piece_width))
        approx_y[indices] = approx_y_piece
    return x, y, approx_y

m=cst.electron_mass
f = lambda x: lennardJones(x)
x_range = (0.01, 1.3)
num_pieces = 4
x, y, approx_y = piecewise_approximation(f, x_range, num_pieces)


def k(E, approx_y):
    E0 = np.array([])
    kcar = np.array([], dtype=complex)
    b = np.array([])#position changement
    E0 = np.append(E0, approx_y[0])
    j = 1
    a=0
    for i in range(len(approx_y)-1):
        a+=1
        if approx_y[i] != approx_y[i+1]:
            j+=1
            E0 = np.append(E0, approx_y[i+1])
            kt = 2*m*(-E0[j-1]+E)/cst.hbar**2
            kcar= np.append(kcar, kt)
            b= np.append(b, i)
            a=0
    return kcar, b
def Mx(k, b):
    M = np.array([[math.exp(k*b), math.exp(-k*b)],[k*math.exp(k*b), -k*math.exp(k*b)]], dtype = complex)
    return M
def M(n, E):
    f = lambda x: lennardJones(x)
    x_range = (0.01, 1.3)
    x, y, approx_y = piecewise_approximation(f, x_range, n)
    K, c = k(E, approx_y)
    M = np.zeros((2,2))
    for i in range(n):
        Mi = Mx(K[n-i], c[n-i])
        if (i == 0):
            M = Mi
        else:
            M = np.linalg.inv(Mi)*M
    return M   
