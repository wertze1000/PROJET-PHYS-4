import numpy as np
import math

from lennardJones import lennardJones
from scipy.signal import argrelextrema
import scipy.constants as cst

def piecewise_approximation(f, x_range, num_pieces):
    x_min, x_max = x_range
    x = np.linspace(x_min, x_max, 1000)
    y = f(x)
    approx_y = np.zeros_like(y)
    piece_width = (x_max - x_min) / num_pieces
    
    start = 0
    for i in range(num_pieces):
        x_piece = np.linspace(x_min + i*piece_width, x_min + (i+1)*piece_width, 1000//num_pieces)
        y_piece = f(x_piece)
        approx_y_piece = np.full_like(y_piece, np.mean(y_piece))
        indices = np.where((x >= x_min + i*piece_width) & (x <= x_min + (i+1)*piece_width))
        approx_y[indices] = approx_y_piece
    return x, y, approx_y, start
def k(E, approx_y):
    m=cst.electron_mass
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
def M(n, E, Eo, sigma):
    f = lambda x: lennardJones(x, Eo, sigma)
    x_range = (0.01, 1.3)
    x, y, approx_y, a = piecewise_approximation(f, x_range, n)
    K, c = k(E, approx_y)
    M = np.zeros((2,2))
    for i in range(n):
        Mi = Mx(K[n-i], c[n-i])
        if (i == 0):
            k1 = K[i]
            M = Mi
        else:
            M = np.linalg.inv(Mi)*M
    return M, k1, a


def nbener(n, Eo, sigma):
    Etest = np.linspace(-Eo, 0, 1000)
    s = 0
    Mt = np.array([])
    for e in Etest:
        m, k1, a = M((n, e, Eo, sigma))
        Mt = np.append(Mt, m[0][0]+m[0][1]*(1/math.tan(k1*a)))
    e = Mt[argrelextrema(Mt, np.less)[0]]
    s = len(e)
    return s, e
def q3():
    n=1
    Eo = 0.45
    sigma = 0.45
    return nbener(n, Eo, sigma)
def q4():
    n = 1
    sigma = np.arrange(0.1, 1, 0.05)
    E0 = np.arrage(1, 10, 0.5)
    soln = np.zeros((len(sigma), len(E0)))#matrice taille (i,j)
    sole = np.zeros((len(sigma), len(E0)))#matrice taille (i,j)
    a=b=0
    for i in sigma:
        for j in E0:
            s, e = nbener(n, i, j) #plus bas e
            soln[b] [a] = s
            sole[b] [a] = np.min(e)
            a+=1
        b+=1         
    #reponse sous plot(brouillon), heat map une n et autre e
