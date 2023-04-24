import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative

def piecewise_approximation(f, x_range, num_pieces):
    x_min, x_max = x_range
    x = np.linspace(x_min, x_max, 1000)
    y = f(x)
    approx_y = np.zeros_like(y)
    derivative_f = np.abs(derivative(f, x))

    for i in range(num_pieces):
        piece_width = derivative_f[int(i*len(x)/num_pieces):int((i+1)*len(x)/num_pieces)].max()
        x_piece = np.linspace(x_min + i*piece_width, x_min + (i+1)*piece_width, 1000//num_pieces + 1)[:-1]
        y_piece = f(x_piece)
        approx_y_piece = np.full_like(y_piece, np.mean(y_piece))
        indices = np.where((x >= x_min + i*piece_width) & (x <= x_min + (i+1)*piece_width))
        approx_y[indices] = approx_y_piece

    x = np.linspace(x_min, x_max, len(approx_y))

    return x, y, approx_y

# Setting lennardJones as the function to be piecewised
def lennardJones(x):
    Eo = 0.45
    sigma = 0.45
    return (4*Eo*((sigma/x)**12 - (sigma/x)**6))

x_range = (0.01, 1.3)
num_pieces = 2

x, y, approx_y = piecewise_approximation(lennardJones, x_range, num_pieces)

plt.plot(x, y, label='True function')
plt.plot(x, approx_y, label='Piecewise approximation')
plt.xlim(left = 0, right = 1.5)
plt.ylim(bottom = -0.6, top = 0.6)
plt.legend()
plt.show()
