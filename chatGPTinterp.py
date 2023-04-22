import numpy as np
import matplotlib.pyplot as plt
from lenardJones import lennardJones

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


#Setting lennardJones as the function to be piecewised
f = lambda x: lennardJones(x)
x_range = (0.01, 1.3)
num_pieces = 25
x, y, approx_y = piecewise_approximation(f, x_range, num_pieces)

plt.plot(x, y, label='true function')
plt.plot(x, approx_y, label=f'{num_pieces} piece approximation')
plt.ylim(bottom = -1, top = 1)
plt.xlim(left = 0, right = 2)
plt.legend()
plt.show()
