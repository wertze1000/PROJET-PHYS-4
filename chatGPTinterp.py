import numpy as np
import matplotlib.pyplot as plt
from lennardJones import lennardJones

def piecewise_approximation(f, x_range, num_pieces):
    x_min, x_max = x_range
    x = np.linspace(x_min, x_max, 1000)
    y = f(x)
    approx_y = np.zeros_like(y)
    piece_widths = np.full(num_pieces, (x_max - x_min) / num_pieces)

    for iteration in range(5):
        # Adjust piece widths based on local curvature
        for i in range(1, num_pieces):
            x_start = x_min + np.sum(piece_widths[:i])
            x_end = x_start + piece_widths[i]
            y_mid = f((x_start + x_end) / 2)
            dy_dx = (f(x_end) - y_mid) / (x_end - x_start)
            curvature = abs(dy_dx) / (1 + dy_dx**2)**1.5
            piece_widths[i] = (x_max - x_min) / num_pieces / (1 + curvature)

        # Update approx_y based on new piece widths
        approx_y = np.zeros_like(y)
        for i in range(num_pieces):
            x_start = x_min + np.sum(piece_widths[:i])
            x_end = x_start + piece_widths[i]
            indices = np.where((x >= x_start) & (x < x_end))
            x_piece = x[indices]
            y_piece = y[indices]
            approx_y_piece = np.full_like(y_piece, np.mean(y_piece))
            approx_y[indices] = approx_y_piece

    return x, y, approx_y


#Setting lennardJones as the function to be piecewised
f = lambda x: lennardJones(x)
x_range = (0.01, 1.3)
num_pieces = 10
x, y, approx_y = piecewise_approximation(f, x_range, num_pieces)

plt.plot(x, y, label='true function')
plt.plot(x, approx_y, label=f'{num_pieces} piece approximation')
plt.ylim(bottom = -0.75, top = 0.75)
plt.xlim(left = 0, right = 1.3)
plt.legend()
plt.show()
