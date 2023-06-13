import matplotlib.pyplot as plt
import numpy as np

def generate_heatmap(E, sigma, nbEntangled):

    E_grid, sigma_grid = np.meshgrid(E, sigma)

    fig,ax = plt.subplots()

    # Plot the heatmap
    heatmap = ax.pcolormesh(E_grid, sigma_grid, nbEntangled, cmap='viridis')

    # Set the colorbar
    cbar = plt.colorbar(heatmap)

    # Set the labels and title
    ax.set_xlabel(r'$E_{0} [eV]$')
    ax.set_ylabel(r'$ \sigma$ [nm]')

    # Show the plot
    plt.show()

E = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
sigma = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
nbEntangled = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
               [10, 9, 8, 7, 6, 5, 4, 3, 2, 1],
               [5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
               [3, 6, 9, 2, 5, 8, 1, 4, 7, 10],
               [7, 4, 1, 8, 5, 2, 9, 6, 3, 10],
               [2, 4, 6, 8, 10, 1, 3, 5, 7, 9],
               [9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
               [10, 5, 2, 8, 1, 3, 9, 7, 6, 4],
               [6, 3, 7, 10, 2, 9, 1, 4, 8, 5]]

generate_heatmap(E, sigma, nbEntangled)
