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

