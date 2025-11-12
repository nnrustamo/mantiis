# /*
#  * Copyright (c) January 2024
#  *
#  * Author: Nijat Rustamov
#  * Organization: University of Wyoming
#  * Email: nrustamo@uwyo.edu
#  *
#  * Academic Supervisor: Saman Aryana
#  * Email: saryana@uwyo.edu
#  * Organization: University of Wyoming
#  *
#  * This file is a part of Lattice Boltzmann Simulation Software
#  * Proprietary Software - All Rights Reserved
#  *
#  * Unauthorized copying, modification, or distribution of this software,
#  * or any portion of it is prohibited
#  */


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys

def plot_imgs(array_2d, color_map='seismic', cbar_title = None, name=''):
    """
    Plot an image from a 2D numpy array and add a colorbar with the specified colormap.

    Parameters:
    - array_2d: 2D numpy array representing the image data.
    - color_map: Colormap to use for the colorbar. Default is 'seismic'.
    - name: Name of the png file to save to

    Returns:
    - None
    """
    if array_2d.ndim != 2:
        raise ValueError("Input array must be 2D.")

    plt.figure(figsize=(8, 6))
    img = plt.imshow(array_2d, cmap=color_map, origin='lower')
    cbar = plt.colorbar(img, orientation='vertical')
    cbar.set_label(cbar_title)

    # cbar_ticks = np.linspace(array_2d.min(), array_2d.max(), 10)
    # rounded_ticks = np.round(cbar_ticks, 2)   
    # cbar.set_ticks(rounded_ticks)

    plt.ylim(0, array_2d.shape[0]) 
    plt.xlim(0, array_2d.shape[0]) 
    tick_positions = np.linspace(0, array_2d.shape[0], 9)
    # plt.xticks(ticks=tick_positions, labels=[f'{int(x)}' for x in tick_positions])
    plt.yticks(ticks=tick_positions, labels=[f'{int(y)}' for y in tick_positions])
    plt.xticks(ticks=tick_positions, labels=[f'{int(y)}' for y in tick_positions])
    plt.savefig(name + ".svg", bbox_inches='tight', pad_inches=0.1)
    plt.close()

def plot_contour(array_2d, color_map='seismic', cbar_title=None, name=''):
    """
    Plot a contour plot from a 2D numpy array and add a colorbar with the specified colormap.

    Parameters:
    - array_2d: 2D numpy array representing the image data.
    - color_map: Colormap to use for the contour plot. Default is 'seismic'.
    - cbar_title: Title for the colorbar. Default is None.
    - name: Name of the PNG file to save to.

    Returns:
    - None
    """
    if array_2d.ndim != 2:
        raise ValueError("Input array must be 2D.")

    plt.figure(figsize=(8, 6))

    # Create the contour plot
    contour = plt.contourf(array_2d, cmap=color_map, levels=20)

    # Add colorbar
    cbar = plt.colorbar(contour, orientation='vertical')
    cbar.set_label(cbar_title)

    # Set axis limits and labels
    plt.xlim(0, array_2d.shape[1])
    plt.ylim(0, array_2d.shape[0])
    
    tick_positions = np.linspace(0, array_2d.shape[0], 9)
    plt.yticks(ticks=tick_positions, labels=[f'{int(y)}' for y in tick_positions])
    
    plt.savefig(name + ".svg", bbox_inches='tight', pad_inches=0.1)
    plt.close()

if __name__ == "__main__":
    # Ensure enough arguments are provided
    if len(sys.argv) < 4:
        print("Usage: python3 plot_heatmaps.py <file_path> <ny> <nx>")
        print("Optionally: python3 create_vtk.py <file_path> <ny> <nx> <output_directory>")
        sys.exit(1)

    file_path = sys.argv[1]
    ny = int(sys.argv[2])
    nx = int(sys.argv[3])

    output_path = file_path
    if len(sys.argv) == 5:
        output_path = sys.argv[4]

    """
     Each file contains linearized 2D array in column first order
    """
    ux = np.loadtxt(file_path + 'ux.txt').reshape((ny, nx))
    print(f"Maximum of ux = {np.max(ux, axis=None)}, minimum of ux = {np.min(ux, axis=None)}")
    plot_imgs(ux, "seismic", cbar_title="Streamwise velocity (m/s)", name = output_path + 'ux')

    uy = np.loadtxt(file_path + 'uy.txt').reshape((ny, nx))
    print(f"Maximum of uy = {np.max(uy, axis=None)}, minimum of uy = {np.min(uy, axis=None)}")
    plot_imgs(uy, "seismic", "Y-dir velocity (m/s)", output_path + 'uy')

    speed = np.sqrt(ux**2 + uy**2)
    plot_imgs(speed, "seismic", "Velocity magnitude (m/s)", output_path + 'speed')

    rho = np.loadtxt(file_path + 'rho.txt').reshape((ny, nx))
    plot_imgs(rho, "seismic", "Density (km/m3)", output_path + 'rho')

    kn = np.loadtxt(file_path + 'Kn.dat').reshape((ny, nx)).transpose()
    plot_imgs(kn, "seismic", "Knudsen number", output_path + 'kn')

    poresize = np.loadtxt(file_path + 'localporesize.dat').reshape((ny, nx)).transpose()
    plot_imgs(poresize, "seismic", "Local characteristic length", output_path + 'localpore')
