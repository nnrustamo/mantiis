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

def plot_imgs(array_2d, color_map='turbo', cbar_title = None, name=''):
    """
    Plot an image from a 2D numpy array and add a colorbar with the specified colormap.

    Parameters:
    - array_2d: 2D numpy array representing the image data.
    - color_map: Colormap to use for the colorbar. Default is 'turbo'.
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
    plt.xlim(0, array_2d.shape[1]) 
    tick_positions = np.linspace(0, array_2d.shape[0], 9)
    # plt.xticks(ticks=tick_positions, labels=[f'{int(x)}' for x in tick_positions])
    plt.yticks(ticks=tick_positions, labels=[f'{int(y)}' for y in tick_positions])
    plt.xticks(ticks=tick_positions, labels=[f'{int(y)}' for y in tick_positions])
    plt.savefig(name + ".svg", bbox_inches='tight', pad_inches=0.1)
    plt.close()


if __name__ == "__main__":
    # Ensure enough arguments are provided
    if len(sys.argv) < 4:
        print("Usage: python3 plot_heatmaps.py <file_path> <ny> <nx>")
        print("Optionally: python3 create_vtk.py <file_path> <ny> <nx> <output_directory>")
        sys.exit(1)

    file_path_multi = sys.argv[1]
    file_path_single = sys.argv[2]
    ny = int(sys.argv[3])
    nx = int(sys.argv[4])

    output_path = file_path_single
    if len(sys.argv) == 5:
        output_path = sys.argv[4]

    """
     Each file contains linearized 2D array in column first order
    """
    ux_multi = np.loadtxt(file_path_multi + 'ux.txt').reshape((ny, nx))
    # plot_imgs(ux_multi, "turbo", cbar_title="Streamwise velocity (m/s)", name = output_path + 'ux')

    ux_single = np.loadtxt(file_path_single + 'ux.txt').reshape((ny, nx))
    # plot_imgs(ux_single, "turbo", "Y-dir velocity (m/s)", output_path + 'uy')

    rel_diff = np.abs(ux_multi - ux_single)/ux_single * 100
    # rel_diff = np.where(ux_single == 0, 0, rel_diff)
    plot_imgs(rel_diff, "turbo", cbar_title="Relative error (%)", name = output_path + 'error')
