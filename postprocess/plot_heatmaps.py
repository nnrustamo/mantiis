import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plot_imgs(array_2d, color_map='turbo', name='turbo'):
    """
    Plot an image from a 2D numpy array and add a colorbar with the specified colormap.

    Parameters:
    - array_2d: 2D numpy array representing the image data.
    - color_map: Colormap to use for the colorbar. Default is 'turbo'.

    Returns:
    - None
    """
    if array_2d.ndim != 2:
        raise ValueError("Input array must be 2D.")

    plt.figure(figsize=(8, 6))
    img = plt.imshow(array_2d, cmap=color_map, origin='lower')
    cbar = plt.colorbar(img, orientation='vertical', shrink=0.8, aspect=10)
    
    cbar_ticks = np.linspace(array_2d.min(), array_2d.max(), 10)
    rounded_ticks = np.round(cbar_ticks, 2)   
    cbar.set_ticks(rounded_ticks)

    plt.xlim(0, array_2d.shape[0]) 
    plt.ylim(0, array_2d.shape[1]) 
    tick_positions = np.linspace(0, array_2d.shape[0], 9)
    plt.xticks(ticks=tick_positions, labels=[f'{int(x)}' for x in tick_positions])
    plt.yticks(ticks=tick_positions, labels=[f'{int(y)}' for y in tick_positions])
    plt.savefig(name + ".png", bbox_inches='tight', pad_inches=0.1)
    plt.close()


if __name__ == "__main__":
    ny = 512
    nx = ny
    pore = np.loadtxt('../input_output/pore.txt').reshape((ny, nx))
    plot_imgs(pore, "gray", "pore")
    kn = np.loadtxt('../input_output/Kn.txt').reshape((ny, nx))
    plot_imgs(kn, "turbo", "kn")
    localporesize = np.loadtxt('../input_output/localporesize.txt').reshape((ny, nx))
    plot_imgs(localporesize, "turbo", "localporesize")
    
    ux = np.loadtxt('../input_output/ux.txt').reshape((ny, nx))
    plot_imgs(ux, "turbo", "ux")
    uy = np.loadtxt('../input_output/uy.txt').reshape((ny, nx))
    plot_imgs(uy, "turbo", "uy")
    rho = np.loadtxt('../input_output/rho.txt').reshape((ny, nx))
    plot_imgs(rho, "turbo", "rho")
