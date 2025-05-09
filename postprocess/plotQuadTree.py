import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import math

def plotQuadTree(quadTreeImage, maxLevel, visibility='off'):
    """
    Plot quadtree image using rectangles
    Parameters:
    - quadTreeImage: 2D numpy array of integers.
    - maxLevel: Integer, the maximum level to display (0, 1, ..., maxLevel).
    - visibility: String, 'on' or 'off' to control figure visibility.

    Returns:
    - fig: Matplotlib figure object.
    """
    
    # Define the facecolors for each level, we only use the levels up to maxLevel
    facecolors = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0.647, 0]]

    facecolors = facecolors[:maxLevel + 1]

    unique_cells = np.unique(quadTreeImage)
    
    fig, ax = plt.subplots()
    
    ax.set_xlim([0, quadTreeImage.shape[1]])
    ax.set_ylim([0, quadTreeImage.shape[0]])
    ax.invert_yaxis()
    
    # Plot each cell as a rectangle
    for cell_value in unique_cells:
        coords = np.argwhere(quadTreeImage == cell_value)
        if coords.size == 0:
            continue
        
        level = int(math.log2(math.sqrt(coords.size / 2)))

        min_row, min_col = coords.min(axis=0)
        max_row, max_col = coords.max(axis=0)
        
        if cell_value < 0:
            color = [0, 0.1, 0.1]
            edgecolor_val = [0, 0, 0]
        else:
            color = facecolors[level]
            edgecolor_val = facecolors[level]
        
        rect = patches.Rectangle((min_col, min_row), max_col - min_col + 1, max_row - min_row + 1,
                                 linewidth=1, edgecolor=edgecolor_val, facecolor=color)
        ax.add_patch(rect)
    
    ax.set_aspect('equal')

    if visibility == 'off':
        plt.axis('off')
    
    cmap = mcolors.ListedColormap(facecolors) 
    bounds = np.arange(maxLevel + 2)  
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    dynamic_ticks = np.arange(maxLevel+1) + 0.5  
    cbar = fig.colorbar(plt.imshow(np.zeros_like(quadTreeImage), cmap=cmap, norm=norm),
                        ticks=dynamic_ticks, ax=ax)
    cbar.ax.set_yticklabels([f'level {i}' for i in range(maxLevel+1)])

    if visibility == 'on':
        plt.show()
    
    return fig


if __name__ == "__main__":
    # Ensure enough arguments are provided
    if len(sys.argv) < 2:
        print("Usage: python3 plotQuadTree.py <ny>")
        sys.exit(1)

    ny = int(sys.argv[1])
    nx = ny 
    maxLevel = int(sys.argv[2])
    file_name = sys.argv[3]
    """
     Each file contains linearized 2D array in column first order
    """
    # Note that this is not recommonded for large images
    # Rendering will take too long and not done properly
    # For large images, use create_vtk script and visualize it in Paraview

    reconstructedImage = np.loadtxt('../reconstructed.dat').reshape((ny, nx))
    print(reconstructedImage)
    fig = plotQuadTree(reconstructedImage, maxLevel, visibility='off')
    fig.savefig(file_name)
