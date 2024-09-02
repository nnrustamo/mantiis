import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plotQuadTree(quadTreeImage, facecolor=None, visibility='off'):
    """
    Plot quadtree image using rectangles.
    
    Parameters:
    - quadTreeImage: 2D numpy array of integers.
    - facecolor: Optional 3D numpy array of shape (rows, cols, 3) specifying the face color.
    - visibility: String, 'on' or 'off' to control figure visibility.
    
    Returns:
    - fig: Matplotlib figure object.
    """
    if facecolor is None:
        facecolor = np.ones((quadTreeImage.shape[0], quadTreeImage.shape[1], 3))
    
    unique_cells = np.unique(quadTreeImage)
    
    fig, ax = plt.subplots()
    
    ax.set_xlim([0, quadTreeImage.shape[1]])
    ax.set_ylim([0, quadTreeImage.shape[0]])
    ax.invert_yaxis()
    
    for cell_value in unique_cells:
        coords = np.argwhere(quadTreeImage == cell_value)
        if coords.size == 0:
            continue
        
        min_row, min_col = coords.min(axis=0)
        max_row, max_col = coords.max(axis=0)
        
        if cell_value < 0:
            color = [0.5, 0.5, 0.5]
            edgecolor_val = [0, 0, 0]
        else:
            color = facecolor[min_row, min_col]
            edgecolor_val = [0, 0, 0]
        
        rect = patches.Rectangle((min_col, min_row), max_col - min_col + 1, max_row - min_row + 1,
                                 linewidth=1, edgecolor=edgecolor_val, facecolor=color)
        ax.add_patch(rect)
    
    ax.set_aspect('equal')
    if visibility == 'off':
        plt.axis('off')
    
    if visibility == 'on':
        plt.show()
    
    return fig

# Note that this is not recommonded for large images
# Rendering will take too long and not done properly
# For large images, use create_vtk script and visualize it in Paraview
ny = 32
nx = ny
reconstructedImage = np.loadtxt('../reconstructed.dat').reshape((ny, nx))
print(reconstructedImage)
fig = plotQuadTree(reconstructedImage, visibility='on')
fig.savefig('quadtree.svg')
