import numpy as np
import vtk

def create_vtk_unstructured_grid_from_quadtree(quadTreeImage):
    """
    Create a VTK Unstructured Grid from a 2D numpy array representing quadtree data with rectangles.

    Parameters:
    - quadTreeImage: 2D numpy array of integers representing the quadtree.

    Returns:
    - vtkUnstructuredGrid object
    """
    rows, cols = quadTreeImage.shape

    # Create VTK structures
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()

    unique_cells = np.unique(quadTreeImage)

    # Set default colors based on cell values
    vtk_colors = vtk.vtkFloatArray()
    vtk_colors.SetNumberOfComponents(1)
    vtk_colors.SetName("Face Colors")

    for cell_value in unique_cells:
        coords = np.argwhere(quadTreeImage == cell_value)
        if coords.size == 0:
            continue

        min_row, min_col = coords.min(axis=0)
        max_row, max_col = coords.max(axis=0)

        if quadTreeImage[min_row, min_col] < 0:
            # Solid nodes (dark gray)
            color = 0.6
        else:
            # Fluid nodes (white)
            color = 0.3
        
        # Create points for the rectangle
        p0 = [min_col, min_row, 0]
        p1 = [max_col + 1, min_row, 0]
        p2 = [max_col + 1, max_row + 1, 0]
        p3 = [min_col, max_row + 1, 0]

        point_ids = []
        for p in [p0, p1, p2, p3]:
            point_id = points.InsertNextPoint(p)
            point_ids.append(point_id)

        cells.InsertNextCell(4, point_ids)
        vtk_colors.InsertNextValue(color)

    # Create an unstructured grid
    unstructured_grid = vtk.vtkUnstructuredGrid()
    unstructured_grid.SetPoints(points)
    unstructured_grid.SetCells(vtk.VTK_QUAD, cells)
    unstructured_grid.GetCellData().SetScalars(vtk_colors)

    return unstructured_grid

def write_vtk_unstructured_grid(filename, unstructured_grid):
    """
    Write a VTK Unstructured Grid to a file.
    
    Parameters:
    - filename: Output VTK file name.
    - unstructured_grid: vtkUnstructuredGrid object
    """
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(unstructured_grid)
    writer.Write()

# Example usage
ny = 64
nx = ny
reconstructedImage = np.loadtxt('../reconstructed.dat').reshape((ny, nx))
unstructured_grid = create_vtk_unstructured_grid_from_quadtree(reconstructedImage)
write_vtk_unstructured_grid('quadtree.vtu', unstructured_grid)
