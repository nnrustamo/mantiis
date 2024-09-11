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

# def create_vtk_structured_grid_from_matrices(matrices, file_name):
#     """
#     Create a VTK Structured Grid from a series of 2D numpy matrices and write it to a VTK file.

#     Parameters:
#     - matrices: List of 2D numpy arrays (all the same shape) representing different fields.
#     - file_name: Output VTK file name.

#     Notes:
#     - The matrices should be in the order [rho, ux, uy, uz, Knudsen number, Localporesize].
#     """
#     if len(matrices) != 6:
#         raise ValueError("There should be exactly 6 matrices (rho, ux, uy, uz, Knudsen number, Localporesize).")

#     num_slices = len(matrices)
#     rows, cols = matrices[0].shape

#     # Create VTK structures
#     structured_grid = vtk.vtkStructuredGrid()

#     # Define the dimensions of the grid
#     structured_grid.SetDimensions(cols, rows, num_slices)

#     # Create points for the structured grid
#     points = vtk.vtkPoints()
#     for z in range(num_slices):
#         for y in range(rows):
#             for x in range(cols):
#                 points.InsertNextPoint(x, y, z)

#     structured_grid.SetPoints(points)

#     # Add scalars for each matrix
#     for i, matrix in enumerate(matrices):
#         scalars = vtk.vtkFloatArray()
#         scalars.SetName(f"Field_{i}")
#         scalars.SetNumberOfComponents(1)
        
#         for z in range(num_slices):
#             for y in range(rows):
#                 for x in range(cols):
#                     scalars.InsertNextValue(matrix[y, x])
        
#         structured_grid.GetCellData().AddArray(scalars)

#     # Write the structured grid to a VTK file
#     writer = vtk.vtkStructuredGridWriter()
#     writer.SetFileName(file_name)
#     writer.SetInputData(structured_grid)
#     writer.Write()

# if __name__ == "__main__":
#     ny = 1024
#     nx = 1024
#     num_slices = 5
#     rows, cols = 10, 10
#     matrices = []

#     reconstructedImage = np.loadtxt('../reconstructed.dat').reshape((ny, nx))
#     matrices.append(reconstructedImage)
#     rho = np.loadtxt('../input_output/rho.txt').reshape((ny, nx))
#     matrices.append(rho)
#     ux = np.loadtxt('../input_output/ux.txt').reshape((ny, nx))
#     matrices.append(ux)
#     uy = np.loadtxt('../input_output/ux.txt').reshape((ny, nx))
#     matrices.append(uy)
#     Kn = np.loadtxt('../input_output/Kn.txt').reshape((ny, nx))
#     matrices.append(uy)
    
#     create_vtk_structured_grid_from_matrices(matrices, "structured_grid.vtk")

#     diff_over_time = np.loadtxt('../input_output/convergence.txt').reshape((ny, nx))

ny = 512
nx = 512
# create_vtk_structured_grid_from_matrices()
reconstructedImage = np.loadtxt('../reconstructed.dat').reshape((ny, nx))
unstructured_grid = create_vtk_unstructured_grid_from_quadtree(reconstructedImage)
write_vtk_unstructured_grid('quadtree.vtu', unstructured_grid)
