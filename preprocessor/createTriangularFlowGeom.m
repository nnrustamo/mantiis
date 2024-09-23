function binaryMatrix = createTriangularFlowGeom(H)
    % Initialize an H x H matrix of ones (reversed color scheme)
    binaryMatrix = ones(H, H);  
    
    % Define the coordinates of the triangle's vertices
    triangleHeight = 0.28 * H;   % Height of the triangle
    triangleWidth = H / 3;       % Width of the base
    
    xCenter = H / 2;             % Horizontal center
    yBase = 2 * H / 3;           % Y-coordinate of the triangle's base
    
    % Define the three vertices of the triangle
    xVertices = [xCenter - triangleWidth / 2, xCenter + triangleWidth / 2, xCenter];
    yVertices = [yBase, yBase, yBase - triangleHeight];
    
    % Create the triangle using poly2mask
    mask = poly2mask(xVertices, yVertices, H, H);
    
    % Set the triangle area to 0 and the rest to 1 (inverted)
    binaryMatrix(mask) = 0;
    
    % Set the top and bottom walls
    binaryMatrix(1, :) = 0;   % Top wall
    binaryMatrix(end, :) = 0; % Bottom wall
end
