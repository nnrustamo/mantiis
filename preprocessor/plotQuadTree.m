function fig = plotQuadTree(quadTreeImage, facecolor, visibility)
%% Plot quadtree image using rectangles
sz = [size(quadTreeImage, 1), size(quadTreeImage, 2), 3];
if nargin<2
    facecolor = ones(sz);
    visibility = 'off';
end

unique_cells = unique(quadTreeImage);

fig = figure('Visible', visibility); hold on;
ylim([0, size(quadTreeImage, 1)]);
xlim([0, size(quadTreeImage, 1)]);
set ( gca, 'ydir', 'reverse' )
for cell_elemenet = 1:length(unique_cells)
    % get cell value
    cell_value = unique_cells(cell_elemenet);
    % get cell coordinates: row = y, col = x
    [row, col] = ind2sub(size(quadTreeImage), find(quadTreeImage == cell_value));
    c1 = facecolor(min(row), min(col), 1);
    c2 = facecolor(min(row), min(col), 2);
    c3 = facecolor(min(row), min(col), 3);
    % plot as rectangles
    if cell_value<0     % solid nodes
        %                       % start pos. (x, y),         width(x-dir),     height (y-dir)
        rectangle("Position",  [min(col), min(row), max(col) - min(col)+1, max(row) - min(row)+1], "FaceColor",[.5 .5 .5]);
    else    % pore nodes
        %                       % start pos. (x, y),         width(x-dir),     height (y-dir)
        rectangle("Position",  [min(col), min(row), max(col) - min(col)+1, max(row) - min(row)+1], 'FaceColor', [c1, c2, c3]);
    end

    % facecolor 
end

end