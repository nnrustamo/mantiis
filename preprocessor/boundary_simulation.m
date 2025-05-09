close all; clear; clc;
%%
% Initialize
N = 256;
ifStream = zeros(N, 9);
xy_norm = zeros(N, 2);
boundary_sim = zeros(3, 3, N);
ifStream(: ,1) = 1;

count  = 0;
count_id = [];
% Center coordinates
i_center = 2; j_center = 2;
% D2Q9 config
Cx = [0 1 0 -1 0 1 -1 -1 1];
Cy = [0 0 1 0 -1 1 1 -1 -1];
remove_based_on_direction_count = [];
remove_based_on_shape = [];
% Iterate through all possible binary numbers for the remaining 8 cells
for i = 0:2^9-1
    % Convert the decimal number to binary and convert it into a 3x3 matrix
    binaryStr = dec2bin(i, 9);
    binaryMatrix = reshape(binaryStr, [3, 3]) - '0';
    % Check if the center cell contains a 1
    if binaryMatrix(2, 2) == 1 
        count = count +1;
        count_id = [count_id, i];
        boundary_sim(:, :, count) = binaryMatrix;
        % categorize ifStream
        if binaryMatrix(8) == 1
            ifStream(count, 2) = 1;
        end
        if binaryMatrix(6) == 1
            ifStream(count, 3) = 1;
        end
        if binaryMatrix(2) == 1
            ifStream(count, 4) = 1;
        end
        if binaryMatrix(4) == 1
            ifStream(count, 5) = 1;
        end
        if binaryMatrix(6) == 1 && binaryMatrix(8) == 1 && binaryMatrix(9) == 1
            ifStream(count, 6) = 1;
        end
        if binaryMatrix(2) == 1 && binaryMatrix(6) == 1 && binaryMatrix(3) == 1
            ifStream(count, 7) = 1;
        end
        if binaryMatrix(2) == 1 && binaryMatrix(4) == 1 && binaryMatrix(1) == 1
            ifStream(count, 8) = 1;
        end
        if binaryMatrix(4) == 1 && binaryMatrix(8) == 1 && binaryMatrix(7) == 1
            ifStream(count, 9) = 1;
        end
        % normal vectors
        direction_counter = 0;
        for ic = 2:9
            if binaryMatrix(i_center + Cy(ic), j_center + Cx(ic)) && ...
                    ifStream(count, ic) == 1
                xy_norm(count, 1) = xy_norm(count, 1) + Cx(ic);
                xy_norm(count, 2) = xy_norm(count, 2) + Cy(ic);
                direction_counter = direction_counter + 1;
            end
        end
        % cleaning 1
        if direction_counter < 3
            remove_based_on_direction_count = [remove_based_on_direction_count, count];
        end
        % cleaning 2
        for row = 1:3   % traverse each row
            rowVec = binaryMatrix(row, :);
            if rowVec(1) == 0 && rowVec(2) == 1 && rowVec(3) == 0
                remove_based_on_shape = [remove_based_on_shape , count];
                break;
            end
        end
        for col = 1:3   % traverse each column
            colVec = binaryMatrix(:, col);
            if colVec(1) == 0 && colVec(2) == 1 && colVec(3) == 0
                remove_based_on_shape = [remove_based_on_shape , count];
                break;
            end
        end
        diags = [diag(binaryMatrix), diag(flipud(binaryMatrix))];
        for d = 1:size(diags, 2)  % traverse each diagonal
            diagVec = diags(:, d);
            if diagVec(1) == 0 && diagVec(2) == 1 && diagVec(3) == 0
                remove_based_on_shape = [remove_based_on_shape , count];
                break;
            end
        end
        
    end
end

% normalize
xy_norm_magnitude = sqrt(xy_norm(: ,1).^2 + xy_norm(: ,2).^2);
xy_norm(:, 1) = xy_norm(:, 1) ./ xy_norm_magnitude;
xy_norm(:, 2) = xy_norm(:, 2) ./ xy_norm_magnitude;
remaining_boundary_ids = [];
nonRemovedTypes = [];
for i = 1:N
    if  ~isnan(xy_norm(i, 1)) && ~isnan(xy_norm(i, 2)) && isempty(find(remove_based_on_direction_count == i)) ...
            && isempty(find(remove_based_on_shape == i))
        nonRemovedTypes = [nonRemovedTypes, i];
        remaining_boundary_ids = [remaining_boundary_ids, count_id(i)];
    end
end

% filter garbage data
xy_norm_new = zeros(length(nonRemovedTypes), 2);
ifStream_new = zeros(length(nonRemovedTypes), 9);
icsr = zeros(length(nonRemovedTypes), 9);

% sample types for plotting
bnd_to_be_plotted = [127, 445, 473, 383];
figure;
subplot_counter = 1;
for i = 1:length(nonRemovedTypes)
    binaryMatrix = boundary_sim(:, :, nonRemovedTypes(i));
    bndmat = reshape(binaryMatrix, [1, 9]);
    decimals = bin2dec(num2str(bndmat));
    if ~isempty(find(bnd_to_be_plotted == decimals))
        subplot(2, 2, subplot_counter);
        imagesc(binaryMatrix);
        title(['Boundary Number: ', num2str(decimals)]);
        axis off;
        colormap(gray);
        subplot_counter = subplot_counter + 1;
    end
end


% filter
vid = figure;
writerObj = VideoWriter('boundary_types.avi');
open(writerObj);
for i = 1:length(nonRemovedTypes)
    ifStream_new(i, :) = ifStream(nonRemovedTypes(i), :);
    xy_norm_new(i, :) = xy_norm(nonRemovedTypes(i), :);
    % calcualte specular reflection direction
    icsr(i, 1) = 1;
    for ic = 2:9
        if ~ifStream_new(i, ic)
            ei_x = round(Cx(ic) - 2 * (Cx(ic) * xy_norm_new(i, 1) + Cy(ic) * xy_norm_new(i, 2)) * xy_norm_new(i, 1));
            ei_y = round(Cy(ic) - 2 * (Cx(ic) * xy_norm_new(i, 1) + Cy(ic) * xy_norm_new(i, 2)) * xy_norm_new(i, 2));
            ei_SR = find(Cx == ei_x & Cy == ei_y);
            if ~isempty(ei_SR)
                icsr(i, ic) = ei_SR;
            else
                icsr(i, ic) = ic;
            end
        else
            icsr(i, ic) = ic;
        end 
    end

    % Write the current frame to the video
    imagesc(binaryMatrix);
    colormap(gray);
    title(['Boundary Number: ', num2str(remaining_boundary_ids(i))]);
    axis off;
    % drawnow;
    pause(0.3);
    writeVideo(writerObj, getframe(vid));

    binaryMatrix = boundary_sim(:, :, nonRemovedTypes(i));
    bndmat = reshape(binaryMatrix, [1, 9]);
    decimals = bin2dec(num2str(bndmat));
end
remaining_boundary_ids = remaining_boundary_ids';
% %% Save
% % Note these files are the same as in ../src/system_files
% close(writerObj);
% close all;
% figure;
% imagesc(ifStream);
% colormap gray;
% colorbar;
% savetoFile(ifStream_new, "ifstream.txt");
% savetoFile(icsr, "icsr.txt");
% savetoFile(xy_norm_new, "xy_norm.txt");
% savetoFile(remaining_boundary_ids, "bnd_types.txt");
