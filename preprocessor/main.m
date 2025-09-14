close all; clear; clc
%% Define domain parameters
% dimensions
ny = 10;
nx = 10;
% write to this folder
% output_folder = '/home/nrustamo/projects/HNO/image_generation/images/nx1024_ny1024_scaling0.8_sigma50_octave40/';
output_folder = '../large_domain/'
% physical conditions
Pressure = 2.0e6; % Pa
Temperature = 300; % K
Resolution = 1.0e-9; % m
[~, ~, mean_free_path] = thermodynamicProperties(Pressure, Temperature);

% create a domain or load existing one
% pore = readBinaryMatrix(output_folder + "binary_image.txt");
% disp(pore)
% pore = reshape(pore, [ny, nx]);
% pore(:, 1:50) = 1;
% pore(:, end-51:end) = 1;
% pore(1, :) = 0;
% pore(end, :) = 0;
% disp(pore)

% pore(:, 1:25) = 1;
% pore(:, end-24:end) = 1;


% pore = randomGaussianDomain(ny, nx);

% pore = randomDomain(ny, nx);
pore = simplePore(ny, nx);
% pore = narrowing_tube(nx, 500, 10);
% pore = createRectangleObstacleFlowGeom(nx);
% pore = createTriangularFlowGeom(nx);
% pore = porous_media_circles(nx, 10, [10, 20]);

% figure(); imagesc(pore); colormap gray;
% saveas(gcf, strcat(output_folder, 'domain.png'))

% clean boundaries
pore = clean_boundaries(pore);
% figure; imagesc(pore); colormap gray;

% local pore size
[kn, localpore, pore, medial_axis] = poreProperties(pore, mean_free_path, Resolution);

% manual treatment for some validation cases only.
% kn_man = 0.15;
% locpore_man = mean_free_path/kn_man/Resolution;
% kn(kn~=0) = kn_man;
% localpore(localpore~=0) = locpore_man;

%% Save
figure(); imagesc(pore); colormap gray;
saveas(gcf, strcat(output_folder, 'domain.png'))
figure(); imagesc(kn); colormap jet; colorbar;
saveas(gcf, strcat(output_folder, 'kn.png'))
figure(); imagesc(localpore); colormap jet; colorbar;
saveas(gcf, strcat(output_folder, 'localpore.png'))

savetoFile(kn, strcat(output_folder ,'Kn.dat'))
savetoFile(localpore, strcat(output_folder ,'localporesize.dat'))
savetoFile(pore, strcat(output_folder ,'pore.dat'))
