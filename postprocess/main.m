close all; clear;
%% Define domain parameters
% dimensions
ny  = 512;
nx = 512;
% write to this folder
output_folder = '../in_512/single/';

% load files
ux = loadTxtFile(strcat(output_folder, 'ux.txt'));
uy = loadTxtFile(strcat(output_folder, 'uy.txt'));
rho = loadTxtFile(strcat(output_folder, 'rho.txt'));
convergence = loadTxtFile(strcat(output_folder, 'convergence.txt'));

% reshape
ux = reshape(ux, [ny, nx])';
uy = reshape(uy, [ny, nx])';
rho = reshape(rho, [ny, nx])';

% visualize
u = sqrt(ux.^2 + uy.^2);
figure; imagesc(u); colormap turbo; colorbar;
figure; plot(convergence);
set(gca, 'Yscale', 'log');

% if multi-block
% reconstructedImage = loadTxtFile(strcat(output_folder, 'recon.dat'));
% reconstructedImage = reshape(reconstructedImage, [ny, nx])';
% fig = plotQuadTree(reconstructedImage);
% shg
% saveas(fig, 'quadtree.svg'); % Save the figure as PNG

