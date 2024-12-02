close all; clear;
%% Define domain parameters
% dimensions
ny  = 500;
nx = 500;
% write to this folder
output_folder = '../slit_pore/';

% physical conditions
Pressure = 2.0e6; % Pa
Temperature = 300; % K
Resolution = 1.0e-9; % m
[~, ~, mean_free_path] = thermodynamicProperties(Pressure, Temperature);

% create a domain or load existing one
% pore = loadTxtFile("../pore_raw.txt");
% pore = reshape(pore, [ny, nx]);
pore = randomGaussianDomain(ny, nx);
% pore = randomDomain(ny, nx);
% pore = simplePore(ny, nx);
% pore = narrowing_tube(nx, 500, 10);
% pore = createRectangleObstacleFlowGeom(nx);
% pore = createTriangularFlowGeom(nx);
% pore = porous_media_circles(nx, 50, [10, 100]);

% clean boundaries
pore = clean_boundaries(pore);
% figure; imagesc(pore); colormap gray;

% local pore size
[kn, localpore, pore] = poreProperties(pore, mean_free_path, Resolution);

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

savetoFile(kn, strcat(output_folder ,'Kn.txt'))
savetoFile(localpore, strcat(output_folder ,'localporesize.txt'))
savetoFile(pore, strcat(output_folder ,'pore.txt'))
