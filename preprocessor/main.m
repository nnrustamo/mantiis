close all; clear;
%% Define domain parameters
% dimensions
ny  = 512;
nx = 512;
% write to this folder
output_folder = '../input_output/';

% physical conditions
Pressure = 2.0e6; % Pa
Temperature = 300; % K
Resolution = 1e-9; % m
[~, ~, mean_free_path] = thermodynamicProperties(Pressure, Temperature);

% create a domain or load existing onee 
% pore = randomGaussianDomain(ny, nx);
% pore = simplePore(ny, nx);
pore = narrowing_tube(nx, 500, 10);

% clean boundaries
% pore = clean_boundaries(pore);
% figure; imagesc(pore); colormap gray;

% local pore size
[kn, localpore, pore] = poreProperties(pore, mean_free_path, Resolution);

%% Save
figure('visible','off'); imagesc(pore); colormap gray;
saveas(gcf, strcat(output_folder, 'domain.svg'))
figure('visible','off'); imagesc(kn); colormap jet; colorbar;
saveas(gcf, strcat(output_folder, 'kn.svg'))
figure('visible','off'); imagesc(localpore); colormap jet; colorbar;
saveas(gcf, strcat(output_folder, 'localpore.svg'))

savetoFile(kn, strcat(output_folder ,'Kn.txt'))
savetoFile(localpore, strcat(output_folder ,'localporesize.txt'))
savetoFile(pore, strcat(output_folder ,'pore.txt'))
