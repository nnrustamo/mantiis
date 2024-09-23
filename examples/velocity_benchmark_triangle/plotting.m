clear; close all;
%% Load Data
nx = 512; ny = 512;
file_md = "md/";
% Load MD data
% md_0 = readmatrix(strcat(file_md, "md_0.csv"));
% md_025 = readmatrix(strcat(file_md, "md_025.csv"));
% md_05 = readmatrix(strcat(file_md, "md_05.csv"));
% md_075 = readmatrix(strcat(file_md, "md_075.csv"));
% save('md.mat', 'md_0', 'md_025', 'md_05', 'md_075');

% Load LBM Data
folder_sm = ["single_grid_results/", "multi_grid_results/"];
x_norm = linspace(0, 1, ny - 2)';
for im = 1:length(folder_sm)
    filename = strcat(folder_sm(im) + "ux.txt");
    ux = loadTxtFile(filename);
    ux = reshape(ux, [ny, nx])';
    ux_0 = flipud(ux(2:end-1, 1)); % reverse the order of values
    ux_0_norm = ux_0./mean(ux_0);
end
