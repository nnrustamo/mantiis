% load('data.mat')
%% Load Data
ny = 64; nx = 64;
folder_sm = ["single_grid_1024", "multi_grid"];
folder_kn = ["kn = 0.01", "kn = 0.1", "kn = 1", "kn = 10"];

lbm = zeros(ny -2, length(folder_kn) * length(folder_sm) + 1);
lbm(:, 1) = linspace(0, 1, ny-2);
count = 2;
for im = 1:length(folder_sm)
    for ik = 1:length(folder_kn)
        filename = strcat(folder_sm(im), "/", folder_kn(ik), "/ux.txt");
        ux = loadTxtFile(filename);
        ux = reshape(ux, [ny, nx])';
        ux_line = ux(2:end-1, 1);
        ux_norm = ux_line./mean(ux_line);
        lbm(:, count) = ux_norm;
        count = count + 1;
    end
end

% experimental data
load('dsmc.mat');

%% 
% Interpolate LB data to match experimental x-values
x_exp = dsmc(:, 1);
y_exp = dsmc(:, 2);
x_lbm = lbm(:, 1);
y_lbm = lbm(:, 2);
y_exp_interpolated = interp1(x_exp, y_exp, x_lbm, 'linear');
percentage_error = abs(y_lbm - y_exp_interpolated) ./ abs(y_exp_interpolated) * 100;


fig = figure;
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 12, 8]);
% fig.Visible = 'off';
subplot(2, 2, 1)
yyaxis left
scatter(lbm(:, 1), lbm(:, 6), 100, 'Marker', 'square', 'MarkerEdgeColor', 'b', ...
    'MarkerFaceColor', 'b',  'DisplayName', 'Multi-block LB');
hold on;
scatter(lbm(:, 1), lbm(:, 2), 20, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 'Single-block LB');
plot(dsmc(:, 1), dsmc(:, 2), 'Color', 'k', 'DisplayName', 'DSMC');

yyaxis right
plot(x_lbm(3:end-3), percentage_error(3:end-3), 'ok--', 'DisplayName', 'Percentage Error');
ylabel('Percentage Error (%)');
ylim([0.9*min(percentage_error), 1.1*max(percentage_error)]);


yyaxis left
set(gca,'FontSize',14);
lgd = legend('Location','south', 'Box', 'off', 'FontSize', 12);
legend 'FontName' 'Times New Roman'
pos = lgd.Position;
pos(2) = pos(2) + 0.1;
lgd.Position = pos;
xlabel("Normalized Position", 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel("Normalized Velocity", 'FontSize', 14, 'FontName', 'Times New Roman');

% set(gca,'YAxisLocation', 'right')
% title('P = 0.5 MPa', 'FontSize', 14, 'FontName', 'Times New Roman')
ax=gca;
ax.FontSize = 14;
set(gca,'FontName','Times New Roman')
set(groot,'defaultAxesTickLabelInterpreter','latex');
axis tight
ylim([0.9*min(min(dsmc(:, [2, 4, 6, 8]))), 1.1*max(max(dsmc(:, [2, 4, 6, 8])))]);
xlim([0, 1.0]);
box on;
text(0.73, 1.5 ,'Kn = 0.01', 'FontSize', 12)

%%
% Interpolate LB data to match experimental x-values
x_exp = dsmc(:, 3);
y_exp = dsmc(:, 4);
x_lbm = lbm(:, 1);
y_lbm = lbm(:, 3);
y_exp_interpolated = interp1(x_exp, y_exp, x_lbm, 'linear');
percentage_error = abs(y_lbm - y_exp_interpolated) ./ abs(y_exp_interpolated) * 100;

subplot(2, 2, 2)
yyaxis left
scatter(lbm(:, 1), lbm(:, 7), 100, 'Marker', 'square', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'DisplayName', 'Multi-block LB');
hold on
scatter(lbm(:, 1), lbm(:, 3), 20, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 'Single-block LB');
plot(dsmc(:, 3), dsmc(:, 4), 'Color', 'k', 'DisplayName', 'DSMC');

yyaxis right
plot(x_lbm, percentage_error, 'ok--', 'DisplayName', 'Percentage Error');
ylabel('Percentage Error (%)');
ylim([0.9*min(percentage_error), 1.1*max(percentage_error)]);


yyaxis left
set(gca,'FontSize',14);
legend('Location','south', 'Box', 'off', 'FontSize', 12);
legend 'FontName' 'Times New Roman'
xlabel("Normalized Position", 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel("Normalized Velocity", 'FontSize', 14, 'FontName', 'Times New Roman');
% set(gca,'YAxisLocation', 'right')
% title('P = 2.0 MPa', 'FontSize', 14, 'FontName', 'Times New Roman')
ax=gca;
ax.FontSize = 14;
set(gca,'FontName','Times New Roman')
set(groot,'defaultAxesTickLabelInterpreter','latex');
axis tight
ylim([0.9*min(min(dsmc(:, [2, 4, 6, 8]))), 1.1*max(max(dsmc(:, [2, 4, 6, 8])))]);
xlim([0, 1.0]);
box on;
text(0.73, 1.5 ,'Kn = 0.1', 'FontSize', 12)


%%
% Interpolate LB data to match experimental x-values
x_exp = dsmc(:, 5);
y_exp = dsmc(:, 6);
x_lbm = lbm(:, 1);
y_lbm = lbm(:, 4);
y_exp_interpolated = interp1(x_exp, y_exp, x_lbm, 'linear');
percentage_error = abs(y_lbm - y_exp_interpolated) ./ abs(y_exp_interpolated) * 100;


subplot(2, 2, 3)
yyaxis left
scatter(lbm(:, 1), lbm(:, 8), 100, 'Marker', 'square', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'DisplayName', 'Multi-block LB');
hold on
scatter(lbm(:, 1), lbm(:, 4), 20, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 'Single-block LB');
plot(dsmc(:, 5), dsmc(:, 6), 'Color', 'k', 'DisplayName', 'DSMC');

yyaxis right
plot(x_lbm, percentage_error, 'ok--', 'DisplayName', 'Percentage Error');
ylabel('Percentage Error (%)');
ylim([0.9*min(percentage_error), 1.1*max(percentage_error)]);


yyaxis left
set(gca,'FontSize',14);
legend('Location','south', 'Box', 'off', 'FontSize', 12);
legend 'FontName' 'Times New Roman'
xlabel("Normalized Position", 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel("Normalized Velocity", 'FontSize', 14, 'FontName', 'Times New Roman');
% set(gca,'YAxisLocation', 'right')
% title('P = 5.0 MPa', 'FontSize', 14, 'FontName', 'Times New Roman')
ax=gca;
ax.FontSize = 14;
ax.FontName = 'Times New Roman';
set(gca,'FontName','Times New Roman')
set(groot,'defaultAxesTickLabelInterpreter','latex');
axis tight
ylim([0.9*min(min(dsmc(:, [2, 4, 6, 8]))), 1.1*max(max(dsmc(:, [2, 4, 6, 8])))]);
xlim([0, 1.0]);
box on;
text(0.73, 1.5 ,'Kn = 1.0', 'FontSize', 12)


%%
% Interpolate LB data to match experimental x-values
x_exp = dsmc(:, 7);
y_exp = dsmc(:, 8);
x_lbm = lbm(:, 1);
y_lbm = lbm(:, 5);
y_exp_interpolated = interp1(x_exp, y_exp, x_lbm, 'linear');
percentage_error = abs(y_lbm - y_exp_interpolated) ./ abs(y_exp_interpolated) * 100;

subplot(2, 2, 4)
yyaxis left
scatter(lbm(:, 1), lbm(:, 9), 100, 'Marker', 'square', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'DisplayName', 'Multi-block LB');
hold on
scatter(lbm(:, 1), lbm(:, 5), 20, 'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 'Single-block LB');
plot(dsmc(:, 7), dsmc(:, 8),  'Color', 'k', 'DisplayName', 'DSMC');

yyaxis right
plot(x_lbm, percentage_error, 'ok--', 'DisplayName', 'Percentage Error');
ylabel('Percentage Error (%)');
ylim([0.9*min(percentage_error), 1.1*max(percentage_error)]);

yyaxis left
set(gca,'FontSize',14);
lgd = legend('Location','southeast', 'Box', 'off', 'FontSize', 12);
legend 'FontName' 'Times New Roman'
pos = lgd.Position;
pos(1) = pos(1) - 0.08;
lgd.Position = pos;
xlabel("Normalized Position", 'FontSize', 14, 'FontName', 'Times New Roman')
ylabel("Normalized Velocity", 'FontSize', 14, 'FontName', 'Times New Roman');
% set(gca,'YAxisLocation', 'right')
% title('P = 10.0 MPa', 'FontSize', 14,'FontName', 'Times New Roman')
ax=gca;
ax.FontSize = 14;
set(gca,'FontName','Times New Roman')
set(groot,'defaultAxesTickLabelInterpreter','latex');
axis tight
ylim([0.9*min(min(dsmc(:, [2, 4, 6, 8]))), 1.1*max(max(dsmc(:, [2, 4, 6, 8])))]);
xlim([0, 1.0]);
box on;
text(0.75, 1.5 ,'Kn = 10.0', 'FontSize', 12)
%% 

saveas(gcf, "velocity_results.png")
