import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams["font.family"] = "Times New Roman"

# Load MD data
file_md = "md/"
md_0 = pd.read_csv(file_md + "md_0.csv", header=None).values
md_025 = pd.read_csv(file_md + "md_025.csv", header=None).values
md_05 = pd.read_csv(file_md + "md_05.csv", header=None).values
md_075 = pd.read_csv(file_md + "md_075.csv", header=None).values

# Load LBM Data
folder_sm = ["single_grid_results/", "multi_grid_results/"]

legends = ["single-grid", "multi-grid"]

# Load single-grid LBM data
nx, ny = 1024, 1024
x_norm_single = np.linspace(0, 1, ny - 2)
ux_single = np.loadtxt(folder_sm[0] + "ux.txt").reshape((ny, nx))

ux_single_0 = np.flipud(ux_single[1:-1, 0])  # Reverse the order of values
ux_single_0_norm = ux_single_0 / np.mean(ux_single_0)

ux_single_05 = np.flipud(ux_single[1:-1, nx//2])  # Reverse the order of values
ux_single_05_norm = ux_single_05 / np.mean(ux_single_05)

# Load multi-grid LBM data
nx, ny = 1024, 1024
x_norm_multi = np.linspace(0, 1, ny - 2)
ux_multi = np.loadtxt(folder_sm[1] + "ux.txt").reshape((ny, nx))

ux_multi_0 = np.flipud(ux_multi[1:-1, 0])  # Reverse the order of values
ux_multi_0_norm = ux_multi_0 / np.mean(ux_multi_0)

ux_multi_05 = np.flipud(ux_multi[1:-1, nx//2])  # Reverse the order of values
ux_multi_05_norm = ux_multi_05 / np.mean(ux_multi_05)

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# First subplot (x/H = 0)
ax1.plot(md_0[:, 0], md_0[:, 1], 'o', label='MD', linewidth=1.5, markersize=7)
ax1.plot(x_norm_single, ux_single_0_norm, label=legends[0], linewidth=1.5, color='r')
ax1.plot(x_norm_multi, ux_multi_0_norm, label=legends[1], linewidth=1.5, color='k')
ax1.set_xlabel("Normalized Position", fontsize=12)
ax1.set_ylabel("Normalized Velocity", fontsize=12)
ax1.set_ylim([0, 2])
ax1.set_xlim([0, 1])
ax1.legend(loc="upper right", frameon=False, fontsize=12)
ax1.set_title('x/H = 0', fontsize=14)
ax1.tick_params(axis='both', labelsize=12)

# Second subplot (x/H = 0.5)
ax2.plot(md_05[:, 0], md_05[:, 1], 'o', label='MD', linewidth=1.5, markersize=7)
ax2.plot(x_norm_single, ux_single_05_norm, label=legends[0], linewidth=1.5, color='r')
ax2.plot(x_norm_multi, ux_multi_05_norm, label=legends[1], linewidth=1.5, color='k')
ax2.set_xlabel("Normalized Position", fontsize=12)
ax2.set_ylabel("Normalized Velocity", fontsize=12)
ax2.set_ylim([0, 2])
ax2.set_xlim([0, 1])
ax2.legend(loc="upper right", frameon=False, fontsize=12)
ax2.set_title('x/H = 0.5', fontsize=14)
ax2.tick_params(axis='both', labelsize=12)

plt.tight_layout()

# Single plot
plt.figure()
plt.plot(md_0[:, 0], md_0[:, 1], 'o', label='MD', linewidth=1.5, markersize=7)
plt.plot(x_norm_single, ux_single_0_norm, label=legends[0], linewidth=1.5, color='r')
plt.plot(x_norm_multi, ux_multi_0_norm, label=legends[1], linewidth=1.5, color='k')
plt.xlabel("Normalized Position", fontsize=12)
plt.ylabel("Normalized Velocity", fontsize=12)
plt.ylim([0, 2])
plt.xlim([0, 1])
plt.legend(loc="upper right", frameon=False, fontsize=12)
plt.tick_params(axis='both', labelsize=12)

plt.show()
