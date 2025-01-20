import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================== SLIT CHANNEL ============================== #
# DSMC
excel_file = 'DSMC_slit_channel.xlsx'
df_dsmcs_simple = pd.read_excel(excel_file, header=2)
columns_mapping = [
    "kn=0.01-x_norm",
    "kn=0.01-ux_norm",
    "kn=0.1-x_norm",
    "kn=0.1-ux_norm",
    "kn=1.0-x_norm",
    "kn=1.0-ux_norm"
]
df_dsmcs_simple = df_dsmcs_simple.iloc[:, :6] 
df_dsmcs_simple.columns = columns_mapping
print("DSMC SIMPLE CHANNEL DATA")
print (df_dsmcs_simple.head())

# LBM
base_folder = os.path.dirname(os.path.abspath(__file__))
folders = ['dbb', 'srbb']
df_lbm_simple = pd.DataFrame()

for folder in folders:
    folder_path = os.path.join(base_folder, folder, 'slit_channel')
    if os.path.exists(folder_path):
        for kn_folder in ['kn=1.0', 'kn=0.1', 'kn=0.01']:
            kn_path = os.path.join(folder_path, kn_folder)
            if os.path.exists(kn_path):
                ux_path = os.path.join(kn_path, 'ux.txt')
                if os.path.exists(ux_path):
                    ux = np.loadtxt(ux_path).reshape(100, 100)
                    slice_column = ux[:, 0][1:-1]
                    mean_value = np.mean(slice_column)
                    normalized_slice = slice_column / mean_value
                    x = np.linspace(0, 1, len(normalized_slice))
                    col_prefix = f"{folder}-{kn_folder}"
                    df_lbm_simple[f"{col_prefix}-x_norm"] = x
                    df_lbm_simple[f"{col_prefix}-ux_norm"] = normalized_slice

print("LBM SIMPLE CHANNEL DATA")
print(df_lbm_simple.head())

# PLOTTING
plt.rcParams.update({
    'font.family': 'Times New Roman',
    'font.size': 14,
    'axes.facecolor': '#f8f8f8',  # Light gray background
    'axes.edgecolor': 'black',
    'grid.color': '#cccccc',      # Subtle gray grid lines
    'grid.alpha': 0.8,
})


fig, axes = plt.subplots(1, 3, figsize=(18, 6))
kn_values = ["kn=1.0", "kn=0.1", "kn=0.01"]

for i, kn in enumerate(kn_values):
    x_dsmc = df_dsmcs_simple[f"{kn}-x_norm"]
    ux_dsmc = df_dsmcs_simple[f"{kn}-ux_norm"]
    axes[i].scatter(
        x_dsmc, ux_dsmc, label="DSMC", edgecolor='blue', facecolor='none', marker='o', s=70, linewidth=1.5
    )

    x_lbm_dbb = df_lbm_simple[f"dbb-{kn}-x_norm"]
    ux_lbm_dbb = df_lbm_simple[f"dbb-{kn}-ux_norm"]
    x_lbm_srbb = df_lbm_simple[f"srbb-{kn}-x_norm"]
    ux_lbm_srbb = df_lbm_simple[f"srbb-{kn}-ux_norm"]
    axes[i].plot(x_lbm_dbb, ux_lbm_dbb, label="LBM (DBB)", color='black', linestyle='-', linewidth=2)
    axes[i].plot(x_lbm_srbb, ux_lbm_srbb, label="LBM (SRBB)", color='darkred', linestyle='--', linewidth=2)

    axes[i].set_title(f"Kn = {kn.split('=')[1]}", fontsize=16)
    axes[i].set_xlabel("Normalized Position", fontsize=16)
    axes[i].set_ylabel("Normalized Velocity", fontsize=16)
    # axes[i].set_xlim(0, 1)
    axes[i].set_ylim(0, 1.5)
    axes[i].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    axes[i].legend(loc='center', fontsize=12, frameon=False)
    axes[i].tick_params(labelsize=14)

plt.tight_layout()
plt.savefig("simple_channel_kn_plot.png", dpi=300)

# ============================== SQUARE OBSTACLE ============================== #

# MD
excel_file = 'MD_square_obstacle.xlsx'
df_md_square = pd.read_excel(excel_file, header=2)
columns_mapping = [
    "x=0-x_norm",
    "x=0-ux_norm",
    "x=0.25-x_norm",
    "x=0.25-ux_norm",
    "x=0.5-x_norm",
    "x=0.5-ux_norm",
    "x=0.75-x_norm",
    "x=0.75-ux_norm"
]

df_md_square.columns = columns_mapping
print("MD SUARE OBSTACLE DATA")
print (df_md_square.head())

# LBM
base_folder = os.path.dirname(os.path.abspath(__file__))
folders = ['dbb', 'srbb']
df_lbm_square = pd.DataFrame()

for folder in folders:
    folder_path = os.path.join(base_folder, folder, 'square_obstacle')
    if os.path.exists(folder_path):
        if os.path.exists(folder_path):
            ux_path = os.path.join(folder_path, 'ux.txt')
            if os.path.exists(ux_path):
                ux = np.loadtxt(ux_path).reshape(500, 500)
                x = np.linspace(0, 1, 500-2)

                ux_0 = ux[:, 0][1:-1]
                col_prefix = f"{folder}-0"
                mean_value = np.mean(ux_0)
                normalized_slice = ux_0 / mean_value
                df_lbm_square[f"{col_prefix}-x_norm"] = x
                df_lbm_square[f"{col_prefix}-ux_norm"] = normalized_slice

                ux_025 = ux[:, 124][1:-1]
                col_prefix = f"{folder}-0.25"
                mean_value = np.mean(ux_025)
                normalized_slice = ux_025 / mean_value
                df_lbm_square[f"{col_prefix}-x_norm"] = x
                df_lbm_square[f"{col_prefix}-ux_norm"] = normalized_slice
                
                ux_05 = ux[:, 249][1:-1].astype(float) 
                mean_value = np.mean(ux_05)
                ux_05[ux_05 == 0] = np.nan
                col_prefix = f"{folder}-0.5"
                normalized_slice = ux_05 / mean_value
                df_lbm_square[f"{col_prefix}-x_norm"] = x
                df_lbm_square[f"{col_prefix}-ux_norm"] = normalized_slice

                x = np.linspace(0, 1, 500-2)
                ux_075 = ux[:, 374][1:-1]
                col_prefix = f"{folder}-0.75"
                mean_value = np.mean(ux_075)
                normalized_slice = ux_075 / mean_value
                df_lbm_square[f"{col_prefix}-x_norm"] = x
                df_lbm_square[f"{col_prefix}-ux_norm"] = normalized_slice

print("LBM SQUARE OBSTACLE DATA")
print(df_lbm_square.head())

# # PLOTTING
plt.rcParams.update({
    'font.family': 'Times New Roman',
    'font.size': 14,
    'axes.facecolor': '#f8f8f8',  # Light gray background
    'axes.edgecolor': 'black',
    'grid.color': '#cccccc',      # Subtle gray grid lines
    'grid.alpha': 0.8,
})


fig, axes = plt.subplots(1, 4, figsize=(18, 6))
xn_values = ["0", "0.25", "0.5", "0.75"]

for i, xn in enumerate(xn_values):
    x_md = df_md_square[f"x={xn}-x_norm"]
    ux_md = df_md_square[f"x={xn}-ux_norm"]
    axes[i].scatter(
        ux_md, x_md, label="MD", edgecolor='blue', facecolor='none', marker='o', s=70, linewidth=1.5
    )

    x_lbm_dbb = df_lbm_square[f"dbb-{xn}-x_norm"]
    ux_lbm_dbb = df_lbm_square[f"dbb-{xn}-ux_norm"]
    x_lbm_srbb = df_lbm_square[f"srbb-{xn}-x_norm"]
    ux_lbm_srbb = df_lbm_square[f"srbb-{xn}-ux_norm"]
    axes[i].plot(ux_lbm_dbb, x_lbm_dbb, label="LBM (DBB)", color='black', linestyle='-', linewidth=2)
    axes[i].plot(ux_lbm_srbb, x_lbm_srbb, label="LBM (SRBB)", color='darkred', linestyle='--', linewidth=2)

    axes[i].set_title(f"x/H = {xn}", fontsize=16)
    axes[i].set_xlabel("Normalized Velocity", fontsize=16)
    axes[i].set_ylabel("Normalized Position", fontsize=16)
    # axes[i].set_xlim(0, 1)
    # axes[i].set_ylim(0, 2.0)
    axes[i].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    if i < 1:
        axes[i].legend(loc='best', fontsize=12, frameon=False)
    else:
        axes[i].legend(loc='right', fontsize=12, frameon=False)
    axes[i].tick_params(labelsize=14)

plt.tight_layout()
plt.savefig("square_obstacle.png", dpi=300)

# ============================== TRIANGULAR OBSTACLE ============================== #
# MD
excel_file = 'MD_triangular_obstacle.xlsx'
df_md_triangles = pd.read_excel(excel_file, header=2)
columns_mapping = [
    "x=0-x_norm",
    "x=0-ux_norm",
    "x=0.25-x_norm",
    "x=0.25-ux_norm",
    "x=0.5-x_norm",
    "x=0.5-ux_norm",
    "x=0.75-x_norm",
    "x=0.75-ux_norm"
]

df_md_triangles.columns = columns_mapping
print("MD TRIANGULAR OBSTACLE DATA")
print (df_md_triangles.head())

# LBM
base_folder = os.path.dirname(os.path.abspath(__file__))
folders = ['dbb', 'srbb']
df_lbm_triangle = pd.DataFrame()

for folder in folders:
    folder_path = os.path.join(base_folder, folder, 'triangle_obstacle')
    if os.path.exists(folder_path):
        if os.path.exists(folder_path):
            ux_path = os.path.join(folder_path, 'ux.txt')
            if os.path.exists(ux_path):
                ux = np.loadtxt(ux_path).reshape(500, 500)
                x = np.linspace(0, 1, 500-2)

                ux_0 = ux[:, 0][1:-1]
                col_prefix = f"{folder}-0"
                mean_value = np.mean(ux_0)
                normalized_slice = ux_0 / mean_value
                df_lbm_triangle[f"{col_prefix}-x_norm"] = x
                df_lbm_triangle[f"{col_prefix}-ux_norm"] = np.flip(normalized_slice)

                ux_025 = ux[:, 124][1:-1]
                col_prefix = f"{folder}-0.25"
                mean_value = np.mean(ux_025)
                normalized_slice = ux_025 / mean_value
                df_lbm_triangle[f"{col_prefix}-x_norm"] = x
                df_lbm_triangle[f"{col_prefix}-ux_norm"] = np.flip(normalized_slice)
                
                ux_05 = ux[:, 249][1:-1].astype(float) 
                mean_value = np.mean(ux_05)
                ux_05[ux_05 == 0] = np.nan
                col_prefix = f"{folder}-0.5"
                normalized_slice = ux_05 / mean_value
                df_lbm_triangle[f"{col_prefix}-x_norm"] = x
                df_lbm_triangle[f"{col_prefix}-ux_norm"] = np.flip(normalized_slice)

                x = np.linspace(0, 1, 500-2)
                ux_075 = ux[:, 374][1:-1]
                col_prefix = f"{folder}-0.75"
                mean_value = np.mean(ux_075)
                normalized_slice = ux_075 / mean_value
                df_lbm_triangle[f"{col_prefix}-x_norm"] = x
                df_lbm_triangle[f"{col_prefix}-ux_norm"] = np.flip(normalized_slice)

print("LBM TRIANGULAR OBSTACLE DATA")
print(df_lbm_triangle.head())

# # PLOTTING
plt.rcParams.update({
    'font.family': 'Times New Roman',
    'font.size': 14,
    'axes.facecolor': '#f8f8f8',  # Light gray background
    'axes.edgecolor': 'black',
    'grid.color': '#cccccc',      # Subtle gray grid lines
    'grid.alpha': 0.8,
})


fig, axes = plt.subplots(1, 4, figsize=(18, 6))
xn_values = ["0", "0.25", "0.5", "0.75"]

for i, xn in enumerate(xn_values):
    x_md = df_md_triangles[f"x={xn}-x_norm"]
    ux_md = df_md_triangles[f"x={xn}-ux_norm"]
    axes[i].scatter(
        ux_md, x_md, label="MD", edgecolor='blue', facecolor='none', marker='o', s=70, linewidth=1.5
    )

    x_lbm_dbb = df_lbm_triangle[f"dbb-{xn}-x_norm"]
    ux_lbm_dbb = df_lbm_triangle[f"dbb-{xn}-ux_norm"]
    x_lbm_srbb = df_lbm_triangle[f"srbb-{xn}-x_norm"]
    ux_lbm_srbb = df_lbm_triangle[f"srbb-{xn}-ux_norm"]
    axes[i].plot(ux_lbm_dbb, x_lbm_dbb, label="LBM (DBB)", color='black', linestyle='-', linewidth=2)
    axes[i].plot(ux_lbm_srbb, x_lbm_srbb, label="LBM (SRBB)", color='darkred', linestyle='--', linewidth=2)

    axes[i].set_title(f"x/H = {xn}", fontsize=16)
    axes[i].set_xlabel("Normalized Velocity", fontsize=16)
    axes[i].set_ylabel("Normalized Position", fontsize=16)
    # axes[i].set_xlim(0, 1)
    # axes[i].set_ylim(0, 2.0)
    axes[i].grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    if i < 1:
        axes[i].legend(loc='best', fontsize=12, frameon=False)
    else:
        axes[i].legend(loc='right', fontsize=12, frameon=False)
    axes[i].tick_params(labelsize=14)

plt.tight_layout()
plt.savefig("triangular_obstacle.png", dpi=300)
