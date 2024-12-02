import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker

plt.rc('font', family='serif')

# Dimensions
nx = 1024
ny = 1024

# ========================================================= Load and plot velocity profile at 2 MPa
ux_mg_2mpa = np.loadtxt("mg/large_media_2MPa/ux.txt")
ux_mg_2mpa =ux_mg_2mpa.reshape(nx, ny)
uy_mg_2mpa = np.loadtxt("mg/large_media_2MPa/uy.txt")
uy_mg_2mpa = uy_mg_2mpa.reshape(nx, ny)

ux_sg_2mpa = np.loadtxt("sg/large_media_2MPa/ux.txt")
ux_sg_2mpa =ux_sg_2mpa.reshape(nx, ny)
uy_sg_2mpa = np.loadtxt("sg/large_media_2MPa/uy.txt")
uy_sg_2mpa = uy_sg_2mpa.reshape(nx, ny)

# u_sg_2mpa_safe = np.where(u_sg_2mpa == 0, np.finfo(float).eps, u_sg_2mpa)
# u_mg_2mpa_safe = np.where(u_mg_2mpa == 0, np.finfo(float).eps, u_mg_2mpa)
# diff = abs(u_mg_2mpa - u_sg_2mpa)
# diff = abs(u_mg_2mpa_safe - u_sg_2mpa_safe)/u_sg_2mpa_safe * 100
# diff = np.where(diff > 100, 0, diff)

# plt.figure()
# plt.imshow(u_sg_2mpa, cmap='jet', origin='lower')
# cbar = plt.colorbar()
# cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# plt.tight_layout()
# # plt.show()
# plt.savefig("u_sg_2MPa.png")

# plt.figure()
# plt.imshow(u_mg_2mpa, cmap='jet', origin='lower')
# cbar = plt.colorbar()
# cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# plt.tight_layout()
# # plt.show()
# plt.savefig("u_mg_2MPa.png")

# plt.figure()
# plt.imshow(diff, cmap='jet', origin='lower')
# cbar = plt.colorbar()
# cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
# plt.tight_layout()
# # plt.show()
# plt.savefig("u_diff_2MPa.png")

# Contour plot for multigrid
x = np.linspace(0, 1024, 1024)
y = np.linspace(0, 1024, 1024)
X, Y = np.meshgrid(x, y)
U = ux_mg_2mpa
V = uy_mg_2mpa

# Calculate the magnitude of the velocity
velocity_magnitude = np.sqrt(U**2 + V**2)
velocity_magnitude_filtered = np.where(velocity_magnitude < 6.0e-6, np.nan, velocity_magnitude)

plt.figure()
contour = plt.contourf(X, Y, velocity_magnitude, 1000, cmap='jet') # this is practically heatmap now
cbar = plt.colorbar(contour, label='Magnitude of velocity (m/s)')
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
plt.savefig("ux_mg_2mpa_contour.png")

plt.figure()
contour = plt.contourf(X, Y, velocity_magnitude_filtered, 1000, cmap='jet')
cbar = plt.colorbar(contour, label='Magnitude of velocity (m/s)')
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
plt.savefig("ux_mg_2mpa_contour_filtered.png")

# Error Distribution
u_mg_2mpa = np.sqrt(np.power(ux_mg_2mpa, 2) + np.power(uy_mg_2mpa, 2))
u_sg_2mpa = np.sqrt(np.power(ux_sg_2mpa, 2) + np.power(uy_sg_2mpa, 2))
u_mg_2mpa_safe = np.where(u_mg_2mpa == 0, np.finfo(float).eps, u_mg_2mpa)
u_sg_2mpa_safe = np.where(u_sg_2mpa == 0, np.finfo(float).eps, u_sg_2mpa)
diff = abs(u_mg_2mpa_safe - u_sg_2mpa_safe)/u_sg_2mpa_safe * 100
diff = diff.reshape(nx*ny, 1)
# diff = diff[np.where(diff>1.0e-1)]

plt.figure()
plt.hist(diff, bins=1000, color='steelblue', edgecolor='black', alpha=0.75)
plt.xlabel('Error (%)', fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Frequency', fontsize=14)
# plt.grid(True, linestyle='--', alpha=0.6)
mean_error = np.mean(diff)
plt.axvline(mean_error, color='red', linestyle='dashed', linewidth=2, label=f'Mean = {round(mean_error, 2)} %')
plt.legend()
plt.tight_layout()
plt.savefig("Errors.png")


# ========================================================= Load and plot convergence profiles
conv_mg_2mpa = np.loadtxt("mg/large_media_2MPa/convergence.txt")
conv_sg_2mpa = np.loadtxt("sg/large_media_2MPa/convergence.txt")

plt.figure()
plt.plot(np.arange(len(conv_sg_2mpa)), conv_sg_2mpa,  color='red', label='single-grid')
plt.plot(np.arange(len(conv_mg_2mpa)), conv_mg_2mpa,  color='black', label='multi-grid')
plt.legend()
plt.yscale('log')
plt.xlabel("Iterations")
plt.ylabel("Error")
plt.savefig("conv_2MPa.png")

# ========================================================= Calculate permeability
T = 300 # K
L = 1024e-9 # meters
NA = 6.02214076e23
R = 8.31446; 
T_c = 190.564  # K, critical temperature
P_c = 4.5992 * 1.0e6  # Pa, critical pressure
M = 16.043 * 1.0e-3  # methane molecular weight, kg/mol
a = 0.45724 * (R ** 2 * T_c ** 2) / P_c 
b = 0.07780 * (R * T_c) / P_c

def peng_robinson_pressure_methane(rho, T):
    Vm = M / rho
    pressure = (R * T) / (Vm - b) - (a / (Vm * (Vm + b) + b * (Vm - b)))
    return pressure

# Simulation mean pressure
P = [2, 3, 4, 5, 6]    # MPa
P_str = ['2', '3', '4', '5', '6']    # MPa
inv_P = []
# Dynamic visc of methane at P: http://www.peacesoftware.de/einigewerte/methan_e.html

mu = [11.5e-6, 11.7e-6, 11.9e-6, 12.2e-6, 12.5e-6]    # Pa s
k_mg = []  # m2
k_sg = []  # m2
diff = []
Cf = 7.72066669163199e+17   # force onversion factor 
F = 1.0e-10  # dimensionless force

for i in range(5):
    # Load single grid
    sg_ux = "sg/large_media_" + str(P_str[i]) + "MPa/" + "ux.txt"
    ux_sg = np.loadtxt(sg_ux)
    ux_sg = ux_sg.reshape(nx, ny)

    sg_rho = "sg/large_media_" + str(P_str[i]) + "MPa/" + "rho.txt"
    rho_sg = np.loadtxt(sg_rho)
    rho_sg = rho_sg.reshape(nx, ny)
    
    # Caluclate single grid
    ux_avg_o = np.mean(ux_sg[1:-1, -1]); # kg/m3

    avg_inlet_rho = np.mean(rho_sg[1:-1, 0]); # kg/m3
    inlet_pressure = peng_robinson_pressure_methane(avg_inlet_rho, T)   # Pa

    avg_outlet_rho = np.mean(rho_sg[1:-1, -1]); # kg/m3
    outlet_pressure = peng_robinson_pressure_methane(avg_outlet_rho, T) # Pa

    # Permeability 
    kg_sg = mu[i] * ux_avg_o / (Cf * F)   # 2*mu[i] * L * ux_avg_o * outlet_pressure / (inlet_pressure**2 - outlet_pressure**2)
    #
    k_sg.append(kg_sg)

    # Multigrid
    mg_ux = "mg/large_media_" + str(P[i]) + "MPa/" + "ux.txt"
    ux_mg = np.loadtxt(mg_ux)
    ux_mg = ux_mg.reshape(nx, ny)

    mg_rho = "mg/large_media_" + str(P[i]) + "MPa/" + "rho.txt"
    rho_mg = np.loadtxt(mg_rho)
    rho_mg = rho_mg.reshape(nx, ny)
    
    # Caluclate single grid
    ux_avg_o = np.mean(ux_mg[1:-1, -1]); # kg/m3

    avg_inlet_rho = np.mean(rho_mg[1:-1, 0]); # kg/m3
    inlet_pressure = peng_robinson_pressure_methane(avg_inlet_rho, T)   # Pa

    avg_outlet_rho = np.mean(rho_mg[1:-1, -1]); # kg/m3
    outlet_pressure = peng_robinson_pressure_methane(avg_outlet_rho, T) # Pa

    # Permeability
    kg_mg = mu[i] * ux_avg_o / (Cf * F) # 2*mu[i] * L * ux_avg_o * outlet_pressure / (inlet_pressure**2 - outlet_pressure**2)
    # mu[i] * ux_avg_o / (Cf * F) #
    k_mg.append(kg_mg)

    err = abs(kg_mg - kg_sg) / kg_sg * 100
    diff.append(err)
    inv_P.append(1/(P[i] * 1.0e6))
    print(f"Pressure: {P[i]}, k_sg: {kg_sg}, k_mg: {kg_mg}")


k_sg = np.array(k_sg)
k_mg = np.array(k_mg)
inv_P = np.array(inv_P)
diff = np.array(diff)

fig, ax1 = plt.subplots()
line1, = ax1.plot(inv_P, k_sg, marker='o', color='red', label='single-grid')
line2, = ax1.plot(inv_P, k_mg, marker='s', color='black', label='multi-grid')

ax1.set_yscale('log') 
ax1.set_xscale('log')
ax1.set_xlabel(r"$1/P  \, (1/\mathrm{Pa})$", fontsize=12)
ax1.set_ylabel(r"$k  \, (\mathrm{m^2})$", fontsize=12)

ax1.set_xlim([min(inv_P)* 0.95, max(inv_P)* 1.05])
ax1.set_ylim([min(min(k_sg), min(k_mg)) * 0.95, max(max(k_sg), max(k_mg)) * 1.05])
ax1.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
ax1.yaxis.set_major_formatter(ticker.LogFormatterSciNotation())

lines = [line1, line2]
labels = [line.get_label() for line in lines]
plt.legend(lines, labels, loc='upper center')
# ax1.grid(True, which='both', axis='both')
plt.tight_layout()
plt.rc('font', family='serif')
plt.savefig("permeability.png")
