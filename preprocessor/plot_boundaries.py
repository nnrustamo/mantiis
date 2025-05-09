import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Times New Roman'
font_size = 14

bnd_to_be_plotted = [127, 445, 473, 383]
bnd1 = np.array([[0, 1, 1], [0, 1, 1], [1, 1, 1]])
bnd2 = np.array([[1, 1, 1], [0, 1, 1], [1, 1, 1]])
bnd3 = np.array([[1, 1, 1], [1, 1, 0], [0, 1, 1]])
bnd4 = np.array([[1, 0, 0], [1, 1, 0], [1, 1, 1]])

bnd1 = np.logical_not(bnd1).astype(int)
bnd2 = np.logical_not(bnd2).astype(int)
bnd3 = np.logical_not(bnd3).astype(int)
bnd4 = np.logical_not(bnd4).astype(int)

fig, axs = plt.subplots(2, 2, figsize=(6, 6))

axs[0, 0].imshow(bnd1, cmap='binary', aspect='equal', interpolation="nearest")
axs[0, 0].scatter(1, 1, color="black", s=50)
axs[0, 0].set_title(f"Boundary number: {bnd_to_be_plotted[0]}", fontsize = font_size)
axs[0, 0].axis('off')
for i in range(bnd1.shape[0] + 1):
    axs[0, 0].hlines(i - 0.5, -0.5, bnd1.shape[1] - 0.5, color="gray", linewidth=0.5)
for j in range(bnd1.shape[1] + 1):
    axs[0, 0].vlines(j - 0.5, -0.5, bnd1.shape[0] - 0.5, color="gray", linewidth=0.5)

axs[0, 1].imshow(bnd2, cmap='binary', aspect='equal', interpolation="nearest")
axs[0, 1].scatter(1, 1, color="black", s=50)
axs[0, 1].set_title(f"Boundary number: {bnd_to_be_plotted[1]}", fontsize = font_size)
axs[0, 1].axis('off')
for i in range(bnd2.shape[0] + 1):
    axs[0, 1].hlines(i - 0.5, -0.5, bnd2.shape[1] - 0.5, color="gray", linewidth=0.5)
for j in range(bnd2.shape[1] + 1):
    axs[0, 1].vlines(j - 0.5, -0.5, bnd2.shape[0] - 0.5, color="gray", linewidth=0.5)

axs[1, 0].imshow(bnd3, cmap='binary', aspect='equal', interpolation="nearest")
axs[1, 0].scatter(1, 1, color="black", s=50)
axs[1, 0].set_title(f"Boundary number: {bnd_to_be_plotted[2]}", fontsize = font_size)
axs[1, 0].axis('off')
for i in range(bnd3.shape[0] + 1):
    axs[1, 0].hlines(i - 0.5, -0.5, bnd3.shape[1] - 0.5, color="gray", linewidth=0.5)
for j in range(bnd3.shape[1] + 1):
    axs[1, 0].vlines(j - 0.5, -0.5, bnd3.shape[0] - 0.5, color="gray", linewidth=0.5)

axs[1, 1].imshow(bnd4, cmap='binary', aspect='equal', interpolation="nearest")
axs[1, 1].scatter(1, 1, color="black", s=50)
axs[1, 1].set_title(f"Boundary number: {bnd_to_be_plotted[3]}", fontsize = font_size)
axs[1, 1].axis('off')
for i in range(bnd4.shape[0] + 1):
    axs[1, 1].hlines(i - 0.5, -0.5, bnd4.shape[1] - 0.5, color="gray", linewidth=0.5)
for j in range(bnd4.shape[1] + 1):
    axs[1, 1].vlines(j - 0.5, -0.5, bnd4.shape[0] - 0.5, color="gray", linewidth=0.5)

plt.tight_layout()
# plt.savefig("bnd.png")
plt.savefig("bnd.svg")
plt.show()
