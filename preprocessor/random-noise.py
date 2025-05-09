import matplotlib.pyplot as plt
from perlin_noise import PerlinNoise
import numpy as np
from scipy.ndimage import gaussian_filter, laplace

noise = PerlinNoise(octaves=2, seed=42)
xpix, ypix = 512, 512
lim_x, lim_y = 10, 10

# Create a grid of coordinates scaled by the limits
x = np.linspace(0, lim_x, xpix)
y = np.linspace(0, lim_y, ypix)
X, Y = np.meshgrid(x, y)

# Generate noise values for the entire grid using list comprehension
pic = np.vectorize(lambda xi, yi: noise([xi, yi], tile_sizes=[20, 5]))(X, Y)
pic = np.transpose(pic)
smoothed_noise = gaussian_filter(pic, sigma=5.0)
edges = np.abs(laplace(smoothed_noise))
threshold = 1.7*np.mean(edges)
binary_structure = (edges > threshold).astype(float)


plt.imshow(binary_structure, cmap="gray")
plt.show()

