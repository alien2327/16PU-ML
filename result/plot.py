import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

import cv2
from skimage.metrics import structural_similarity
from skimage.metrics import mean_squared_error

fig = plt.figure(figsize=(6,6))

data_set = np.loadtxt("true.dat", dtype="float", delimiter=",")

x, y = [], []
for data in data_set:
    x.append(data[0])
    y.append(data[1])

ax1 = fig.add_subplot(1,1,1)
ax1.hist2d(x, y, bins=100, cmap=cm.inferno)
ax1.set_xlabel("x[mm]")
ax1.set_ylabel("y[mm]")
ax1.set_xlim([-82.5, 82.5])
ax1.set_ylim([-82.5, 82.5])

ax1.axis("off")

plt.savefig("true.png", bbox_inches="tight", pad_inches=0.0)

fig = plt.figure(figsize=(6,6))
data_set = np.loadtxt(f"test.dat", dtype="float", delimiter=",")
x, y = [], []
for data in data_set:
    x.append(data[0])
    y.append(data[1])
ax1 = fig.add_subplot(1,1,1)
ax1.hist2d(x, y, bins=100, cmap=cm.inferno)
ax1.set_xlabel("x[mm]")
ax1.set_ylabel("y[mm]")
ax1.set_xlim([-82.5, 82.5])
ax1.set_ylim([-82.5, 82.5])
ax1.axis("off")
plt.savefig(f"test.png", bbox_inches="tight", pad_inches=0.0)

# Read images from file.
im1 = cv2.imread('./true.png')
im2 = cv2.imread(f'./test.png')

#tempDiff = cv2.subtract(im1, im2)
gray1 = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
gray2 = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)
mse = mean_squared_error(im1, im2)
print(f"MSE : {mse: .5f}")
plt.close()



"""
fig = plt.figure(figsize=(12,6))

data_set = np.loadtxt("true.dat", dtype="float", delimiter=",")

x, y = [], []
for data in data_set:
    x.append(data[0])
    y.append(data[1])

ax1 = fig.add_subplot(1,2,1)
H = ax1.hist2d(x, y, bins=100, cmap=cm.inferno)
ax1.set_xlabel("x[mm]")
ax1.set_ylabel("y[mm]")
ax1.set_xlim([-60, 60])
ax1.set_ylim([-60, 60])

driver1 = make_axes_locatable(ax1)
axHistx1 = driver1.append_axes("top", 1.2, pad=0.1, sharex=ax1)
axHisty1 = driver1.append_axes("right", 1.2, pad=0.1, sharey=ax1)

axHistx1.xaxis.set_tick_params(labelbottom=False)
axHisty1.yaxis.set_tick_params(labelleft=False)
axHistx1.hist(x, bins=100)
axHisty1.hist(y, bins=100, orientation='horizontal')

data_set = np.loadtxt("test.dat", dtype="float", delimiter=",")
x, y = [], []
for data in data_set:
    x.append(data[0])
    y.append(data[1])

ax2 = fig.add_subplot(1,2,2)
H = ax2.hist2d(x, y, bins=100, cmap=cm.inferno)
ax2.set_xlabel("x[mm]")
ax2.set_ylabel("y[mm]")
ax2.set_xlim([-60, 60])
ax2.set_ylim([-60, 60])

driver2 = make_axes_locatable(ax2)
axHistx2 = driver2.append_axes("top", 1.2, pad=0.1, sharex=ax2)
axHisty2 = driver2.append_axes("right", 1.2, pad=0.1, sharey=ax2)

axHistx2.xaxis.set_tick_params(labelbottom=False)
axHisty2.yaxis.set_tick_params(labelleft=False)
axHistx2.hist(x, bins=100)
axHisty2.hist(y, bins=100, orientation='horizontal')

"""