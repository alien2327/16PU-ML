import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

b = 100

data_set_init = np.loadtxt("./init.csv", dtype="float", delimiter=",")
data_set_fin = np.loadtxt("./final.csv", dtype="float", delimiter=",")

fig = plt.figure(figsize=(12,6))

x_init, y_init = [], []
x_fin, y_fin = [], []
for data in data_set_fin:
    x_fin.append(data[0])
    y_fin.append(data[1])
for data in data_set_init:
    x_init.append(data[0])
    y_init.append(data[1])

ax1 = fig.add_subplot(1,2,1)
ax1.hist2d(x_init, y_init, bins=b, cmap=cm.plasma)
ax1.set_xlabel("x[mm]")
ax1.set_ylabel("y[mm]")
ax1.set_xlim([-82.5, 82.5])
ax1.set_ylim([-82.5, 82.5])

driver = make_axes_locatable(ax1)
axHistx = driver.append_axes("top", 1.2, pad=0.1, sharex=ax1)
axHisty = driver.append_axes("right", 1.2, pad=0.1, sharey=ax1)

axHistx.xaxis.set_tick_params(labelbottom=False)
axHisty.yaxis.set_tick_params(labelleft=False)
axHistx.hist(x_init, bins=b)
axHisty.hist(y_init, bins=b, orientation='horizontal')

ax1 = fig.add_subplot(1,2,2)
ax1.hist2d(x_fin, y_fin, bins=b, cmap=cm.plasma)
ax1.set_xlabel("x[mm]")
ax1.set_ylabel("y[mm]")
ax1.set_xlim([-82.5, 82.5])
ax1.set_ylim([-82.5, 82.5])

driver = make_axes_locatable(ax1)
axHistx = driver.append_axes("top", 1.2, pad=0.1, sharex=ax1)
axHisty = driver.append_axes("right", 1.2, pad=0.1, sharey=ax1)

axHistx.xaxis.set_tick_params(labelbottom=False)
axHisty.yaxis.set_tick_params(labelleft=False)
axHistx.hist(x_fin, bins=b)
axHisty.hist(y_fin, bins=b, orientation='horizontal')

plt.show()