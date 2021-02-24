import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

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

plt.show()