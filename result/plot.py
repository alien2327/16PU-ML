import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

data_set = np.loadtxt("./final.dat", dtype="float", delimiter=",")

fig = plt.figure(figsize=(6,6))

x, y = [], []
for data in data_set:
    x.append(data[0])
    y.append(data[1])

ax1 = fig.add_subplot(1,1,1)
H = ax1.hist2d(x, y, bins=(75,75), cmap=cm.inferno)
ax1.set_xlabel("x[mm]")
ax1.set_ylabel("y[mm]")
#ax1.set_xlim([-82.5, 82.5])
#ax1.set_ylim([-82.5, 82.5])

driver = make_axes_locatable(ax1)
axHistx = driver.append_axes("top", 1.2, pad=0.1, sharex=ax1)
axHisty = driver.append_axes("right", 1.2, pad=0.1, sharey=ax1)

axHistx.xaxis.set_tick_params(labelbottom=False)
axHisty.yaxis.set_tick_params(labelleft=False)
axHistx.hist(x, bins=75)
axHisty.hist(y, bins=75, orientation='horizontal')
plt.show()