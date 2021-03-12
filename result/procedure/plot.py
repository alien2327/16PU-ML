import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


for i in range(50):
    fig = plt.figure(figsize=(6,6))
    data_set = np.loadtxt(f"test_{i}.dat", dtype="float", delimiter=",")
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
    plt.savefig(f"test_{i}.png", bbox_inches="tight", pad_inches=0.0)
    plt.close()