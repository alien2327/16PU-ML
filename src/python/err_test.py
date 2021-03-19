import numpy as np
from scipy.stats import norm
import copy
import math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm

def plot_bunch(true_bunch, try_bunch):
    fig = plt.figure(figsize=(12,6))
    x, y = true_bunch
    x[0] = -82.5
    x[1] = 82.5
    y[0] = 82.5
    y[1] = -82.5
    ax1 = fig.add_subplot(1,2,1)
    H = ax1.hist2d(x, y, bins=100, cmap=cm.inferno)
    ax1.set_xlabel("x[mm]")
    ax1.set_ylabel("y[mm]")
    ax1.set_xlim([-82.5, 82.5])
    ax1.set_ylim([-82.5, 82.5])
    driver1 = make_axes_locatable(ax1)
    axHistx1 = driver1.append_axes("top", 1.2, pad=0.1, sharex=ax1)
    axHisty1 = driver1.append_axes("right", 1.2, pad=0.1, sharey=ax1)
    axHistx1.xaxis.set_tick_params(labelbottom=False)
    axHisty1.yaxis.set_tick_params(labelleft=False)
    axHistx1.hist(x, bins=100)
    axHisty1.hist(y, bins=100, orientation='horizontal')
    x, y = try_bunch
    x[0] = -82.5
    x[1] = 82.5
    y[0] = 82.5
    y[1] = -82.5
    ax2 = fig.add_subplot(1,2,2)
    H = ax2.hist2d(x, y, bins=100, cmap=cm.inferno)
    ax2.set_xlabel("x[mm]")
    ax2.set_ylabel("y[mm]")
    ax2.set_xlim([-82.5, 82.5])
    ax2.set_ylim([-82.5, 82.5])
    driver2 = make_axes_locatable(ax2)
    axHistx2 = driver2.append_axes("top", 1.2, pad=0.1, sharex=ax2)
    axHisty2 = driver2.append_axes("right", 1.2, pad=0.1, sharey=ax2)
    axHistx2.xaxis.set_tick_params(labelbottom=False)
    axHisty2.yaxis.set_tick_params(labelleft=False)
    axHistx2.hist(x, bins=100)
    axHisty2.hist(y, bins=100, orientation='horizontal')
    plt.show()
    return

def gen_bunch(x_mu, x_sigma, y_mu, y_sigma):
    x = np.random.normal(x_mu, x_sigma, 1000)
    y = np.random.normal(y_mu, y_sigma, 1000)
    bunch = np.concatenate((x, y))
    bunch = np.reshape(bunch, (2, 1000))
    return bunch

def trans_method(p, c, dx, dy):
    for i in range(len(p[0])):
        x = copy.copy(p[0][i])
        y = copy.copy(p[1][i])
        if c == 0:
            p[0][i] += dx
            p[1][i] += dy
        elif c == 1:
            p[0][i] *= (1 + dx)
            p[1][i] *= (1 + dy)
        elif c == 2:
            p[0][i] *= (1 + y*dx)
            p[1][i] *= (1 + x*dy)
        elif c == 3:
            p[0][i] *= (1 + abs(y)*dx)
            p[1][i] *= (1 + abs(x)*dy)
        elif c == 4:
            theta = math.atan(dy/dx)
            p[0][i] = x*math.cos(theta) - y*math.sin(theta)
            p[1][i] = x*math.sin(theta) + y*math.cos(theta)
        elif c == 5:
            p[0][i] = x*dx + y*dy
            p[1][i] = x*dy + y*dx
    return p

def get_vol(bunch):
    mat = np.loadtxt('./calibration/wire/A_3-4MHz_13.txt')
    gain = np.loadtxt('./calibration/BBGC/gain_13.txt')
    buf = np.zeros(16)
    vol = np.zeros(16)
    x, y = bunch
    for i in range(len(x)):
        position = x[i] + y[i]*1j
        for j in range(16):
            if (j%2 == 0):
                buf[j] += pow(position, (j+1)/2).imag
            else:
                buf[j] += pow(position, (j+1)/2).real
    for i in range(16):
        buf[i] /= 1000
    for i in range(16):
        for j in range(16):
            vol[i] += mat[i][j] * buf[j]
    vol = np.linalg.norm(vol * gain)
    return vol

def loss(true_bunch, try_bunch):
    w = 0.01
    mse = np.square(np.subtract(get_vol(true_bunch), get_vol(try_bunch))).mean()
    x_true, y_true = true_bunch
    x_try, y_try = try_bunch
    sig_x_true, sig_y_true = np.std(x_true), np.std(y_true)
    sig_x_try, sig_y_try = np.std(x_try), np.std(y_try)
    ex = np.exp(w*pow(sig_x_true - sig_x_try, 2)) + np.exp(w*pow(sig_y_true - sig_y_try, 2))
    print("::")
    print(f":: MSE: {mse: 1.8f}\tEP-2: {ex-2: 1.8f}\tRES: {mse + ex - 2: 1.8f}")
    return mse + ex - 2

if __name__ == "__main__":
    true_bunch = gen_bunch(0, 4, 0, 10)
    trans_method(true_bunch, 4, 1, -1)
    try_bunch = gen_bunch(0, 1, 0, 10)
    trans_method(try_bunch, 4, 1, -1)
    loss(true_bunch, try_bunch)
    print("::")
    print(":: Position X: ", np.mean(true_bunch[0]), "[mm]")
    print(":: Position Y: ", np.mean(true_bunch[1]), "[mm]")
    print(":: Size X: ", np.std(true_bunch[0]), "[mm^2]")
    print(":: Size Y: ", np.std(true_bunch[1]), "[mm^2]")
    print("::")
    print(":: Position X: ", np.mean(try_bunch[0]), "[mm]")
    print(":: Position Y: ", np.mean(try_bunch[1]), "[mm]")
    print(":: Size X: ", np.std(try_bunch[0]), "[mm^2]")
    print(":: Size Y: ", np.std(try_bunch[1]), "[mm^2]")
    print("::")
    plot_bunch(true_bunch, try_bunch)
    print(":: Done")