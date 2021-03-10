import matplotlib.pyplot as plt
import multiprocessing as mp
import matplotlib.cm as cm
import numpy as np
import cv2
from skimage.metrics import structural_similarity
from skimage.metrics import mean_squared_error

MAX_GEN = 20

def plot_bunch(true_bunch, try_bunch):
    fig = plt.figure(figsize=(6,6))
    x, y = true_bunch
    x[0] = -82.5
    x[1] = 82.5
    y[0] = 82.5
    y[1] = -82.5
    ax1 = fig.add_subplot(1,1,1)
    ax1.hist2d(x, y, bins=100, cmap=cm.inferno)
    ax1.set_xlabel("x[mm]")
    ax1.set_ylabel("y[mm]")
    ax1.set_xlim([-82.5, 82.5])
    ax1.set_ylim([-82.5, 82.5])
    ax1.axis("off")
    plt.savefig("true.png", bbox_inches="tight", pad_inches=0.0)
    fig = plt.figure(figsize=(6,6))
    x, y = try_bunch
    x[0] = -82.5
    x[1] = 82.5
    y[0] = 82.5
    y[1] = -82.5
    ax1 = fig.add_subplot(1,1,1)
    ax1.hist2d(x, y, bins=100, cmap=cm.inferno)
    ax1.set_xlabel("x[mm]")
    ax1.set_ylabel("y[mm]")
    ax1.set_xlim([-82.5, 82.5])
    ax1.set_ylim([-82.5, 82.5])
    ax1.axis("off")
    plt.savefig("test.png", bbox_inches="tight", pad_inches=0.0)
    im1 = cv2.imread('./true.png')
    im2 = cv2.imread('./test.png')
    gray1 = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
    gray2 = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)
    mse = mean_squared_error(gray1, gray2)
    score = structural_similarity(gray1, gray2)
    print("::")
    print(f":: MSE : {mse: .5f}")
    print(f":: SSIM: {score: .5f}")

def gen_bunch():
    maxR = 50
    minR = -50
    maxR2 = 10
    minR2 = 1
    x_mu = minR + (maxR - minR) * np.random.rand()
    x_sigma = minR2 + (maxR2 - minR2) * np.random.rand()
    y_mu = minR + (maxR - minR) * np.random.rand()
    y_sigma = minR2 + (maxR2 - minR2) * np.random.rand()
    x = np.random.normal(x_mu, x_sigma, 1000)
    y = np.random.normal(y_mu, y_sigma, 1000)
    bunch = np.concatenate((x, y))
    bunch = np.reshape(bunch, (2, 1000))
    return bunch

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
    return mse + ex - 2

def trans_method(p, c):
    dx = 0.001 + (0.1 - 0.001) * np.random.rand()
    dy = 0.001 + (0.1 - 0.001) * np.random.rand()
    for i in range(len(p[0])):
        x, y = p[0][i], p[1][i]
        if c == 0:
            p[0][i] += dx
            p[1][i] += dy
        elif c == 1:
            p[0][i] -= dx
            p[1][i] += dy
        elif c == 2:
            p[0][i] += dx
            p[1][i] -= dy
        elif c == 3:
            p[0][i] -= dx
            p[1][i] -= dy
        elif c == 4:
            p[0][i] *= (1 + dx)
            p[1][i] *= (1 + dy)
        elif c == 5:
            p[0][i] *= (1 - dx)
            p[1][i] *= (1 + dy)
        elif c == 6:
            p[0][i] *= (1 + dx)
            p[1][i] *= (1 - dy)
        elif c == 7:
            p[0][i] *= (1 - dx)
            p[1][i] *= (1 - dy)
        elif c == 8:
            p[0][i] *= (1 + y*dx)
            p[1][i] *= (1 + x*dy)
        elif c == 9:
            p[0][i] *= (1 - y*dx)
            p[1][i] *= (1 + x*dy)
        elif c == 10:
            p[0][i] *= (1 + y*dx)
            p[1][i] *= (1 - x*dy)
        elif c == 11:
            p[0][i] *= (1 - y*dx)
            p[1][i] *= (1 - x*dy)
        elif c == 12:
            p[0][i] *= (1 + abs(y)*dx)
            p[1][i] *= (1 + abs(x)*dy)
        elif c == 13:
            p[0][i] *= (1 - abs(y)*dx)
            p[1][i] *= (1 + abs(x)*dy)
        elif c == 14:
            p[0][i] *= (1 + abs(y)*dx)
            p[1][i] *= (1 - abs(x)*dy)
        elif c == 15:
            p[0][i] *= (1 - abs(y)*dx)
            p[1][i] *= (1 - abs(x)*dy)
    return p

def trans_beam(true_bunch, bunch):
    loss_value = []
    res_bunch = [bunch]
    for i in range(16):
        res_bunch.append(trans_method(bunch, i))
    for i, try_bunch in enumerate(res_bunch):
        loss_value.append((i, loss(true_bunch, try_bunch)))
    loss_value = sorted(loss_value, key = lambda x: x[1])
    return res_bunch[loss_value[0][0]]

def initial_seq(true_bunch):
    candidate_bunch = []
    candidate_bunch_selected = []
    loss_value = []
    for i in range(100):
        candidate_bunch.append(gen_bunch())
        loss_value.append((i, loss(true_bunch, candidate_bunch[i])))
    loss_value = sorted(loss_value, key=lambda x: x[1])
    for i in range(100):
        candidate_bunch_selected.append(candidate_bunch[loss_value[int(i/10)][0]])
    candidate_bunch_selected = mutate(candidate_bunch_selected)
    return candidate_bunch_selected

def mutate(bunch_list):
    mrate = 0.01
    for i in range(len(bunch_list)):
        if mrate > np.random.rand():
            bunch_list[i] = gen_bunch()
    return bunch_list

def generation(true_bunch, try_bunch, gen):
    candidate_bunch = []
    candidate_bunch_selected = []
    loss_value = []
    p = mp.Pool()
    if gen == 0:
        for i in range(100):
            candidate_bunch.append(gen_bunch())
            loss_value.append((i, loss(true_bunch, candidate_bunch[i])))
        loss_value = sorted(loss_value, key=lambda x: x[1])
        for i in range(len(candidate_bunch)):
            candidate_bunch_selected.append(candidate_bunch[loss_value[int(i/10)][0]])
        candidate_bunch_selected = mutate(candidate_bunch_selected)
        return candidate_bunch
    elif gen == 1:
        candidate_bunch = generation(true_bunch, candidate_bunch, 0)
        async_res = [p.apply_async(trans_beam, args=(true_bunch, candidate_bunch[i],)) for i in range(len(candidate_bunch))]
        res_bunch = [ar.get() for ar in async_res]
        for i in range(len(res_bunch)):
            loss_value.append((i, loss(true_bunch, res_bunch[i])))
        loss_value = sorted(loss_value, key=lambda x: x[1])
        for i in range(len(res_bunch)):
            candidate_bunch_selected.append(res_bunch[loss_value[int(i/10)][0]])
        candidate_bunch_selected = mutate(candidate_bunch_selected)
        return candidate_bunch_selected
    else:
        candidate_bunch = try_bunch
        async_res = [p.apply_async(trans_beam, args=(true_bunch, candidate_bunch[i],)) for i in range(len(candidate_bunch))]
        res_bunch = [ar.get() for ar in async_res]
        for i in range(len(res_bunch)):
            loss_value.append((i, loss(true_bunch, res_bunch[i])))
        loss_value = sorted(loss_value, key=lambda x: x[1])
        for i in range(len(res_bunch)):
            candidate_bunch_selected.append(res_bunch[loss_value[int(i/10)][0]])
        if gen != MAX_GEN - 1:
            candidate_bunch_selected = mutate(candidate_bunch_selected)
        return candidate_bunch_selected

if __name__ == "__main__":
    true_bunch = gen_bunch()
    try_bunch = []
    loss_value = []
    for i in range(1, MAX_GEN):
        print(":: Generation #", i)
        try_bunch = generation(true_bunch, try_bunch, i)
    for i in range(len(try_bunch)):
        try_bunch[i] = trans_beam(true_bunch, try_bunch[i])
        loss_value.append((i, loss(true_bunch, try_bunch[i])))
    loss_value = sorted(loss_value, key=lambda x: x[1])
    plot_bunch(true_bunch, try_bunch[loss_value[0][0]])
    print("::")
    print(":: Position X: ", np.mean(true_bunch[0]), "[mm]")
    print(":: Position Y: ", np.mean(true_bunch[1]), "[mm]")
    print(":: Size X: ", np.std(true_bunch[0]), "[mm^2]")
    print(":: Size Y: ", np.std(true_bunch[1]), "[mm^2]")
    print("::")
    print(":: Position X: ", np.mean(try_bunch[loss_value[0][0]][0]), "[mm]")
    print(":: Position Y: ", np.mean(try_bunch[loss_value[0][0]][1]), "[mm]")
    print(":: Size X: ", np.std(try_bunch[loss_value[0][0]][0]), "[mm^2]")
    print(":: Size Y: ", np.std(try_bunch[loss_value[0][0]][1]), "[mm^2]")
    print("::")
    print(":: Done")

"""
def generation(true_bunch, try_bunch, gen):
    candidate_bunch = []
    candidate_bunch_selected = []
    loss_value = []
    p = mp.Pool()
    if gen == 0:
        for i in range(100):
            candidate_bunch.append(gen_bunch())
        return candidate_bunch
    elif gen == 1:
        async_res = [p.apply_async(trans_method, args=(x, i,)) for i, x in enumerate(bunch_list)]
        res_bunch = [ar.get() for ar in async_res]
        candidate_bunch = generation(true_bunch, candidate_bunch, 0)
        for i in range(len(candidate_bunch)):
            candidate_bunch[i] = trans_beam(true_bunch, candidate_bunch[i])
            loss_value.append((i, loss(true_bunch, candidate_bunch[i])))
        loss_value = sorted(loss_value, key=lambda x: x[1])
        for i in range(len(candidate_bunch)):
            candidate_bunch_selected.append(candidate_bunch[loss_value[int(i/10)][0]])
        candidate_bunch_selected = mutate(candidate_bunch_selected)
        return candidate_bunch_selected
    else:
        candidate_bunch = try_bunch
        for i in range(len(candidate_bunch)):
            candidate_bunch[i] = trans_beam(true_bunch, candidate_bunch[i])
            loss_value.append((i, loss(true_bunch, candidate_bunch[i])))
        loss_value = sorted(loss_value, key=lambda x: x[1])
        for i in range(len(candidate_bunch)):
            candidate_bunch_selected.append(candidate_bunch[loss_value[int(i/10)][0]])
        if gen != MAX_GEN - 1:
            candidate_bunch_selected = mutate(candidate_bunch_selected)
        return candidate_bunch_selected
"""