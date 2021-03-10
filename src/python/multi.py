from concurrent.futures import ProcessPoolExecutor as Pool
import numpy as np

def gen_bunch():
    maxR = 60
    minR = -60
    maxR2 = 20
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

def trans_method(x, c):
    d = 0.01 * np.random.rand()
    for i in range(len(x[0])):
        if c == 0:
            x[0][i] += d
        elif c == 1:
            x[0][i] -= d
        elif c == 2:
            x[0][i] *= (1 + d)
        elif c == 3:
            x[0][i] *= (1 - d)
        else:
            x[0][i] = x[0][i]
    d = 0.01 * np.random.rand()
    for i in range(len(x[1])):
        if c == 0:
            x[1][i] += d
        elif c == 1:
            x[1][i] -= d
        elif c == 2:
            x[1][i] *= (1 + d)
        elif c == 3:
            x[1][i] *= (1 - d)
        else:
            x[1][i] = x[1][i]
    return x

def trans_beam(true_bunch, bunch):
    loss_value = []
    x, y = bunch
    for k in range(len(x)):
        x[k] = trans_method(x[k], i)
    try_bunch = np.concatenate((x, y))
    try_bunch = np.reshape(try_bunch, (2, 1000))
    loss_value.append((loss(true_bunch, try_bunch)))
    return bunch

if __name__ == '__main__':
    true = gen_bunch()
    x_list = []
    procs = []
    x = gen_bunch()
    for i in range(5):
        x_list.append(x)
    print(x_list[0][0][0])
    p = Pool()
    async_res = [p.apply_async(trans_method, args=(x, i,)) for i, x in enumerate(x_list)]
    res = [ar.get() for ar in async_res]
    print(res[0][0][0])