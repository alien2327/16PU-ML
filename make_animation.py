import struct, os, shutil, sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import multiprocessing
import numpy as np
import cv2
import datetime
from tqdm import tqdm
import requests

now = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

class Grid(object):
    def __init__(self, n, addr="#13"):
        self.n = n
        self.x = 80
        self.y = 80
        self.dx = self.x / self.n
        self.dy = self.y / self.n
        self.grid_coor = np.zeros((self.n, self.n, 2))
        self.grid = np.zeros((self.n, self.n))
        self.alloc_coor()

    def alloc_coor(self):
        for x in range(self.n):
            for y in range(self.n):
                self.grid_coor[y, x] = ((-(self.x / 2) + self.dx * x)+self.dx/2, ((self.y / 2) - self.dy * y)-self.dy/2)

    def calHist(self, bunch):
        self.grid = np.zeros((self.n, self.n))
        for i in range(len(bunch)):
            a, b = bunch[i][0], bunch[i][1]
            try:
                a_p, b_p = int((a + (self.x - self.dx)/2)//self.dx), self.n - 1 - int((b + (self.y - self.dy)/2)//self.dy)
                self.grid[b_p, a_p] += (abs(self.grid_coor[b_p - 1, a_p + 1][0]-a)*abs(self.grid_coor[b_p - 1, a_p + 1][1]-b)) / (self.dx*self.dy)
                self.grid[b_p - 1, a_p] += (abs(self.grid_coor[b_p, a_p + 1][0]-a)*abs(self.grid_coor[b_p, a_p + 1][1]-b)) / (self.dx*self.dy)
                self.grid[b_p, a_p + 1] += (abs(self.grid_coor[b_p - 1, a_p][0]-a)*abs(self.grid_coor[b_p - 1, a_p][1]-b)) / (self.dx*self.dy)
                self.grid[b_p - 1, a_p + 1] += (abs(self.grid_coor[b_p, a_p][0]-a)*abs(self.grid_coor[b_p, a_p][1]-b)) / (self.dx*self.dy)
            except:
                pass

def worker(g, gen, parent_dir):
    offset = np.full((10000, 2), [0.7799954861648153, -0.17114181002093012])
    chr_fname = f"/individual_{gen}.dat"
    print(gen, end=", ")
    f = open(parent_dir+chr_fname, "rb")
    dlen = 10000
    fig, axes = plt.subplots(6, 6, sharex='all', sharey='all', figsize=(40, 40))
    for i in range(6):
        for j in range(6):
            if i == 0 and j == 0:
                data_raw = struct.unpack("d" * 2*dlen, f.read(8*2*dlen))
                data_raw = list(data_raw)
                x, y = data_raw[:dlen], data_raw[dlen:]
                pos = list(zip(x, y))
                pos = np.array(pos)#*offset
                g.calHist(pos)
                x, y, z = g.grid_coor[:, :, 0], g.grid_coor[:, :, 1], g.grid
                axes[i][j].pcolormesh(x, y, z, shading="gouraud", cmap=cm.inferno)
                axes[i][j].pcolormesh(x, y, z, shading="gouraud", cmap=cm.inferno)
                axes[i][j].spines["top"].set_color("blue")
                axes[i][j].spines["left"].set_color("blue")
                axes[i][j].spines["right"].set_color("blue")
                axes[i][j].spines["bottom"].set_color("blue")
                axes[i][j].spines["top"].set_linewidth(3)
                axes[i][j].spines["left"].set_linewidth(3)
                axes[i][j].spines["right"].set_linewidth(3)
                axes[i][j].spines["bottom"].set_linewidth(3)
            else:
                data_raw = struct.unpack("d" * 2*dlen, f.read(8*2*dlen))
                data_raw = list(data_raw)
                x, y = data_raw[:dlen], data_raw[dlen:]
                pos = list(zip(x, y))
                pos = np.array(pos)#*offset
                g.calHist(pos)
                x, y, z = g.grid_coor[:, :, 0], g.grid_coor[:, :, 1], g.grid
                axes[i][j].pcolormesh(x, y, z, shading="gouraud", cmap=cm.inferno)
    f.close()
    os.remove(parent_dir+chr_fname)
    plt.savefig(f"./png_temp/individual_{gen}.png", bbox_inches='tight')
    plt.close()

def MSE(origin, compare, grid_num):
    num = grid_num**2
    origin = origin.reshape(num)
    compare = compare.reshape(num)
    return np.sum((origin-compare)**2)/num

def plot_best(g, gen):
    global now
    fig = plt.figure(figsize=(10, 10))
    parent_dir = "result"
    dlen = 10000
    chr_fname = f"/individual_{gen}.dat"
    fb = open(parent_dir+chr_fname, "rb")
    data_raw_b = struct.unpack("d"*2*dlen, fb.read(8*2*dlen))
    x, y = data_raw_b[:dlen], data_raw_b[dlen:]

    mx = np.mean(x)
    my = np.mean(y)
    sx = np.sqrt(np.mean(abs(np.array(x) - mx)**2))
    sy = np.sqrt(np.mean(abs(np.array(y) - my)**2))

    pos = list(zip(x, y))
    pos = np.array(pos)
    g.calHist(pos)
    x_b, y_b, z_b = g.grid_coor[:, :, 0], g.grid_coor[:, :, 1], g.grid

    #err = MSE(ref, z_b, 100)

    axes = fig.add_subplot(1,1,1)
    axes.pcolormesh(x_b, y_b, z_b, shading="gouraud", cmap=cm.inferno)
    axes.set_xlabel("x [mm]")
    axes.set_ylabel("y [mm]")
    axes.grid(color="gray")

    axes.text(15, 22, f"$\langle x\\rangle$: {mx:1.5e} mm\n$\langle y\\rangle$: {my:1.5e} mm\n$\sigma_x$: {sx:1.5e}\n$\sigma_y$: {sy:1.5e}", bbox={'facecolor': 'white', 'pad': 5})
    hist_x = np.zeros(200)
    hist_y = np.zeros(200)

    z_b_t = g.grid.T

    for i in range(200):
        hist_x[i] = sum(z_b_t[i])
        hist_y[i] = sum(z_b[i])

    driver = make_axes_locatable(axes)
    axHistx = driver.append_axes("top", 1.2, pad=0.1, sharex=axes)
    axHistx.set_title("Best Fit")
    axHisty = driver.append_axes("right", 1.2, pad=0.1, sharey=axes)

    axHistx.xaxis.set_tick_params(labelbottom=False)
    axHisty.yaxis.set_tick_params(labelleft=False)

    line = np.linspace(-38, 38, 200)

    axHistx.plot(line, hist_x)
    axHisty.plot(hist_y, line)#, orientation='horizontal')

    plt.savefig(f"./{now}_best_scaled.png", bbox_inches='tight')
    plt.close()

def plot_ref(g):
    global now
    dlen = 10000
    ref_fname = "ref.dat"
    fb = open(ref_fname, "rb")

    data_raw_b = struct.unpack("d"*2*dlen, fb.read(8*2*dlen))

    fig = plt.figure(figsize=(10, 10))
    x, y = data_raw_b[:dlen], data_raw_b[dlen:]

    mx = np.mean(x)
    my = np.mean(y)
    sx = np.sqrt(np.mean(abs(np.array(x) - mx)**2))
    sy = np.sqrt(np.mean(abs(np.array(y) - my)**2))

    pos = list(zip(x, y))
    pos = np.array(pos)
    g.calHist(pos)
    x_b, y_b, z_b = g.grid_coor[:, :, 0], g.grid_coor[:, :, 1], g.grid

    axes = fig.add_subplot(1,1,1)
    axes.pcolormesh(x_b, y_b, z_b, shading="gouraud", cmap=cm.inferno)
    axes.set_xlabel("x [mm]")
    axes.set_ylabel("y [mm]")
    axes.grid(color="gray")

    axes.text(13, 20, f"$\langle x\\rangle$: {mx:1.5e} mm\n$\langle y\\rangle$: {my:1.5e} mm\n$\sigma_x$: {sx:1.5e}\n$\sigma_y$: {sy:1.5e}", bbox={'facecolor': 'white', 'pad': 5})
    hist_x = np.zeros(200)
    hist_y = np.zeros(200)

    z_b_t = g.grid.T

    for i in range(200):
        hist_x[i] = max(z_b_t[i])
        hist_y[i] = max(z_b[i])

    driver = make_axes_locatable(axes)
    axHistx = driver.append_axes("top", 1.2, pad=0.1, sharex=axes)
    axHistx.set_title("Reference")
    axHisty = driver.append_axes("right", 1.2, pad=0.1, sharey=axes)

    axHistx.xaxis.set_tick_params(labelbottom=False)
    axHisty.yaxis.set_tick_params(labelleft=False)

    line = np.linspace(-30, 30, 200)

    axHistx.plot(line, hist_x)
    axHisty.plot(hist_y, line)#, orientation='horizontal')

    plt.savefig(f"./{now}_ref_scaled.png", bbox_inches='tight')
    plt.close()

    return z_b

def send_line():
    global now
    url = "https://notify-api.line.me/api/notify"
    acc_tok = "bh2SkHk5HJgAsrrrzttUz6OQiEr4Ao9kdeZfSyDPFiH"
    headers = {'Authorization': 'Bearer ' + acc_tok}
    parent_dir = "result"
    mes = f"{parent_dir} {now} data is processed"
    payload = {'message': mes}
    files_1 = {'imageFile': open(f"./{now}_best_scaled.png", 'rb')}
    return requests.post(url, headers=headers, params=payload, files=files_1,)

if __name__ == "__main__":

    parent_dir = "result"
    g = Grid(200)

    turn = sys.argv[1]
    #ref = plot_ref(g)
    plot_best(g, turn)
    #send_line()

    shutil.copyfile(parent_dir+"/data_result.csv", f"./{now}_data_result.csv")
    shutil.copyfile(parent_dir+f"/individual_{turn}.dat", f"./{now}_individual_{turn}.dat")

    print("\ndone")
"""
    print("\nGenerating png files")
    for i in range(50):
        processes = []
        print("Processing chr: ", end="")
        for j in range(10):
            processes.append(multiprocessing.Process(target=worker, args=(Grid(100), 50*i+j, parent_dir, )))
        for p in processes:
            p.start()
        for p in processes:
            p.join()
        print()

    print("Make animated mp4 file")
    fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
    video = cv2.VideoWriter(f'{now}_result.mp4',fourcc, 5.0, (1500, 1500))
    for i in tqdm(range(turn)):
        file = f"./png_temp/individual_{i}.png"
        img = cv2.imread(file)
        img = cv2.resize(img, (1500, 1500))
        video.write(img)

    video.release()
"""