from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
#import matplotlib
# matplotlib.use('TkAgg')
import matplotlib as mpl

mpl.rcParams["font.size"] = 14

plt.ion()

ncfile = "/home1/datawork/groullet/data/Nyles/tgv/tgv_{rank:02}_hist.nc"
ncfile = "/home1/datawork/groullet/data/Nyles/tank_128/tank_128_{rank:02}_hist.nc"
ncfile = "/home1/datawork/groullet/data/Nyles/tank_256/tank_256_{rank:02}_hist.nc"
#ncfile = "/home1/datawork/groullet/data/Nyles/tank_128_bis/tank_128_bis_{rank:02}_hist.nc"

#ncfile = "/home1/datawork/groullet/data/Nyles/tgv_256/tgv_256_{rank:02}_hist.nc"
#ncfile = "/home/roullet/data/Nyles/tgv_new/tgv_new_{rank:02}_hist.nc"
#ncfile = "/home/roullet/data/Nyles/tgv_old/tgv_old_{rank:02}_hist.nc"
#ncfile = "/home1/datawork/groullet/data/Nyles/tgv/tgv_{rank:02}_hist.nc"

nprocs = 64

with Dataset(ncfile.format(rank=0)) as nc:
    nx = len(nc.dimensions["x_u"])
    nt = len(nc.dimensions["t"])
    time = nc.variables["t"][:]


npx = npy = 8
nx = ny = nz = 256
nx2 = nx//npx
ny2 = ny//npy


kt = 15
j = ny2//2

j0 = 0


def get_b(kt, j0):
    b = np.zeros((nz, nx))
    for rank in range(j0, j0+npy):
        i0 = rank % npx
        idx = slice(nx2*i0, nx2*(i0+1))
        with Dataset(ncfile.format(rank=rank)) as nc:
            b[:, idx] = nc.variables["b"][kt][:, j]
    return b


vmin = -1.
vmax = 1.
j0 = 8


def plot_sequence_snapshots():
    fig, axs = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(10.5, 7))
    kts = [0, 5, 10, 15, 20, 30]
    for j, i in np.ndindex((2, 3)):
        k = (i+j*3)
        kt = 2*kts[k]
        b = get_b(kt, j0)

        axs[j, i].imshow(b, origin="lower", extent=[0, 1, 0, 1],
                         vmin=vmin, vmax=vmax, cmap="RdBu_r")
        axs[j, i].text(0.1, 0.1, f"t={time[kt]:.1f}", color="white")

    axs[0, 0].set_ylabel("z")
    axs[1, 0].set_ylabel("z")
    axs[1, 0].set_xlabel("x")
    axs[1, 1].set_xlabel("x")
    axs[1, 2].set_xlabel("x")

    plt.tight_layout()


def produce_frames_for_movie():
    """
    ffmpeg -f image2 -framerate 12 -i snap_%03d.png -s 800x600 -vcodec libx264   -qscale 200 movie.mp4

    """
    fig, ax = plt.subplots(figsize=(8, 6.), dpi=100)
    first = True
    for kt in range(nt):
        b = get_b(kt, j0)
        title = f"t={time[kt]:.1f}"
        if first:
            im = ax.imshow(b, origin="lower",
                           extent=[0, 1, 0, 1],
                           vmin=vmin, vmax=vmax,
                           cmap="RdBu_r")
            plt.colorbar(im)
            ti = ax.set_title(title)
            ax.set_xlabel("x")
            ax.set_ylabel("z")
            first = False
        else:
            im.set_data(b)
            ti.set_text(title)
        # plt.pause(1e-4)
        plt.savefig(f"snap_{kt:03}.png")
