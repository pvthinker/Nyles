from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
import multiprocessing as mp

#import matplotlib
# matplotlib.use('TkAgg')


plt.ion()


ncfile = "/home/roullet/data/Nyles/tgv/tgv_{rank:02}_hist.nc"
# ncfile = "/home/roullet/data/Nyles/tgv_new/tgv_new_{rank:02}_hist.nc"
# ncfile = "/home/roullet/data/Nyles/tgv_old/tgv_old_{rank:02}_hist.nc"
#ncfile = "/home/roullet/data/Nyles/tgv_2x2/tgv_2x2_{rank:02}_hist.nc"
ncfile = "/home1/datawork/groullet/data/Nyles/tank_256/tank_256_{rank:02}_hist.nc"

nprocs = 64

with Dataset(ncfile.format(rank=0)) as nc:
    nx = len(nc.dimensions["x_u"])
    nt = len(nc.dimensions["t"])
    time = nc.variables["t"][:]

dx = 2*np.pi/nx

print(f"#snapshots: {nt}")


def compute_ke(kt):
    ke = 0
    for rank in range(nprocs):
        with Dataset(ncfile.format(rank=rank)) as nc:
            u = nc.variables["u"][kt, :, :, :]
            v = nc.variables["v"][kt, :, :, :]
            w = nc.variables["w"][kt, :, :, :]
        ke += (0.5/dx**2)*(np.mean(u**2)+np.mean(v**2)+np.mean(w**2))
    return ke


def compute_ke(nt):
    ke = np.zeros((nt,))
    for rank in range(nprocs):
        print(f"{rank}/{nprocs}")
        with Dataset(ncfile.format(rank=rank)) as nc:
            for kt in range(nt):
                u = nc.variables["u"][kt, :, :, :]
                v = nc.variables["v"][kt, :, :, :]
                w = nc.variables["w"][kt, :, :, :]
                ke[kt] += (0.5/dx**2)*(np.mean(u**2) +
                                       np.mean(v**2)+np.mean(w**2))
    return ke


def compute_kerank(rank):
    ke = np.zeros((nt,))
    print(f"{rank}/{nprocs}")
    with Dataset(ncfile.format(rank=rank)) as nc:
        for kt in range(nt):
            u = nc.variables["u"][kt, :, :, :]
            v = nc.variables["v"][kt, :, :, :]
            w = nc.variables["w"][kt, :, :, :]
            ke[kt] = (0.5/dx**2)*(np.mean(u**2)+np.mean(v**2)+np.mean(w**2))
    return ke


def compute_enstrophy(kt):
    with Dataset(ncfile) as nc:
        u = nc.variables["vor_i"][kt, :, :, :]
        v = nc.variables["vor_j"][kt, :, :, :]
        w = nc.variables["vor_k"][kt, :, :, :]
    ens = (0.5/dx**4)*(np.mean(u**2)+np.mean(v**2)+np.mean(w**2))
    return ens


#ke = compute_ke(nt)
with mp.Pool(processes=16) as pool:
    kes = pool.starmap(compute_kerank, range(nprocs))
ke = np.sum(np.asarray(kes), axis=0)

#ke = np.asarray([compute_ke(kt) for kt in range(nt)])
#ens = np.asarray([compute_enstrophy(kt) for kt in range(nt)])

dKdt = -np.diff(ke)/np.diff(time)

fig, ax = plt.subplots(figsize=(6.55, 4.92))
# plt.clf()
ax.plot(time[1:], dKdt)
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$-dK/dt$")
ax.set_xlim([0, 20])
ax.set_ylim([0, 0.018])
plt.grid()
