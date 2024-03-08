from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
#import matplotlib
#matplotlib.use('TkAgg')
import matplotlib as mpl
#import multiprocessing as mp
from scipy import signal

mpl.rcParams["font.size"]=14

plt.ion()

ncfile = "/home1/datawork/groullet/data/Nyles/tgv_256/tgv_256_{rank:02}_hist.nc"


nprocs = 64

with Dataset(ncfile.format(rank=0)) as nc:
    nx = len(nc.dimensions["x_u"])
    nt = len(nc.dimensions["t"])
    time =  nc.variables["t"][:]


npx = npy = 8
nx = ny =  nz = 256
nx2 = nx//npx
ny2 = ny//npy

kt = 15
j = ny2//2


def get_velcomp(kt, j0=0, comp="u"):
    assert comp in ["u","v","w"]
    u = np.zeros((nz,nx))
    for rank in range(j0,j0+npy):
        i0 = rank % npx
        idx = slice(nx2*i0, nx2*(i0+1))
        with Dataset(ncfile.format(rank=rank)) as nc:
            u[:,idx] = nc.variables[comp][kt][:,j]
    return u

def get_spec(kt):
    for comp in "uvw":
        u=get_velcomp(60,comp=comp)
        kx,pxxs=signal.periodogram(u)
        if comp=="u":
            pxx=pxxs.mean(axis=0)
        else:
            pxx+= pxxs.mean(axis=0)
    pxx*=0.5
    return kx, pxx

kt = 40
kx, pxx = get_spec(kt)

fig, ax = plt.subplots()
ax.loglog(kx,pxx,lw=2)
ax.loglog(kx,10e-4*(kx/1e-2)**(-5/3),"k",label=r"$k^{-5/3}$")
ax.set_xlabel("k")
ax.set_ylabel("KE spectrum")
ax.set_ylim([1e-8,1e-2])
ax.legend()
ax.set_title(f"t={time[kt]:.1f}")
plt.tight_layout()

