from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
#import matplotlib
# matplotlib.use('TkAgg')
import matplotlib as mpl
import multiprocessing as mp

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


nbins = 110
e = 1e-3
bins = np.linspace(-1.1-e, 1.1+e, nbins+1)


rank = 0


def get_pdf_rank(rank, bins):
    nbins = len(bins)-1
    print(rank)
    with Dataset(ncfile.format(rank=rank)) as nc:
        nt = len(nc.dimensions["t"])
        pdf = np.zeros((nt, nbins))
        for kt in range(nt):
            b = nc.variables["b"][kt]
            p, _ = np.histogram(b.ravel(), bins, normed=True)
            pdf[kt] = p

    return pdf

# pdfs = [get_pdf_rank(rank, bins)
#         for rank in range(nprocs)]


with mp.Pool(processes=16) as pool:
    args = [(rank, bins) for rank in range(nprocs)]
    pdfs = pool.starmap(get_pdf_rank, args)

pdf = sum(pdfs)/nprocs

fig, ax = plt.subplots()
b = 0.5*(bins[1:]+bins[:-1])
kts = [5, 10, 15, 20, 25, 30]
for k in range(6):  # range(10,nt,10):
    kt = 2*kts[k]
    ax.semilogy(b, pdf[kt], label=f"t={time[kt]:.1f}", lw=2)

ax.legend()
ax.set_xlabel("b")
ax.set_ylabel("count per bin")
ax.set_ylim([1e-2, 1e2])
plt.tight_layout()
