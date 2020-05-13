from matplotlib import pyplot as plt
from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('TkAgg')


plt.ion()


ncfile = "/home/roullet/data/Nyles/tank_toy_1/tank_toy_1_00_hist.nc"

with Dataset(ncfile) as nc:
    b = nc.variables["b"][:, :, :, :]
    nt = len(nc.dimensions["t"])

kt = 20
j = 2

maxi = 1e-8

fig = plt.figure(1)
plt.clf()
first = True
for kt in range(nt):
    b2 = b[kt, :, j, :]
    z2d = b2 + np.flip(np.flip(b2, axis=1), axis=0)
    if first:
        first = False
        im = plt.imshow(z2d, origin="xy", cmap=plt.get_cmap("RdBu_r"),
                        vmin=-maxi, vmax=maxi)
        plt.colorbar(im)
    else:
        im.set_data(z2d)
    fig.canvas.draw()
