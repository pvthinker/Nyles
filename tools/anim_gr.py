import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import movietools

plt.ion()


fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
movie = movietools.Movie(fig, "movie")

ncfile = "out.nc"
j = 8
maxi = 1.2

template = "time = %4.2f"

with nc.Dataset(ncfile) as fid:
    nt = len(fid.dimensions["t"])

    for kt in range(nt):
        print("\rframe=%3i/%i" % (kt, nt-1), end="")
        b = fid["b"][kt, :, j, :]
        t = fid["t"][kt]
        t_str = template % t
        if kt == 0:
            im = plt.imshow(b, origin="xy", cmap=plt.get_cmap('RdBu_r'), vmin=-maxi, vmax=maxi)
            ti = plt.title(t_str)
            plt.colorbar(im)
        else:
            im.set_data(b)
            ti.set_text(t_str)
        fig.canvas.draw()
        movie.addframe()
    print()

movie.finalize()
