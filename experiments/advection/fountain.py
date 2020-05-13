import numpy as np
import matplotlib.pyplot as plt

import nyles as nyles_module
import parameters


nh = 3

nxglo = 32
nyglo = 32
nzglo = 32

npx = 2
npy = 2
npz = 1

nx = nxglo//npx
ny = nyglo//npy
nz = nzglo//npz

Lx = 1.
Ly = 1.
Lz = 1.


# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "advection"
param.model["geometry"] = "closed"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "fountain1"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['b', 'u']
param.IO["timestep_history"] = 0.1

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 10.0
param.time["auto_dt"] = False
# parameter if auto_dt is False
param.time["dt"] = .01
# parameters if auto_dt is True
param.time["cfl"] = 0.4
param.time["dt_max"] = 1.0

param.discretization["global_nx"] = nxglo
param.discretization["global_ny"] = nyglo
param.discretization["global_nz"] = nzglo
param.discretization["orderVF"] = 3  # upwind-order for vortex-force term
param.discretization["orderA"] = 3  # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = npx
param.MPI["npy"] = npy
param.MPI["npz"] = npz


nyles = nyles_module.Nyles(param)

b = nyles.model.state.b.view('i')
w = nyles.model.state.u['k'].view('i')
z = nyles.grid.z_b.view('i')/Lz
y = nyles.grid.y_b.view('i')/Ly
x = nyles.grid.x_b.view('i')/Lx

y1d = nyles.grid.y_b_1D/Ly
x1d = nyles.grid.x_b_1D/Lx

def idx2slice(idx):
    if len(idx) == 0:
        s = slice(0,0)
    else:
        s= slice(idx[0], idx[-1]+1)
    return s

i1 = nx//3
i2 = 2*nx//3
j1 = ny//3
j2 = 2*ny//3
ii = np.where((x1d>1./3) & (x1d<2./3))[0]
jj = np.where((y1d>1./3) & (y1d<2./3))[0]
#w[:-1, j1:j2, i1:i2] = 1.
dz = Lz/nzglo
w[:-1, idx2slice(jj), idx2slice(ii)] = 1.*dz
print('and make it divergent-free')
nyles.model.make_u_divergentfree()

b[:] = np.round(x*6) % 2 + np.round(y*6) % 2 + np.round(z*6) % 2


nyles.run()
