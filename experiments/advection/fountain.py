import numpy as np
import matplotlib.pyplot as plt

import nyles as nyles_module
import parameters


nh = 2
nz = 64
ny = 32
nx = 32
Lz = nz
Ly = ny
Lx = nx


# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "advection"
param.model["geometry"] = "closed"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "fountain0"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['u', 'b']
param.IO["timestep_history"] = 0.2

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 120.0
param.time["auto_dt"] = False
# parameter if auto_dt is False
param.time["dt"] = .5
# parameters if auto_dt is True
param.time["cfl"] = 0.4
param.time["dt_max"] = 1.0

param.discretization["global_nx"] = nx
param.discretization["global_ny"] = ny
param.discretization["global_nz"] = nz
param.discretization["orderVF"] = 5  # upwind-order for vortex-force term
param.discretization["orderA"] = 3  # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = 1
param.MPI["npy"] = 1
param.MPI["npz"] = 1


nyles = nyles_module.Nyles(param)

b = nyles.model.state.b.view('i')
w = nyles.model.state.u['k'].view('i')
z = nyles.grid.z_b.view('i')/Lz
y = nyles.grid.y_b.view('i')/Ly
x = nyles.grid.x_b.view('i')/Lx

i1 = nx//3
i2 = 2*nx//3
j1 = ny//3
j2 = 2*ny//3
w[:-1, j1:j2, i1:i2] = 1.
print('and make it divergent-free')
nyles.model.make_u_divergentfree()

b[:] = np.round(x*6) % 2 + np.round(y*6) % 2 + np.round(z*6) % 2


nyles.run()
