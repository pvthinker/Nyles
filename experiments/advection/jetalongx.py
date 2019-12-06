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
param.model["geometry"] = "perio_x"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "jetx"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['u', 'b']
param.IO["timestep_history"] = 0.2

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 12.0
param.time["auto_dt"] = False
# parameter if auto_dt is False
param.time["dt"] = .05
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

state = nyles.model.state
grid = nyles.grid

b = state.b.view('i')
u = state.u['i'].view('i')
z = grid.z_b.view('i')/Lz
y = grid.y_b.view('i')/Ly
x = grid.x_b.view('i')/Lx

i1 = nx//3
i2 = 2*nx//3
j1 = ny//3
j2 = 2*ny//3

u[...] = z-0.5

b[:] = np.round(x*6) % 2 + np.round(y*6) % 2 + np.round(z*6) % 2
b *= 0.1

nyles.model.diagnose_var(state)
nyles.run()
