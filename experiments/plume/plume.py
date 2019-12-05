import numpy as np

import nyles as nyles_module
import parameters


nh = 2
nz = 32
ny = 32
nx = 32
# It must be possible to set the following lengths to arbitrary values,
# but currently, due to a problem in the handling of the metric in the
# calculation of p, it is necessary to ensure dx = dy = dz = 1. #TODO
Lz = 1.0 * nz
Ly = 1.0 * ny
Lx = 1.0 * nx


# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "LES"
param.model["geometry"] = "closed"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "plume"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = "all"
param.IO["timestep_history"] = 0

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 80.0
param.time["auto_dt"] = False
# parameter if auto_dt is False
param.time["dt"] = 0.02
# parameters if auto_dt is True
param.time["cfl"] = 0.4
param.time["dt_max"] = 1.0

param.discretization["global_nx"] = nx
param.discretization["global_ny"] = ny
param.discretization["global_nz"] = nz
param.discretization["orderVF"] = 5  # upwind-order for vortex-force term
param.discretization["orderA"] = 5  # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = 1
param.MPI["npy"] = 1
param.MPI["npz"] = 1


nyles = nyles_module.Nyles(param)

b = nyles.model.state.b.view('i')
x = nyles.grid.x_b.view('i') / Lx
y = nyles.grid.y_b.view('i') / Ly
z = nyles.grid.z_b.view('i') / Lz

b[:] = 10.*z

b[(x > 0.4) & (x < 0.6) & (y > 0.4) & (y < 0.6) & (z < .2)] = 8


# Add noise, uniformly distributed from -1 to +1 (times a factor)
#np.random.seed(1)
#noise = np.random.uniform(size = np.shape(b))
#noise = noise*2-1

#b += noise*.01

nyles.run()
