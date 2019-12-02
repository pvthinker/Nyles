import numpy as np

import nyles as nyles_module
import parameters


nh = 2
nz = 32
ny = 32
nx = 128
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
param.IO["expname"] = "khi_0"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['b']
param.IO["timestep_history"] = 1.  # 0.0 saves every frame
param.IO["disk_space_warning"] = 0.5  # in GB

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 100.0
param.time["auto_dt"] = False
# parameter if auto_dt is False
param.time["dt"] = 0.1
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
u = nyles.model.state.u['i'].view('i')
x = nyles.grid.x_b.view('i') / Lx
y = nyles.grid.y_b.view('i') / Ly
z = nyles.grid.z_b.view('i') / Lz

b[:] = np.tanh(((x-0.5)+(z-0.5))/.02)


# Another initial buoyancy profile (as in Fluid2D):
# def sigmoid(x, delta):
#     return 1 / (1 + np.exp(-(x-0.5)/delta))
# def stratif():
#     sigma = nyles.grid.dx/2  # width of the interface
#     return sigmoid(z/nyles.grid.Lz, sigma/nyles.grid.Lz)
# b[:] = (1 - stratif() - 0.5)

# Add noise, uniformly distributed from -1 to +1 (times a factor)
noise = np.random.uniform(size = np.shape(b)) * 2 - 1
#b += noise * 1e-2

nyles.model.diagnose_var(nyles.model.state)

nyles.run()