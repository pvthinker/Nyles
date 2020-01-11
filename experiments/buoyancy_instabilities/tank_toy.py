import numpy as np

import nyles as nyles_module
import parameters


nh = 3
nxglo = 64
nyglo = 8
nzglo = 64

npx = 2
npy = 1
npz = 2

nx = nxglo//npx
ny = nyglo//npy
nz = nzglo//npz

Lx = 8.
Ly = 1.
Lz = 8.


# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "LES"
param.model["geometry"] = "closed"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "tank_toy_1"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['b', 'u', 'vor', 'div']
param.IO["timestep_history"] = 1.  # 0.0 saves every frame
param.IO["disk_space_warning"] = 0.5  # in GB

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 10.0
param.time["auto_dt"] = True
# parameter if auto_dt is False
param.time["dt"] = 0.2
# parameters if auto_dt is True
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.2

param.discretization["global_nx"] = nxglo
param.discretization["global_ny"] = nyglo
param.discretization["global_nz"] = nzglo
param.discretization["orderVF"] = 5  # upwind-order for vortex-force term
param.discretization["orderA"] = 5  # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = npx
param.MPI["npy"] = npy
param.MPI["npz"] = npz

#param.multigrid["nglue"] = 1

# param.physics["diff_coef"] = {"u": 1e-1}

nyles = nyles_module.Nyles(param)

b = nyles.model.state.b.view('i')
u = nyles.model.state.u['i'].view('i')
x = nyles.grid.x_b.view('i') / Lx
y = nyles.grid.y_b.view('i') / Ly
z = nyles.grid.z_b.view('i') / Lz

b[:] = np.tanh(((x-0.5)+(z-0.5))/.05)


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
