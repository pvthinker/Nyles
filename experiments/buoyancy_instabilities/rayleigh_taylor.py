import numpy as np

import nyles as nyles_module
import parameters


nh = 3

nz = 64
ny = 32
nx = 32

npx = 1
npy = 1
npz = 1

Lz = 1.0
Ly = 1.0
Lx = 1.0


# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "LES"
param.model["geometry"] = "closed"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "test"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ["b", "u", "vor"]
param.IO["timestep_history"] = 0.1

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 2.0
param.time["auto_dt"] = True
# parameter if auto_dt is False
param.time["dt"] = 0.01
# parameters if auto_dt is True
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.01

param.discretization["global_nx"] = nx*npx
param.discretization["global_ny"] = ny*npy
param.discretization["global_nz"] = nz*npz
param.discretization["orderVF"] = 5  # upwind-order for vortex-force term
param.discretization["orderA"] = 5  # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = npx
param.MPI["npy"] = npy
param.MPI["npz"] = npz


nyles = nyles_module.Nyles(param)

b = nyles.model.state.b.view('i')
x = nyles.grid.x_b.view('i') / Lx
y = nyles.grid.y_b.view('i') / Ly
z = nyles.grid.z_b.view('i') / Lz

b[:] = 15 - 5*np.tanh((np.cos(np.pi*x)*0.+z-0.5)/.05)


# Another initial buoyancy profile (as in Fluid2D):
# def sigmoid(x, delta):
#     return 1 / (1 + np.exp(-(x-0.5)/delta))
# def stratif():
#     sigma = nyles.grid.dx/2  # width of the interface
#     return sigmoid(z/nyles.grid.Lz, sigma/nyles.grid.Lz)
# b[:] = (1 - stratif() - 0.5)

# Add noise, uniformly distributed from -1 to +1 (times a factor)
np.random.seed(1)
noise = np.random.uniform(size = np.shape(b))
noise = noise*2-1

b += noise*.01

nyles.run()
