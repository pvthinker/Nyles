import numpy as np

import nyles as nyles_module
import parameters


nh = 3

factor = 1

nzglo = 32*factor
nyglo = 32*factor
nxglo = 32*factor

npx = 2
npy = 2
npz = 1

nx = nxglo // npx
ny = nyglo // npy
nz = nzglo // npz

Lz = 2*np.pi
Ly = 2*np.pi
Lx = 2*np.pi


# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "Euler3d"
param.model["geometry"] = "perio_xyz"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

#param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "tgv"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ["u", "vor"]
param.IO["timestep_history"] = 0.25

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 10
param.time["auto_dt"] = True
# parameter if auto_dt is False
param.time["dt"] = 0.01/factor
# parameters if auto_dt is True
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.4/factor

param.discretization["global_nx"] = nx*npx
param.discretization["global_ny"] = ny*npy
param.discretization["global_nz"] = nz*npz
param.discretization["orderVF"] = 5  # upwind-order for vortex-force term
param.discretization["orderA"] = 5  # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = npx
param.MPI["npy"] = npy
param.MPI["npz"] = npz

param.multigrid["nglue"] = 16

nyles = nyles_module.Nyles(param)

#b = nyles.model.state.b.view('i')
u = nyles.model.state.u['i'].view('i')
v = nyles.model.state.u['j'].view('i')
x = nyles.grid.x_b.view('i')
y = nyles.grid.y_b.view('i')
z = nyles.grid.z_b.view('i')

#b[:] = 0.

a=1.2
b=1.8
c=0.5
u[:] = + np.sin(x+a) * np.cos(y+b) * np.cos(z+c)
v[:] = - np.cos(x+a) * np.sin(y+b) * np.cos(z+c)

dx = nyles.grid.dx

u *= dx
v *= dx

# Another initial buoyancy profile (as in Fluid2D):
# def sigmoid(x, delta):
#     return 1 / (1 + np.exp(-(x-0.5)/delta))
# def stratif():
#     sigma = nyles.grid.dx/2  # width of the interface
#     return sigmoid(z/nyles.grid.Lz, sigma/nyles.grid.Lz)
# b[:] = (1 - stratif() - 0.5)

# Add noise, uniformly distributed from -1 to +1 (times a factor)
# np.random.seed(1)
# noise = np.random.uniform(size = np.shape(b))
# noise = noise*2-1

#b += noise*1e-16



nyles.run()
