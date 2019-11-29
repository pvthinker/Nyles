import numpy as np
import matplotlib.pyplot as plt

import nyles as nyles_module
import parameters


nh = 2
nz = 32
ny = 32
nx = 64
Lz = nz
Ly = ny
Lx = nx


# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "LES"
param.model["geometry"] = "closed"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "wave0"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = "all"
param.IO["timestep_history"] = 0.0

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 20.0
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
z = nyles.grid.z_b.view('i')
x = nyles.grid.x_b.view('i')

b[:] = 15 + 5*np.tanh((np.cos(np.pi*x/Lx)*0.2+z/Lz-0.5)/.05)

# Another initial buoyancy profile (as in Fluid2D):
# def sigmoid(x, delta):
#     return 1 / (1 + np.exp(-(x-0.5)/delta))
# def stratif():
#     sigma = nyles.grid.dx/2  # width of the interface
#     return sigmoid(z/nyles.grid.Lz, sigma/nyles.grid.Lz)
# b[:] = (1 - stratif() - 0.5)

noise = np.random.uniform(size = np.shape(b))
noise = noise*2-1

#b += noise*.1

nyles.run()

plt.figure()
plt.pcolor(b[:,10,:])
plt.colorbar()

plt.show()
