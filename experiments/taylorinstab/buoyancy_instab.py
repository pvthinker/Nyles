import numpy as np
import matplotlib.pyplot as plt

import nyles as nyles_module
import parameters

nh = 2
nz = 32
ny = 32
nx = 48
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

param.IO["expname"] = "buoyancy instab"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = "all"
param.IO["timestep_history"] = 0.0

param.time["timestepping"] = "EF"
param.time["tend"] = 1.0
param.time["auto_dt"] = False
param.time["dt"] = 0.1

param.discretization["global_nx"] = nx
param.discretization["global_ny"] = ny
param.discretization["global_nz"] = nz
param.discretization["orderVF"] = 1  # upwind-order for vortex-force term
param.discretization["orderA"] = 1  # upwind-order for advection term

param.MPI["nh"] = nh


nyles = nyles_module.Nyles(param)

b = nyles.model.state.b.view('i')
z = nyles.grid.z_b.view('i')
x = nyles.grid.x_b.view('i')

b[:] = 15+ 5*np.tanh((np.cos(np.pi*x)*0.1+z-0.5)/.02)

noise = np.random.uniform(size = np.shape(b))
noise = noise*2-1

b += noise*.1

nyles.run()

plt.figure()
plt.pcolor(b[:,10,:])
plt.colorbar()

plt.show()
