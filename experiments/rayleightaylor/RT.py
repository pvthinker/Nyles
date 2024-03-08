import numpy as np

import nyles as nyles_module
import parameters


nh = 3

factor = 1

npx = 1
npy = 1
npz = 1

nz = 64*factor
ny = 32*factor//npy
nx = 32*factor//npx

Lx = 8.0
Ly = 8.0
Lz = 16.0


# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "LES"
param.model["geometry"] = "closed"#"perio_xy"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

#param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "RT"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ["b"]
param.IO["timestep_history"] = 0.5

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 15.0
param.time["auto_dt"] = True
# parameter if auto_dt is False
param.time["dt"] = 0.01/factor
# parameters if auto_dt is True
param.time["cfl"] = 0.8
param.time["dt_max"] = 0.1/factor

param.discretization["global_nx"] = nx*npx
param.discretization["global_ny"] = ny*npy
param.discretization["global_nz"] = nz*npz
param.discretization["orderVF"] = 5  # upwind-order for vortex-force term
param.discretization["orderA"] = 5  # upwind-order for advection term

param.MPI["nh"] = nh
param.MPI["npx"] = npx
param.MPI["npy"] = npy
param.MPI["npz"] = npz

param.multigrid["nglue"] = 1

nyles = nyles_module.Nyles(param)

b = nyles.model.state.b.view('i')
x = nyles.grid.x_b.view('i')
y = nyles.grid.y_b.view('i')
z = nyles.grid.z_b.view('i')

dx = nyles.grid.dx
noise = 1e-2*np.random.normal(size=z.shape)

b[:] = np.tanh((0.5*Lz-z+noise)/(dx))

#nyles.model.halo.fill(b)

nyles.run()
