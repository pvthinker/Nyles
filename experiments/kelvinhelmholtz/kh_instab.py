import numpy as np

import nyles as nyles_module
import parameters


nh = 2
nz = 32
ny = 32
nx = 64
# It must be possible to set the following lengths to arbitrary values,
# but currently, due to a problem in the handling of the metric in the
# calculation of p, it is necessary to ensure dx = dy = dz = 1. #TODO
Lz = 1.0
Ly = 1.0
Lx = 2.0


# Get the default parameters, then modify them as needed
param = parameters.UserParameters()

param.model["modelname"] = "LES"
param.model["geometry"] = "perio_x"
param.model["Lx"] = Lx
param.model["Ly"] = Ly
param.model["Lz"] = Lz

param.IO["datadir"] = "~/data/Nyles"
param.IO["expname"] = "khi_0"
param.IO["mode"] = "overwrite"
param.IO["variables_in_history"] = ['b', 'u']
param.IO["timestep_history"] = .2  # 0.0 saves every frame
param.IO["disk_space_warning"] = 0.5  # in GB

param.time["timestepping"] = "LFAM3"
param.time["tend"] = 40.0
param.time["auto_dt"] = True
# parameter if auto_dt is False
param.time["dt"] = 0.2
# parameters if auto_dt is True
param.time["cfl"] = 0.8
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
grid = nyles.grid

b = nyles.model.state.b.view('i')
u = nyles.model.state.u['i'].view('i')
x = nyles.grid.x_b.view('i') / Lx
y = nyles.grid.y_b.view('i') / Ly
z = nyles.grid.z_b.view('i') / Lz

perturb = np.sin(x*8*np.pi)

sigma = .1
Deltab = .1
DeltaU = 1.
Ri = Deltab/(sigma*Lz) / (DeltaU/(sigma*Lz))**2
print("Richardson number: %.2g" % Ri)

b[...] = Deltab*0.5*np.tanh((1e-2*perturb+(z-0.5))/sigma)
u[...] = DeltaU*0.5*np.tanh((1e-2*perturb+(z-0.5))/sigma)
u *= grid.dx  # to make it covariant!!!

nyles.model.diagnose_var(nyles.model.state)

nyles.run()
