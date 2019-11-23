import model_les_AL as LES
import topology as topo
import numpy as np
import matplotlib.pyplot as plt
import grid

import nyles

procs = [1, 1, 1]
topo.topology = 'closed'
myrank = 0
nh = 2
nz = 32
ny = 32
nx = 48
Lz = 1.0
Ly = 1.0
Lx = 1.0

loc = topo.rank2loc(myrank, procs)
neighbours = topo.get_neighbours(loc, procs)

param = {
    # I/O options
    "datadir": "~/data/Nyles",
    "expname": "buoyancy instab",
    "timestep_history": 0.0,
    "variables_in_history": "all",
    "mode": "overwrite",
    # General model options
    "modelname": "LES",
    "geometry": "closed",
    # Grid options
    "Lx": Lx,
    "Ly": Ly,
    "Lz": Lz,
    "nx": nx,
    "ny": ny,
    "nz": nz,
"nh": nh,
    # Timestepping options
    "timestepping": "EF",
    "tend": 1.0,
    "auto_dt": False,
    "cfl": 0.5,
    "dt": 0.1,
    # Spatial discretization options
    "orderVF": 1,  # upwind-order for vortex-force term
    "orderA": 1,  # upwind-order for advection term

    "neighbours": neighbours,
    "procs": procs, "topology": topo.topology,
    "npre": 3, "npost": 3, "omega": 0.8, "ndeepest": 20, "maxite": 20, "tol": 1e-6
}


t = 0.
dx = 1.
cfl = 0.85
dt = cfl*dx

nyles = nyles.Nyles(param)
grd = grid.Grid(param)

b = nyles.model.state.b.view('i')
z = grd.z_b.view('i')
x = grd.x_b.view('i')

b[:] = 15+ 5*np.tanh((np.cos(np.pi*x)*0.1+z-0.5)/.02)

noise = np.random.uniform(size = np.shape(b))
noise = noise*2-1

b += noise*.1

nyles.run()

plt.figure()
plt.pcolor(b[:,10,:])
plt.colorbar()

plt.show()
