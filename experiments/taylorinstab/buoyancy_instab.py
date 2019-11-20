import model_les_AL as LES
import topology as topo
import numpy as np
import matplotlib.pyplot as plt

import nyles

procs = [1, 1, 1]
topo.topology = 'closed'
myrank = 0
nh = 2
nz = 32
ny = 32
nx = 128

loc = topo.rank2loc(myrank, procs)
neighbours = topo.get_neighbours(loc, procs)

param = {
    # I/O options
    "datadir": "~/data/Nyles",
    "expname": "buoyancy instab",
    "timestep_history": 1.0,
    "variables_in_history": ["vor","b"],
    "mode": "overwrite",
    # General model options
    "modelname": "LES",
    "geometry": "closed",
    # Grid options
    "nx": nx,
    "ny": ny,
    "nz": nz,
    "nh": nh,
    # Timestepping options
    "timestepping": "LFAM3",
    "tend": 10.0,
    "cfl": 1.5,
    # Spatial discretization options
    "orderVF": 5,  # upwind-order for vortex-force term
    "orderA": 5,  # upwind-order for advection term

    "neighbours": neighbours,
    "procs": procs, "topology": topo.topology,
    "npre": 3, "npost": 3, "omega": 0.8, "ndeepest": 20, "maxite": 20, "tol": 1e-6
}


t = 0.
dx = 1.
cfl = 0.85
dt = cfl*dx

nyles = nyles.Nyles(param)

b = nyles.model.state.b.view('k')
b[:,:,:nz//2] = 20
b[:,:,nz//2:] = 10

noise = np.random.normal(size = np.shape(b))
noise = noise/np.max(noise)

b = b + noise

nyles.run()

plt.figure()
plt.pcolor(b[0,:,:])
plt.colorbar()

plt.show()
