"""

provides all the tools to advect a tracer with
an imposed velocity field

"""

import numpy as np
import variables as var
import tracer as tracer
import timescheme as ts
from matplotlib import pyplot as plt
import topology as topo

class Advection(object):
    """
    Pure advection model

    It advects the buoyancy with a prescribed
    contravariant velocity

    """
    def __init__(self, param, grid):
        self.grid = grid
        self.state = var.get_state(param)
        self.traclist = ['b']
        self.order = param['orderA']
        self.timescheme = ts.Timescheme(param, self.state)
        self.timescheme.set(self.rhs)

    def forward(self, t, dt):
        self.timescheme.forward(self.state, t, dt)

    def rhs(self, state, t, dstate):
        tracer.rhstrac(state, dstate, self.grid, self.traclist, self.order)


if __name__ == '__main__':
    from grid import Grid


    procs = [1, 1, 1]
    topo.topology = 'closed'
    myrank = 0
    nh = 2

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    nz, ny, nx = 32, 32, 128
    param = {
        'nx': nx, 'ny': ny, 'nz': nz, 'nh': nh,
        'Lx': 1.0, 'Ly': 1.0, 'Lz': 1.0,
        'neighbours': neighbours,
        # Choose a timestepping method
        'timestepping': 'LFAM3',
        #'timestepping': 'EF',
        # Set the order of the upwind scheme
        'orderA': 5,
        #'orderA': 3,
        #'orderA': 1,
    }

    grid = Grid(param)

    model = Advection(param, grid)

    # set up a uniform velocity along i
    # It does not enforce the no-flow BC,
    # but who cares as long as we don't
    # integrate too long
    Ui = model.state.U['i'].view()
    Ui[:, :, :] = 1.

    # minimalist grid (dx should be equal to 1 to be consistent
    # with the hard-coded dx in tracer.py)
    dx = 1.
    Lx = nx*dx

    # set up a gaussian 'b' along i
    b = model.state.b.view()
    for i in range(nx):
        x = (i+0.5)*dx-0.2*Lx
        b[:, :, i] = np.exp(-0.5*(x/(dx*5))**2)

    plt.figure()
    t = 0.
    cfl = 0.85
    dt = cfl*dx
    for kt in range(41):
        model.forward(t, dt)
        t += dt
        if kt % 5 == 0:
            plt.plot(b[0, 0, :], label='%i' % kt)
    plt.legend()
    plt.show()
