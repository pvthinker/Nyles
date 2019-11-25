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
    def __init__(self, param):
        self.state = var.get_state(param)
        self.traclist = ['b']
        self.order = param['orderA']
        self.timescheme = ts.Timescheme(param, self.state)
        self.timescheme.set(self.rhs)

    def forward(self, t, dt):
        self.timescheme.forward(self.state, t, dt)

    def rhs(self, state, t, dstate):
        tracer.rhstrac(state, dstate, self.traclist, self.order)


if __name__ == '__main__':
    procs = [1, 1, 1]
    topo.topology = 'closed'
    myrank = 0
    nh = 2

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs)

    nz, ny, nx = 32, 32, 128
    Lx, Ly, Lz = 1.0, 1.0, 1.0
    dx, dy, dz = Lx/nx, Ly/ny, Lz/nz
    param = {
        'nx': nx, 'ny': ny, 'nz': nz, 'nh': nh,
        'Lx': Lx, 'Ly': Ly, 'Lz': Lz,
        'neighbours': neighbours,
        # Choose a timestepping method
        'timestepping': 'LFAM3',
        #'timestepping': 'EF',
        # Set the order of the upwind scheme
        'orderA': 5,
        #'orderA': 3,
        #'orderA': 1,
    }

    model = Advection(param)

    # set up a uniform velocity along i
    # It does not enforce the no-flow BC,
    # but who cares as long as we don't
    # integrate too long
    U0 = 1.0
    Ui = model.state.U['i'].view()
    Ui[...] = U0

    # set up a gaussian 'b' along i
    b = model.state.b.view()
    for i in range(nx):
        x = (i+0.5)*dx-0.2*Lx
        b[:, :, i] = np.exp(-0.5*(x/(dx*5))**2)

    plt.figure()
    t = 0.
    cfl = 0.5
    dt = cfl * U0
    for kt in range(41):
        model.forward(t, dt)
        t += dt
        if kt % 5 == 0:
            plt.plot(b[0, 0, :], label='%i' % kt)
    plt.legend()
    plt.show()
