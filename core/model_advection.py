"""

provides all the tools to advect a tracer with
an imposed velocity field

"""

import numpy as np
import variables as var
import tracer as tracer
import timescheme as ts
from matplotlib import pyplot as plt


class Advection(object):
    """
    Pure advection model

    It advects the buoyancy with a prescribed
    contravariant velocity

    """
    def __init__(self, param):
        self.state = var.get_state(param)
        self.traclist = ['b']
        self.timescheme = ts.Timescheme(param, self.state)
        self.timescheme.set(self.rhs)

    def forward(self, t, dt):
        self.timescheme.forward(self.state, t, dt)

    def rhs(self, state, t, dstate):
        tracer.rhstrac(state, dstate, self.traclist)


if __name__ == '__main__':

    nz, ny, nx = 32, 32, 128
    param = {
        'nx': nx, 'ny': ny, 'nz': nz, 'nh': 2,
        'timestepping': 'LFAM3',
        # For test purposes, you can also try the following:
        # 'timestepping': 'EF',
    }

    model = Advection(param)

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
