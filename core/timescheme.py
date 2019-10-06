"""

timescheme.py provides the various time schemes

for the forward integration in time

it relies on the abstract 'rhs' function that is
defined outside of this module. The header for 'rhs'
is

    rhs(state, t, dstate)

where 'state' and 'dstate' are 'State' instances
and 't' is the model time

To allow for any time step, a forward() function is
defined. This function simply points to one of the
timeschemes. See explanations in the code below.

To push the model state one time step ahead we simply
do

    model.forward(state, t, dt)


Note that dstate is actually allocated *within* this
module during the __init__. The reason is that for multistages
timeschemes, we need to store more than just one dstate

"""


class Timescheme(object):
    """

    Catalog of time schemes

    It provides a forward() method that uses a generic 'rhs' function

    """

    def __init__(self, param, state):
        # list of prognostic *scalars* that should be
        # integrated in time
        # *scalars* => vectors should be passed as
        # list of components
        self.prognostic_variables = param['prognostic_variables']
        self.timestepping = param['timestepping']

        # dictionary of *functions*
        self.timeschemelist = {'EF': self.EulerForward,
                               # 'LF': self.LeapFrog,
                               # 'Heun': self.Heun,
                               # 'AB2': self.AB2,
                               # 'AB3': self.AB3,
                               'LFAM3': self.LFAM3}
        # 'RK3_SSP': self.RK3_SSP,
        # 'RK3': self.RK3,
        # 'RK4_LS': self.RK4_LS}

        if self.timestepping in self.timeschemelist.keys():
            pass
        else:
            raise ValueError('Unknown time scheme')

        # internal arrays
        self.dstate = state.duplicate()

        if self.timestepping == 'LFAM3':
            self.stateb = state.duplicate()
            self.state = state.duplicate()

        # for schemes with Leap-Frog time step
        self.first = True

    def set(self, rhs):
        """

        assign the 'rhs' and
        the 'timestepping' to the 'forward' method

        """
        self.rhs = rhs
        # HUGE trick: forward is a ... function
        # it points to one of the functions defined below
        # its header is
        # forward(state, t, dt, **kwargs)
        # like all functions below!
        self.forward = self.timeschemelist[self.timestepping]

    # ----------------------------------------
    def EulerForward(self, state, t, dt, **kwargs):
        self.rhs(state, t, self.dstate)
        # MR: the following loop can probably be optimized by making use of the
        #     new attribute "prognostic" of every variable.
        for v in self.prognostic_variables:
            # Get a view on the data without changing its orientation
            s = getattr(state, v).view()
            # Get a view on dstate in the same orientation as state
            ds = getattr(self.dstate, v).viewlike(getattr(state, v))
            s += dt * ds

    # ----------------------------------------
    def LFAM3(self, state, t, dt, **kwargs):
        # Predictor
        self.rhs(state, t, self.dstate)

        if self.first:
            # Euler Forward if very first time step
            # MR: the following loop can probably be optimized by making use of the
            #     new attribute "prognostic" of every variable.
            for v in self.prognostic_variables:
                s = getattr(state, v).view()
                ds = getattr(self.dstate, v).viewlike(getattr(state, v))
                sb = getattr(self.stateb, v).viewlike(getattr(state, v))
                sn = getattr(self.state, v).viewlike(getattr(state, v))
                sn[:] = s
                sb[:] = s
                s += dt * ds
            self.first = False

        else:
            # MR: the following loop can probably be optimized by making use of the
            #     new attribute "prognostic" of every variable.
            for v in self.prognostic_variables:
                s = getattr(state, v).view()
                ds = getattr(self.dstate, v).viewlike(getattr(state, v))
                sb = getattr(self.stateb, v).viewlike(getattr(state, v))
                sn = getattr(self.state, v).viewlike(getattr(state, v))
                # backup state into 'now' state
                sn[:] = s

                # LF s is at n+1
                s[:] = sb + (2.*dt) * ds
                # AM3 gives s at n+1/2
                s[:] = (1./12.)*(5.*s+8.*sn-sb)
                #
                # previous two lines combined in one => SLOWER
                # s[:] = (1./12.)*(5.*sb + (10*dt) * ds + 8.*sn-sb)

                # and backup former 'now' state into 'before' state
                sb[:] = sn

            # Corrector step at n+1/2
            self.rhs(state, t+dt*.5, self.dstate)

            # move from n to n+1
            # MR: the following loop can probably be optimized by making use of the
            #     new attribute "prognostic" of every variable.
            for v in self.prognostic_variables:
                s = getattr(state, v).view()
                ds = getattr(self.dstate, v).viewlike(getattr(state, v))
                sn = getattr(self.state, v).viewlike(getattr(state, v))
                s[:] = sn + dt*ds


# ----------------------------------------------------------------------
if __name__ == '__main__':

    from mpi import topology as topo
    import variables as var

    procs = [4, 2, 1]
    topology = 'closed'
    myrank = 3
    nh = 3

    loc = topo.rank2loc(myrank, procs)
    neighbours = topo.get_neighbours(loc, procs, topology)

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': nh,
             'neighbours': neighbours,
             'timestepping': 'LFAM3',
             'prognostic_variables': ['b', 'u_i', 'u_j', 'u_k']}

    state = var.get_state(param)

    timescheme = Timescheme(param, state)

    def rhs(state, t, dstate):
        # this routine should clearly update dstate using state
        # but to test this module, we need nothing more than
        # this empty function
        print('I am pretending to compute a rhs but actually I do nothing')

    timescheme.set(rhs)

    t = 0.
    dt = 0.1
    for kt in range(10):
        timescheme.forward(state, t, dt)
        t += dt
