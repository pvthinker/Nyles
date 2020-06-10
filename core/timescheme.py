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
    """Catalogue and handler of timeschemes.

    The user mainly interacts with an object of this class through its
    "forward" method, which uses an "rhs" function to calculate the
    "right hand side" of the prognostic model equation.  This "rhs"
    function must be set by calling the "set" method.

    Attributes:
     - prognostic_scalars: list with the names of prognostic scalars
        in the model, which can be Scalar variables or the components
        of Vector variables
     - dstate: internal instance of the State class used to store
        intermediate steps in the calculation; this State object
        contains only prognostic variables
     - other attributes may exist depending on the chosen timescheme

    Methods:
     - set
     - rhs: this method is not defined within this class definition,
        but must be implemented in the model and given to Timescheme
        by calling the "set" method
     - forward: this method is defined in the constructor as an alias
        for one of the following time stepping methods
     - EulerForward
     - LFAM3

    """

    def __init__(self, param, state):
        """Initialise the time stepping method specified in param."""
        self.prognostic_scalars = state.get_prognostic_scalars()

        # Activate the timescheme given in param
        timestepping = param['timestepping']
        # ... which must be in the following dictionary of *functions*
        timeschemelist = {
            'EF': self.EulerForward,
            # 'LF': self.LeapFrog,
            # 'Heun': self.Heun,
            # 'AB2': self.AB2,
            # 'AB3': self.AB3,
            'LFAM3': self.LFAM3,
            'RK3_SSP': self.RK3_SSP,
            # 'RK3': self.RK3,
            # 'RK4_LS': self.RK4_LS,
        }
        try:
            self.forward = timeschemelist[timestepping]
            # This is a HUGE trick: forward is a function (and not a
            # variable).  It points to one of the timestepping
            # functions defined below.  They all have the same header:
            # forward(state, t, dt, **kwargs)
        except KeyError:
            raise ValueError('unknown time scheme ' + repr(timestepping))

        # Create internal arrays
        self.dstate = state.duplicate_prognostic_variables()
        if timestepping == 'LFAM3':
            self.stateb = state.duplicate_prognostic_variables()
            self.state = state.duplicate_prognostic_variables()
            self.first = True
        if timestepping == 'RK3_SSP':
            self.ds0 = self.dstate
            self.ds1 = state.duplicate_prognostic_variables()
            self.ds2 = state.duplicate_prognostic_variables()


    def set(self, rhs, diagnose_var):
        """Assign the right hand side of the model.

        The argument 'rhs' of this method is a function with the
        signature rhs(state, t, dstate), where 'state' is the current
        state of the model and 't' is the current timestep.  The
        result is written to 'dstate'.
        """
        self.rhs = rhs
        self.diagnose_var = diagnose_var

    # ----------------------------------------
    def EulerForward(self, state, t, dt, **kwargs):
        self.rhs(state, t, self.dstate, last=True)
        for scalar_name in self.prognostic_scalars:
            scalar = state.get(scalar_name)
            # Get a view on the data without changing its orientation
            s = scalar.view()
            # Get a view on dstate in the same orientation as state
            ds = self.dstate.get(scalar_name).viewlike(scalar)
            s += dt * ds
        self.diagnose_var(state)

    # ----------------------------------------
    def LFAM3(self, state, t, dt, **kwargs):

        if self.first:
            self.rhs(state, t, self.dstate, last=True)
            # Euler Forward if very first time step
            for scalar_name in self.prognostic_scalars:
                scalar = state.get(scalar_name)
                s = scalar.view()
                ds = self.dstate.get(scalar_name).viewlike(scalar)
                sb = self.stateb.get(scalar_name).viewlike(scalar)
                sn = self.state.get(scalar_name).viewlike(scalar)
                sn[:] = s
                sb[:] = s
                s += dt * ds
            self.first = False
            self.diagnose_var(state)

        else:
            # Predictor
            self.rhs(state, t, self.dstate)
            for scalar_name in self.prognostic_scalars:
                scalar = state.get(scalar_name)
                s = scalar.view()
                ds = self.dstate.get(scalar_name).viewlike(scalar)
                sb = self.stateb.get(scalar_name).viewlike(scalar)
                sn = self.state.get(scalar_name).viewlike(scalar)
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

            self.diagnose_var(state)

            # Corrector step at n+1/2
            self.rhs(state, t+dt*.5, self.dstate, last=True)

            # move from n to n+1
            for scalar_name in self.prognostic_scalars:
                scalar = state.get(scalar_name)
                s = scalar.view()
                ds = self.dstate.get(scalar_name).viewlike(scalar)
                sn = self.state.get(scalar_name).viewlike(scalar)
                s[:] = sn + dt*ds

            self.diagnose_var(state)
    # ----------------------------------------

    def RK3_SSP(self, state, t, dt, **kwargs):
        """ RK3 SSP

        The three stages are

        s1 = s^n + dt*L(s^n)
        s2 = s^n + (dt/4)*( L(s^n)+L(s1) )
        s^n+1 =  s^n + (dt/6)*( L(s^n)+L(s1)+4*L(s2) )

        or equivalently

        ds0 = L(s)
        s = s+dt*ds0
        ds1 = L(s)
        s = s+(dt/4)*(ds1-3*ds0)
        ds2 = L(s)
        s = s+(dt/12)*(8*ds2-ds0-ds1)

        """
        self.rhs(state, t, self.ds0, last=False)
        for scalar_name in self.prognostic_scalars:
            s = state.get(scalar_name).view("i")
            ds = self.ds0.get(scalar_name).view("i")
            s += dt * ds
        self.diagnose_var(state)

        self.rhs(state, t+dt, self.ds1, last=False)
        for scalar_name in self.prognostic_scalars:
            s = state.get(scalar_name).view("i")
            ds0 = self.ds0.get(scalar_name).view("i")
            ds1 = self.ds1.get(scalar_name).view("i")
            s += (dt/4.) * (ds1-3*ds0)
        self.diagnose_var(state)

        self.rhs(state, t+dt*0.5, self.ds2, last=True)
        for scalar_name in self.prognostic_scalars:
            s = state.get(scalar_name).view("i")
            ds0 = self.ds0.get(scalar_name).view("i")
            ds1 = self.ds1.get(scalar_name).view("i")
            ds2 = self.ds2.get(scalar_name).view("i")
            s += (dt/12.) * (8*ds2-ds0-ds1)
        self.diagnose_var(state)


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
             }

    state = var.get_state(param)

    timescheme = Timescheme(param, state)

    # A counter to see how often the dummy "rhs" function defined below is invoked.
    i = 0

    def rhs(state, t, dstate):
        # this routine should clearly update dstate using state
        # but to test this module, we need nothing more than
        # this empty function
        global i
        print(i, 'I am pretending to compute a rhs but actually I am just counting.')
        i += 1

    timescheme.set(rhs)

    t = 0.
    dt = 0.1
    for kt in range(10):
        timescheme.forward(state, t, dt)
        t += dt
