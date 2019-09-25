"""

timescheme.py provides the various time schemes

for the forward integration in time

it relies on the abstract 'rhs' function that is
defined outside of this module. The header for 'rhs'
is

    rhs(state, t, dstate)

where 'state' and 'dstate' are 'State' instances
and t is the model time

To allow for any time step, a forward() function is
defined. This function simply points to one of the
timescheme. See explanations in the code.

To push the model state one time step ahead we simply
do

model.forward(state, t, dt)


note that dstate is actually allocated *within* this
module during the __init__. The reason is that for multistages
timeschemes, we need to store more than just one dstate

"""


class Timescheme(object):
    """

    Catalog of time schemes

    It provides a forward() method that uses a generic rhs function

    """

    def __init__(self, param, state):
        self.prognostic_variables = ['b', 'u']
        timestepping = param['timestepping']

        # dictionnary of *functions*

        self.timeschemelist = {'EF': self.EulerForward,
                               # 'LF': self.LeapFrog,
                               # 'Heun': self.Heun,
                               # 'AB2': self.AB2,
                               # 'AB3': self.AB3,
                               'LFAM3': self.LFAM3}
        # 'RK3_SSP': self.RK3_SSP,
        # 'RK3': self.RK3,
        # 'RK4_LS': self.RK4_LS}
        # internal arrays
        self.dstate = state.duplicate()
        if timestepping == 'LFAM3':
            self.stateb = state.duplicate()
            self.staten = state.duplicate()

        # for schemes with Leap-Frog time step
        self.first = True

    def set(self, rhs, timestepping):
        """

        assign the 'rhs' and
        the 'timestepping' to the 'forward' method

        """
        self.rhs = rhs
        # HUGE trick: forward is a ... function
        # that is one of the function defined below
        # its header is
        # forward(state, t, dt, **kwargs)
        # like all functions below
        self.forward = self.timeschemelist[timestepping]

    # ----------------------------------------
    def EulerForward(self, state, t, dt, **kwargs):
        self.rhs(state, t, self.dstate)
        for var in self.prognostic_variables:
            s = state.get(var).view()
            # make sure that state and dstate have the same convention
            ds = self.dstate.get(var).viewlike(s)
            s += dt * ds

    # ----------------------------------------
    def LFAM3(self, state, t, dt, **kwargs):
        # Predictor
        self.rhs(state, t, self.dstate)

        if self.first:
            # Euler Forward
            for var in self.prognostic_variables:
                s = state.get(var).view()
                ds = self.dstate.get(var).viewlike(s)
                s += dt * ds
            self.first = False

        else:
            for var in self.prognostic_variables:
                s = state.get(var).view()
                ds = self.dstate.get(var).viewlike(s)
                sb = self.stateb.get(var).viewlike(s)
                sn = self.state.get(var).viewlike(s)
                # LF s is at n+1
                s[:] = sb + (2*dt) * ds
                # AM3 gives s at n+1/2
                s[:] = (1./12.)*(5.*s + 8.*sn-sb)

            # Corrector step at n+1/2
            self.rhs(state, t+dt*.5, self.dstate)

            # move from n to n+1
            for var in self.prognostic_variables:
                s = state.get(var).view()
                ds = self.dstate.get(var).viewlike(s)
                s[:] = sn + dt*ds

        # backup former 'now' state into 'before' state
        for var in self.prognostic_variables:
            sn = self.state.get(var).view()
            sb = self.stateb.get(var).viewlike(s)
            sb[:] = sn
            # self.stateb.copyfrom(self.state, self.prognostic_variables)
