"""Implementation of the LES model proposed by Markus REINERT, 25 Oct 2019.

The update of the state variables is implemented like this:
 1. Calculate delta_b, delta_u, delta_U and update p in compute_rhs.
 2. Increment b, u, U by their respective deltas in the timescheme.
 3. Calculate omega and kinetic energy in compute_diagnostic_variables.

Advantages:
 - The velocity has directly the correct value instead of an
   intermediate value.
 - No need to handle p separately.

Disadvantage:
 - Non-divergence is not imposed explicitly, but only implicitly.

Alternatively, it is possible to move the calculation of U into step 3.
This increases the amount of calculations necessary by one matrix
multiplication, but restores the distinction between prognostic and
diagnostic variables.  Also it prevents the accumulation of errors in
U with respect to u.

Conversely, it is possible to handle omega like a prognostic variable.
This, together with removing KE from the model variables, allows to get
rid of step 3, in which case no modification of the time-stepping is
necessary.  However, this has the potential risk that small errors in
the calculation of delta_omega accumulate and thus the relation between
u and omega is violated over time.  Also it removes the distinction
between prognostic and diagnostic variables almost completely, except
for p.
"""

from variables import get_state
from timescheme import Timescheme
from vortex_force import vortex_force
from bernoulli import bernoulli
from kinenergy import kinenergy
from vorticity import vorticity_all_comp as vorticity
from timing import timing
import fortran_upwind


class LES(object):
    def __init__(self, param):
        # Create a State object for the LES model
        self.state = get_state(param)
        # Make U and all its components prognostic.  This should
        # actually be done before the state is created, but I don't want
        # to mess with variables.py at the moment.  This must be done
        # before the timescheme is created.
        self.state.U.prognostic = True
        for i in "ijk":
            self.state.U[i].prognostic = True
        # Create and set-up a handler for the timescheme
        self.timescheme = Timescheme(param, self.state)
        self.timescheme.set(compute_rhs)
        # TODO: also inform the timescheme about the function to compute diagnostic variables

    def forward(self, t, dt):
        self.timescheme.forward(self.state, t, dt)


def compute_rhs(state, t, dstate):
    """Compute deltas for b, u, U and update p.

    The deltas for b, u, U are saved in dstate, while the updated p is
    saved in state.  The time t has currently no influence.  The
    contravariant velocity U is handled like a prognostic variable,
    because its delta is a by-product of the pressure-calculation, even
    though it does not appear in a time derivative of the LES equations.
    """
    # Advect buoyancy as in rhstrac
    # TOOO: these should not be hardcoded
    ds2 = 1.  # 1/dx**2
    vol = 1.  # vol=dx*dy*dz
    cff = vol/ds2  # for z coordinates
    order = 5
    for i in 'ijk':
        # Overwrite the value in dstate if i == 'i', otherwise add the new value to dstate
        fortran_upwind.upwind(
            state.b.view(i), state.U[i].view(i), dstate.b.view(i), cff, int(i == 'i')
        )
    # Calculate time-derivative of u except for the p-term
    vortex_force(state, dstate, order)
    bernoulli(state, dstate)  # This assumes that KE was already computed!
    # Calculate the p-term, update the delta_u and calculate delta_U
    calculate_U_from_u(dstate)  # This calculates here actually dU
                                # (dstate.U) from du (dstate.u), but the
                                # function is more general, as its name
                                # suggests
    calculate_p_from_dU(state, dstate)  # This solves the poisson
                                        # equation with dU (dstate.U),
                                        # stores the result in p
                                        # (state.p) and updates dU
    calculate_u_from_U(dstate)  # This calculates here actually du
                                # (dstate.u) from dU (dstate.U), but the
                                # function is more general, as its name
                                # suggests


def compute_diagnostic_variables(state):
    """Update vorticity and kinetic energy in state."""
    # Calculate the new vorticity
    vorticity(state)
    # Calculate the new kinetic energy
    kinenergy(state)


@timing
def calculate_U_from_u(state):
    # This is copied from lotsofstuff.py with slight modifications
    # TODO: replace these dummy values by the actual ones
    idx2, idy2, idz2 = 1., 1., 1.

    # Several options for the metric of which currently only one is implemented
    # TODO: remove this when the other coordinate systems are implemented
    # metric = 'sigma2D'  # like ROMS/CROCO, h(x,y)
    # metric = 'sigma1D'  # h(x)
    metric = 'cartesian'  # dx, dy, dz are uniform in space, but not necessarily equal

    if metric == 'cartesian':
        u = state.u['i'].view('k')
        v = state.u['j'].view('k')
        w = state.u['k'].view('k')

        U = state.U['i'].view('k')
        V = state.U['j'].view('k')
        W = state.U['k'].view('k')

        U[:] = u * idx2
        V[:] = v * idy2
        W[:] = w * idz2

    elif metric == 'sigma1D':
        # TODO
        raise NotImplementedError(
            'metric "sigma1D" requires slope and gamma, which are not yet defined'
        )
        def sxp(array2d):
            return np.roll(array2d, 1, axis=1)
        def syp(array2d):
            return np.roll(array2d, 1, axis=0)
        def sxm(array2d):
            return np.roll(array2d, -1, axis=1)
        def sym(array2d):
            return np.roll(array2d, -1, axis=0)
        # cf Roullet et al, OM2017
        # slope = dz/dx, at cell center
        u = state.u['i'].view('j')
        v = state.u['j'].view('j')
        w = state.u['k'].view('j')

        U = state.U['i'].view('j')
        V = state.U['j'].view('j')
        W = state.U['k'].view('j')

        V[:] = v * idy2
        for j in range(ny):
            U[j, :, :] = u[j, :, :] - sxp(slope[j][:, :]*sym(w[j, :, :]))
            W[j, :, :] = gamma[j][:, :]*w[j][:, :] - syp(slope[j][:, :]*sxm(u[j, :, :]))

    else:
        # TODO: possibly implement sigma2D
        raise ValueError("unknown metric:", repr(metric), '(or not yet defined)')


@timing
def calculate_u_from_U(state):
    """Make the inverse operation of calculate_U_from_u."""
    # TODO: replace these dummy values by the actual ones
    idx2, idy2, idz2 = 1., 1., 1.
    metric = 'cartesian'

    if metric == 'cartesian':
        u = state.u['i'].view('k')
        v = state.u['j'].view('k')
        w = state.u['k'].view('k')

        U = state.U['i'].view('k')
        V = state.U['j'].view('k')
        W = state.U['k'].view('k')

        u[:] = U / idx2
        v[:] = V / idy2
        w[:] = W / idz2
    else:
        # TODO: implement other coordinate systems
        raise ValueError("unknown metric:", repr(metric), '(or not yet implemented)')


@timing
def calculate_p_from_dU(state, dstate):
    # TODO: implement
    if __name__ == "__main__":
        print("!!! calculate_p_from_dU is not yet implemented.")
    else:
        raise NotImplementedError("calculate_p_from_dU is not yet implemented.")


if __name__ == "__main__":
    param = {
        'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2,
        'timestepping': 'LFAM3',
    }
    l = LES(param)
    print("step 1")
    l.forward(t=0, dt=0.1)
    print("step 2")
    l.forward(t=0.1, dt=0.1)
