"""Module to compute -div( tracer * U ) for Cartesian coordinates."""

import fortran_upwind as fortran
import variables as var
from timing import timing


@timing
def rhstrac(state, rhs, traclist, order):
    """Compute negative advection of a tracer in Cartesian coordinates.

    Arguments:
     - state: State instance containing the current value of the tracer
       and the contravariant velocity
     - rhs: State instance containing prognostic variables to save the
       result within
     - traclist: list of prognostic variables that are advected, e.g.,
       traclist = ['b']
     - order: can be 1, 3, or 5; sets the order of the upwind scheme
       used to calculate the advection of the tracer
    """

    assert(order in {1, 3, 5})

    for tracname in traclist:
        trac = state.get(tracname)  # trac is a 'Scalar' instance
        dtrac = rhs.get(tracname)

        for direction in 'ijk':
            velocity = state.U[direction].view(direction)
            field = trac.view(direction)
            dfield = dtrac.view(direction)

            if direction == 'i':
                # overwrite rhs
                fortran.upwind(field, velocity, dfield, order, 1)
            else:
                # add to rhs
                fortran.upwind(field, velocity, dfield, order, 0)


# ----------------------------------------------------------------------
if __name__ == '__main__':
    param = {
        'nx': 64, 'ny': 32, 'nz': 128,
        'Lx': 1.0, 'Ly': 1.0, 'Lz': 1.0,
        'nh': 2, 'neighbours': {},
    }
    state = var.get_state(param)
    ds = state.duplicate_prognostic_variables()

    rhstrac(state, ds, ['b'], order=3)
