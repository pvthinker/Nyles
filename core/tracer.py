"""

module to compute -div( vol * trac )

"""
import fortran_upwind as fortran
import variables as var
from timing import timing

@timing
def rhstrac(state, rhs, grid, traclist, order):
    """

    traclist = ['b']

    or a longer list, each variable is present in both state and rhs,
    i.e., each variable in traclist must be prognostic

    """

    upw_orders = {1,3,5}

    assert(order in upw_orders)

    # in sigma coordinates use the contravariant velocity state.U
    # in z coordinates use the covariant velocity state.u
    # and account for the metric term by tweaking the volume vol=>vol/ds**2
    for tracname in traclist:
        trac = state.get(tracname)  # trac is a 'Scalar' instance
        dtrac = rhs.get(tracname)

        for direction in 'ijk':
            cff = grid.vol_per_ds2[direction]  # for z coordinates

            velocity = state.U[direction].view(direction)
            field = trac.view(direction)
            dfield = dtrac.view(direction)

            if direction == 'i':
                # overwrite rhs
                fortran.upwind(field, velocity, dfield, cff, order, 1)
            else:
                # add to rhs
                fortran.upwind(field, velocity, dfield, cff, order, 0)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    from grid import Grid


    param = {
        'nx': 64, 'ny': 32, 'nz': 128,
        'Lx': 1.0, 'Ly': 1.0, 'Lz': 1.0,
        'nh': 2, 'neighbours': {},
    }
    state = var.get_state(param)
    grid = Grid(param)

    ds = state.duplicate_prognostic_variables()

    rhstrac(state, ds, grid, ['b'], order=3)
