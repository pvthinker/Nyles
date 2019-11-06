"""

module to compute -div( vol * trac )

"""
import fortran_upwind as fortran
import variables as var
from timing import timing

@timing
def rhstrac(state, rhs, traclist, order):
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
            ds2 = 1.  # 1/dx**2
            vol = 1.  # vol=dx*dy*dz
            cff = vol/ds2  # for z coordinates

            component = state.U[direction]
            field = trac.view(direction)
            dfield = dtrac.view(direction)
            velocity = component.view(direction)

            if direction == 'i':
                # overwrite rhs
                fortran.upwind(field, velocity, dfield, cff, order, 1)
            else:
                # add to rhs
                fortran.upwind(field, velocity, dfield, cff, order, 0)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}
    state = var.get_state(param)
    ds = state.duplicate_prognostic_variables()

    rhstrac(state, ds, ['b'], 3)
