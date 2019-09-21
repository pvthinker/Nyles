"""

module to compute -div( vol * trac )

"""
import fortran_upwind as fortran
import variables as var
from timing import timing


@timing
def rhstrac(state, rhs, traclist):
    """

    traclist = ['b']

    or a longer list, each variable is present in both state and rhs

    """
    # in sigma coordinates use 'U' the contravariant velocity
    # U = state.get('U') # U is a 'Velocity' instance
    # in z coordinates use 'u' the covariant
    # and account for the metric term by tweaking the volume vol=>vol/ds**2
    U = state.get('u')

    for tracname in traclist:
        trac = state.get(tracname)  # trac is a 'Scalar' instance
        dtrac = rhs.get(tracname)

        for k, direction in enumerate('ijk'):
            ds2 = 1.  # 1/dx**2
            vol = 1.  # vol=dx*dy*dz
            cff = vol/ds2  # for z coordinates

            component = getattr(U, direction)

            field = trac.view(direction)
            dfield = dtrac.view(direction)
            velocity = component.view(direction)

            if k == 0:
                # overwrite rhs
                fortran.upwind(field, velocity, dfield, cff, 1)
            else:
                # add to rhs
                fortran.upwind(field, velocity, dfield, cff, 0)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}
    state = var.get_state(param)
    ds = state.duplicate()

    rhstrac(state, ds, ['b'])
