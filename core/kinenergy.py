import fortran_kinenergy as fortran
import variables as var
from timing import timing


@timing
def kinenergy(state):
    """
    compute the kinetic energy function from the model state

    """
    # loop over direction
    for direction in 'ijk':
        # Get array of covariant velocity u
        u = state.u[direction].view(direction)
        # In sigma-coordinates, use the contravariant velocity U
        # U = state.U[direction].view(direction)

        # Get array of kinetic energy ke
        ke = state.ke.view(direction)
        # In z-coordinates, the kinetic energy is ke = u**2 / ds**2,
        # where ds2 is the diagonal term of the inverse metric tensor
        ds2 = 1.  # TODO: use dx**-2, dy**-2, dz**-2

        if direction == 'i':
            fortran.kin(u, u, ke, ds2, 1)  # overwrite ke
        else:
            fortran.kin(u, u, ke, ds2, 0)  # add to ke


# ----------------------------------------------------------------------
if __name__ == '__main__':

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}
    state = var.get_state(param)

    kinenergy(state)
