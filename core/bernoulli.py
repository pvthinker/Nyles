"""

add b*grad(z) + grad(ke) to rhs

TODO:
  to increase computational intensity, compute delta[b*grad(z) + grad(ke)]
  and add it to the rhs of the multigrid

"""

import fortran_bernoulli as fortran
from timing import timing


@timing
def bernoulli(state, rhs):

    dz = 1. # the vertical grid size

    for direction in 'ijk':
        u_component = rhs.u[direction]
        if direction in 'ij':
            fortran.gradke(state.ke.view(direction), u_component.view(direction))
        else:
            fortran.gradkeandb(state.ke.view(direction), state.b.view(direction),
                               u_component.view(direction), dz)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    import variables as var
    
    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}
    state = var.get_state(param)
    ds = state.duplicate()

    bernoulli(state, ds)
