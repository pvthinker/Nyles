"""

add b*grad(z) + grad(ke) to rhs

TODO:
  to increase computational intensity, compute delta[b*grad(z) + grad(ke)]
  and add it to the rhs of the multigrid

"""

import fortran_bernoulli as fortran
import variables as var
from timing import timing


@timing
def bernoulli(state, rhs):

    du = rhs.get('u')
    b = state.get('b')
    ke = state.get('ke')

    dz = 1. # the vertical grid size

    for d in 'ijk':
        comp = getattr(du, d)
        if d in 'ij':
            fortran.gradke(ke.view(d), comp.view(d))

        else:
            fortran.gradkeandb(ke.view(d), b.view(d), comp.view(d), dz)


# ----------------------------------------------------------------------
if __name__ == '__main__':

    param = {'nx': 40, 'ny': 50, 'nz': 60, 'nh': 2}
    state = var.get_state(param)
    ds = state.duplicate()

    bernoulli(state, ds)
