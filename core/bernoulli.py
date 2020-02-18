"""Module to compute the bernoulli term without pressure.

TODO:
  to increase computational intensity, compute delta[b*grad(z) - grad(ke)]
  and add it to the rhs of the multigrid

"""

import fortran_bernoulli as fortran
from timing import timing


@timing
def bernoulli(state, rhs, grid):
    """Add b*grad(z)-grad(ke) to the rhs."""
    for i in 'ijk':
        du_i = rhs.u[i].view(i)
        ke = state.ke.view(i)
        if i in 'ij':
            fortran.gradke(ke, du_i)
        else:
            b = state.b.view(i)
            fortran.gradkeandb(ke, b, du_i, grid.dz)

@timing
def bernoulli2D(state, rhs, grid):
    """Add -grad(ke+g*h) to the rhs."""

    g = 1.
    for i in 'ij':
        du_i = rhs.u[i].view(i)
        ke = state.ke.view(i)
        h = state.h.view(i)
        fortran.gradber(ke, h, du_i, g)


# ----------------------------------------------------------------------
if __name__ == '__main__':
    import variables as var
    from grid import Grid

    param = {
        'nx': 40, 'ny': 50, 'nz': 60,
        'Lx': 1.0, 'Ly': 1.0, 'Lz': 1.0,
        'neighbours': {}, 'nh': 3,
    }

    grid = Grid(param)
    state = var.get_state(param)
    ds = state.duplicate_prognostic_variables()

    bernoulli(state, ds, grid)
