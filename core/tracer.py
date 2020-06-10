"""Module to compute -div( tracer * U ) for Cartesian coordinates."""

import fortran_upwind as fortran
import fortran_dissipation as dissip
import variables as var
from timing import timing


class Tracer_numerics(object):
    def __init__(self, param, grid, traclist, order, diff_coef=[]):
        """
        - grid: nyles grid object
        - traclist: list of prognostic variables that are advected, e.g.,
        traclist = ['b']
        - order: can be 1, 3, or 5; sets the order of the upwind scheme
        used to calculate the advection of the tracer
        - diffcoef is a dictionnary with diffusion coefficient for each variable
        """
        self.traclist = traclist
        assert(order in {1, 2, 3, 4, 5})
        self.order = order
        diffusion = len(diff_coef) > 0
        self.diffusion = diffusion
        if diffusion:
            self.ids2 = grid.ids2
            self.diff_coef = diff_coef
        self.i0 = {}
        self.i1 = {}
        ngbs = param["neighbours"]
        for d in 'ijk':
            i0 = 0
            i1 = 0
            if d == 'i' and not((0, 0, -1) in ngbs.keys()): i0 = 1
            if d == 'i' and not((0, 0, +1) in ngbs.keys()): i1 = 1
            if d == 'j' and not((0, -1, 0) in ngbs.keys()): i0 = 1
            if d == 'j' and not((0, +1, 0) in ngbs.keys()): i1 = 1
            if d == 'k' and not((-1, 0, 0) in ngbs.keys()): i0 = 1
            if d == 'k' and not((+1, 0, 0) in ngbs.keys()): i1 = 1

            self.i0[d] = i0
            self.i1[d] = i1

    @timing
    def rhstrac(self, state, rhs, last=False):
        """Compute negative advection of a tracer in Cartesian coordinates.

        Note: diffusion term should be integrated in time with a Euler forward
        For multi-stages timeschemes handling diffusion, set add_diffusion = True
        only for the last stage.


        Arguments:
         - state: State instance containing the current value of the tracer
           and the contravariant velocity
         - rhs: State instance containing prognostic variables to save the
           result within
        """

        for tracname in self.traclist:
            trac = state.get(tracname)  # trac is a 'Scalar' instance
            dtrac = rhs.get(tracname)

            for direction in 'ijk':
                velocity = state.U[direction].view(direction)
                field = trac.view(direction)
                dfield = dtrac.view(direction)
                i0 = self.i0[direction]
                i1 = self.i1[direction]

                fortran.upwind(field, velocity, dfield, self.order)

                if self.diffusion and last:
                    if (tracname in self.diff_coef.keys()):
                        coef = self.diff_coef[tracname]*self.ids2[direction]
                        dissip.add_laplacian(field, dfield, coef)


# ----------------------------------------------------------------------
if __name__ == '__main__':
    param = {
        'nx': 64, 'ny': 32, 'nz': 128,
        'Lx': 1.0, 'Ly': 1.0, 'Lz': 1.0,
        'nh': 2, 'neighbours': {},
    }
    state = var.get_state(param)
    ds = state.duplicate_prognostic_variables()

    grid = None
    tracnum = Tracer_numerics(grid, ['b'], order=3)
    tracnum.rhstrac(state, ds)
