import vortex_force as vortf


class VFwork:
    """ diagnose the point wise vortex-force work"""

    def __init__(self, model, grid):
        self.model = model
        self.grid = grid

    def compute(self):
        state = self.model.state
        dstate = self.model.timescheme.dstate

        self.compute_vortexforce(state, dstate)
        self.innerproduct(state.U, dstate.u, state.work)
        k0, k1, j0, j1, i0, i1 = self.grid.domainindices
        self.worksum = state.work.view()[k0:k1, j0:j1, i0:i1].sum()

    def compute_vortexforce(self, state, dstate):
        for i in "ijk":
            var = dstate.get("u")[i].view()
            var[:] = 0.0

        vortf.vortex_force(state, dstate, self.model.orderVF)

    def innerproduct(self, vec1, vec2, result):
        for direction in 'ijk':
            u1 = vec1[direction].view(direction)
            u2 = vec2[direction].view(direction)
            res = result.view(direction)
            if direction == "i":
                res[:] = 0.
            product = u1*u2*0.5
            res[:, :, 1:] += product[:, :, 1:]+product[:, :, :-1]
