
class Snapshots:
    def __init__(self, param, state, varname, slices):
        self.param = param
        self.state = state
        self.varname = varname
        self.slices = slices

    def set_specs(self):
        specs = [npx, npy, nx, ny, nz, i0, j0, k0]
        self.specs = np.asarray(specs, dtype="i")

    def create_file(self):
        self.filename = f"{out_dir}/{expname}/snap_{self.varname}.dat"

        if self.param["myrank"] > 0:
            return

        with open(self.filename, "br") as fid:
            fid.write(self.specs.tobytes())

    def save(self, time):
        var = self.state.get(self.varname)
        z3d = var.view("i")
        z2d = z3d[self.slices]
