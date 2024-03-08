import build
import numpy as np
import mgmod

class MG:
    def __init__(self, npx, npy, nx, ny, nz, nh, topology=1):
        self.nh = nh

        isfloat32 = False
        vertices = False
        short = False
        is3d = True

        get_mg = mgmod.get("get_ptrmg", isfloat32)
        get_mg.f.restype = build.POINTER(build.c_int64)

        print_mg = mgmod.get("print_mginfos", isfloat32)

        self.libsolve = mgmod.get("solve", isfloat32)
        self.get_pyshape = mgmod.get("get_pyshape", isfloat32)
        self.set_pyarray = mgmod.get("set_pyarray", isfloat32)
        self.get_pyarray = mgmod.get("get_pyarray", isfloat32)

        self.mg = get_mg(npx, npy, nx, ny, nz, vertices, short, is3d, topology)

        print_mg.f(self.mg)

        self.shape = self.get_arrayshape()
        print("shape=", self.shape)
        print("should be=", (nx+2*nh, ny+2*nh, nz+2*nh))

        self.stats = {"normb": 0, "res": [0], "blowup": False}

        #assert np.allclose(self.shape, (nx+2*nh,ny+2*nh,nz+2*nh))

    def get_idx_from_neighbours(self, neighbours):

        nh = self.nh

        i0 = 0 if (0, 0, -1) in neighbours else nh
        i1 = None if (0, 0, 1) in neighbours else -nh

        j0 = 0 if (0, -1, 0,) in neighbours else nh
        j1 = None if (0, 1, 0) in neighbours else -nh

        k0 = 0 if (-1, 0, 0) in neighbours else nh
        k1 = None if (1, 0, 0) in neighbours else -nh
        return (slice(k0, k1), slice(j0, j1), slice(i0, i1))

    def preallocate_for_nyles(self, dx, neighbours, halo):

        self.dx = dx
        self.idx = self.get_idx_from_neighbours(neighbours)

        self.x = np.zeros(self.shape)
        self.b = np.zeros(self.shape)

        self.halo = halo

        n1, n2, n3 = self.shape
        lev = 1

        self.bptr = (self.mg,) + build.ptr((lev, 2, n3, n2, n1, self.b))
        self.xptr = (self.mg,) + build.ptr((lev, 1, n3, n2, n1, self.x))

    def solve_directly(self, p, div):
        nh = self.nh
        #print(self.idx, self.b.shape, div.shape)

        self.halo.fill(div)

        self.b[self.idx] = div

        self.set_pyarray.f(*self.bptr)
        self.libsolve.f(self.mg)
        self.get_pyarray.f(*self.xptr)

        p[:, :, :] = self.x[self.idx]*self.dx**2

    def get_arrayshape(self):
        arrayshape = np.zeros((3,), dtype="i")
        lev = 1
        self.get_pyshape.f(self.mg, build.ptr(lev), build.ptr(arrayshape))
        return arrayshape

    def set_array(self, array, ivar=1):
        n1, n2, n3 = array.shape
        lev = 1
        # print(build.ptr(array))
        pargs = (self.mg,) + build.ptr((lev, ivar, n3, n2, n1, array))
        self.set_pyarray.f(*pargs)

    def get_array(self, array, ivar=1):
        n1, n2, n3 = array.shape
        lev = 1
        # print(build.ptr(array))
        pargs = (self.mg,) + build.ptr((lev, ivar, n3, n2, n1, array))
        self.get_pyarray.f(*pargs)

    def solve(self, x, b):
        #assert np.allclose(self.shape, x.shape)

        self.set_array(b, ivar=2)
        self.libsolve.f(self.mg)
        self.get_array(x, ivar=1)


def testmpi():
    npx = 1
    npy = 1
    nx, ny, nz = 64, 64, 64
    nh = 3

    mg = MG(npx, npy, nx, ny, nz, nh)

    print(f"array shape={mg.shape}")

    b = np.zeros(mg.shape)
    x = np.zeros(mg.shape)

    b[4, 15, 15] = -1
    b[45, 45, 45] = 1

    mg.solve(x, b)


if __name__ == "__main__":

    testmpi()
