import numpy as np

def substractmean(mg, lev, b):
    g = mg.grid[lev]
    k0, k1, j0, j1, i0, i1 = g.domainindices
    bb = np.reshape(b, g.size)
    mb = np.mean(bb[k0:k1, j0:j1, i0:i1])
    b -= mb


def test_smooth(mg, lev, nite=10):
    b = np.reshape(mg.grid[lev].b, mg.grid[lev].size)
    x = mg.grid[lev].x
    x *= 0
    b *= 0
    b[1, 1, 1] = 1
    b[-2, -2, 2] = -1
    #substractmean(mg, lev, b)
    mg.grid[lev].halofill('b')
    for k in range(nite):
        mg.funcs[lev].smooth(1)
        mg.funcs[lev].residual()
        res = mg.funcs[lev].norm()
        print(lev, k, res)


def test_Vcycle(mg, lev):
    b = mg.grid[lev].b
    b[1, 1, 1] = 1
    b[-2, -2, 2] = -1
    b -= np.mean(b)

    mg.grid[lev].x[:] = 0.
    for k in range(10):
        mg.Vcycle(lev)
        mg.funcs[lev].residual()
        res = mg.funcs[lev].norm()
        print(lev, k, res)


def test_twolevels(mg, lev0):
    mg.grid[lev0].b[500] = 1.
    mg.grid[lev0].b[2500] = -1.
    substractmean(mg, lev0, b)
    mg.grid[lev0].halofill('b')
    mg.grid[lev0].x[:] = 0.
    for k in range(10):
        for lev in range(lev0, lev0+1):
            mg.funcs[lev].smooth(mg.npre)
            mg.funcs[lev].residual()
            mg.inter[lev].fine2coarse()

        mg.funcs[lev0+1].smooth(20)

        for lev in range(lev0, lev0-1, -1):
            mg.inter[lev].coarse2fine()
            mg.funcs[lev].smooth(mg.npost)

        mg.funcs[lev0].residual()
        res = mg.funcs[lev0].norm()
        print(lev0, k, res)
