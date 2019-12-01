import numpy as np
import grids as grd
import intergrids as intergrd
import mpitools as mpi
import subdomains as subdom
import laplacian as laplac


class Multigrid(object):
    def __init__(self, param, modelgrid):

        for key in ['npre', 'npost', 'maxite', 'tol', 'procs']:
            setattr(self, key, param[key])
        self.verbose = False

        myrank = mpi.get_myrank(self.procs)

        # 1/ determine the hierarchy: grids and subdomains partitions
        allgrids = grd.define_grids(param)
        if myrank == 0 and self.verbose:
            grd.print_grids(allgrids)

        nlevs = len(allgrids)
        self.nlevs = nlevs
        subdomains = subdom.set_subdomains(allgrids)
        subdom.attach_subdomain_to_grids(allgrids, subdomains, myrank)

        # 2/ grid instances: contains the arrays x, b and r, the halofill and the information to glue
        self.grid = [1] * nlevs
        for lev, g in enumerate(allgrids):
            self.grid[lev] = grd.Grid(g, param)

        # 3/ intergrid functions
        self.inter = [1] * (nlevs-1)
        for lev in range(nlevs-1):
            coarse = self.grid[lev+1]
            fine = self.grid[lev]
            self.inter[lev] = intergrd.Intergrids(fine, coarse)

        # 4/ sparse matrices: laplacian and smoothing
        self.As = [1] * nlevs
        for lev, g in enumerate(self.grid):
            if lev == 0:
                Afinest = laplac.set_finest(g, modelgrid)
                g.set_ADS(Afinest)
            else:
                Afine = self.grid[lev-1].A
                intergrid = self.inter[lev-1]
                Acoarse = intergrid.Restrict*Afine*intergrid.Interpol
                g.set_ADS(Acoarse)

        # store the statistics of the last 'solve_directly'
        self.stats = {}
        
    def solve(self, x, b, cycle='V'):
        """

        driver to solve A*x=b using the Vcycle
        x is the first guess, b the right hand side

        """
        g = self.grid[0]

        # Copy x (first guess) and b (rhs) to solve directly on local attributes
        g.x[:] = x
        g.b[:] = b

        nite, res = self.solve_directly(cycle='V')

        # Copy x (solution) back to the variable which is held by the caller
        x[:] = g.x

        return nite, res

    def solve_directly(self, cycle='V'):
        """

        driver to solve A*x=b using the Vcycle

        the first guess and the rhs have been assigned outside of the routine

        likewise, the solution is to be copied from outside

        """
        g = self.grid[0]

        # No need to copy x and b, since they have already been assigned before

        g.residual()

        normb = g.norm(which='b')
        if self.verbose:
            print('||b|| = ', normb)

        if normb > 0:
            res0 = g.norm()/normb
        else:
            return 0, 0.

        if normb > 1e6:
            raise ValueError('blowup')

        res = res0
        reslist = [res]
        nite = 0
        nite_diverge = 0
        # improve the solution until one of this condition is wrong
        ok = True
        while (nite < self.maxite) and (res0 > self.tol) and ok:
            if cycle == 'V':
                self.Vcycle(0)

            elif cycle == 'F':
                self.Fcycle()

            else:
                raise ValueError('use cycle V or F')

            g.residual()

            res = g.norm() / normb
            conv = res0 / res

            res0 = res
            reslist += [res]
            nite += 1
            if self.verbose:
                template = ' ite = {} / res = {:.2e} / conv = {:8.4f}'
                print(template.format(nite, res, conv))

            if (conv < 1):
                nite_diverge += 1

            if (nite_diverge > 4):
                ok = False
                print('solver is not converging')
                print('Abort!')
                #raise ValueError('solver is not converging')

        # No need to copy x back to an external variable

        # store the statistics
        self.stats = {'normb': normb, 'res': reslist}

        return nite, res

    def Vcycle(self, lev0):
        for lev in range(lev0, self.nlevs-1):
            g = self.grid[lev]
            g.smooth(self.npre)
            g.residual()
            self.inter[lev].fine2coarse()

        self.grid[-1].solveexactly()

        for lev in range(self.nlevs-2, lev0-1, -1):
            self.inter[lev].coarse2fine()
            self.grid[lev].smooth(self.npost)

    def Fcycle(self):
        for lev in range(self.nlevs-1):
            self.inter[lev].fine2coarse(which='b')

        self.grid[-1].solveexactly()

        for lev in range(self.nlevs-2, -1, -1):
            self.inter[lev].coarse2fine()
            self.Vcycle(lev)


if __name__ == '__main__':

    import topology as topo
    import matplotlib.pyplot as plt
    import test_mg as tmg

    procs = [1, 1, 1]
    topology = 'perio_xy'
    topo.topology = topology
    nh = 1

    param = {'nx': 128, 'ny': 32, 'nz': 32, 'nh': nh, 'procs': procs, 'topology': topology,
             'npre': 3, 'npost': 3, 'omega': 0.8, 'ndeepest': 20, 'maxite': 10, 'tol': 1e-6}

    mg = Multigrid(param)

    lev = 0
    g = mg.grid[lev]
    #tmg.test_smooth(mg, 0, nite=3)
    y = np.reshape(g.x, g.size)
    plt.imshow(y[1, :, :])

    lev = 1
    g = mg.grid[lev]
    #tmg.test_smooth(mg, lev, nite=4)
    y = np.reshape(g.x, g.size)
    plt.imshow(y[1, :, :])

    g = mg.grid[0]
    b = g.b
    bb = np.reshape(b, g.size)
    bb *= 0
    bb[1, 1, 1] = 1
    bb[-2, -2, -2] = -1
    #tmg.substractmean(mg, 0, b)
    x = np.zeros_like(b)
    mg.solve(x, b, cycle='F')
    y = np.reshape(g.x, g.size)
    plt.imshow(y[1, :, :])
