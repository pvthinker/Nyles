import numpy as np
import mpitools as mpi


class BackgroundProfile:
    def __init__(self, zmin, zmax, nz, nh=0, nbins=500):
        self.zmin = zmin
        self.zmax = zmax
        self.nbins = nbins
        self.nz = nz
        self.nh = nh
        if nh == 0:
            self.interior = tuple([slice(None), slice(None), slice(None)])
        else:
            self.interior = tuple(
                [slice(None), slice(nh, -nh), slice(nh, -nh)])

        #self.z1d = zmin+(np.arange(nz+2*nh)+0.5-nh)*(zmax-zmin)/nz
        self.z1d = np.zeros((nz+2*nh,))
        self.z1d[nh:-nh] = zmin+(np.arange(nz)+0.5)*(zmax-zmin)/nz

    def compute(self, b):
        self.bmin = mpi.Min(np.min(b))
        self.bmax = mpi.Max(np.max(b))
        self.deltab = (1+1e-6)*(self.bmax-self.bmin)
        deltaz = self.zmax-self.zmin

        bins = np.linspace(self.bmin, self.bmax, self.nbins+1)
        h, bins = np.histogram(b[self.interior], bins)

        self.dz = h*deltaz/sum(h)
        self.bins = bins
        self.db = self.deltab/self.nbins
        self.b = 0.5*(bins[1:]+bins[:-1])
        self.z = np.zeros((self.nbins+1,))

        self.z[1:] = np.cumsum(self.dz)
        self.z += self.zmin

        self.b1d = np.interp(self.z1d, self.z, self.bins)

    def zr(self, binput):
        b = binput[self.interior].copy()
        N = b.size
        z = np.zeros((N,))
        bb = b.ravel()
        #bb = np.round(bb, decimals=3)
        idx = np.argsort(bb, kind="stable")
        z[idx] = np.linspace(self.zmin, self.zmax, N)
        z.shape = b.shape

        zoutput = np.zeros_like(binput)
        zoutput[self.interior] = z
        return zoutput

    def zr_v0(self, b):
        idx, frac = np.divmod(self.nbins*(b-self.bmin)/self.deltab, 1)
        idx = np.asarray(idx, dtype="i")
        return (self.z[idx]+frac*self.dz[idx])

    def dzrdb(self, b):
        idx = np.floor(self.nbins*(b-self.bmin)/self.deltab)
        idx = np.asarray(idx, dtype="i")
        return self.dz[idx]/self.db.data

    def Deltaz(self, b):
        return self.zr(b)-self.z1d[np.newaxis, :, np.newaxis]

    def Deltab(self, b):
        return b-self.b1d[np.newaxis, :, np.newaxis]

    def rhs(self, b):
        return 0.5*(self.Deltab(b)*self.dzrdb(b)+self.Deltaz(b))
