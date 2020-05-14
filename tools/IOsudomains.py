"""
Tools to read variables from native, unjoined, history files
"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import os

class Variable(object):
    def __init__(self, nctemplate, varname, debug=False):
        """

        the name of the template should be passed here,

        template is a string such that ncfile = template % rank,
        returns the history file of rank

        varname is the variable to be read

        to read the buoyancy "b", first create the
        Variable object, then read any slice like
        you would do with a single file

        for instance, if the one netCDF file case is
        >>> with Dataset(ncfile) as nc:
        >>>     b = nc.variables["b"][:, 4, :, -1]

        the multi subdomains case will be
        >>> v = Variable(nctemplate,"b")
        >>> b = v[:, 4, :, -1]

        """
        nt, shape, procs = self.get_infos(nctemplate % 0)
        self.nt = nt
        self.shape = shape # space dimension
        self.procs = procs
        self.varname = varname
        self.nctemplate = nctemplate
        self.debug = debug

    def get_infos(self, ncfile):
        paramlist = ["nx", "ny", "nz",
                     "npx", "npy", "npz", "include_halo"]
        param = {}
        with Dataset(ncfile) as nc:            
            for p in paramlist:
                param[p] = nc.getncattr(p)
            nt = len(nc.dimensions["t"])
        print("include_halo: ", param["include_halo"])
        shape = (param["nz"], param["ny"], param["nx"])
        procs = (param["npz"], param["npy"], param["npx"])
        print("shape of local subdomain: ", shape)
        print("number of frames        : ", nt)
        print("domain partition        : ", procs)
        return nt, shape, procs

    def __repr__(self):
        msg = ["subdomain dimensions: %r" % (self.dimensions,)]
        msg += ["variable name: %s" % self.name]
        return "\n".join(msg)

    def __getitem__(self, elem):
        """ __getitem__ is a special Python method that allows to access an
object with bracket and indexing. See

        https://docs.python.org/2.7/reference/datamodel.html#special-method-names

        and

        https://docs.scipy.org/doc/numpy/reference/arrays.indexing.html#advanced-indexing

        The 'elem' argument is directly inspired from class MFDataset(Dataset) found in

        https://github.com/Unidata/netcdf4-python/blob/master/netCDF4/_netCDF4.pyx

        """
        gloshape = []
        shape = self.shape
        procs = self.procs
        nprocs = np.prod(procs)

        if self.debug:
            print("*"*60)
            print("elem:", elem)

        s = slice(None)
        slices = [s, s, s, s]
        for k, e in enumerate(elem):
            if k == 0: # time
                if type(e) is slice:
                    gloshape += [self.nt]
                else:
                    slices[k] = e
            else:
                kk = k-1
                if type(e) is slice:
                    gloshape += [shape[kk]*procs[kk]]
                elif type(e) is int:
                    if self.debug:
                        print("cut in direction %i at %i" % (kk, e))
                    slices[k] = e
                    k0 = kk
                    e0 = e//shape[kk]

        if self.debug:
            print("global returned variable slice is %r" % (gloshape,))
            print("slice is :", slices)

        # allocate the returned array
        data = np.zeros(gloshape)

        # loop trough all ranks
        nfreads = 0
        for rank in range(nprocs):
            loc = set_loc(rank, procs)

            # is there any data to read from this rank ?
            # yes, if the hyperslab runs through that rank
            if True:
                time = slices[0]
                localelem = [time]
                if type(time) is slice:
                    gloelem = [time]
                else:
                    gloelem = []
                skip = False
                for k, e in enumerate(elem[1:]):
                    i0 = loc[k]*shape[k]
                    if type(e) is slice:
                        localelem += [e]
                        gloelem += [slice(i0, i0+shape[k])]
                    elif type(e) is int:
                        if e >= 0:
                            i1 = e-loc[k]*shape[k]
                        else:
                            i1 = (shape[k]*procs[k]+e)-loc[k]*shape[k]
                        if (i1 >= 0) and (i1 < shape[k]):
                            localelem += [i1]
                        else:
                            skip = True

                #print("rank %i / local elem is %r" % (rank, gloelem,))
                if skip:
                    pass
                else:
                    ncfile = self.nctemplate % rank
                    if self.debug:
                        print("gloelem  :", gloelem)
                        print("localelem:", localelem, "rank=",rank)
                    with Dataset(ncfile, "r") as nc:
                        b = nc["b"][tuple(localelem)]
                        data[tuple(gloelem)] = b
                        nfreads += 1
            else:
                # no data to read
                pass

        if self.debug:
            print("to complete this slice, I read %i netcdf files" % nfreads)

        return data

def set_rank(loc):
    """
    convert subdomain location (k,j,i) into rank
    """
    k, j, i = loc
    rank = i+procs[0]*(j+procs[1]*k)
    return rank


def set_loc(rank, procs):
    loc = [rank // (procs[2]*procs[1]),
           (rank // procs[2]) % procs[1],
           rank % procs[2]]
    return loc
