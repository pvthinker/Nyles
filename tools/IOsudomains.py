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

        The method also allows for b=v[-1] or b=b[-1,0]
        that return respectively a 3D and a 2D array

        """
        self.nctemplate = nctemplate
        self.varname = varname
        self.debug = debug

        self.check_sanity(varname)

        dimensions, shape, localshape, procs, nprocs = self.get_dimensions(
            nctemplate % 0)

        self.dimensions = dimensions
        self.shape = shape
        self.localshape = localshape
        self.procs = procs
        self.nprocs = nprocs

    def check_sanity(self, varname):
        with Dataset(self.nctemplate % 0) as nc:
            var = nc.variables
        varlist = ""
        for v in var:
            varlist += str(v)+", "
        msg = "varname should be in ["+varlist[:-1]+"]"
        assert varname in var, msg

    def get_dimensions(self, ncfile):
        with Dataset(ncfile) as nc:
            var = nc.variables[self.varname]
            dimensions = var.dimensions

            procs = {}
            nprocs = 1
            for d in "xyz":
                np = nc.getncattr("np"+d)
                procs[d] = np
                nprocs *= np

            localshape = []
            shape = []

            for d in dimensions:
                n = len(nc.dimensions[d])
                localshape += [n]
                if d in "xyz":
                    n *= procs[d]
                shape += [n]
        return dimensions, tuple(shape), tuple(localshape), procs, nprocs

    def __repr__(self):
        msg = ["template    : %s" % self.nctemplate]
        msg += ["variable    : %s" % self.varname]
        msg += ["dimensions  : " + str(self.dimensions)]
        msg += ["global shape: " + str(self.shape)]
        msg += ["partition   : " + str(self.procs)]
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
        dimensions = self.dimensions
        shape = self.shape
        localshape = self.localshape
        procs = self.procs
        nprocs = self.nprocs
        ndims = len(self.dimensions)

        if type(elem) in [int, slice]:
            elem = [elem]
        else:
            elem = list(elem)

        nelem = len(elem)

        # complete missing dimensions with full slices
        # this means that b[-1] stands for b[-1, :, :, :]
        if nelem < ndims:
            elem += [slice(None)]*(ndims-nelem)

        if self.debug:
            print("*"*60)
            print("elem:", elem)

        s = slice(None)
        slices = [s]*ndims
        for k, e in enumerate(elem):
            d = dimensions[k]
            if d == "t":  # time
                if type(e) is slice:
                    # <= time is necessary the first dimension
                    gloshape += [shape[0]]
                else:
                    slices[k] = e
            else:
                if type(e) is slice:
                    gloshape += [shape[k]]
                elif type(e) is int:
                    if self.debug:
                        print("cut in direction %s at %i" % (d, e))
                    slices[k] = e
                    k0 = k
                    e0 = e//shape[k]

        if ndims == 1:
            d = dimensions[0]
            npx = procs["x"]
            npy = procs["y"]
            npz = procs["z"]
            if d == "t":
                procs_to_scan = [0]
            elif d == "x":
                procs_to_scan = range(0, npx)
            elif d == "y":
                procs_to_scan = range(0, npx*npy, npx)
            elif d == "z":
                procs_to_scan = range(0, npx*npy*npz, npx*npy)

        else:
            procs_to_scan = range(nprocs)

        if self.debug:
            print("global returned variable slice is %r" % (gloshape,))
            print("slice is :", slices)

        # allocate the returned array
        data = np.zeros(gloshape)

        # loop trough subdomains
        nfreads = 0
        for rank in procs_to_scan:
            loc = set_loc(rank, procs)
            # is there any data to read from this rank ?
            # yes, if the hyperslab runs through that rank

            skip = False
            gloelem = []
            localelem = []
            for k, e in enumerate(elem):
                d = dimensions[k]
                if d == "t":  # time
                    if type(e) is slice:
                        gloelem += [s]
                        localelem += [s]
                    else:
                        localelem += [e]

                else:
                    i0 = loc[d]*localshape[k]
                    if type(e) is slice:
                        localelem += [e]
                        gloelem += [slice(i0, i0+localshape[k])]
                    elif type(e) is int:
                        if e >= 0:
                            i1 = e-loc[d]*localshape[k]
                        else:
                            i1 = (shape[k]+e)-loc[d]*localshape[k]
                        if (i1 >= 0) and (i1 < localshape[k]):
                            localelem += [i1]
                        else:
                            skip = True

            if skip:
                # no data to read
                pass

            else:
                ncfile = self.nctemplate % rank
                if self.debug:
                    print("gloelem  :", gloelem)
                    print("localelem:", localelem, "rank=", rank)

                # this is where we read the subdomain file and
                # store the result in data, at the good place
                # the good place is controlled by "gloelem"
                with Dataset(ncfile, "r") as nc:
                    data[tuple(gloelem)] = nc[self.varname][tuple(localelem)]
                    nfreads += 1

        if self.debug:
            print("to complete this slice, I read %i netcdf files" % nfreads)

        return data


def set_loc(rank, procs):
    """how subdomains are organized

    loc[d] is the subdomain coordinate in direction d

    """
    loc = {"z": rank // (procs["x"]*procs["y"]),
           "y": (rank // procs["x"]) % procs["y"],
           "x": rank % procs["x"]}
    return loc
