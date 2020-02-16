""" 
    join multiple history files
    do it only on 'b' (so far)
"""

import numpy as np
import netCDF4 as nc
import topology as topo

integers = "0123456789"

varname = "b"

def read_param(ncfile):
    param = {}
    with nc.Dataset(ncfile, "r") as fid:
        param_list = fid.ncattrs()
        # print(param_list)
        for p in param_list:
            val = fid.getncattr(p)
            if type(val) is str:
                if val in ["False", "True"]:
                    val = (val == "True")
                elif "class 'list" in val:
                    val = val.split('>:')[-1].strip()
                    val = val[1:-1].split(', ')
                    if val[0][0] in integers:
                        val = [int(e) for e in val if e[0] in integers]
                    elif val[0][0] is "'":
                        val = [e.strip("'") for e in val]
            param[p] = val
    return param


def join(param):
    outfile="%s.nc" % varname
    procs = param["procs"]
    nh = param["nh"]
    nx = param["nx"]
    ny = param["ny"]
    nz = param["nz"]

    topo.topology = param["geometry"]

    global_nx = param["global_nx"]
    global_ny = param["global_ny"]
    global_nz = param["global_nz"]
    size = [nz, ny, nx]

    template = param["expname"]+"_%02i_hist.nc"
    with nc.Dataset(outfile, "w") as fid, nc.Dataset(template % 0, "r") as fin:
        nt = len(fin.dimensions["t"])
        dims = fin.dimensions
        var = fin.variables
        assert varname in var
        
        fid.createDimension("t", nt)
        fid.createDimension("x", global_nx)
        fid.createDimension("y", global_ny)
        fid.createDimension("z", global_nz)

        v = fid.createVariable(varname, "f", ("t", "z", "y", "x"))
        v.long_name = var[varname].long_name
        v.units = var[varname].units
        
        for d in dims:
            v = fid.createVariable(d, "f", (d, ))
            v.long_name = var[d].long_name
            v.units = var[d].units

    nranks = np.prod(procs)
    print("subdomains partition :", procs)
    print("found %i snapshots" % nt)
    print("start joining '%s'" % varname)

    with nc.Dataset(outfile, "r+") as fid:
        with nc.Dataset(template % 0, "r") as fin:
            for kt in range(nt):
                fid["t"][kt] = fin["t"][kt]
                
            for k in range(procs[0]):
                j, i = 0, 0
                loc = [k, j, i]
                ngs = topo.get_neighbours(loc, procs)
                shape, domainindices = topo.get_variable_shape(
                    size, ngs, nh)
                k0, k1, j0, j1, i0, i1 = domainindices
                ka, kb = k*nz, (k+1)*nz                
                fid["z"][ka:kb] = fin["z"][k0:k1]

            for j in range(procs[1]):
                k, i = 0, 0
                loc = [k, j, i]
                ngs = topo.get_neighbours(loc, procs)
                shape, domainindices = topo.get_variable_shape(
                    size, ngs, nh)
                k0, k1, j0, j1, i0, i1 = domainindices
                ja, jb = j*ny, (j+1)*ny
                fid["y"][ja:jb] = fin["y"][j0:j1]

            for i in range(procs[0]):
                k, j = 0, 0
                loc = [k, j, i]
                ngs = topo.get_neighbours(loc, procs)
                shape, domainindices = topo.get_variable_shape(
                    size, ngs, nh)
                k0, k1, j0, j1, i0, i1 = domainindices
                ia, ib = i*nx, (i+1)*nx
                fid["x"][ia:ib] = fin["x"][i0:i1]
                
        for k in range(procs[0]):
            for j in range(procs[1]):
                for i in range(procs[2]):
                    loc = [k, j, i]
                    rank = topo.loc2rank(loc, procs)
                    ncfile = template % rank
                    ngs = topo.get_neighbours(loc, procs)
                    shape, domainindices = topo.get_variable_shape(
                        size, ngs, nh)
                    k0, k1, j0, j1, i0, i1 = domainindices
                    ka, kb = k*nz, (k+1)*nz
                    ja, jb = j*ny, (j+1)*ny
                    ia, ib = i*nx, (i+1)*nx
                    #print(rank, i0, i1, ia, ib, ngs)
                    with nc.Dataset(ncfile, "r") as fin:
                        for kt in range(nt):
                            print("\r %3i/%3i - %3i/%3i"
                                 % (rank, nranks-1, kt, nt-1), end="")
                            var = fin[varname][kt][:, :, :]
                            z2d = var[k0:k1, j0:j1, i0:i1]
                            fid[varname][kt, ka:kb, ja:jb, ia:ib] = z2d
    print()
    print("%s has been joined into '%s'" % (varname, outfile))


if __name__ == '__main__':
    import sys
    import glob

    ncfile = sorted(glob.glob("*hist.nc") )[0]
    varname = sys.argv[-1]

    param = read_param(ncfile)

    msg = "call join.py with one of the following variables: " +  ", ".join(param["variables_in_history"])
    assert varname in param["variables_in_history"], msg
    
    join(param)
