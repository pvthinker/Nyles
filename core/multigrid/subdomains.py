import numpy as np
import topology as topo
import DOESNTEXISTMODULE  # to force an error, subdomains has been superseded


def set_subdomains(allgrids):
    grid = allgrids[0]
    procs = grid['procs'].copy()
    procs0 = procs.copy()
    incr = [1, 1, 1]
    nprocs = np.prod(procs)
    ranks = np.arange(nprocs)
    loc = topo.rank2loc(ranks, procs)
    loc0 = loc.copy()
    k0, j0, i0 = loc0
    gathers = []
    doms = []
    first = True
    tags = ranks*0
    subdom_partition = []

    for lev, grid in enumerate(allgrids):
        c = grid['coarsen']
        keys = c.keys()
        if any(['g' in k for k in keys]) or lev == 0:
            g = []
            levs = [lev]
            for n, direc in enumerate('kji'):
                if 'g'+direc in keys:
                    factor = 2*c['g'+direc]
                    incr[n] *= factor
                    procs[n] //= factor
                    loc[n] //= factor
                    g += [direc]
            if g != '':
                gathers += [g]

            else:
                pass
            k, j, i = loc

            nz, ny, nx = procs0
            incz, incy, incx = incr
            dom = [0]*nprocs
            for r in ranks:
                l = [k[r], j[r], i[r]]
                # d is the index of the subdomain
                # rank with same d do the same computation
                d = topo.myloc(l, incr, [0, 0, 0], procs0)
                dom[r] = d

            doms += [dom]
            dd = set(dom)
            neighbours = [{}]*len(ranks)
            domloc = {}
            glued = {}
            for d in dd:
                family = [r for r in ranks if dom[r] == d]
                # in a family, all ranks have the same loc
                # because they all compute the same subdomain
                locs = [(k0[r], j0[r], i0[r]) for r in family]

                # tags is defined from the previous domain decomposition
                fam = family[0]
                glued[d] = [[r for r in family if tags[r] == t]
                            for t in set(tags[family])]
                if len(glued[d]) == 1:
                    glued[d] = glued[d][0]

                # print(' - subdomain : ', set([d]))
                # print('     located at ', locs[0])
                # print('     computed by ' % d, family)
                # if not(first):
                #     print('     obtained by gluing ranks:', glued)
                # here is the list of neighbours
                #
                # the print has been discarded to lighten the output
                # but clearly this information is very important
                # note that the neighbours relationships for a given
                # rank depends on the subdomain decomposition
                # the neighbours get further and further as
                # subdomains are glued. This is accounted by
                # 'incr'.  incr == 1 on the finest grid, then
                # whenever there is gluing in one direction
                # incr is doubled in that direction, until
                # the point where the is only one subdomain is that
                # direction. Each rank is then neighbour with itself
                #
                # print('     neighbours:')
                for r in family:
                    l = [k0[r], j0[r], i0[r]]
                    nghbs = topo.get_neighbours(
                        l, procs0, incr=incr, extension=26)
                    ngbs = fixneighbours(nghbs, r, d)
                    # print('rank : %i / ' % r, ngbs)
                    neighbours[r] = ngbs

                r0 = d  # family[0]
                domloc[d] = [loc0[0][r0], loc0[1][r0], loc0[2][r0]]
#                print('domain: %i' % d, domloc[d])
            localpartition = {'procs': procs.copy(),
                              'incr': incr.copy(),
                              'gluing': g.copy(),
                              'glued': glued.copy(),
                              'domloc': domloc.copy(),
                              'dom': dom.copy(),
                              'levels': levs,
                              'allneighbours': neighbours.copy()}
#            print(levs, domloc)
            #            if not(first):
            localpartition['tags'] = tags.copy()
            subdom_partition += [localpartition]

            #
            # tags is used to identify which ranks are glued together
            #
            # glued ranks had identical tag in the previous
            # subdomain decomposition
            #
            for d in dd:
                family = [r for r in ranks if dom[r] == d]
                for k, r in enumerate(family):
                    tags[r] = k
            first = False
        else:
            levs += [lev]

    return subdom_partition


def attach_subdomain_to_grids(allgrids, subdomains, myrank):
    for isub, subd in enumerate(subdomains):
        ngbs = subd['allneighbours'][myrank]
        for count, lev in enumerate(subd['levels']):
            g = allgrids[lev]
            shape = g['shape']
            nh = g['nh']
            size, domainindices = topo.get_variable_shape(shape, ngbs, nh)
            g['subd'] = subd
            g['subdomain'] = isub
            g['neighbours'] = ngbs
            g['N'] = np.prod(size)
            g['size'] = size
            g['domainindices'] = domainindices
            if (count == 0) and (len(subd["gluing"]) > 0):
                g["glueflag"] = True
                g["gluing"] = subd["gluing"]
                dom = subd["dom"]
                g["glued"] = subd["glued"]
            else:
                g["glueflag"] = False


def fixneighbours(neighbours, myrank, masterrank):
    """

    for subdomains computed by several ranks
    the neighbours are correct only for the master rank
    (smallest rank of the family)
    for the others, the neighbours are recovered by
    shifting the value by (myrank-masterrank)

    this is what is doing this fix

    """
    for key, val in neighbours.items():
        neighbours[key] += (myrank - masterrank)
    # note that a dictionnary is mutable so there is
    # need to return something, neighbours is change in-place
    return neighbours


def print_subdomains(subdomains):
    print('-'*80)
    for isub, subdom in enumerate(subdomains):
        dom = subdom['dom']
        procs = subdom['procs']
        incr = subdom['incr']
        tags = subdom['tags']
        levs = subdom['levels']
        glued = subdom['glued']
        domloc = subdom['domloc']
        ranks = np.arange(len(dom))
        k, j, i = topo.rank2loc(ranks, procs)
        if isub == 0:
            k0 = k.copy()
            j0 = j.copy()
            i0 = i.copy()
            procs0 = procs.copy()

        print(' Partition %1i' % isub)
        print('=============')

        print('Subdomains are : ', dom)
        print('corresponding to domain decomposition: ', procs)
        print('increments in each directions: ', incr)
        print('grid levels having this partition: ', levs)
        print('Location of each subdomain')
        print('k: ', k)
        print('j: ', j)
        print('i: ', i)
#        for idx in gr.ranktoloc(ranks, procs):
#            print('x: ', idx)
        if 'gluing' in subdom.keys():
            print('gluing was in directions ', subdom['gluing'])
        for d in set(dom):
            family = [r for r in ranks if dom[r] == d]
            # in a family, all ranks have the same loc
            # because they all compute the same subdomain
            locs = [(k0[r], j0[r], i0[r]) for r in family]

            # tags is defined from the previous domain decomposition
            # glued = [[r for r in family if tags[r] == t]
            #          for t in set(tags[family])]

            print(' - subdomain : ', set([d]))
            print('     located at ', domloc[d])  # s[0])
            print('     computed by ' % d, family)
            if 'gluing' in subdom.keys():
                print('     obtained by gluing ranks:', glued)
            print('     neighbours:')
            for r in family:
                print('rank : %i / ' % r, subdom['allneighbours'][r])


if __name__ == '__main__':

    import grids as grd

    #procs = [4, 8, 4]
    procs = [2, 4, 2]
    topo.topology = 'closed'
    nh = 2
    myrank = 3
    param = {'nx': 128, 'ny': 128, 'nz': 256, 'nh': nh, 'procs': procs,
             "nglue": 32, "ncellscoarsest": 16}

    allgrids = grd.define_grids(param)

    subdomains = set_subdomains(allgrids)

    for s in subdomains:
        print(s["glued"])
        print("*"*80)

    print_subdomains(subdomains)
    exit(0)
    attach_subdomain_to_grids(allgrids, subdomains, myrank)
    #print(isolate(subdomains, 3, 3))

    # print(len(subdomains))
    for lev, g in enumerate(allgrids):
        isubd = g["subdomain"]
        subd = subdomains[isubd]
        levs = subd["levels"]
        gl = subd["gluing"]
        glueflag = g["glueflag"]
        if glueflag:
            print(lev, isubd, subd["procs"], levs, glueflag, gl, g["glued"])
        else:
            print(lev, isubd, subd["procs"], levs, glueflag)
