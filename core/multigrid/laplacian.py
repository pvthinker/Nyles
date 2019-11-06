import numpy as np
from scipy import sparse


def set_finest(fine):
    """

    Return the Laplacian matrix on the finest grid 
    for the Cartesian case, 7-points uniform Laplacian

    """
    Nf = fine.N
    nmax = Nf*7
    data = np.zeros((nmax,))
    row = np.zeros((nmax,), dtype=int)
    col = np.zeros((nmax,), dtype=int)
    counter = 0

    # todo: set the correct coefficients here
    Az, Ay, Ax = 1., 1., 1.
    
    nk, nj, ni = fine.size
    Gf = np.arange(fine.N).reshape(fine.size)

    k0, k1, j0, j1, i0, i1 = fine.domainindices
    for k in range(k0, k1):
        if k==0:
            dk, ck = [1], [Az]
        elif k==nk-1:
            dk, ck = [-1], [Az]
        else:
            dk, ck = [-1, 1], [Az, Az]
        for j in range(j0, j1):
            if j==0:
                dj, cj = [1], [Ay]
            elif j==nj-1:
                dj, cj = [-1], [Ay]
            else:
                dj, cj = [-1, 1], [Ay, Ay]
            for i in range(i0, i1):
                if i==0:
                    di, ci = [1], [Ax]
                elif i==ni-1:
                    di, ci = [-1], [Ax]
                else:
                    di, ci = [-1, 1], [Ax, Ax]

                Js = []
                coefs = []
                for d in dk:
                    Js += [Gf[k+d, j, i]]
                    coefs += [Az]
                for d in dj:
                    Js += [Gf[k, j+d, i]]
                    coefs += [Ay]
                for d in di:
                    Js += [Gf[k, j, i+d]]
                    coefs += [Ax]
                # diagonal term => Neumann => -sum(extra diag terms)
                Js += [Gf[k, j, i]]
                coefs += [-np.sum(coefs)]
                
                nblock = len(Js)
                Is = [Gf[k, j, i]]*nblock
 
                data[counter:counter+nblock] = coefs
                row[counter:counter+nblock] = Is
                col[counter:counter+nblock] = Js
                counter += nblock
    A = sparse.coo_matrix((data[:counter], (row[:counter], col[:counter])),
                          shape=(Nf, Nf))
    
    return A.tocsr()


