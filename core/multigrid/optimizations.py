from ctypes import POINTER, c_void_p, c_int, c_char, c_double, byref, cdll

mkl = cdll.LoadLibrary("libmkl_rt.so")
SpMV = mkl.mkl_cspblas_dcsrgemv

mkl_matrix_size = 32**3

def MatVecmult(A, x, y):
    (m, n) = A.shape
    if A.nnz > mkl_matrix_size:
        MKL_mult(A, x, y)
    else:
        classical_mult(A, x, y)


def classical_mult(A, x, y):
    """ scipy default sparse matrix-vector multiplication
    """
    y[:] = A.dot(x)


def MKL_mult(A, x, y):
    """ MKL sparse matrix-vector multiplication
    """
    # https://stackoverflow.com/questions/17158893/does-scipy-support-multithreading-for-sparse-matrix-multiplication-when-using-mk#23294826
    (m, n) = A.shape
    # The data of the matrix
    data = A.data.ctypes.data_as(POINTER(c_double))
    indptr = A.indptr.ctypes.data_as(POINTER(c_int))
    indices = A.indices.ctypes.data_as(POINTER(c_int))

    np_x = x.ctypes.data_as(POINTER(c_double))
    np_y = y.ctypes.data_as(POINTER(c_double))
    SpMV(byref(c_char(b"N")), byref(c_int(m)),
         data, indptr, indices, np_x, np_y)

class MatVecObj(object):
    def __init__(self, A, x, y):
        (m, n) = A.shape
        # The data of the matrix
        self.data = A.data.ctypes.data_as(POINTER(c_double))
        self.indptr = A.indptr.ctypes.data_as(POINTER(c_int))
        self.indices = A.indices.ctypes.data_as(POINTER(c_int))

        self.np_x = x.ctypes.data_as(POINTER(c_double))
        self.np_y = y.ctypes.data_as(POINTER(c_double))
        self.m = m

    def do(self):
        SpMV(byref(c_char(b"N")), byref(c_int(self.m)),
             self.data, self.indptr, self.indices, self.np_x, self.np_y)
