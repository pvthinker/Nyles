import tests as mgmod
import build
import matplotlib.pyplot as plt
import numpy as np
import ctypes

isfloat32 = False

get_mg = mgmod.get("get_ptrmg", isfloat32)
get_mg.f.restype = build.POINTER(build.c_int64)
print_mg = mgmod.get("print_mginfos", isfloat32)

npy, npx = 1, 1
vertices = False
short = False
is3d = True
topology = 1

mg = get_mg(npx, npy, 32, 32, 32, vertices, short, is3d, topology)

print_mg.f(mg)


class MG_Param(ctypes.Structure):
    """ creates a struct to match MG_Param of mg_types.f90 """

    _fields_ = [('npre', ctypes.c_int32),
                ('npost', ctypes.c_int32),
                ('nexact', ctypes.c_int32),
                ('nh', ctypes.c_int32),
                ('maxite', ctypes.c_int32),
                ('tol', ctypes.c_double),
                ('omega', ctypes.c_double),
                ('debug', ctypes.c_int32),
                ('vertices', ctypes.c_bool),
                ('is2d', ctypes.c_bool)]

    def __repr__(self):
        string = []
        for key in self.__dir__():
            if key[0] != "_":
                value = getattr(self, key)
                string += [f"{key} = {value}"]
        return "\n".join(string)

    @property
    def _ptr(self):
        return ctypes.cast(ctypes.pointer(self), ctypes.c_void_p)


p = MG_Param()

get_param = mgmod.get("get_param", isfloat32)

args = (mg, p._ptr)
get_param.f(*args)

print(p)
