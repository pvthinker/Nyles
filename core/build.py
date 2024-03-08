"""
Why use Makefiles when Python does it so easily?

Compile the Fortran subroutines and generate the dynamic libraries

"""
import os
import time
import numpy as np
import subprocess
from ctypes import POINTER, c_char, c_double, c_float, c_int64, c_int32, c_int8, byref, cdll, CDLL, c_int, c_bool, c_wchar_p, create_string_buffer, c_char_p


def bold(string):
    return f"\033[1m{string}\033[0m"


MAP_TYPES = {
    "float64": c_double,
    "float32": c_float,
    "int8": c_int8,
    "int32": c_int32,
    "int64": c_int64,
    "bool": c_bool
}


def localpath(pathtofile):
    return os.path.abspath(os.path.dirname(pathtofile))


def setup_fortran_func(fortranfunc, common_args):
    tags = {}
    com_args = ptr(common_args)

    def func(*args, tag=None):
        if tag not in tags:
            tags[tag] = ptr(args)
        local_args = tags[tag]
        fortranfunc(*local_args, *com_args)
    return func


def ptr(thing):
    """ return pointer(s)"""
    if isinstance(thing, tuple) or isinstance(thing, list):
        args = [ptr(t) for t in thing]
        return tuple(args)

    elif isinstance(thing, np.ndarray):
        # return np.ctypeslib.as_ctypes(thing)
        return thing.ctypes.data_as(POINTER(MAP_TYPES[str(thing.dtype)]))

    elif isinstance(thing, int):
        return byref(c_int32(thing))

    elif isinstance(thing, np.int32):
        return byref(c_int32(thing))

    elif isinstance(thing, np.int64):
        return byref(c_int64(thing))

    elif isinstance(thing, np.float32):
        return byref(c_float(thing))

    elif isinstance(thing, float):
        return byref(c_double(thing))

    elif isinstance(thing, bool):
        return byref(c_bool(thing))

    elif isinstance(thing, bytes):
        return byref(c_char(thing))

    elif isinstance(thing, str):
        # return byref(c_wchar_p(thing))
        # return byref(c_char_p(bytes(thing, "utf-8")))
        return byref(create_string_buffer(bytes(thing, "utf-8")))

    else:
        raise ValueError(
            f"problem with argument {thing} of type {type(thing)}")


def isfile1newerthanfile2(file1, file2, path):
    t1 = os.path.getmtime(f"{path}/{file1}")
    if os.path.isfile(f"{path}/{file2}"):
        t2 = os.path.getmtime(f"{path}/{file2}")
        return t1 > t2
    else:
        return True


def srcs_newer_than_lib(srcs, lib, path):
    if isinstance(srcs, list):
        return any([isfile1newerthanfile2(src, lib, path) for src in srcs])
    elif isinstance(srcs, str):
        return isfile1newerthanfile2(srcs, lib, path)
    else:
        raise ValueError


def load(library):
    return CDLL(library)


def import_from_library(library, function, path="."):
    lines = subprocess.check_output(
        ["nm", f"{path}/{library}"]).decode().split("\n")
    line = [l for l in lines if "MOD_"+function in l]
    if len(line) == 1:
        name = line[0].split()[-1]
        lib = load(f"{path}/{library}")
        return getattr(lib, name)

    elif len(line) > 1:
        msg = f"found multiple occurences of {function} in {library}"
        msg += "\n"+str(line)
    else:
        msg = f"did not found {function} in {library}"
    raise ValueError(msg)


def make_executable(srcs, binary, isfloat32=False):
    compiler, flags = get_compiler(isfloat32)
    objs = " ".join(["".join(src.rsplit(".")[:-1]+[".o"])
                     for src in srcs])
    #os.system("rm *.mod")
    for src in srcs:
        command = f"{compiler} {flags} -c {src}"
        print(command)
        os.system(command)

    command = f"{compiler} {flags} {objs} -o {binary}"
    print(command)
    os.system(command)


def compile_fortran(compiler, flags, srcs, library, path):
    if isinstance(srcs, list):
        srcs = " ".join(srcs)

    #currentpath = os.path.abspath(os.path.curdir)
    command = f"cd {path};{compiler} {flags} -fPIC -shared  {srcs} -o {library}"
    tic = time.time()
    os.system(command)
    toc = time.time()
    elapsed = toc-tic
    #print(f"{elapsed:8.3}s : {command}")
    return elapsed, command


def get_compiler(isfloat32, ismpi=True):
    nodename = os.uname().nodename
    if "irene" in nodename:
        compiler = "ifort"
        flags = "-cpp -O3 -mavx2 -axCORE-AVX512,MIC-AVX512"
        if not isfloat32:
            flags += " -r8"

    else:
        if ismpi:
            compiler = "mpif90.mpich"
        else:
            compiler = "gfortran"

        flags = "-cpp -Ofast -march=x86-64 -mtune=native -Wno-overflow -ffree-line-length-none -Wno-all"
        flags = "-cpp -Ofast"
        #flags = "-cpp -g"
        #flags += " -floop-nest-optimize -floop-parallelize-all -fdump-fortran-original"
        #flags += " -Wvector-operation-performance"
        #flags = "-cpp  -Ofast -march=native -ffast-math  -ffree-line-length-none"
        if not isfloat32:
            flags += " -freal-4-real-8"

    return compiler, flags


def build(modules, path=None, isfloat32=False, skip=True):
    """Build dynamic libraries from Fortran files

    Parameters
    ----------

    modules: dict

    whose keys are the library name (without lib and .so), and values
    are, either the name of the Fortran filename, or the list of
    Fortran filenames

    """
    assert isinstance(modules, dict)

    assert path is not None, "path is required when calling build"

    compiler, flags = get_compiler(isfloat32)

    for library, srcs in modules.items():
        if srcs_newer_than_lib(srcs, library, path) or not skip:
            elapsed, command = compile_fortran(
                compiler, flags, srcs, library, path)
            print(bold(f"{library:>20}: [compiled in {elapsed:4.2} s]"))
            for line in command.split(";"):
                print(" "*4, line)
        else:
            print(bold(f"{library:>20}: [ok]"))


def buildmodule(fortranfile, path=".", skip=False, modulename=None):
    """ transform a fortran script into a python module"""

    if isinstance(fortranfile, str):
        functions = retrieve_functions(fortranfile, path=path)

        modulename = os.path.splitext(fortranfile)[0]

    elif isinstance(fortranfile, list):
        functions = {}
        for f in fortranfile:
            functions.update(retrieve_functions(f, path=path))

        assert isinstance(modulename, str), "module name should be specified"

    else:
        assert False

    #lib = f"lib{modulename}.so"
    for isfloat32 in [False]:
        lib = get_libname(f"lib{modulename}.so", isfloat32)
        build({lib: fortranfile}, path=path, isfloat32=isfloat32, skip=skip)

    write_module(modulename, functions, path)
    print(f"create Python module: {bold(modulename)}")
    print("that contains:")
    for name, signature in functions.items():
        print(f"  - {name}: {signature}")


def retrieve_functions(fortranfile, path="."):
    with open(path+"/"+fortranfile) as fid:
        fortran_code = fid.readlines()

    func_lines = [" ".join(l.strip().split()[1:])
                  for l in fortran_code
                  if l.strip().startswith("function")
                  or l.strip().startswith("subroutine")]
    functions = {}
    for func in func_lines:
        istart = func.find("(")
        name = func[:istart]
        signature = func[istart:]
        functions[name] = signature
    return functions


def get_libname(lib, isfloat32):
    digit = "32" if isfloat32 else "64"
    return lib[:-3]+digit+".so"


class Function:
    def __init__(self, lib, name, signature, isfloat32=False,path="."):
        self.name = name
        self.signature = signature
        self.lib = get_libname(lib, isfloat32)
        self.isfloat32 = isfloat32
        self.f = import_from_library(self.lib, name,path=path)
        self.f.restype = c_float if isfloat32 else c_double
        self.tags = {}

    def __repr__(self):
        dtype = "float32" if self.isfloat32 else "float64"
        return f"<Fortran function> [from {self.lib}] {self.name}{self.signature} "

    def __call__(self, *args, tag=None):
        if tag is None:
            return self.f(*ptr(args))
        else:
            if tag in self.tags:
                a = self.tags[tag]
            else:
                a = ptr(args)
                self.tags[tag] = a
            return self.f(*a)


template = """from build import Function

functions = {strfunctions}

def get(name, isfloat32=False):
    if name not in functions:
        print(f"{{name}} not in {{list(functions.keys())}}")
        return None
    signature = functions[name]
    path = '{path}'
    return Function('{lib}', name, signature, isfloat32, path=path)

"""


def write_module(modulename, functions, path):
    absolute_path = os.getcwd()+"/"+path
    lib = f"lib{modulename}.so"
    strfunctions = ",\n".join(
        [f"   '{name}': '{signature}'"
         for name, signature in functions.items()])
    strfunctions = "{\n"+strfunctions+"}\n"
    infos = {"lib": lib,
             "strfunctions": strfunctions,
             "path": absolute_path}
    with open(f"{modulename}.py", "w") as fid:
        for line in template.split("\n"):
            fid.write(line.format(**infos)+"\n")


# def reload(module):
#     pythonname = os.path.basename(module.__file__)
#     name = os.path.splitext(pythonname)[0]
#     lib = CDLL(f"./lib{name}.so")
#     h = lib._handle
#     del lib
#     dlclose = CDLL(None).dlclose
#     dlclose(h)

def maketestgluesplit(isfloat32):
    srcs = ["types.f90",
            "mg_enums.f90",
            "mod_tiles.f90",
            "mod_io.f90",
            "mod_gluesplit.f90",
            "testgluesplit.f90"
            ]
    make_executable(srcs, "testgluesplit", isfloat32)


def makedemo(isfloat32):
    srcs = ["types.f90",
            "mg_enums.f90",
            "mod_io.f90",
            "mod_tiles.f90",
            "mod_halo.f90",
            "mod_gluesplit.f90",
            "mg_types.f90",
            "mg_log.f90",
            "basicoperators.f90",
            "operators.f90",
            "mg_setup.f90",
            "solvers.f90",
            "tuning.f90",
            "tests.f90",
            "demo.f90",]
    make_executable(srcs, "demo", isfloat32)


def makedemohalo(isfloat32):
    srcs = ["types.f90",
            "mg_enums.f90",
            "mg_types.f90",
            "mg_log.f90",
            "basicoperators.f90",
            "operators.f90",
            "mg_setup.f90",
            "mod_tiles.f90",
            "mod_halo.f90",
            "demo_halo.f90"]
    make_executable(srcs, "demo_halo", isfloat32)


if __name__ == "__main__":
    import sys
    argv = sys.argv[1:]

    opts = [a for a in argv
            if a[0] == "-"]

    isfloat32 = "-32" in opts
    print(isfloat32, opts)
    makedemo(isfloat32)
    # maketestgluesplit(isfloat32)
    # makedemohalo(isfloat32)
