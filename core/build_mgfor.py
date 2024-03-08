import build


def build_mgfor():
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
            "pytools.f90",
            "solvers.f90",
            "tuning.f90",
            "tests.f90"]

    build.buildmodule(srcs, path="mgfor", modulename="mgmod")


if __name__ == "__main__":
    build_mgfor()
