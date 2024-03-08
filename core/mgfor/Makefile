FC := mpif90.mpich
FFLAGS := -cpp -Ofast -freal-4-real-8

SRCS := types.f90 \
       mg_enums.f90 \
       mod_io.f90 \
       mod_tiles.f90 \
       mod_halo.f90 \
       mod_gluesplit.f90 \
       mg_types.f90 \
       mg_log.f90 \
       basicoperators.f90 \
       operators.f90 \
       mg_setup.f90 \
       solvers.f90 \
       tuning.f90 \
       tests.f90 \
       demo.f90

OBJS := $(SRCS:%:=%.o)

%.o: %.f90
	$(FC) -c $(FFLAGS) $(CPPFLAGS) $< -o $@


all: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o demo
