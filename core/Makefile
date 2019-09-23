%.so : %.f90
	$(F2PY) -m $* -c $< --opt='-O3 -fPIC'
	echo $*


SRCS = $(wildcard *.f90)
OBJS = $(patsubst %.f90, %.so, $(SRCS))

all:	$(OBJS)

