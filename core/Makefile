%.so : %.f90 ; $(F2PY) -m $* -c $< --opt='-O3 -fPIC' 


all:	fortran_interp.so fortran_kinenergy.so fortran_upwind.so fortran_bernoulli.so

