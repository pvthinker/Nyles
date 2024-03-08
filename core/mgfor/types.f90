module types
!#define DUAL

#ifdef DUAL
#define REAL type(dual)
#else
#define REAL real
#endif

  type dual
     real:: m,e
  end type dual

#ifdef DUAL
  type(dual)::zero=dual(0.,0.), one=dual(1.,0.)
#else
  real::zero=0., one=1.
#endif


  interface operator(+)
     module procedure add_dual
  end interface operator(+)

  interface operator(-)
     module procedure sub_dual
  end interface operator(-)

  interface operator(*)
     module procedure mult_dual, mult_int, mult_coef
  end interface operator(*)

contains

  elemental function add_dual(x,y) result(z)
    type(dual),intent(in):: x,y
    type(dual):: z

    z%m  = x%m + y%m
    z%e  = x%e + y%e
  end function add_dual

  elemental function sub_dual(x,y) result(z)
    type(dual),intent(in):: x,y
    type(dual):: z

    z%m  = x%m - y%m
    z%e  = x%e - y%e
  end function sub_dual

  elemental function mult_dual(x,y) result(z)
    type(dual),intent(in):: x,y
    type(dual):: z

    z%m  = x%m * y%m
    z%e  = x%m*y%e + y%m*x%e
  end function mult_dual

  elemental function mult_int(x,y) result(z)
    integer,intent(in):: x
    type(dual),intent(in):: y
    type(dual):: z

    z%m  = x * y%m
    z%e  = x * y%e
  end function mult_int

  elemental function mult_coef(x,y) result(z)
    real,intent(in):: x
    type(dual),intent(in):: y
    type(dual):: z

    z%m  = x * y%m
    z%e  = x * y%e
  end function mult_coef

end module types
