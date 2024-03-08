module solvers
  use operators
  implicit none

  integer:: ncalls=0
contains

  subroutine solve(mg)
    type(MGinfos) :: mg
    integer:: nite
    REAL:: normb, res
    logical:: verbose

    verbose = ((mg%verbose>0).and.(mod(ncalls,20)==0))

    nite = 0
    normb = norm(mg, 1, bnorm)
    if(verbose) print '("||b||= ",(e9.3))',normb
    res = normresidual(mg, normb)
    if(verbose) print '("  ",(i2), ". ||res||= ",(e9.3))',nite,res
    do
       if(res<mg%param%tol) exit
       call vcycle(mg, 1)
       !call smooth(mg, 1, 5000)
       nite = nite+1
       if(nite>=mg%param%maxite) exit
       res = normresidual(mg, normb)
       if(verbose)print '("  ",(i2), ". ||res||= ",(e9.3))',nite,res
    enddo
    mg%stats%nite = nite
    mg%stats%res = res
    ncalls = ncalls + 1
  end subroutine solve

  subroutine vcycle(mg, lev0)
    type(MGinfos) :: mg
    integer:: lev, lev0, lev1

    lev1 = mg%nlevels-1
    do lev = lev0, lev1
       call smooth(mg, lev, mg%param%npre)
       call residual(mg, lev)
       call restriction(mg, lev, rb)
       if(mg%param%debug>5) call check_l2norms(mg, lev)
    enddo

    call smooth(mg, lev1+1, mg%param%nexact)
    if(mg%param%debug>5) call check_l2norms(mg, lev1+1)

    do lev = lev1, lev0, -1
       call prolongation(mg, lev)
       call smooth(mg, lev, mg%param%npost)
       if(mg%param%debug>5) call check_l2norms(mg, lev)
    enddo
  end subroutine vcycle


  subroutine fcycle(mg, lev0)
    type(MGinfos) :: mg
    integer:: lev, lev0, lev1

    lev1 = mg%nlevels-1
    do lev = 1, lev1
       if(lev>1) then
          call restriction(mg, lev, bb)
       else
          call restriction(mg, lev, rb)
       endif
    enddo

    call smooth(mg, lev1+1, mg%param%nexact)

    do lev = lev1, 1, -1
       call prolongation(mg, lev)
       call Vcycle(mg, lev)
    enddo
  end subroutine fcycle


  subroutine check_l2norms(mg, lev)
    type(MGinfos) :: mg
    integer:: lev

    call residual(mg, lev)
    print*,"lev",lev,&
         sum(mg%data(lev)%x**2),&
         sum(mg%data(lev)%b**2),&
         sum(mg%data(lev)%r**2)
  end subroutine check_l2norms

end module solvers
