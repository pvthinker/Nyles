module tuning

  use solvers
  implicit none

  integer:: debug_tuning = 0

contains

  subroutine setup_initial_state(mg)
    type(MGinfos) :: mg
    integer:: i,j,k, lev
    type(MG_grid):: grid

    i = mg%hierarchy%grid(1)%nx/2
    j = mg%hierarchy%grid(1)%ny/2
    mg%data(1)%x(:,:,:) = zero
    mg%data(1)%b(:,:,:) = zero
    if(mg%hierarchy%tile(1)%rank == 0) then
       if(mg%param%is2d) then
          mg%data(1)%b(i-5,j,1) = one
          mg%data(1)%b(i,j+5,1) = -one
       else
          mg%data(1)%b(i-5,j,5) = one
          mg%data(1)%b(i,j+5,5) = -one
       endif
    endif
    lev = 1
    grid = mg%hierarchy%grid(lev)

    call fill(mg%halo(lev), mg%data(lev)%b, grid%nx, grid%ny, grid%nz, grid%nh)

  end subroutine setup_initial_state

  subroutine tune(mg)
    type(MGinfos) :: mg
    REAL:: omega, omega0
    REAL:: domega=0.2, res, res0
    integer:: i, nite=15
    integer:: verbose
    type(MG_Param):: param
    character*40::msg

    param = mg%param
    verbose = mg%verbose

    mg%verbose = 0
    mg%param%maxite = 7

    omega = 0.7

    call setup_initial_state(mg)
    call solve(mg)
    res0 = mg%stats%res
    omega = omega + domega
    do i = 1, nite
       omega0 = omega
       mg%param%omega=omega
       call setup_initial_state(mg)
       call solve(mg)
       res = mg%stats%res
       if(res>res0) then
          omega = omega+domega
          domega = -domega/1.5
       else
          omega = omega-domega
          domega = domega/1.5
       endif
       res0 = res

       if(debug_tuning>0) print*,res,omega
    enddo
    write(msg,'("tuning: best omega = ",(F8.5))') omega
    call logprint(mg, msg)!print*,"tuning: best omega = ", omega
    mg%param = param
    mg%param%omega = omega
    mg%verbose = verbose
    ncalls = 0
    mg%data(1)%b(:,:,:) = zero
    mg%data(1)%x(:,:,:) = zero
  end subroutine tune

end module tuning
