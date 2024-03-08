module operators
  use mg_types
  use mg_log
  !use mod_tiles
  !use mg_setup
  use basicoperators
  use mod_halo

  implicit none

  enum, bind(c)
     enumerator :: rnorm=1, bnorm, rb, bb
  endenum

  pointer :: fsmoother, fresidual, frestrict,fprolongation, fnorm

  interface
     subroutine fsmoother(x,b,y,idiag,omega,nx,ny,nz,nh)
       integer::nh
       integer:: nx,ny,nz
       REAL,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::x,y,b,idiag
       real:: omega
       intent(inout)::x,y
       intent(in)::b,idiag
     end subroutine fsmoother

     subroutine fresidual(x,b,r,msk,diag,nx,ny,nz,nh)
       integer::nh
       integer:: nx,ny,nz
       REAL,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz):: b,x,r,diag
       integer,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz):: msk
       intent(out):: r
       intent(in)::x,b,msk,diag
     end subroutine fresidual

     subroutine frestrict(xf, xc, coef, nxf,nyf,nzf,nxc,nyc,nzc,nh)
       integer::nh
       integer:: nxf,nyf,nzf
       integer:: nxc,nyc,nzc
       REAL,dimension(1-nh:nxf+nh, 1-nh:nyf+nh, 1-nh:nzf+nh):: xf
       REAL,dimension(1-nh:nxc+nh, 1-nh:nyc+nh, 1-nh:nzc+nh):: xc, coef
       intent(in) :: xf,coef
       intent(out):: xc
     end subroutine frestrict

     subroutine fprolongation(xf, xc, coef, nxf,nyf,nzf,nxc,nyc,nzc,nh)
       integer::nh
       integer:: nxf,nyf,nzf
       integer:: nxc,nyc,nzc
       REAL,dimension(1-nh:nxf+nh, 1-nh:nyf+nh, 1-nh:nzf+nh):: xf, coef
       REAL,dimension(1-nh:nxc+nh, 1-nh:nyc+nh, 1-nh:nzc+nh):: xc
       !REAL,dimension(:,:,:):: xf, coef
       !REAL,dimension(:,:,:):: xc
       intent(in) :: xc,coef
       intent(out):: xf
     end subroutine fprolongation

     function fnorm(msk, x, nx,ny,nz,nh) result(sum)
       integer:: nx,ny,nz,nh
       REAL,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::x
       integer,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::msk
       REAL:: sum
     end function fnorm
  end interface

contains

  subroutine fillmsk(mg,lev)
    type(MGinfos):: mg
    type(MG_Grid):: g
    integer:: lev

    g = mg%hierarchy%grid(lev)
    mg%data(lev)%y = mg%oper(lev)%msk
    call fill(mg%halo(lev), mg%data(lev)%y, g%nx, g%ny, g%nz, g%nh)
    mg%oper(lev)%msk = mg%data(lev)%y

  end subroutine fillmsk


  function normresidual(mg, normb)
    type(MGinfos) :: mg
    REAL:: normb, res
    REAL:: normresidual
    integer:: lev

    lev = 1

    if(normb>0.)then
       call residual(mg, lev)
       res = norm(mg, lev, rnorm)/normb
    else
       res = 0.
    end if
    normresidual = res
  end function normresidual

  function norm(mg, lev, which)
    type(MGinfos) :: mg
    type(MG_Grid):: g
    integer:: lev, which, n
    REAL:: norm
    REAL,dimension(1):: send,recv
    integer:: ierr

    !n = mg%hierarchy%grid(lev)%n
    n = size(mg%data(lev)%r)

    g = mg%hierarchy%grid(lev)

    if(which==rnorm) then
       norm = fnorm(mg%oper(lev)%msk, mg%data(lev)%r, g%nx, g%ny, g%nz, g%nh)

    elseif (which==bnorm) then
       norm = fnorm(mg%oper(lev)%msk, mg%data(lev)%b, g%nx, g%ny, g%nz, g%nh)

    else
       stop
    endif

    send(1) = norm
    call MPI_Allreduce(send, recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
    norm = recv(1)

  end function norm

  subroutine residual(mg, lev)
    type(MGinfos) :: mg
    type(MG_Grid):: g
    integer:: lev

    g = mg%hierarchy%grid(lev)

    call fresidual(&
         mg%data(lev)%x,&
         mg%data(lev)%b, &
         mg%data(lev)%r, &
         mg%oper(lev)%msk, &
         mg%oper(lev)%diag, &
         g%nx, g%ny, g%nz, &
         mg%param%nh)

!!$    print*, "----------------------------------------"
!!$    print*,shape(mg%data(lev)%r)
!!$    print*, g%nx, g%ny, g%nz, g%nh
!!$    print*, "----------------------------------------"
    call fill(mg%halo(lev), mg%data(lev)%r, g%nx, g%ny, g%nz, g%nh)

  end subroutine residual

  subroutine smooth(mg, lev, nite)
    type(MGinfos) :: mg
    type(MG_Grid):: g
    integer:: lev, nite

    g = mg%hierarchy%grid(lev)

    call fsmoother(&
         mg%data(lev)%x,&
         mg%data(lev)%b, &
         mg%data(lev)%y, &
         !mg%oper(lev)%msk, &
         mg%oper(lev)%idiag, &
         mg%param%omega, &
         g%nx, g%ny, g%nz, &
         mg%param%nh)

    call fill(mg%halo(lev), mg%data(lev)%x, g%nx, g%ny, g%nz, g%nh)

  end subroutine smooth

  subroutine restriction(mg, lev, which)
    type(MGinfos), target :: mg
    type(MG_Grid):: g,gc
    integer:: lev, which

    real,dimension(:,:,:), pointer:: coarse, fine!, Rcoef


    g = mg%hierarchy%grid(lev)

    if(which == rb)then
       fine => mg%data(lev)%r
    else
       fine => mg%data(lev)%b
    endif

    if(mg%hierarchy%tile(lev+1)%glue)then
       coarse => mg%data(lev+1)%small
       gc = mg%hierarchy%smallgrid(lev+1)
    else
       coarse => mg%data(lev+1)%b
       gc = mg%hierarchy%grid(lev+1)
    endif

    call frestrict(&
         fine,&
         coarse, &
         mg%oper(lev+1)%Rcoef, &
         g%nx, g%ny, g%nz,&
         gc%nx, gc%ny, gc%nz, mg%param%nh)

    if(mg%hierarchy%tile(lev+1)%glue)then
       !print*,"CALL GLUE"
       call glue(mg%data(lev+1)%b, coarse, mg%hierarchy%infos(lev+1))
       !print*,"DONE GLUE"
    endif

    mg%data(lev+1)%x = zero
    gc = mg%hierarchy%grid(lev+1)
    call fill(mg%halo(lev+1), mg%data(lev+1)%b, gc%nx, gc%ny, gc%nz, gc%nh)
    !print*,"DONE FILL"
  end subroutine restriction

  subroutine prolongation(mg, lev)
    type(MGinfos), target :: mg
    type(MG_Grid):: g,gc
    integer:: lev

    real,dimension(:,:,:), pointer:: coarse

    g = mg%hierarchy%grid(lev)


    if(mg%hierarchy%tile(lev+1)%glue)then
       coarse => mg%data(lev+1)%small
       call split(mg%data(lev+1)%x, coarse, mg%hierarchy%infos(lev+1))
       gc = mg%hierarchy%smallgrid(lev+1)
    else
       coarse => mg%data(lev+1)%x
       gc = mg%hierarchy%grid(lev+1)
    endif


    call fprolongation(&
            mg%data(lev)%x,&
            coarse, &
            mg%oper(lev)%Pcoef, &
            g%nx, g%ny, g%nz,&
            gc%nx, gc%ny, gc%nz, mg%param%nh)

    call fill(mg%halo(lev), mg%data(lev)%x, g%nx, g%ny, g%nz, g%nh)

  end subroutine prolongation

  subroutine apply_default_msk(mg, lev)
    type(MGinfos):: mg
    integer:: lev, nx,ny,nz,i0,j0,npx,npy
    type(MPI_tile):: g
    type(MPI_neighbours) :: nghb
    integer:: east,west,north,south

    nx = mg%hierarchy%grid(lev)%nx
    ny = mg%hierarchy%grid(lev)%ny
    nz = mg%hierarchy%grid(lev)%nz

    g = mg%halo(lev)%tile
    nghb = mg%halo(lev)%nghb
    south = nghb%infos(3,2)
    west = nghb%infos(3,4)
    east = nghb%infos(3,6)
    north = nghb%infos(3,8)

    npx = mg%halo(lev)%tile%npx
    npy = mg%halo(lev)%tile%npy

    if(mg%param%vertices)then
       i0 = 2
       j0 = 2
    else
       i0 = 1
       j0 = 1
    endif
    ! msk should be completely initialized when used with the
    ! dynamic library

    !print*,"****************************************"
    !print*,">>>",npx,npy,g%i,g%j,g%rank

    mg%data(lev)%y = 1
    if(west==-1)then
       mg%data(lev)%y(:i0-1,:,:) = 0
    endif
    if(south==-1)then
       mg%data(lev)%y(:,:j0-1,:) = 0
    endif
    if(east==-1)then
       mg%data(lev)%y(nx+1:,:,:) = 0
    endif
    if(north==-1)then
       mg%data(lev)%y(:,ny+1:,:) = 0
    endif
    !mg%data(lev)%y(i0:nx,j0:ny,1:nz) = 1

    mg%oper(lev)%msk = mg%oper(lev)%msk*mg%data(lev)%y
    mg%data(lev)%y = 0
  end subroutine apply_default_msk

  subroutine compute_msk(mg, lev)
    type(MGinfos) :: mg
    integer:: lev
    integer:: threshold, ierr

    if(mg%param%vertices)then
       threshold = 2
    else
       threshold = 0
    endif

    !if(mg%param%debug>0) print*,"compute msk at level=",lev+1
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    mg%oper(lev+1)%msk = 1
    !call apply_default_msk(mg,lev+1)
    mg%oper(lev+1)%Rcoef = 1!mg%oper(lev+1)%msk
    mg%data(lev)%b = mg%oper(lev)%msk
    call restriction(mg, lev, bb)

    mg%oper(lev+1)%msk = zero
    where (mg%data(lev+1)%b>threshold)
       mg%oper(lev+1)%msk = one
    end where

    !call fillmsk(mg, lev+1)
    call apply_default_msk(mg,lev+1)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !print*,"DONE FILLMSK"
    if(mg%hierarchy%tile(lev+1)%glue)then
       mg%data(lev+1)%y = mg%oper(lev+1)%msk
       call split(mg%data(lev+1)%y, mg%data(lev+1)%small, mg%hierarchy%infos(lev+1))
       mg%oper(lev+1)%smallmsk = mg%data(lev+1)%small
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !print*,"DONE SPLIT"

  end subroutine compute_msk

  subroutine compute_Rcoef(mg, lev)
    type(MGinfos), target :: mg
    type(MG_grid):: g
    integer:: lev
    real,dimension(:,:,:), pointer:: indicatrix

    if(mg%param%debug>0) print*,"compute Rcoef at level=",lev+1

    if(mg%hierarchy%tile(lev+1)%glue)then

       mg%oper(lev+1)%Rcoef =  mg%oper(lev+1)%smallmsk

       indicatrix => mg%data(lev+1)%small

    else

       mg%oper(lev+1)%Rcoef = mg%oper(lev+1)%msk

       indicatrix => mg%data(lev+1)%b

    endif

    if(mg%param%vertices) then
       mg%data(lev)%b = mg%oper(lev)%msk
    else
       mg%data(lev)%b = one
    endif

    call restriction(mg, lev, bb)
    !print*,"RESTRICTION OK"
    mg%oper(lev+1)%Rcoef = zero

    indicatrix => mg%data(lev+1)%b

    if(mg%param%vertices) then
       where (indicatrix>0.)
          mg%data(lev+1)%y = 0.25
       end where
    else
       where (indicatrix>0.)
          mg%data(lev+1)%y = 4./indicatrix
       end where
    endif

    mg%data(lev+1)%y =mg%oper(lev+1)%msk*mg%data(lev+1)%y

    if(mg%hierarchy%tile(lev+1)%glue)then
       call split(mg%data(lev+1)%y, mg%oper(lev+1)%Rcoef, mg%hierarchy%infos(lev+1))
       !mg%oper(lev+1)%Rcoef = mg%oper(lev+1)%smallmsk*mg%oper(lev+1)%Rcoef
    else
       mg%oper(lev+1)%Rcoef = mg%data(lev+1)%y
    endif

    !g = mg%hierarchy%grid(lev+1)
    !call fill(mg%halo(lev+1), mg%oper(lev+1)%Rcoef, g%nx, g%ny, g%nz, g%nh)

  end subroutine compute_Rcoef

  subroutine compute_Pcoef(mg, lev)
    type(MGinfos) :: mg
    type(MG_grid):: g
    integer:: lev

    if(mg%param%debug>0) print*,"compute Pcoef at level=",lev

    mg%data(lev)%x = zero
    if(mg%param%vertices) then
       mg%data(lev+1)%x = one
    else
       mg%data(lev+1)%x = mg%oper(lev+1)%msk
    endif

    mg%oper(lev)%Pcoef = mg%oper(lev)%msk!one

    call prolongation(mg, lev)

    mg%oper(lev)%Pcoef = zero
    where (mg%data(lev)%x>0.)
       mg%oper(lev)%Pcoef = 1./mg%data(lev)%x
    end where

    mg%oper(lev)%Pcoef = mg%oper(lev)%msk*mg%oper(lev)%Pcoef

    g = mg%hierarchy%grid(lev)
    !call fill(mg%halo(lev), mg%oper(lev)%Pcoef, g%nx, g%ny, g%nz, g%nh)


  end subroutine compute_Pcoef

  subroutine compute_diag(mg, lev)
    type(MGinfos) :: mg
    type(MG_grid):: g
    integer:: lev

    !mg%data(lev)%y = mg%oper(lev)%msk

    if(mg%param%vertices) then
       mg%data(lev)%x = one
       !mg%oper(lev)%msk = 1 !<- DESTROY the mask !!! ERROR
    else
       mg%data(lev)%x = mg%oper(lev)%msk
    endif

    mg%data(lev)%b = zero
    mg%oper(lev)%diag = zero
    call residual(mg, lev)

    !mg%oper(lev)%msk = mg%data(lev)%y !<- DESTROY the mask !!! ERROR
    mg%oper(lev)%diag = -mg%data(lev)%r

    !g = mg%hierarchy%grid(lev)
    !call fill(mg%halo(lev), mg%oper(lev)%diag, g%nx, g%ny, g%nz, g%nh)

    mg%oper(lev)%idiag = zero
    where (mg%oper(lev)%diag>0.)
       mg%oper(lev)%idiag = 1./mg%oper(lev)%diag
    end where
    mg%data(lev)%x = zero


  end subroutine compute_diag



  subroutine setup_operators(mg)
    type(MGinfos) :: mg
    integer:: lev

    if(mg%param%is2d) then
       fsmoother => fsmoother2d
       fresidual => fresidual2d
       fnorm => fnorm2d
       if(mg%param%vertices) then
          call logprint(mg, "solve on vertices")
          frestrict => frestrict_vertices2d
          fprolongation => fprolongation_vertices2d
       else
          call logprint(mg, "solve on centers")
          frestrict => frestrict_centers2d
          fprolongation => fprolongation_centers2d
       endif
    else
       fsmoother => fsmoother3d
       fresidual => fresidual3d
       fnorm => fnorm3d
       if(mg%param%vertices) then
          call logprint(mg, "solve on vertices")
          !frestrict => frestrict_vertices3d
          !fprolongation => fprolongation_vertices3d
       else
          !call logprint(mg, "3D case not yet implemented")
          call logprint(mg, "solve on centers")
          frestrict => frestrict_centers3d
          fprolongation => fprolongation_centers3d
       endif

    endif

    do lev = 1, mg%nlevels-1
       call compute_msk(mg, lev)
    enddo
    do lev = 1, mg%nlevels-1
       call compute_Rcoef(mg, lev)
       call compute_Pcoef(mg, lev)
    enddo
    do lev = 1, mg%nlevels
       call compute_diag(mg, lev)
    enddo
  end subroutine setup_operators

end module operators
