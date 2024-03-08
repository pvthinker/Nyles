module tests

  use mg_setup
  use solvers
  use tuning
  use mg_log
  use mod_halo
  use mod_io

  implicit none

contains

  subroutine test3d()
    type(MGinfos),pointer:: mg
    type(MPI_tile):: g, tile
    type(MG_grid):: grid
    integer:: nxglo, nyglo, nx, ny, nz, nh
    integer:: nprocs,ierr, lev,nt=1,kt
    real:: t0,t1

    nxglo = 64*2
    nyglo = 64*2
    nz = 64*2

    g%npx = 1
    g%npy = 1

    nx = nxglo/g%npx
    ny = nyglo/g%npy

    nprocs = g%npx*g%npy
    mg => get_ptrmg(g%npx,g%npy,nx,ny,nz,&
         vertices=.false.,short=.false.,is3d=.true.,topology=closed)

    !call tune(mg)
    mg%param%omega = 0.85152

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, g%rank, ierr)


    if(g%rank==0)then
       print*,"****************************************"
       call print_mginfos(mg)
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)


    mg%param%maxite=10
    mg%param%tol = 1e-9
    mg%data(1)%x= 0.
    mg%data(1)%b = 0
    if(g%rank==0)then
       mg%data(1)%b(1,1,4) = 1.
       mg%data(1)%b(4,4,4) = -1.
    endif

    lev = 1
    grid = mg%hierarchy%grid(lev)

    call fill(mg%halo(lev), mg%data(lev)%b, grid%nx, grid%ny, grid%nz, grid%nh)


    call cpu_time(t0)
    do kt = 1, nt
    !mg%data(lev)%x = 0.
    call solve(mg)
    enddo
    call cpu_time(t1)

    if(g%rank==0)then
       print*,"time to solve", (t1-t0)/nt
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    lev = 2
    call compute_Pcoef(mg, lev)
    !lev=lev+1
    lev = 2
    grid = mg%hierarchy%grid(lev)
    call residual(mg, lev)
    !mg%data(lev)%x = mg%oper(lev)%Pcoef
    mg%data(lev)%r = mg%oper(lev)%msk
    call write_array(mg%data(lev)%r, size(mg%data(lev)%r), grid%nx, grid%ny, grid%nz, grid%nh, g%npx, g%npy, g%rank)

  end subroutine test3d

  subroutine test2d()
    type(MGinfos),pointer:: mg
    type(MPI_tile):: g, tile
    type(MG_grid):: grid
    integer:: nxglo, nyglo, nx, ny, nz, nh
    integer:: nprocs,ierr, lev,nt=1,kt
    real:: t0,t1

    nxglo = 64*2
    nyglo = 64*2
    nz = 1

    g%npx = 1
    g%npy = 1

    nx = nxglo/g%npx
    ny = nyglo/g%npy

    nprocs = g%npx*g%npy
    mg => get_ptrmg(g%npx,g%npy,nx,ny,nz,&
         vertices=.true.,short=.false.,is3d=.false.,topology=closed)

    !call tune(mg)
    mg%param%omega = 0.85152

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, g%rank, ierr)


    if(g%rank==0)then
       print*,"****************************************"
       call print_mginfos(mg)
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)


    mg%param%maxite=10
    mg%param%tol = 1e-9
    mg%data(1)%x= 0.
    mg%data(1)%b = 0
    if(g%rank==0)then
       mg%data(1)%b(1,1,1) = 1.
       mg%data(1)%b(4,4,1) = -1.
    endif

    lev = 1
    grid = mg%hierarchy%grid(lev)

    call fill(mg%halo(lev), mg%data(lev)%b, grid%nx, grid%ny, grid%nz, grid%nh)


    call cpu_time(t0)
    do kt = 1, nt
    !mg%data(lev)%x = 0.
    call solve(mg)
    enddo
    call cpu_time(t1)

    if(g%rank==0)then
       print*,"time to solve", (t1-t0)/nt
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)


  end subroutine test2d

  subroutine testmpi()
    type(MGinfos),pointer:: mg
    type(MPI_tile):: g, tile
    type(MG_grid):: grid
    type(MPI_Halo):: halo

    integer:: nprocs, ierr, lev
    logical:: isinit

    integer:: r1, r2
    integer:: nxglo, nyglo, nx, ny, nh

    integer:: icase, i, j

    real::t0,t1

    nxglo = 256
    nyglo = 256

    do icase = 2,2
       if(icase==1) then
          g%npx = 1
          g%npy = 1
       elseif(icase==2) then
          g%npx = 2
          g%npy = 2
       elseif(icase==3) then
          g%npx = 4
          g%npy = 4
       else
          print*,"not implemented"
          stop
       endif

       nx = nxglo/g%npx
       ny = nyglo/g%npy

       nprocs = g%npx*g%npy
       mg => get_ptrmg(g%npx,g%npy,nx,ny,3,short=.true.)

       call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, g%rank, ierr)


       if(g%rank==0)then
          print*,"****************************************"
          call print_mginfos(mg)
       endif
       call MPI_Barrier(MPI_COMM_WORLD, ierr)

       if(g%rank==2)then
          mg%oper(1)%msk(8:16,8:16,:) = 0
       endif

       call setup_fine_msk(mg)
       call setup_operators(mg)

       if(.false.) then
          lev = 3
          grid = mg%hierarchy%grid(lev)
          tile = mg%hierarchy%tile(lev)
          nh = 0!grid%nh
          do j = 1-nh,grid%ny+nh
             do i = 1-nh,grid%nx+nh
                mg%data(lev)%x(i,j,:) = i*j!+(tile%i*grid%nx) + (tile%j*grid%ny)
                !mg%data(lev)%x(i,j,:) = 3!tile%i + tile%j
             enddo
          enddo
          !print*,tile%rank,tile%i, tile%j
          print*,"lbound 00", lbound(mg%data(lev)%x),lbound(mg%data(lev)%small)
          call split(mg%data(lev)%x, mg%data(lev)%small, mg%hierarchy%infos(lev))
          !call write_array(mg%data(lev)%x, size(mg%data(lev)%r), grid%nx, grid%ny, grid%nz, grid%nh, g%npx, g%npy, g%rank)
          call MPI_Barrier(MPI_COMM_WORLD, ierr)

          grid = mg%hierarchy%smallgrid(lev)
          print*,grid
          call write_array(mg%data(lev)%small, size(mg%data(lev)%small), grid%nx, grid%ny, grid%nz, grid%nh, g%npx, g%npy, g%rank)

          call MPI_Barrier(MPI_COMM_WORLD, ierr)
          stop
       endif

       mg%param%maxite=10
       mg%data(1)%x= 0.
       mg%data(1)%b = 0
       if(g%rank==0)then
          mg%data(1)%b(1,1,1) = 1.
          mg%data(1)%b(40,40,1) = -1.
       endif
       !stop
       !lev = 1
       !grid = mg%hierarchy%grid(lev)
       !mg%data(lev)%x = 0
       !mg%data(lev)%x(1:grid%nx,1:grid%ny,:) = g%rank
       !call fill(mg%halo(lev), mg%data(lev)%x, grid%nx, grid%ny, grid%nz, grid%nh)
       !call write_array(mg%data(lev)%x, size(mg%data(lev)%r), grid%nx, grid%ny, grid%nz, grid%nh, g%npx, g%npy, g%rank)

       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       !call MPI_Barrier(MPI_COMM_WORLD, ierr)
       !call mpi_finalize(ierr)
       !stop


       call cpu_time(t0)
       call solve(mg)
       call cpu_time(t1)

       if(g%rank==0)then
          print*,"time to solve", (t1-t0)
       endif

       call MPI_Barrier(MPI_COMM_WORLD, ierr)

       lev = 2
       call compute_Pcoef(mg, lev)
       !lev=lev+1
       lev = 1
       grid = mg%hierarchy%grid(lev)
       call residual(mg, lev)
       !mg%data(lev)%r = mg%oper(lev)%Pcoef
       !mg%data(lev)%r = mg%oper(lev)%msk
       call write_array(mg%data(lev)%x, size(mg%data(lev)%r), grid%nx, grid%ny, grid%nz, grid%nh, g%npx, g%npy, g%rank)

       grid = mg%hierarchy%smallgrid(lev+1)
       !call write_array(mg%data(lev+1)%small, size(mg%data(lev+1)%small), grid%nx, grid%ny, grid%nz, grid%nh, g%npx, g%npy, g%rank)
       !mg%data(lev+1)%small = mg%oper(lev+1)%smallmsk
       !call write_array(mg%data(lev+1)%small, size(mg%data(lev+1)%small), grid%nx, grid%ny, grid%nz, grid%nh, g%npx, g%npy, g%rank)

    enddo
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    print*,g%rank,mg%halo(1)%nghb%infos(3,:)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call mpi_finalize(ierr)

  end subroutine testmpi

  subroutine testmpi_v0()
    type(MGinfos),pointer:: mg
    type(MPI_tile):: g
    type(MG_grid):: grid
    type(MPI_Halo):: halo

    integer:: nprocs, ierr, lev
    logical:: isinit

    integer:: r1, r2

    g%npx = 2
    g%npy = 2

    nprocs = g%npx*g%npy
    !call mpi_finalize(ierr)
    !call MPI_INIT_thread(4, 1, MPI_THREAD_FUNNELED)
!!$
!!$    if(g%npx*g%npy /= nprocs) then
!!$       print*,"wrong number of mpi cores, should be :", g%npx*g%npy
!!$       stop
!!$    endif
!!$
!!$    call set_loc(g)

    !mg => get_ptrmg(g%npx,g%npy,256,256,3)
    mg => get_ptrmg(g%npx,g%npy,64,64,3)

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, g%rank, ierr)

!!$    mg%verbose=-g%rank+1
!!$
!!$    mg%param%maxite=20
!!$    if(g%rank==0)then
!!$       mg%param%debug=1
!!$    else
!!$       mg%param%debug=0
!!$    endif
!!$
!!$    do lev=1, mg%nlevels
!!$       call set_halo(g,mg%hierarchy%grid(lev),mg%halo(lev))
!!$    enddo
!!$    if(mg%verbose>0)then
!!$       do lev=1, mg%nlevels
!!$          grid = mg%hierarchy%grid(lev)
!!$          print*,lev, grid%nx,grid%ny,grid%nh,grid%topology
!!$       enddo
!!$    endif
!!$    if(mg%verbose>0)then
!!$       do lev=1, mg%nlevels
!!$          halo = mg%halo(lev)
!!$          print*,lev, halo%nghb%rank, halo%nghb%infos
!!$       enddo
!!$    endif
!!$    call logprint(mg, "halos are set")
!!$
!!$    call setup_fine_msk(mg)
!!$    call setup_operators(mg)

!!$    do lev = 1, mg%nlevels
!!$       print*,lev,shape(mg%oper(lev)%msk),mg%oper(lev)%msk(:,4,1)
!!$    enddo


    if(g%rank==0)then
       call print_mginfos(mg)
    endif

    mg%data(1)%x= 0.
    if(g%rank==0)then
       mg%data(1)%b(1,1,1) = 1.
       mg%data(1)%b(40,40,1) = -1.
    endif
    mg%param%maxite=10
    call solve(mg)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

!!$    lev = 2
!!$    call residual(mg, lev)
!!$    do r1 = 0, nprocs-1
!!$       if(g%rank==r1)then
!!$          !print*,"*** rank=",r1,mg%oper(lev)%diag(:,1,1)
!!$          print*,"*** rank=",r1,mg%data(lev)%r(:,1,1)
!!$       endif
!!$       call MPI_Barrier(MPI_COMM_WORLD, ierr)
!!$    enddo
    lev = 1
    grid = mg%hierarchy%grid(lev)
    call write_array(mg%data(lev)%r, size(mg%data(lev)%r), grid%nx, grid%ny, grid%nz, grid%nh, g%npx, g%npy, g%rank)

    call mpi_finalize(ierr)
  end subroutine testmpi_v0


  subroutine testdebug()
    type(MGinfos),pointer:: mg

    mg => get_ptrmg(1,1,128,128,128)
    mg%param%maxite=1
    mg%param%debug=1
    call setup_operators(mg)

    call print_mginfos(mg)

    mg%data(1)%x= 0.
    mg%data(1)%b(1,1,1) = 1.
    mg%data(1)%b(80,80,1) = -1.

    !print*,"shape(b)=",shape(mg%data(1)%b), norm(mg, 1, bnorm),mg%hierarchy%grid(1)%n
    call solve(mg)

  end subroutine testdebug

  subroutine testsolve()
    type(MGinfos),pointer:: mg
    !mg => get_ptrmg(1,1,256,256,1)
    mg => get_ptrmg(1,1,128,128,128)
    !mg => get_ptrmg(1,1,128,128,1)
    !mg%param%vertices = .true.
    call setup_operators(mg)


    call print_mginfos(mg)

    mg%verbose = 0
    call tune(mg)

    mg%data(1)%x= zero
    mg%data(1)%b(20,20,1) = 1
    mg%data(1)%b(80,80,1) = -1

    mg%verbose = 1
    call solve(mg)
  end subroutine testsolve


  subroutine testspeed()
    type(MGinfos),pointer:: mg
    real:: t0, t1
    integer:: i, nt=1000

    !mg => get_ptrmg(1,1,256,256,1)
    mg => get_ptrmg(1,1,128,128,128)

    call print_mginfos(mg)

    call cpu_time(t0)
    do i = 1, nt
       call vcycle(mg, 1)
    enddo
    call cpu_time(t1)


    print*,"time per Vcycle", (t1-t0)/nt
    print*,"time per Vcycle per dof", (t1-t0)/(nt*mg%ndofs)

  end subroutine testspeed

  subroutine testtune()
    type(MGinfos),pointer:: mg
    REAL::omega
    mg => get_ptrmg(1,1,256,256,1)
    debug_tuning = 1
    call tune(mg)
  end subroutine testtune

end module tests
