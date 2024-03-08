module mg_setup

  use mg_types
  use mg_log
  use operators
  implicit none

  integer:: max_predef_mg
  parameter(max_predef_mg=10)
  type(MGinfos),dimension(max_predef_mg),target:: predef_mg
  integer,save:: nb_predef_mg=0

contains

  subroutine set_shape(grid)
    !type(MG_Param):: param
    type(MG_Grid):: grid
    !integer:: i0,i1,j0,j1,k0,k1
    !integer:: nx,ny,nz

    grid%i0 = 1-grid%nh
    grid%i1 = grid%nx+grid%nh
    grid%j0 = 1-grid%nh
    grid%j1 = grid%ny+grid%nh
    grid%k0 = 1
    grid%k1 = grid%nz
  end subroutine set_shape

  subroutine allocate_array(grid, array)
    type(MG_Grid)::grid
    REAL,dimension(:,:,:),allocatable::array
    !integer:: i0,i1,j0,j1,k0,k1

    !call get_shape(grid, i0,i1,j0,j1,k0,k1)
    allocate(array(grid%i0:grid%i1,grid%j0:grid%j1,grid%k0:grid%k1))

  end subroutine allocate_array


  subroutine print_level(mg, lev)
    integer:: lev
    type(MGinfos):: mg

    type(MG_Grid):: grid
    type(MPI_tile):: tile
    character(len=60)::msg

    grid = mg%hierarchy%grid(lev)
    tile = mg%hierarchy%tile(lev)

    if(tile%glue)then
       write(msg, '( (I3)," x",(I3)," x",(I3)," | ",(I3)," x",(I3)," | ",(I2)," x",(I2))') &
            grid%nx,grid%ny,grid%nz,&
            tile%npx/tile%incrx,tile%npy/tile%incry, &
            tile%ni, tile%nj
    else
       write(msg, '( (I3)," x",(I3)," x",(I3)," | ",(I3)," x",(I3) )') &
            grid%nx,grid%ny,grid%nz,&
            tile%npx/tile%incrx,tile%npy/tile%incry
!            tile%incrx, tile%incry
    endif
    call logprint(mg, msg)
  end subroutine print_level

  subroutine print_mginfos(mg)
    type(MGinfos):: mg
    integer:: i,nx,ny,nz
    character(len=40)::msg

    write(msg, '("MG has",(I2), "levels")') mg%nlevels
    call logprint(mg, msg)
    print*, "parameters=", mg%param
    do i = 1, mg%nlevels
       call print_level(mg, i)
!       print*, i,mg%hierarchy%grid(i)
    enddo
!!$    do i = 1, mg%nlevels
!!$       print*, i,mg%hierarchy%tile(i)
!!$    enddo
  end subroutine print_mginfos

  subroutine setup_mg(mg, nx,ny,nz, topology)
    type(MGinfos):: mg
    integer:: nx,ny,nz,nh,npx,npy, rank
    logical:: is3d, iszperio
    integer:: topology

    if(mg%param%debug>0)then
       print*,"setting up mg"
    endif
    npx = mg%npx
    npy = mg%npy
    rank = mg%hierarchy%tile(1)%rank
    is3d = .not.(mg%param%is2d)
    if(topology == zperio .or. topology == xyzperio)then
       iszperio = .true.
    else
       iszperio = .false.
    endif
    nh = mg%param%nh
    call create_hierarchy(mg%hierarchy, nx,ny,nz, nh, npx, npy, rank, is3d, topology)
    call allocate_mg(mg, mg%hierarchy)
    call set_default_msk(mg, is3d, iszperio)

  end subroutine setup_mg


  subroutine allocate_mg(mg, h)
    type(MGinfos):: mg
    type(MG_Hierarchy):: h
    type(MG_Grid)::grid, gsmall

    integer:: i
    integer:: i0,i1,j0,j1,k0,k1
    integer:: nx,ny,nz

    if(mg%param%debug>0)then
       print*,"allocating mg"
    endif
    mg%nlevels = h%nlevels
    mg%hierarchy = h

    allocate(mg%data(h%nlevels))
    allocate(mg%oper(h%nlevels))
    allocate(mg%halo(h%nlevels))

    do i = 1, h%nlevels
       grid = h%grid(i)
       grid%nh=mg%param%nh
       call set_shape(grid)

       call allocate_array(grid, mg%data(i)%x)
       call allocate_array(grid, mg%data(i)%b)
       call allocate_array(grid, mg%data(i)%r)
       call allocate_array(grid, mg%data(i)%y)

       if(h%tile(i)%glue) then
          !call logprint(mg, "level i needs glue")
          !print*, "level needs glue",i
          gsmall = h%smallgrid(i)
          gsmall%nh = mg%param%nh
          call set_shape(gsmall)

          !call allocate_array(gsmall, mg%oper(i)%smallmsk)
          allocate(mg%oper(i)%smallmsk(gsmall%i0:gsmall%i1,gsmall%j0:gsmall%j1,gsmall%k0:gsmall%k1))
          call allocate_array(gsmall, mg%data(i)%small)

          call allocate_array(gsmall, mg%oper(i)%Rcoef)
          !call allocate_array(gsmall, mg%oper(i)%Pcoef)

       else

          call allocate_array(grid, mg%oper(i)%Rcoef)

       endif

       call allocate_array(grid, mg%oper(i)%Pcoef)

       call allocate_array(grid, mg%oper(i)%diag)
       call allocate_array(grid, mg%oper(i)%idiag)
       allocate(mg%oper(i)%msk(grid%i0:grid%i1,grid%j0:grid%j1,grid%k0:grid%k1))
       !call allocate_array(grid, mg%oper(i)%msk)

       nx = grid%nx
       ny = grid%ny
       nz = grid%nz
       !mg%oper(i)%msk = 0
       !mg%oper(i)%msk(1:nx,1:ny,1:nz) = 1

       grid%n = size(mg%data(i)%x)

       mg%hierarchy%grid(i) = grid
       mg%hierarchy%smallgrid(i) = gsmall

    end do
    mg%ndofs = h%grid(1)%nx*h%grid(1)%ny*h%grid(1)%nz

  end subroutine allocate_mg

  subroutine set_default_msk(mg, is3d, iszperio)
    type(MGinfos):: mg
    logical:: is3d, iszperio
    integer:: lev, nx,ny,nz,nh
    type(MG_Grid)::grid

    do lev = 1, mg%hierarchy%nlevels

       grid = mg%hierarchy%grid(lev)
       nx = grid%nx
       ny = grid%ny
       nz = grid%nz
       nh = grid%nh

       mg%oper(lev)%msk = 0

       if(is3d)then


          if(iszperio) then
             mg%oper(lev)%msk(1:nx,1:ny,:) = 1
          else
             mg%oper(lev)%msk(1:nx,1:ny,1+nh:nz-nh) = 1
          endif

       else

          mg%oper(lev)%msk(1:nx,1:ny,1:nz) = 1

       endif
    enddo
  end subroutine set_default_msk

  subroutine setup_fine_msk(mg)
    type(MGinfos):: mg
    type(MG_Grid):: g
    integer:: lev

    lev = 1
    g = mg%hierarchy%grid(lev)
    mg%data(lev)%y = mg%oper(lev)%msk
    call fill(mg%halo(lev), mg%data(lev)%y, g%nx, g%ny, g%nz, g%nh)
    mg%oper(lev)%msk = mg%data(lev)%y
   end subroutine setup_fine_msk

  subroutine create_hierarchy(h, nx,ny,nz, nh, npx, npy, rank, is3d, topology)
    integer:: nx,ny,nz, nh
    type(MG_Hierarchy),intent(inout):: h
    integer:: i,x,y,z
    integer:: npx,npy,ni,nj,incx,incy,npx0,npy0, rank
    logical:: is3d
    integer:: topology

    integer:: ngather=16*16, maxglue=8


    x = nx
    y = ny
    z = nz
    if(is3d) z = nz+2*nh

    i = 1
    h%grid(i)%nx = x
    h%grid(i)%ny = y
    h%grid(i)%nz = z

    h%tile(i)%npx = npx
    h%tile(i)%npy = npy
    h%tile(i)%rank = rank
    h%tile(i)%topology = topology
    !print*,"RANK=",rank

    npx0 = npx
    npy0 = npy

    incx = 1
    incy = 1

    h%tile(i)%incrx = incx
    h%tile(i)%incry = incy
    call set_loc(h%tile(i))

    do
       if((x==2).or.(y==2))exit

       x = x/2
       y = y/2
       if(is3d) z=z/2+nh
       i = i+1

       h%tile(i)%npx = npx
       h%tile(i)%npy = npy
       h%tile(i)%rank = rank
       h%tile(i)%topology = topology

       if(((x*y<=ngather).or.(x<8).or.(y<8)).and.(npx0*npy0>1))then
          ! glue
          h%smallgrid(i)%nx = x
          h%smallgrid(i)%ny = y
          h%smallgrid(i)%nz = z

          ni = min(maxglue, npx0)
          nj = min(maxglue, npy0)
          npx0 = npx0/ni
          npy0 = npy0/nj
          x = x*ni
          y = y*nj
          h%tile(i)%ni = ni
          h%tile(i)%nj = nj
          h%tile(i)%glue = .true.

          incx = incx*ni
          incy = incy*nj
       endif
       h%tile(i)%incrx = incx
       h%tile(i)%incry = incy
       call set_loc(h%tile(i))

       h%grid(i)%nx = x
       h%grid(i)%ny = y
       h%grid(i)%nz = z


    end do

    h%nlevels = i

  end subroutine create_hierarchy

  subroutine get_param(mg, param)
    type(MGinfos),intent(in):: mg
    type(MG_Param),intent(inout):: param

    param = mg%param
    !print*,param
  end subroutine get_param


  function get_ptrmg(npx,npy,nx,ny,nz,vertices,short,is3d, topology) result(mg)
    type(MGinfos),pointer:: mg
    integer::nx,ny,nz,npx,npy
    logical,optional:: vertices, short,is3d
    integer,optional:: topology

    type(MPI_tile):: g
    integer:: nprocs, ierr, lev
    integer:: nxs, nys, z, nh, ni, nj
    character*20::msg

    mg => null()
    !if(nb_predef_mg==0) then
       call MPI_INIT(ierr)
    !endif

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, g%rank, ierr)

    if(npx*npy /= nprocs) then
       print*,"wrong number of mpi cores, should be :", npx*npy
       stop
    endif
    g%npx = npx
    g%npy = npy
    call set_loc(g)


    if(nb_predef_mg<max_predef_mg)then
       nb_predef_mg = nb_predef_mg+1
       mg => predef_mg(nb_predef_mg)
       if(g%rank==0)then
          print '((A)(i1))', "define mg #", nb_predef_mg
       endif

       if(present(vertices))then
          if(vertices) mg%param%vertices = vertices
       endif
       if(present(is3d))then
          if(is3d) mg%param%is2d = .false.
       endif
       if(present(topology))then
          g%topology = topology
       else
          g%topology = closed
       endif
       write(msg, '("topology=",(I1))') g%topology
       call logprint(mg, msg)

       mg%npx = npx
       mg%npy = npy
       mg%hierarchy%tile(1) = g
!!$       mg%hierarchy%tile(1)%rank = g%rank
!!$       mg%hierarchy%tile(1)%i = g%i
!!$       mg%hierarchy%tile(1)%j = g%j

       call setup_mg(mg, nx,ny,nz, g%topology)
       !call setup_operators(mg)

    else
       print*,"too many mg already allocated"
       print*,"STOP"
       stop
    endif

    mg%verbose=-g%rank+1

    if(g%rank==0)then
       mg%param%debug=1
    else
       mg%param%debug=0
    endif

    do lev=1, mg%nlevels
       call set_halo(mg%hierarchy%tile(lev),mg%hierarchy%grid(lev),mg%halo(lev))
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       !print*,"-----------------",lev
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       !call exchange_with_my_neighbours(mg%halo(lev))
       !print "(i3, ': ', 9(i3))",mg%hierarchy%tile(lev)%rank, mg%halo(lev)%nghb%infos(3,:)
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       !print*,mg%hierarchy%tile(lev)
       !call print_neighbours(mg%halo(lev)%nghb)
       !call exchange_with_my_neighbours(mg%halo(lev))
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
    enddo
    !stop
    do lev=1, mg%nlevels
       if(mg%hierarchy%tile(lev)%glue)then
          nxs = mg%hierarchy%smallgrid(lev)%nx
          nys = mg%hierarchy%smallgrid(lev)%ny
          z = mg%hierarchy%smallgrid(lev)%nz
          nh = mg%hierarchy%smallgrid(lev)%nh
          ni = mg%hierarchy%tile(lev)%ni
          nj = mg%hierarchy%tile(lev)%nj

          g = mg%hierarchy%tile(lev)
          g%incrx = g%incrx/ni
          g%incry = g%incry/nj
          mg%hierarchy%infos(lev) = setup_glue(&
               g,nxs, nys, z, nh, ni, nj, mg%data(lev)%small)
          !print*,"ok",nxs,nys!lev,nxs,nys,ni,nj
       endif
    enddo
    !call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !stop


    call logprint(mg, "halos are set")

    if((present(short)).and.(short))then
       ! do nothing
       print*,"[WARNING] operators are not defined yet"
    else
       call setup_fine_msk(mg)
       call setup_operators(mg)
    endif


  end function get_ptrmg

end module mg_setup
