module gluesplit

  use mpi_f08
  use mod_tiles

  implicit none

  type T_glue
     real,dimension(:,:), allocatable:: buffers
     integer, dimension(:), allocatable:: offsets
     integer:: smallsize, buffersize
     integer:: nxs, nys, nz, nh, ic, jc, color, pos
     integer:: ni, nj
     TYPE(MPI_Comm) :: localcomm
     type(MPI_Datatype):: datatype
     type(MPI_Request) :: allscatter_request
  end type T_glue

  type T_colorpos
     integer:: color, pos, jc, ic
  end type T_colorpos

contains

  function get_color1d(npx, incx, ni) result(color)
    integer:: npx, incx, ni
    integer,dimension(0:npx-1):: color

    integer:: i, n2

    n2 = ni*incx
    do i = 0, npx-1
       color(i) = mod(i, incx)+n2*(i/n2)
    end do

  end function get_color1d

  function get_color2d(npx, npy, incx, incy, ni, nj) result(colors)
    integer:: npx, npy, incx, incy, ni, nj
    integer,dimension(0:npx*npy-1) :: colors

    integer,dimension(0:npx-1):: cx
    integer,dimension(0:npy-1):: cy
    integer:: i,j,k

    cx = get_color1d(npx, incx, ni)
    cy = get_color1d(npy, incy, nj)
    !print*,"cxcy,npx,incx,ni=",cx,cy,npx, incx, ni
    k = 0
    do j = 0, npy-1
       do i = 0, npx-1
          colors(k) = cx(i)+cy(j)*npx
          k = k+1
       enddo
    enddo
  end function get_color2d

  function get_colorpos(tile, ni, nj) result(colorpos)
    type(MPI_tile):: tile
    integer:: ni, nj
    type(T_colorpos):: colorpos

    integer,dimension(0:tile%npx*tile%npy-1) :: colors
    integer:: mycolor, pos, r

    colors = get_color2d(tile%npx, tile%npy, &
         tile%incrx, tile%incry, ni, nj)

    mycolor = colors(tile%rank)
    !print*,"ninj, rank=",ni,nj, tile%rank
    pos = -1
    do r = 0, tile%rank
       if(colors(r).eq.mycolor)then
          pos = pos+1
       endif
    enddo
    colorpos%color = mycolor
    colorpos%pos = pos+1
    colorpos%ic = mod(pos, ni)
    colorpos%jc = pos/ni

  end function get_colorpos

  function setup_glue(tile, nxs, nys, nz, nh, ni, nj, small_array) result(infos)
    type(MPI_tile):: tile
    type(T_glue):: infos
    integer:: nxs, nys, nz, nh, ni, nj
    real,dimension(nxs+2*nh,nys+2*nh,nz) :: small_array

    type(T_colorpos):: colorpos
    integer:: nxb, nyb, smallsize, ierr, color, i

    TYPE(MPI_Info) :: info
    !TYPE(MPI_Request), INTENT(OUT) :: request

    nxb = nxs*ni
    nyb = nys*nj

    smallsize = (nxs+2*nh)*(nys+2*nh)*nz

    infos%smallsize = smallsize
    infos%datatype = MPI_DOUBLE

    colorpos = get_colorpos(tile, ni, nj)
    infos%nxs = nxs
    infos%nys = nys
    infos%nz = nz
    infos%nh = nh
    infos%ni = ni
    infos%nj = nj
    infos%color = colorpos%color
    infos%ic = colorpos%ic
    infos%jc = colorpos%jc
    !print*,"icjc=",infos%ic,infos%jc,tile%rank,tile%npx,tile%npy
    !print*,"colorpos=",colorpos
    !call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Comm_Split(MPI_COMM_WORLD, colorpos%color, &
         colorpos%pos, infos%localcomm, ierr)
    !print*,"SPLIT OK", colorpos%color
    !call MPI_Barrier(MPI_COMM_WORLD, ierr)
    allocate(infos%buffers(smallsize,ni*nj))
    !allocate(infos%buffersize(smallsize*ni*nj))
    allocate(infos%offsets(ni*nj))

    infos%buffersize = smallsize*ni*nj
    infos%offsets(1) = 1
    do i = 2, ni*nj
       infos%offsets(i) = infos%offsets(i-1)+smallsize
    enddo

    !call MPI_info_set(info, "host", "localhost",ierr)

!!$    call MPI_Allgather_init(&
!!$         small_array, infos%smallsize, MPI_DOUBLE, &
!!$         infos%buffers, infos%smallsize, MPI_DOUBLE, &
!!$         infos%localcomm, &
!!$         MPI_info_null, infos%allscatter_request, ierr)

  end function setup_glue


  subroutine glue(big_array, small_array, infos)
    real,dimension(:,:,:) :: big_array, small_array
    type(T_glue):: infos

    integer:: nxb, nyb, nxs, nys, nz, nh, ni, nj

    nxs = infos%nxs
    nys = infos%nys
    nz = infos%nz
    nh = infos%nh
    ni = infos%ni
    nj = infos%nj
    nxb = nxs*ni
    nyb = nys*nj
    !print*,"buffer shape=",shape(infos%buffers), nxs+2*nh, nys+2*nh, nz
    call gathersubdomains(small_array, infos, nxs+2*nh, nys+2*nh, nz)

    call gth(big_array, infos%buffers, nxb, nyb, nxs, nys, nz, nh, ni, nj)
  end subroutine glue


  subroutine split(big_array, small_array, infos)
    real,dimension(:,:,:) :: big_array, small_array
    type(T_glue):: infos

    integer:: nxb, nyb, nxs, nys, nz, nh, ni, nj, ic, jc

    nxs = infos%nxs
    nys = infos%nys
    nz = infos%nz
    nh = infos%nh
    ni = infos%ni
    nj = infos%nj
    nxb = nxs*ni
    nyb = nys*nj
    ic = infos%ic
    jc = infos%jc

    call splt(big_array, small_array, nxb, nyb, nxs, nys, nz, nh, ic, jc)
    !print*,"info nh",nh
  end subroutine split


  subroutine gathersubdomains(small_array, infos, n, m, l)
    real,dimension(n,m,l) :: small_array
    type(T_glue):: infos
    integer:: n,m,l
    type(MPI_status):: status

    integer:: ierr

    !call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !print*,"small vs buffer", shape(small_array), shape(infos%buffers)
    !call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Allgather(&
         small_array, infos%smallsize, MPI_DOUBLE, &
         infos%buffers, infos%smallsize, MPI_DOUBLE, &
         infos%localcomm, ierr)

!!$    call MPI_Start(infos%allscatter_request, ierr)
!!$    call MPI_Wait(infos%allscatter_request, status, ierr)

    !print*,"IERR=",ierr
    !call MPI_Barrier(MPI_COMM_WORLD, ierr)
    !print*,"allgather ok"
  end subroutine gathersubdomains

  subroutine splitsubdomains(small_array, infos, n, m, l)
    real,dimension(n,m,l) :: small_array
    type(T_glue):: infos
    integer:: n,m,l

    integer:: ierr

!      int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
!       void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)

    ! call MPI_scatter()

  end subroutine splitsubdomains

  subroutine gth(big, buffers, nxb, nyb, nxs, nys, nz, nh, ni, nj)

    real, dimension(nxs+2*nh, nys+2*nh, nz, ni, nj):: buffers
    integer:: nxb, nyb, nxs, nys, nz, nh
    integer:: ni, nj
    real, dimension(1-nh: nxb+nh, 1-nh: nyb+nh, nz):: big

    integer:: isub, jsub, i0, j0, j, n1, n2, k

    n1 = nxs+2*nh
    n2 = nys+2*nh

    do k = 1, nz
       jsub = 1
       j0 = 1
       do j = 1-nh, nyb+nh
          !print*, j, j0, jsub
          do isub = 1, ni
             !i0 = (isub-1)*nxs+1-nh
             i0 = 1+(isub-1)*nxs
             big(i0: i0+nxs, j, k) = buffers(nh+1: nh+nxs, j0, k, isub, jsub)
          enddo
          j0 = j0+1
          if ((j0.gt.(nys+nh)).and.(jsub.lt.nj))then
             j0 = 1+nh
             jsub = jsub+1
          end if
       enddo
    enddo

  end subroutine gth

  subroutine splt(big, small, nxb, nyb, nxs, nys, nz, nh, ic, jc)

    integer :: nxb, nyb, nxs, nys, nz, nh
    integer :: ic, jc
    real, dimension(1-nh:nxb+nh, 1-nh:nyb+nh, nz):: big
    real, dimension(1-nh:nxs+nh, 1-nh:nys+nh, nz):: small

    integer:: i, j, n1, i0, j0, k

    n1 = nxs+2*nh

    !print*,"shape", shape(big),shape(small), ic, jc, lbound(big),lbound(small)
    i0 = ic*nxs
    j0 = jc*nys
    do k = 1, nz
       do j = 1-nh, nys+nh
          do i = 1-nh, nxs+nh
             small(i , j, k) = big(i+i0, j+j0, k)
          enddo
       enddo
    enddo
  end subroutine splt



end module gluesplit
