module mod_halo

  use mod_tiles

  implicit none

  type MPI_halobuffers
     real,dimension(:,:,:,:),allocatable:: i,j,ij
  end type MPI_halobuffers

  type MPI_halo
     type(MG_grid):: grid
     type(MPI_tile):: tile
     type(MPI_neighbours):: nghb
     type(MPI_halobuffers):: send,recv
     type(MPI_Request), dimension(:), allocatable:: recv_reqs, send_reqs
     type(MPI_Status) ,dimension(:),allocatable::status
     integer:: request_number
     TYPE(MPI_Comm) :: comm
     logical:: iszperio=.false.
  end type MPI_halo

contains

  subroutine set_buffers(halo)
    type(MPI_halo),intent(inout):: halo
    integer::nx,ny,nz,nh

    nx = halo%grid%nx
    ny = halo%grid%ny
    nz = halo%grid%nz
    nh = halo%grid%nh
    !print*,"set buffers", nx,nx,ny,nz

    allocate(halo%send%i (nh,ny,nz,2))
    allocate(halo%recv%i (nh,ny,nz,2))

    allocate(halo%send%j (nx,nh,nz,2))
    allocate(halo%recv%j (nx,nh,nz,2))

    allocate(halo%send%ij(nh,nh,nz,4))
    allocate(halo%recv%ij(nh,nh,nz,4))

  end subroutine set_buffers

  subroutine set_requests(halo)
    type(MPI_halo),intent(inout),target:: halo
    integer:: i,n,di,dj,idx
    integer:: myrank, destrank,tag1, tag2,req_idx
    real,dimension(:,:,:),pointer:: sbuffer,rbuffer
    real::dummy
    integer:: ierr
    type(MPI_datatype):: mpitype

    !print*,"digits=",digits(dummy)
    if(digits(dummy).le.30)then
       mpitype=MPI_REAL
    else
       mpitype=MPI_DOUBLE
    end if

    n = halo%nghb%n

    allocate(halo%recv_reqs(n))
    allocate(halo%send_reqs(n))
    allocate(halo%status(n))

    myrank = halo%nghb%rank
    req_idx = 0
    halo%request_number = 0
    do i=1,9
       destrank = halo%nghb%infos(3,i)
       if((destrank/=-1).and.(destrank/=myrank))then
          di = halo%nghb%infos(1,i)
          dj = halo%nghb%infos(2,i)
          tag1 = 3+di+dj*2
          tag2 = 3-di-dj*2

          if(di/=0 .and. dj/=0)then
             idx = 1+(di+1)/2+(dj+1)
             sbuffer => halo%send%ij(:,:,:,idx)
             rbuffer => halo%recv%ij(:,:,:,idx)
          elseif(dj==0)then
             idx = 1+(di+1)/2
             sbuffer => halo%send%i(:,:,:,idx)
             rbuffer => halo%recv%i(:,:,:,idx)
          else
             idx = 1+(dj+1)/2
             sbuffer => halo%send%j(:,:,:,idx)
             rbuffer => halo%recv%j(:,:,:,idx)
          end if
          req_idx = req_idx+1
          !print'((I2),"<->",(I2)," i=",(I1)," tags=",(I3),",",(I3))',myrank,destrank,i,tag1,tag2
          halo%request_number = halo%request_number + 1
          call MPI_Send_init(sbuffer, size(sbuffer), mpitype, destrank, tag1, &
               halo%comm, halo%send_reqs(req_idx),ierr)

          call MPI_Recv_init(rbuffer, size(rbuffer), mpitype, destrank, tag2, &
               halo%comm, halo%recv_reqs(req_idx),ierr)
       endif
    enddo

  end subroutine set_requests

  subroutine set_halo(g, grid, halo)
    type(MPI_tile) :: g
    type(MG_Grid):: grid
    type(MPI_halo),intent(inout):: halo
    integer:: ierr

    halo%grid = grid
    halo%tile = g
    halo%comm = MPI_COMM_WORLD
    call set_neighbours(g, halo%nghb)
    !call print_neighbours(halo%nghb)
    call set_buffers(halo)
    call set_requests(halo)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    halo%iszperio = .false.
    if(g%topology == zperio .or. g%topology == xyzperio) then
       halo%iszperio = .true.
    else
       halo%iszperio = .false.
    endif
  end subroutine set_halo

  subroutine copytobuffer(halo, array, nx,ny,nz,nh)
    real,dimension(1-nh:nx+nh,1-nh:ny+nh,1:nz)::array
    type(MPI_halo):: halo
    integer,dimension(9):: rank
    integer::nx,ny,nz,nh

    rank = halo%nghb%infos(3,:)
!!$    print*,"shape array", shape(array)
!!$    print*,"shape send%i", shape( halo%send%i   )
!!$    print*,"shape send%j", shape( halo%send%j   )
!!$    print*,"shape send%ij", shape( halo%send%ij   )
    if(rank(4)/=-1) halo%send%i(:,:,:,1)  = array(1:nh,1:ny,1:nz)
    if(rank(6)/=-1) halo%send%i(:,:,:,2)  = array(nx-nh+1:nx,1:ny,1:nz)
    if(rank(2)/=-1) halo%send%j(:,:,:,1)  = array(1:nx,1:nh,1:nz)
    if(rank(8)/=-1) halo%send%j(:,:,:,2)  = array(1:nx,ny-nh+1:ny,1:nz)
    if(rank(1)/=-1) halo%send%ij(:,:,:,1) = array(1:nh,1:nh,1:nz)
    if(rank(3)/=-1) halo%send%ij(:,:,:,2) = array(nx-nh+1:nx,1:nh,1:nz)
    if(rank(7)/=-1) halo%send%ij(:,:,:,3) = array(1:nh,ny-nh+1:ny,1:nz)
    if(rank(9)/=-1) halo%send%ij(:,:,:,4) = array(nx-nh+1:nx,ny-nh+1:ny,1:nz)
  end subroutine copytobuffer

  subroutine copytoarray(halo, array, nx,ny,nz,nh)
    real,dimension(1-nh:nx+nh,1-nh:ny+nh,1:nz)::array
    type(MPI_halo):: halo
    integer,dimension(9):: rank
    integer::nx,ny,nz,nh

    rank = halo%nghb%infos(3,:)
    if(rank(4)/=-1) array(1-nh:0,1:ny,1:nz)           = halo%recv%i(:,:,:,1)
    if(rank(6)/=-1) array(nx+1:nx+nh,1:ny,1:nz)       = halo%recv%i(:,:,:,2)
    if(rank(2)/=-1) array(1:nx,1-nh:0,1:nz)           = halo%recv%j(:,:,:,1)
    if(rank(8)/=-1) array(1:nx,ny+1:ny+nh,1:nz)       = halo%recv%j(:,:,:,2)
    if(rank(1)/=-1) array(1-nh:0,1-nh:0,1:nz)         = halo%recv%ij(:,:,:,1)
    if(rank(3)/=-1) array(nx+1:nx+nh,1-nh:0,1:nz)     = halo%recv%ij(:,:,:,2)
    if(rank(7)/=-1) array(1-nh:0,ny+1:ny+nh,1:nz)     = halo%recv%ij(:,:,:,3)
    if(rank(9)/=-1) array(nx+1:nx+nh,ny+1:ny+nh,1:nz) = halo%recv%ij(:,:,:,4)
  end subroutine copytoarray

  subroutine check_recv_buffers(halo)
    type(MPI_halo):: halo
    integer,dimension(9):: rank
    integer::r

    r = halo%nghb%rank
    rank = halo%nghb%infos(3,:)
    if(rank(4)/=-1) print*,4,r,  halo%recv%i(:,:,:,1)
    if(rank(6)/=-1) print*,6,r,  halo%recv%i(:,:,:,2)
    if(rank(2)/=-1) print*,2,r,  halo%recv%j(:,:,:,1)
    if(rank(8)/=-1) print*,8,r,  halo%recv%j(:,:,:,2)
    if(rank(1)/=-1) print*,1,r,  halo%recv%ij(:,:,:,1)
    if(rank(3)/=-1) print*,3,r,  halo%recv%ij(:,:,:,2)
    if(rank(7)/=-1) print*,7,r,  halo%recv%ij(:,:,:,3)
    if(rank(9)/=-1) print*,9,r, halo%recv%ij(:,:,:,4)

  end subroutine check_recv_buffers

  subroutine check_send_buffers(halo)
    type(MPI_halo):: halo
    integer,dimension(9):: rank
    integer::r

    r = halo%nghb%rank
    rank = halo%nghb%infos(3,:)
    if(rank(4)/=-1) print*,4,r,  halo%send%i(:,:,:,1)
    if(rank(6)/=-1) print*,6,r,  halo%send%i(:,:,:,2)
    if(rank(2)/=-1) print*,2,r,  halo%send%j(:,:,:,1)
    if(rank(8)/=-1) print*,8,r,  halo%send%j(:,:,:,2)
    if(rank(1)/=-1) print*,1,r,  halo%send%ij(:,:,:,1)
    if(rank(3)/=-1) print*,3,r,  halo%send%ij(:,:,:,2)
    if(rank(7)/=-1) print*,7,r,  halo%send%ij(:,:,:,3)
    if(rank(9)/=-1) print*,9,r,  halo%send%ij(:,:,:,3)
  end subroutine check_send_buffers

  subroutine fill(halo, array, nx,ny,nz,nh)
    real,dimension(1-nh:nx+nh,1-nh:ny+nh,1:nz)::array
    type(MPI_halo):: halo
    integer::nx,ny,nz,nh
    integer::count, ierr

    count = halo%request_number

    if(count > 0) then

       call copytobuffer(halo, array, nx,ny,nz,nh)
       !call check_send_buffers(halo)
       !print*,">> halo"
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       call MPI_STARTALL(count, halo%recv_reqs, ierr)
       call MPI_STARTALL(count, halo%send_reqs, ierr)

       call MPI_WAITALL(count, halo%send_reqs, halo%status, ierr)
       call MPI_WAITALL(count, halo%recv_reqs, halo%status, ierr)

       !print*,"=============="
       !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       !call check_recv_buffers(halo)

       !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
       !print*,"<< "
       call copytoarray(halo, array, nx,ny,nz,nh)

    endif

    call exchange_with_myself(halo, array, nx,ny,nz,nh)

  end subroutine fill

  subroutine exchange_with_myself(halo, array, nx,ny,nz,nh)
    real,dimension(1-nh:nx+nh,1-nh:ny+nh,1:nz)::array
    type(MPI_halo):: halo
    integer::nx,ny,nz,nh

    integer:: myrank
    integer,dimension(9):: rank

    rank = halo%nghb%infos(3,:)

    myrank = halo%nghb%rank

    if(rank(6)== myrank) array(nx+1:nx+nh,1:ny,1:nz)       = array(1:nh,1:ny,1:nz)
    if(rank(4)== myrank) array(1-nh:0,1:ny,1:nz)           = array(nx-nh+1:nx,1:ny,1:nz)

    if(rank(8)== myrank) array(1:nx,ny+1:ny+nh,1:nz)       = array(1:nx,1:nh,1:nz)
    if(rank(2)== myrank) array(1:nx,1-nh:0,1:nz)           = array(1:nx,ny-nh+1:ny,1:nz)

    if(rank(9)== myrank) array(nx+1:nx+nh,ny+1:ny+nh,1:nz) = array(1:nh,1:nh,1:nz)
    if(rank(7)== myrank) array(1-nh:0,ny+1:ny+nh,1:nz)     = array(nx-nh+1:nx,1:nh,1:nz)
    if(rank(3)== myrank) array(nx+1:nx+nh,1-nh:0,1:nz)     = array(1:nh,ny-nh+1:ny,1:nz)
    if(rank(1)== myrank) array(1-nh:0,1-nh:0,1:nz)         = array(nx-nh+1:nx,ny-nh+1:ny,1:nz)

    if(halo%iszperio)then
       array(:,:,1:nh) = array(:,:,nz-2*nh+1:nz-nh)
       array(:,:,nz-nh+1:nz) = array(:,:,nh+1:2*nh)
    endif
  end subroutine exchange_with_myself

  subroutine exchange_with_my_neighbours(halo)
    ! demo of peer to peer exchange
    type(MPI_halo):: halo
    integer,dimension(9):: msg
    integer:: i,count,n,di,dj
    integer:: myrank, destrank,tag1, tag2, ierr

    n = halo%nghb%n
    myrank = halo%nghb%rank
    count=0
    msg(:)=-1
    do i=1,9
       di = halo%nghb%infos(1,i)
       dj = halo%nghb%infos(2,i)
       destrank = halo%nghb%infos(3,i)
       tag1 = (2+di)+(2+dj)*10
       tag2 = (2-di)+(2-dj)*10
       if(destrank>-1)then
          call MPI_SEND(myrank, 1, MPI_INT, destrank, tag1, halo%comm, ierr)
          call MPI_RECV(msg(i), 1, MPI_INT, destrank, tag2, halo%comm, MPI_STATUS_IGNORE, ierr)
       endif

    enddo
    print*,"halo test:",myrank,msg

  end subroutine exchange_with_my_neighbours


end module mod_halo
