module mod_io

  use mpi_f08
  !use mod_halo
  implicit none

contains
  function get_pos(length, rank) result(pos)
    integer:: length, rank
    integer:: pos
    integer:: realsize, integersize

    realsize = 8
    integersize = 4

    pos = rank*realsize*length+7*integersize
  end function get_pos

  subroutine write_array(vector, length, nx, ny, nz, nh, npx, npy, rank)
    real,dimension(length):: vector
    integer:: nx, ny, nz, nh
    integer:: length, npx,npy,rank
    character(len=20):: filename
    integer:: fid, pos
    integer:: ierr

    fid = 20
    filename = "out.dat"
    pos = get_pos(length, rank)
    !print*,"rank=",rank,"pos=",pos

    open(fid, access="stream", form="unformatted")

    if ((rank==npx*npy-1).and.(rank>0)) then

       write(fid, pos=pos+1) vector
       call MPI_Barrier(MPI_COMM_WORLD, ierr)

    else

       call MPI_Barrier(MPI_COMM_WORLD, ierr)

       if(rank==0)then
          write(fid) nx, ny, nz, nh
          write(fid) length, npx, npy
       endif

       write(fid, pos=pos+1) vector

    endif

    close(fid)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end subroutine write_array

end module mod_io
