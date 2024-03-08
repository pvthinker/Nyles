module mod_tiles

  use mpi_f08
  use mg_enums

  implicit none
  ! https://curc.readthedocs.io/en/latest/programming/MPI-Fortran.html

  type MPI_neighbours
     integer:: n, rank
     integer,dimension(3,9):: infos
  end type MPI_neighbours

  enum, bind(c)
     enumerator :: allranksworking=1,onlyrootworking,allrankswithprerequest
  endenum

  integer:: gathersplitmethod = allranksworking

  type MPI_tile
     integer:: npx=1,npy=1,incrx=1,incry=1
     integer:: i=0,j=0,rank=0,root=0
     !root is the rank in the local communicator in charge of the gather and scatter
     ! a rank for whose rank != root has all its neighbours deactivated -> no halo filling
     integer:: topology
     integer:: ni, nj
     logical:: glue=.false.
  end type MPI_tile


contains


  subroutine set_loc(g)
    type(MPI_tile):: g
    g%i = mod(g%rank, g%npx)
    g%j = g%rank/g%npx
  end subroutine set_loc

  function get_rank_at(g, i, j) result(rank)
    type(MPI_tile):: g
    integer::i,j,rank

    if(g%topology == xperio .or. g%topology == xyperio .or. g%topology == xyzperio) i = mod(i+g%npx,g%npx)
    if(g%topology == yperio .or. g%topology == xyperio .or. g%topology == xyzperio) j = mod(j+g%npy,g%npy)
    !print*,">>>",g%rank,i,j
    if((0<=i).and.(i<g%npx).and.(0<=j).and.(j<g%npy))then
       rank = j*g%npx+i
    else
       rank = -1
    end if
  end function get_rank_at

  subroutine set_neighbours(g, nghb)
    type(MPI_tile):: g
    type(MPI_neighbours), intent(inout):: nghb
    integer:: k,di,dj,count,rank

    rank = 0
    count = 0
    k = 1
    do dj=-1,1
       do di=-1,1
          if(di==0 .and. dj==0) then
             rank = -1
          else
             rank = get_rank_at(g, g%i+di*g%incrx, g%j+dj*g%incry)
          endif

          if(rank/=-1)count = count+1

          nghb%infos(1,k) = di
          nghb%infos(2,k) = dj
          nghb%infos(3,k) = rank
          k=k+1
       enddo
    enddo
    nghb%rank = g%rank
    nghb%n = count

  end subroutine set_neighbours

  subroutine print_neighbours(nghb)
    type(MPI_neighbours), intent(in):: nghb
    integer:: i

    print'("rank=",(i2)," has ",(i2)," neighbours")',nghb%rank,nghb%n
    do i=1,9
       if(nghb%infos(3,i)/=-1) print*,nghb%rank, nghb%infos(:,i)
    enddo
  end subroutine print_neighbours

end module mod_tiles
