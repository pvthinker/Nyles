module mg_types
  use types
  use mg_enums
  use mod_halo
  use mod_tiles
  use gluesplit
  USE ISO_C_BINDING

  implicit none

  integer :: nmaxlevels
  parameter(nmaxlevels=50)


  type, bind(C) :: MG_Param
     integer:: npre=1
     integer:: npost=1
     integer:: nexact=20
     integer:: nh=3
     integer:: maxite=20
     REAL:: tol=1e-6
     REAL:: omega=0.9
     integer:: debug=0
     logical:: vertices=.false.
     logical:: is2d=.true.
  end type MG_Param

  type ::  MG_Data
     REAL,dimension(:,:,:),allocatable:: x,r,b,y
     real,dimension(:,:,:),allocatable:: small
  end type MG_Data

  type MG_Operators
     REAL,dimension(:,:,:),allocatable:: diag, idiag, Rcoef, Pcoef
     integer,dimension(:,:,:),allocatable::msk, smallmsk
  end type MG_Operators

  type MG_Hierarchy
     integer:: nlevels
     type(MG_Grid),dimension(nmaxlevels):: grid
     type(MG_Grid),dimension(nmaxlevels):: smallgrid
     type(MPI_tile),dimension(nmaxlevels):: tile
     type(T_glue),dimension(nmaxlevels):: infos
  end type MG_Hierarchy

  type MG_Stats
     integer:: nite
     REAL:: res
  end type MG_Stats

  type MGinfos
     integer:: npx, npy
     integer:: nlevels
     integer:: ndofs
     integer:: verbose=1
     type(MG_Param):: param
     type(MG_Stats):: stats
     type(MG_Hierarchy):: hierarchy
     type(MG_Data),dimension(:), allocatable:: data
     type(MG_Operators),dimension(:), allocatable :: oper
     type(MPI_halo), dimension(:), allocatable:: halo
  end type MGinfos

end module mg_types
