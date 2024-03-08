module mg_enums

  implicit none

  enum, bind(c)
     enumerator :: closed=1,xperio,yperio,zperio,xyperio,xyzperio
  endenum

  type MG_Grid
     integer:: nx,ny,nz,nh,n
     integer:: i0,i1,j0,j1,k0,k1
     integer:: topology=closed
  end type MG_Grid

end module mg_enums
