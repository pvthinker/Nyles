module pytools

  use mg_types
  implicit none

contains

  subroutine get_pyshape(mg, lev, gshape)
    type(MGinfos):: mg
    integer:: lev
    integer,intent(out):: gshape(3)
    gshape(3) = mg%hierarchy%grid(lev)%nx+2*mg%param%nh
    gshape(2) = mg%hierarchy%grid(lev)%ny+2*mg%param%nh
    gshape(1) = mg%hierarchy%grid(lev)%nz
  end subroutine get_pyshape

  subroutine get_pyarray(mg, lev, ivar, n1,n2,n3,x)
    type(MGinfos):: mg
    integer:: lev, ivar
    integer::n1,n2,n3
    REAL,dimension(n1,n2,n3),intent(inout):: x
    if(ivar==1)then
       x = mg%data(lev)%x
    elseif(ivar==2)then
       x = mg%data(lev)%b
    elseif(ivar==3)then
       x = mg%data(lev)%r
    elseif(ivar==4)then
       x = mg%data(lev)%y
    elseif(ivar==5)then
       x = mg%oper(lev)%diag
    elseif(ivar==6)then
       x = mg%oper(lev)%idiag
    elseif(ivar==7)then
       x = mg%oper(lev)%msk
    elseif(ivar==8)then
       x = mg%oper(lev)%Rcoef
    elseif(ivar==9)then
       x = mg%oper(lev)%Pcoef
    end if

  end subroutine get_pyarray

  subroutine set_pyarray(mg, lev, ivar, n1,n2,n3,x)
    type(MGinfos):: mg
    integer:: lev, ivar
    integer::n1,n2,n3
    REAL,dimension(n1,n2,n3),intent(in):: x
    if(ivar==1)then
       mg%data(lev)%x = x
    elseif(ivar==2)then
       mg%data(lev)%b = x
    elseif(ivar==3)then
       mg%data(lev)%r = x
    elseif(ivar==4)then
       mg%data(lev)%y = x
    elseif(ivar==7)then
       mg%oper(lev)%msk = x
    end if
  end subroutine set_pyarray

end module pytools
