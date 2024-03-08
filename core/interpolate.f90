include "weno.f90"

subroutine interpolate(vU, qp, qm, order, n)

  implicit none

  integer, intent(in):: n, order

  real*8, dimension(n), intent(in) :: vU
  real*8, dimension(n), intent(out) :: qm
  real*8, dimension(0:n-1), intent(out) :: qp

  !f2py intent(inplace):: vU, qm, qp

  integer:: i

  real*8:: c1,c2,c3, b1, b2, b3, b4, b5, d1, e1,e2

  ! second order
  d1 = 0.5
  ! third order
  c1 = -1./6.
  c2 =  5./6.
  c3 =  2./6.
  ! fourth order
  e1 = -1./12
  e2 = 7./12.
  ! fifth order
  b1 = 2./60.
  b2 = -13./60.
  b3 = 47./60.
  b4 = 27./60.
  b5 = -3./60.

  if (order.eq.5) then
     !
     qp(0) = 0.
     !
     i=1
     qp(i) = vU(i)
     qm(i) = vU(i)
     !
     i=2
     qp(i) = c1*vU(i-1) + c2*vU(i) + c3*vU(i+1)
     qm(i) = c3*vU(i-1) + c2*vU(i) + c1*vU(i+1)
     !
     do i=3,n-3
        qp(i) = b1*vU(i-2) + b2*vU(i-1) + b3*vU(i) + b4*vU(i+1) + b5*vU(i+2)
        qm(i) = b5*vU(i-2) + b4*vU(i-1) + b3*vU(i) + b2*vU(i+1) + b1*vU(i+2)
     enddo
     !
     i=n-2
     qp(i) = c1*vU(i-1) + c2*vU(i) + c3*vU(i+1)
     qm(i) = c3*vU(i-1) + c2*vU(i) + c1*vU(i+1)
     !
     i=n-1
     qp(i) = vU(i)
     qm(i) = vU(i)
     !
     qm(n) = vU(n)
     !
  elseif (order.eq.3) then
     !
     qp(0) = 0.
     !
     i=1
     qp(i) = vU(i)
     qm(i) = vU(i)
     !
     do i=2,n-2
        qp(i) = c1*vU(i-1) + c2*vU(i) + c3*vU(i+1)
        qm(i) = c3*vU(i-1) + c2*vU(i) + c1*vU(i+1)
     enddo
     !
     i=n-1
     qp(i) = vU(i)
     qm(i) = vU(i)
     !
     qm(n) = vu(n)
     !
  !1st order
  elseif (order.eq.1) then
     !
     qp(0) = 0.
     !
     do i=1,n-1
        qp(i) = vU(i)
        qm(i) = vU(i)
     enddo
     !
     qm(n) = vU(n)
     !
  ! 2nd order
  elseif (order.eq.2) then
     qm(1)=0.5*vu(1)
     do i=2,n
        qm(i)=0.5*(vU(i-1)+vU(i))
     enddo
  ! 4th order
  elseif (order.eq.4) then
     i=1
     qm(i)=e2*(vU(i))+e1*(vU(i+1))
     i=2
     qm(i)=e2*(vU(i-1)+vU(i))+e1*(vU(i+1))
     do i=3,n-1
        qm(i)=e2*(vU(i-1)+vU(i))+e1*(vU(i-2)+vU(i+1))
     enddo
     i=n
     qm(i)=e2*(vU(i-1)+vU(i))+e1*(vU(i-2))
  endif

end subroutine interpolate
