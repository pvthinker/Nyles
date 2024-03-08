include "weno.f90"

subroutine interpolate(q, qp, qm, order, n)
  !
  ! Interpolate q from tracer point to U-point for abritary order
  !
  ! - if order%2 == 1 (upwind biased)
  !    qp and qm are the right and left interpolations
  !    qp(i) and qm(i+1) are the value that go with u(i)
  !
  ! - if order%2 == 0 (centered)
  !    only qp is returned, qp takes values from i=1..n-1
  !    qp(i) is the centered value for u(i)
  ! 
  
  implicit none

  integer, intent(in):: n, order

  real*8, dimension(n), intent(in) :: q
  real*8, dimension(n), intent(out) :: qm, qp

  !f2py intent(inplace):: q, qm, qp

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
     i=1
     qp(i) = q(i)
     qm(i) = q(i)
     !
     i=2
     qp(i) = c1*q(i-1) + c2*q(i) + c3*q(i+1)
     qm(i) = c3*q(i-1) + c2*q(i) + c1*q(i+1)
     !
     do i=3,n-2
        qp(i) = b1*q(i-2) + b2*q(i-1) + b3*q(i) + b4*q(i+1) + b5*q(i+2)
        qm(i) = b5*q(i-2) + b4*q(i-1) + b3*q(i) + b2*q(i+1) + b1*q(i+2)
     enddo
     !
     i=n-1
     qp(i) = c1*q(i-1) + c2*q(i) + c3*q(i+1)
     qm(i) = c3*q(i-1) + c2*q(i) + c1*q(i+1)
     !
     i=n
     qp(i) = q(i)
     qm(i) = q(i)
     !
  elseif (order.eq.3) then
     !
     i=1
     qp(i) = q(i)
     qm(i) = q(i)
     !
     do i=2,n-1
        qp(i) = c1*q(i-1) + c2*q(i) + c3*q(i+1)
        qm(i) = c3*q(i-1) + c2*q(i) + c1*q(i+1)
     enddo
     !
     i=n
     qp(i) = q(i)
     qm(i) = q(i)
     !
  !1st order
  elseif (order.eq.1) then
     !
     do i=1,n
        qp(i) = q(i)
        qm(i) = q(i)
     enddo
     !
  ! 2nd order
  elseif (order.eq.2) then
     !
     do i=1,n-1
        qp(i)=0.5*(q(i)+q(i+1))
     enddo
     !
  ! 4th order
  elseif (order.eq.4) then
     !
     i=1
     qp(i)=e2*(q(i)+q(i+1))+e1*(q(i+2))
     !
     do i=2,n-2
        qp(i)=e2*(q(i)+q(i+1))+e1*(q(i-1)+q(i+2))
     enddo
     !
     i=n-1
     qp(i)=e2*(q(i)+q(i+1))+e1*(q(i-1))
     !
  endif

end subroutine interpolate
