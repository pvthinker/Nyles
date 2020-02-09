include "interpolate_tracer.f90"
!----------------------------------------
subroutine upwind(trac, u, dtrac, order, l, m, n)
  !
  ! compute dtrac = -delta[ trac*U ]
  ! where delta[ ] is the finite difference in the i direction
  ! which is the third entry of the 3D array
  ! U is the contravariant component  
  ! trac is upwinded at u-point
  !
  ! note that with arbitrary coordinates we should have
  ! dtrac = -delta[ trac*U*vol ] / vol
  !
  ! in Cartesian coordinates, vol is uniform, so vol
  ! can be drop out
  !
  implicit none

  integer, intent(in):: order, l, m, n
  real*8, dimension(l, m, n), intent(in) :: trac, u
  real*8, dimension(l, m, n), intent(inout) :: dtrac

  !f2py intent(inplace):: trac, dtrac, u

  integer:: i, j, k
  real*8::c1,c2,c3, b1, b2, b3, b4, b5
  real*8:: UU, fx, fxm
  real*8,dimension(n) :: phi, qm, qp, um, up


  do k = 1, l
     do j = 1, m
        do i = 1, n
           UU = abs(u(k,j,i))
           up(i) = 0.5*(u(k,j,i)+UU) ! right-going flux
           um(i) = 0.5*(u(k,j,i)-UU) ! left-going flux
           phi(i) = trac(k,j,i)
        enddo
        !
        call interpolate(phi, qp, qm, order, n)
        !
        if (mod(order, 2).eq.0) then
           fxm = 0.
           do i=1,n-1
              fx = u(k,j,i)*qp(i)
              dtrac(k, j, i) = dtrac(k, j, i) + fxm - fx
              fxm = fx
           enddo
           i=n
           fx = 0.
           dtrac(k, j, i) = dtrac(k, j, i) + fxm - fx
        else
           fxm = 0.
           do i=1,n-1
              fx = up(i)*qp(i) + um(i)*qm(i+1)
              dtrac(k, j, i) = dtrac(k, j, i) + fxm - fx
              fxm = fx
           enddo
           i=n
           fx = 0.
           dtrac(k, j, i) = dtrac(k, j, i) + fxm - fx
        endif
     enddo
  enddo

end subroutine upwind
