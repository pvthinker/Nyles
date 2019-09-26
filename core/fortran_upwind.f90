!----------------------------------------
subroutine upwind(trac, u, dtrac, vol, iflag, l, m, n)
  !
  ! compute dtrac = delta[ vol*trac*u ]
  ! where delta[ ] is the finite difference in the 'u' direction
  ! the 'u' direction is the third entry of the 3D array
  ! trac is upwinded at U point
  !
  implicit none

  integer, intent(in):: iflag, l, m, n
  real*8 :: vol
  real*8, dimension(l, m, n), intent(in) :: trac, u
  real*8, dimension(l, m, n), intent(inout) :: dtrac

  !f2py intent(inplace):: trac, dtrac, u

  integer:: i, j, k
  real*8::c1,c2,c3
  real*8:: UU, up, um, qm, qp, fx, fxm

  c1 = -1./6.
  c2 =  5./6.
  c3 =  2./6.

  do k = 1, l
     do j = 1, m
        fxm = 0.
        do i = 1, n-1
           UU = abs(u(k,j,i))
           up = 0.5*(u(k,j,i)+UU) ! right-going flux
           um = 0.5*(u(k,j,i)-UU) ! left-going flux
           if (i.gt.1) then
              ! 3rd order upwind interpolation at U point              
              qp = c1*trac(k,j,i-1)+c2*trac(k,j,i  )+c3*trac(k,j,i+1)
           else
              ! 1st order
              qp = trac(k,j,i)
           endif
           if (i.le.(n-2))then
              ! 3rd order upwind
              qm = c3*trac(k,j,i  )+c2*trac(k,j,i+1)+c1*trac(k,j,i+2)
           elseif (i.le.(n-1))then
              ! 1st order
              qm = trac(k,j,i+1)
           else
              ! no flux through the boundary
              ! TODO: change this in case of inflow/outflow
              ! cf Winters 2012
              qm = 0.
           endif

           fx = vol*(up*qp + um*qm)
           if (iflag.eq.1) then ! overwrite rhs
              dtrac(k, j, i) = fxm - fx
           else
              dtrac(k, j, i) = dtrac(k, j, i) + fxm - fx
           endif
           fxm = fx
        enddo
     enddo
  enddo
  
end subroutine upwind
