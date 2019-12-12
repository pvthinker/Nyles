!----------------------------------------
subroutine upwind(trac, u, dtrac, order, iflag, i0, l, m, n)
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

  integer, intent(in):: order, iflag, i0, l, m, n
  real*8, dimension(l, m, n), intent(in) :: trac, u
  real*8, dimension(l, m, n), intent(inout) :: dtrac

  !f2py intent(inplace):: trac, dtrac, u

  integer:: i, j, k, ii
  real*8::c1,c2,c3, b1, b2, b3, b4, b5
  real*8:: UU, up, um, qm, qp, fx, fxm

  c1 = -1./6.
  c2 =  5./6.
  c3 =  2./6.

  b1 = 2./60.
  b2 = -13./60.
  b3 = 47./60.
  b4 = 27./60.
  b5 = -3./60.

  do k = 1, l
     do j = 1, m
        fxm = 0.
        if (i0.eq.1)then
           dtrac(k, j, 1) = 0.
        endif
        do i = 1+i0, n
           ii = i-i0
           UU = abs(u(k,j,i))
           up = 0.5*(u(k,j,i)+UU) ! right-going flux
           um = 0.5*(u(k,j,i)-UU) ! left-going flux

           if ((ii.gt.2).and.(i.le.(n-2)).and.(order.eq.5)) then
              qp = b1*trac(k,j,i-2) + b2*trac(k,j,i-1) + b3*trac(k,j,i) + b4*trac(k,j,i+1) + b5*trac(k,j,i+2)
           elseif ((ii.gt.1).and.(i.le.(n-1)).and.(order.ge.3)) then
              ! 3rd order upwind interpolation at U point
              qp = c1*trac(k,j,i-1)+c2*trac(k,j,i  )+c3*trac(k,j,i+1)
           else
              ! 1st order
              qp = trac(k,j,i)
           endif


           if ((ii.gt.1).and.(i.le.(n-3)).and.(order.eq.5)) then
              qm = b5*trac(k,j,i-1) + b4*trac(k,j,i) + b3*trac(k,j,i+1) + b2*trac(k,j,i+2) + b1*trac(k,j,i+3)
           elseif ((i.le.(n-2)).and.(order.ge.3)) then
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

           fx = up*qp + um*qm
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
