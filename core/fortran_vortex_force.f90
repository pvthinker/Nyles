!-----------------------------------------
!
! compute du_i/dt += epsilon * w_j * U^k (epsilon=-1)
! arrays are in [j, i, k] convention
! U^k is averaged at i point
! w_j upwinded in the k direction
! 

subroutine vortex_force_calc(U, vort, res, epsilon, order, m, n, l)

  implicit none

  integer, intent(in):: l, m, n, epsilon, order
  real*8, dimension(m, n, l), intent(in) :: U, vort
  real*8, dimension(m, n, l), intent(inout) :: res

  !f2py intent(inplace):: vort, res, U

  integer:: i, j, k
  real*8:: c1,c2,c3, b1, b2, b3, b4, b5
  real*8:: UU, up, um, qm, qp, U_interp, U_omega

  c1 = -1./6.
  c2 =  5./6.
  c3 =  2./6.

  b1 = 2./60.
  b2 = -13./60.
  b3 = 47./60.
  b4 = 27./60.
  b5 = -3./60.


  do j = 1,m
     do i = 1,n-1 ! because u_i[:,:,-1] = 0. (BC)
        do k = 1,l
          if (k.gt.1) then
             U_interp = 0.25 * ( U(j,i,k  ) + U(j,i+1,k) &
                               + U(j,i,k-1) + U(j,i+1,k-1))
          else
            U_interp = 0.25 * ( U(j,i,k) + U(j,i+1,k)) !no normal velocity
          endif

          UU = abs(U_interp)
          up = 0.5*(U_interp+UU) ! right-going flux
          um = 0.5*(U_interp-UU) ! left-going flux

          if ((k.gt.3).and.(k.lt.(l-1)).and.(order.eq.5)) then
             qp = b1*vort(j,i,k-3) + b2*vort(j,i,k-2) + b3*vort(j,i,k-1) + b4*vort(j,i,k) + b5*vort(j,i,k+1)
          elseif ((k.gt.2).and.(k.lt.(l)).and.(order.ge.3)) then
             ! 3rd order upwind
             qp = c1*vort(j,i,k-2) + c2*vort(j,i,k-1) + c3*vort(j,i,k)
          elseif (k.gt.1) then
             ! 1st order
             qp = vort(j,i,k-1)
          else
             qp = 0.!vort(j,i,k) ! <- very hazardous, CHECK!!!!
          endif


          if ((k.gt.2).and.(k.lt.(l-2)).and.(order.eq.5)) then
             qm = b5*vort(j,i,k-2) + b4*vort(j,i,k-1) + b3*vort(j,i,k) + b2*vort(j,i,k+1) + b1*vort(j,i,k+2)
          elseif ((k.gt.1).and.(k.lt.(l-1)).and.(order.ge.3)) then
             ! 3rd order upwind
             qm = c3*vort(j,i,k-1) + c2*vort(j,i,k) + c1*vort(j,i,k+1)
          elseif (k.le.(l-1)) then
             ! 1st order
             qm = vort(j,i,k)
          else
             ! no flux through the boundary
             ! TODO: change this in case of inflow/outflow
             ! cf Winters 2012
             qm = 0.!vort(j,i,k) ! <- very hazardous, CHECK!!!!
          endif

          ! this term is U_j * omega_k
          U_omega = up*qp + um*qm

          res(j,i,k) = res(j,i,k) - epsilon * U_omega

      end do
    end do
  end do

end subroutine vortex_force_calc
