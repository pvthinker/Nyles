!-----------------------------------------

subroutine vortex_force_calc(U, vorticity, res, epsilon, l, m, n)

  implicit none

  integer, intent(in):: l, m, n, epsilon
  real*8, dimension(l, m, n), intent(in) :: U, vorticity
  real*8, dimension(l, m, n), intent(inout) :: res

  !f2py intent(inplace):: vorticity, res, U

  integer:: i, j, k
  real*8:: c1,c2,c3
  real*8:: UU, up, um, qm, qp, vorticity_interp, U_interp, U_omega

  c1 = -1./6.
  c2 =  5./6.
  c3 =  2./6.

  do k = 1,l
    do j = 1,m
      do i = 1,n

        if (j.gt.1) then

          if (i.lt.n) then
            U_interp = 0.25 * ( U(k, j-1, i) + U(k, j-1, i+1) + U(k, j, i) + U(k, j, i+1) )
          else
            U_interp = 0.25 * ( U(k, j-1, i) + U(k, j, i)) !no normal velocity
          endif

          UU = abs(U_interp)
          up = 0.5*(U_interp+UU) ! right-going flux
          um = 0.5*(U_interp-UU) ! left-going flux

          if (i.gt.1) then
             ! 3rd order upwind interpolation at U point
             qp = c1*vorticity(k,j,i-1) + c2*vorticity(k,j,i) + c3*vorticity(k,j,i+1)
          else
             ! 1st order
             qp = vorticity(k,j,i)
          endif

          if (i.le.(n-2)) then
             ! 3rd order upwind
             qm = c3*vorticity(k,j,i) + c2*vorticity(k,j,i+1) + c1*vorticity(k,j,i+2)
          elseif (i.le.(n-1)) then
             ! 1st order
             qm = vorticity(k,j,i+1)
          else
             ! no flux through the boundary
             ! TODO: change this in case of inflow/outflow
             ! cf Winters 2012
             qm = 0.
          endif

          ! this term is U_j * omega_k
          U_omega = up*qp + um*qm

          res(k,j,i) = res(k,j,i) - epsilon * U_omega

        else

          !no slip BC => do nothing

        endif

      end do
    end do
  end do

end subroutine vortex_force_calc
