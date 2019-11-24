!-----------------------------------------
subroutine vortex_force_calc(U_k, U_j, vort_k, vort_j, res, order, l, m, n)

  implicit none

  integer, intent(in):: l, m, n, order
  real*8, dimension(l, m, n), intent(in) :: U_k, U_j , vort_k, vort_j
  real*8, dimension(l, m, n), intent(inout) :: res

  !f2py intent(inplace):: vort_k, vort_j, res, U_k, U_j

  integer:: i, j, k
  integer, parameter :: out_unit=30
  real*8:: c1,c2,c3, b1, b2, b3, b4, b5
  real*8:: UU_k, up_k, um_k, qm_k, qp_k, U_k_interp, UU_j, up_j, um_j, qm_j, qp_j, U_j_interp, U_omega

  c1 = -1./6.
  c2 =  5./6.
  c3 =  2./6.

  b1 = 2./60.
  b2 = -13./60.
  b3 = 47./60.
  b4 = 27./60.
  b5 = -3./60.


  do k = 1,l
    do j = 1,m
      do i = 1,n

        !w_j U_k
        if ((k.gt.1).and.(j.lt.m)) then
          U_k_interp = 0.25 * ( U_k(k, j, i) + U_k(k-1, j, i) + U_k(k-1, j+1, i) + U_k(k, j+1, i) )
        elseif (j.lt.m) then
          U_k_interp = 0.25 * ( U_k(k, j+1, i) + U_k(k, j, i)) !no normal velocity
        else
        endif

        UU_k = abs(U_k_interp)
        up_k = 0.5*(U_k_interp+UU_k) ! right-going flux
        um_k = 0.5*(U_k_interp-UU_k) ! left-going flux

        !write(out_unit,*) "up_k ", up_k, "um_k", um_k

        if ((i.gt.2).and.(order.eq.5)) then
           qp_j = b1*vort_j(k,j,i-2) + b2*vort_j(k,j,i-1) + b3*vort_j(k,j,i) + b4*vort_j(k,j,i+1) + b5*vort_j(k,j,i+2)
        elseif ((i.gt.1).and.(order.ge.3)) then
           qp_j = c1*vort_j(k,j,i-1) + c2*vort_j(k,j,i) + c3*vort_j(k,j,i+1)
        else
           qp_j = vort_j(k,j,i)
        endif

        if ((i.le.(n-3)).and.(order.eq.5)) then
           qm_j = b5*vort_j(k,j,i-1) + b4*vort_j(k,j,i) + b3*vort_j(k,j,i+1) + b2*vort_j(k,j,i+2) + b1*vort_j(k,j,i+3)
        elseif ((i.le.(n-2)).and.(order.ge.3)) then
           qm_j = c3*vort_j(k,j,i) + c2*vort_j(k,j,i+1) + c1*vort_j(k,j,i+2)
        elseif (i.le.(n-1)) then
           qm_j = vort_j(k,j,i+1)
        else
           ! no flux through the boundary
           ! TODO: change this in case of inflow/outflow
           ! cf Winters 2012
           qm_j = 0.
        endif

        !w_k U_j
        if ((k.gt.1).and.(i.lt.m)) then
          U_j_interp = 0.25 * ( U_j(k, j, i) + U_j(k-1, j, i) + U_j(k-1, j, i+1) + U_j(k, j, i+1) )
        elseif (k.lt.l) then
          U_j_interp = 0.25 * ( U_j(k+1, j, i) + U_j(k, j, i)) !no normal velocity
        else
        endif

        UU_j = abs(U_j_interp)
        up_j = 0.5*(U_j_interp+UU_j) ! right-going flux
        um_j = 0.5*(U_j_interp-UU_j) ! left-going flux

        !write(out_unit,*) "up_j ", up_j, "um_j", um_j

        if ((j.gt.2).and.(order.eq.5)) then
           qp_k = b1*vort_k(k,j-2,i) + b2*vort_k(k,j-1,i) + b3*vort_k(k,j,i) + b4*vort_k(k,j+1,i) + b5*vort_k(k,j+2,i)
        elseif ((j.gt.1).and.(order.ge.3)) then
           qp_k = c1*vort_k(k,j-1,i) + c2*vort_k(k,j,i) + c3*vort_k(k,j+1,i)
        else
           qp_k = vort_k(k,j,i)
        endif

        if ((j.le.(m-3)).and.(order.eq.5)) then
           qm_k = b5*vort_k(k,j-1,i) + b4*vort_k(k,j,i) + b3*vort_k(k,j+1,i) + b2*vort_k(k,j+2,i) + b1*vort_k(k,j+3,i)
        elseif ((j.le.(m-2)).and.(order.ge.3)) then
           qm_k = c3*vort_k(k,j,i) + c2*vort_k(k,j+1,i) + c1*vort_k(k,j+2,i)
        elseif (j.lt.m) then
           qm_k = vort_k(k,j+1,i)
        else
           ! no flux through the boundary
           ! TODO: change this in case of inflow/outflow
           ! cf Winters 2012
           qm_k = 0.
        endif

        !write(out_unit,*), "qm_j", qm_j, "qp_j", qp_j, "qp_k", qp_k, "qm_k", qm_k

        ! this term is U_j * omega_k

        U_omega =  (up_k*qp_j + um_k*qm_j) - ( up_j * qp_k + um_j * qm_k)

        !write(out_unit,*) "U_omega ", U_omega


        res(k,j,i) = res(k,j,i) - U_omega


      end do
    end do
  end do

end subroutine vortex_force_calc
