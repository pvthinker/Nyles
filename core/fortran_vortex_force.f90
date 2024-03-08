include "interpolate.f90"
!-----------------------------------------
!
! compute du_i/dt += epsilon * w_j * U^k (epsilon=-1)
! arrays are in [j, i, k] convention
! U^k is averaged at i point
! w_j upwinded in the k direction
!

subroutine vortex_force_direc(U, vort, res, order, m, n, l)

  implicit none

  integer, intent(in):: l, m, n, order
  real*8, dimension(m, n, l), intent(in) :: U, vort
  real*8, dimension(m, n, l), intent(inout) :: res

  !f2py intent(inplace):: vort, res, U

  integer:: i, j, k
  real*8:: c1,c2,c3, b1, b2, b3, b4, b5
  real*8:: UU, UU_0, UU_1, U_omega
  real*8,dimension(l) :: qm
  real*8,dimension(l) :: q, flux, u1d
  real*8,dimension(0:l-1):: qp
  real*8,dimension(l) :: vU, U_interp, um, up

  logical:: linear=.false.

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
        if (linear) then
           UU_0 = 0.
           do k=1,l
              UU_1 = 0.5 * ( U(j,i,k) + U(j,i+1,k))
              vU(k) =  vort(j,i,k)!*UU_1
              U_interp(k) = 0.5*(UU_0+UU_1)
              UU = abs(U_interp(k))
              up(k) = 0.5*(U_interp(k)+UU) ! right-going flux
              um(k) = 0.5*(U_interp(k)-UU) ! left-going flux
              UU_0 = UU_1
           enddo
           !
           call interpolate(vU, qp, qm, order, l)
           !
           if (mod(order, 2).eq.0) then
              do k=1,l
                 res(j,i,k) = res(j,i,k)-qm(k)
              enddo
           else
              do k=1,l
                 res(j,i,k) = res(j,i,k)-qp(k-1)*up(k)-qm(k)*um(k)
              enddo
           endif
        else
           UU_0 = 0.
           do k=1,l
              UU_1 = 0.5 * ( U(j,i,k) + U(j,i+1,k))
              u1d(k) = 0.5*(UU_0+UU_1)
              UU_0 = UU_1
           enddo
           q(1)=0.
           do k=2,l
              q(k) = vort(j,i,k-1)
           enddo
           call flux1d(u1d, q, flux, l)
           do k=1,l
              res(j,i,k) = res(j,i,k)-flux(k)!qp(k-1)*up(k)-qm(k)*um(k)
           enddo
        endif
    end do
  end do

end subroutine vortex_force_direc

subroutine vortex_force_flip(U, vort, res, order, m, n, l)
  ! directions 1 and 2 are swapped compared to above
  ! the force is counted with a +
  implicit none

  integer, intent(in):: l, m, n, order
  real*8, dimension(m, n, l), intent(in) :: U, vort
  real*8, dimension(m, n, l), intent(inout) :: res

  !f2py intent(inplace):: vort, res, U

  integer:: i, j, k
  real*8:: c1,c2,c3, b1, b2, b3, b4, b5
  real*8:: UU, UU_0, UU_1, U_omega
  real*8,dimension(n) :: q, flux, u1d
  real*8,dimension(n):: qm
  real*8,dimension(0:n-1):: qp
  real*8,dimension(n) :: vU, U_interp, up, um

  logical:: linear=.false.


  c1 = -1./6.
  c2 =  5./6.
  c3 =  2./6.

  b1 = 2./60.
  b2 = -13./60.
  b3 = 47./60.
  b4 = 27./60.
  b5 = -3./60.

  do j = 1,m
     do k=1,l-1
        if (linear) then
           UU_0 = 0.
           do i=1,n
              UU_1 = 0.5 * ( U(j,i,k) + U(j,i,k+1))
              vU(i) =  vort(j,i,k)!*UU_1
              U_interp(i) = 0.5*(UU_0+UU_1)
              UU = abs(U_interp(i))
              up(i) = 0.5*(U_interp(i)+UU) ! right-going flux
              um(i) = 0.5*(U_interp(i)-UU) ! left-going flux
              UU_0 = UU_1
           enddo!
           !
           call interpolate(vU, qp, qm, order, n)

           if (mod(order, 2).eq.0) then
              do i=1,n
                 res(j,i,k) = res(j,i,k)+qm(i)
              enddo
           else
              do i=1,n
                 res(j,i,k) = res(j,i,k)+qp(i-1)*up(i)+qm(i)*um(i)
              enddo
           endif

        else ! WENO
           UU_0 = 0.
           do i=1,n
              UU_1 = 0.5 * ( U(j,i,k) + U(j,i,k+1))
              u1d(i) = 0.5*(UU_0+UU_1)
              UU_0 = UU_1
           enddo
           q(1)=0.
           do i=2,n
              q(i) = vort(j,i-1,k)
           enddo
           call flux1d(u1d, q, flux, n)
           do i=1,n
              res(j,i,k) = res(j,i,k)+flux(i)!qp(i-1)*up(i)+qm(i)*um(i)
           enddo
        endif

     enddo
  enddo

end subroutine vortex_force_flip
