include "interpolate.f90"
!----------------------------------------
subroutine kin(u, v, ke, ds2, order, l, m, n)

  ! compute the kinetic energy of the component associated with
  ! the direction of the last entry (the one contiguous in memory)
  !
  ! algo: "square first, average then"
  ! compute ke at velocity point (say U) and average at T-point (cell center)
  !
  ! convention u[k, j, i] is the velocity on the right face of the cell [k, j, i]
  !
  !
  ! ds2 is the diagonal metric term 1/dx**2 or 1/dy**2 or 1/dz**2
   
  implicit none

  integer, intent(in):: l, m, n, order
  real*8 :: ds2
  real*8, dimension(l, m, n), intent(in) :: u, v
  real*8, dimension(l, m, n), intent(inout) :: ke 

  !f2py intent(inplace):: u, v, ke

  real*8 :: u2, u22, cff, cff2, UU, up, um
  real*8,dimension(n) :: zke, kem
  real*8,dimension(0:n-1) :: kep
  integer:: i, j, k
  real*8::c1,c2,c3, b1, b2, b3, b4, b5, d1, e1,e2

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

  cff2 = 0.5*ds2

  do k = 1, l
     do j = 1, m
        do i = 2, n
           ke(k,j,i) = ke(k,j,i)+cff2*0.5*(u(k,j,i)**2+u(k,j,i-1)**2)
        enddo
     enddo
  enddo
end subroutine kin
