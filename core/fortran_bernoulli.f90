!----------------------------------------
subroutine gradke(ke, du, l, m, n)
  !
  ! du += delta[ ke ]
  ! where delta[ ] is the finite difference in the 'u' direction
  ! the 'u' direction is the third entry of the 3D array
  !
  implicit none

  integer, intent(in):: l, m, n
  real*8, dimension(l, m, n), intent(in) :: ke
  real*8, dimension(l, m, n), intent(inout) :: du

  !f2py intent(inplace):: ke, du

  integer:: i, j, k
  
  do k = 1, l
     do j = 1, m
        do i = 1, n-1
           du(k,j,i) = du(k,j,i) + ke(k,j,i+1)-ke(k,j,i)
        enddo
     enddo
  enddo
  
end subroutine gradke

!----------------------------------------
subroutine gradkeandb(ke, b, du, dz, l, m, n)
  !
  ! du += delta[ ke ] + b*dz
  ! where delta[ ] is the finite difference in the 'u' direction
  ! the 'u' direction is the third entry of the 3D array
  !
  implicit none

  integer, intent(in):: l, m, n
  real*8:: dz
  real*8, dimension(l, m, n), intent(in) :: ke, b
  real*8, dimension(l, m, n), intent(inout) :: du

  !f2py intent(inplace):: ke, b, du

  integer:: i, j, k
  real*8:: cff

  cff = 0.5*dz
  
  do k = 1, l
     do j = 1, m
        do i = 1, n-1
           du(k,j,i) = du(k,j,i) + ke(k,j,i+1)-ke(k,j,i) &
                + 0.5*(b(k,j,i+1)+b(k,j,i))
        enddo
     enddo
  enddo
  
end subroutine gradkeandb
