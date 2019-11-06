!----------------------------------------
subroutine norm(x, res, k0, k1, j0, j1, i0, i1, l, m, n)
  !
  ! compute the sum of x**2 in the interior domain
  ! (excluding the halo)
  !
  implicit none

  integer, intent(in):: l, m, n
  integer, intent(in):: k0, k1, j0, j1, i0, i1
  real*8, dimension(l, m, n), intent(in) :: x
  real*8, intent(out) :: res

  !f2py intent(inplace):: x

  integer:: i, j, k

  res = 0.
  do k = k0+1, k1
     do j = j0+1, j1
        do i = i0+1, i1
           res = res+ x(k,j,i)*x(k,j,i)
        enddo
     enddo
  enddo
  
end subroutine norm
