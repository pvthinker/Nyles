!----------------------------------------
subroutine var2mg(x, y, idx, l, m, n, ll, mm, nn)
  !
  ! copy model variable 'x' to multigrid variable 'y'
  !
  implicit none

  integer, intent(in):: l, m, n,ll, mm, nn
  real*8, dimension(l, m, n), intent(inout) :: x
  real*8, dimension(ll, mm, nn), intent(inout) :: y
  integer, dimension(6) :: idx

  !f2py intent(inplace):: x, y

  integer:: i, j, k, ii, jj, kk
  
  do k = 1, ll
     kk = k+idx(1)
     do j = 1, mm
        jj = j+idx(3)
        do i = 1, nn
           ii = i+idx(5)
           y(k,j,i) = x(kk,jj,ii)
        enddo
     enddo
  enddo
  
end subroutine var2mg

!----------------------------------------
subroutine mg2var(x, y, idx, l, m, n, ll, mm, nn)
  !
  ! copy multigrid variable 'y' to model variable 'x'
  !
  implicit none

  integer, intent(in):: l, m, n,ll, mm, nn
  real*8, dimension(l, m, n), intent(inout) :: x
  real*8, dimension(ll, mm, nn), intent(inout) :: y
  integer, dimension(6) :: idx

  !f2py intent(inplace):: x, y

  integer:: i, j, k, ii, jj, kk
  
  do k = 1, ll
     kk = k+idx(1)
     do j = 1, mm
        jj = j+idx(3)
        do i = 1, nn
           ii = i+idx(5)
           x(kk,jj,ii) = y(k,j,i)
        enddo
     enddo
  enddo
  
end subroutine mg2var
