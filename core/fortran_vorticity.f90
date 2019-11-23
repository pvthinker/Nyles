!----------------------------------------
subroutine vorticity(ui, uj, wk, l, m, n)
  !
  ! omega_k = delta_i[u_j] - delta_j[u_i]
  !
  ! direction i should be the first one,
  ! imposing to view arrays in the j direction
  !
  implicit none

  integer, intent(in):: l, m, n
  real*8, dimension(l, m, n), intent(in) :: ui, uj
  real*8, dimension(l, m, n), intent(inout) :: wk

  !f2py intent(inplace):: ui, uj, wk

  integer:: i, j, k
  
  do k = 1, l
     do j = 1, m-1
        do i = 1, n-1
           wk(k,j,i) = uj(k,j,i+1) - uj(k,j,i) - ui(k,j+1,i) + ui(k,j,i)
        enddo
     enddo
  enddo
  
end subroutine vorticity

!----------------------------------------
subroutine vorticity_all_comp(ui, uj, uk, wi, wj, wk, l, m, n)
  !
  ! compute the three components with all arrays in the
  ! same convention (k, j, i)
  !
  implicit none

  integer, intent(in):: l, m, n
  real*8, dimension(l, m, n), intent(in) :: ui, uj, uk
  real*8, dimension(l, m, n), intent(inout) :: wi, wj, wk

  !f2py intent(inplace):: ui, uj, uk, wi, wj, wk

  integer:: i, j, k
  
  do k = 1, l-1
     do j = 1, m-1
        do i = 1, n
           wi(k,j,i) = uk(k,j+1,i) - uk(k,j,i) - uj(k+1,j,i) + uj(k,j,i)
        enddo
     enddo
  enddo
  do k = 1, l-1
     do j = 1, m
        do i = 1, n-1
           wj(k,j,i) = ui(k+1,j,i) - ui(k,j,i) - uk(k,j,i+1) + uk(k,j,i)
        enddo
     enddo
  enddo
  do k = 1, l
     do j = 1, m-1
        do i = 1, n-1
           wk(k,j,i) = uj(k,j,i+1) - uj(k,j,i) - ui(k,j+1,i) + ui(k,j,i)
        enddo
     enddo
  enddo
  
end subroutine vorticity_all_comp

