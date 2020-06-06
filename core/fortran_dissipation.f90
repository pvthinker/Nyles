!----------------------------------------
subroutine add_laplacian(phi, dphi, coef, l, m, n)
  !
  ! dphi += delta_i[ coef*delta_i[phi] ]
  !
  ! direction i should be the first one
  !
  implicit none

  integer, intent(in):: l, m, n
  real*8, dimension(l, m, n), intent(in) :: phi
  real*8, dimension(l, m, n), intent(inout) :: dphi
  real*8 :: coef, fx, fxm

  !f2py intent(inplace):: phi, dphi

  integer:: i, j, k

  do k = 1, l
     do j = 1, m
        ! assume no flux at left boundary
        fxm = 0.
        do i = 1, n-1
           fx = phi(k,j,i+1)-phi(k,j,i)
           dphi(k,j,i) = dphi(k,j,i) +  coef*(fx-fxm)
           fxm = fx
        enddo
        ! assume no flux at right boundary
        fx = 0.
        i = n
        dphi(k,j,i) = dphi(k,j,i) +  coef*(fx-fxm)
     enddo
  enddo
  
end subroutine add_laplacian

!----------------------------------------
subroutine add_laplacian_uxx(phi, dphi, coef, l, m, n)
  !
  ! dphi += delta_i[ coef*delta_i[phi] ]
  !
  ! direction i should be the first one
  !
  ! special case for velocity 'u' along 'x'
  !
  implicit none

  integer, intent(in):: l, m, n
  real*8, dimension(l, m, n), intent(in) :: phi
  real*8, dimension(l, m, n), intent(inout) :: dphi
  real*8 :: coef, fx, fxm

  !f2py intent(inplace):: phi, dphi

  integer:: i, j, k

  do k = 1, l
     do j = 1, m
        ! this is the exception for u_xx
        fxm = phi(k,j,1)
        do i = 1, n-1
           fx = phi(k,j,i+1)-phi(k,j,i)
           dphi(k,j,i) = dphi(k,j,i) +  coef*(fx-fxm)
           fxm = fx
        enddo
        ! assume no flux at right boundary
        !fx = 0.
        !i = n
        !dphi(k,j,i) = dphi(k,j,i) +  coef*(fx-fxm)
     enddo
  enddo

end subroutine add_laplacian_uxx

