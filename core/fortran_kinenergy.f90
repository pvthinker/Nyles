!----------------------------------------
subroutine kin(u, v, ke, ds2, iflag, l, m, n)

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

  integer, intent(in):: l, m, n, iflag
  real*8 :: ds2
  real*8, dimension(l, m, n), intent(in) :: u, v
  real*8, dimension(l, m, n), intent(inout) :: ke 

  !f2py intent(inplace):: u, v, ke

  real*8 :: u2, u22, cff
  integer:: i, j, k

  cff = 0.25*ds2

  if (iflag.eq.1) then ! overwrite ke
     do k = 1, l
        do j = 1, m
           u2 = cff*u(k, j, 1)*v(k, j, 1)
           ke(k, j, 1) = u2 
           do i = 2, n
              u22 = cff*u(k, j, i)*v(k, j, i)
              ke(k, j, i) = u2 + u22
              u2 = u22
           enddo
        enddo
     enddo
  else ! add to ke
     do k = 1, l
        do j = 1, m
           u2 = cff*u(k, j, 1)*v(k, j, 1)
           ke(k, j, 1) = ke(k, j, 1) + u2
           do i = 2, n
              u22 = cff*u(k, j, i)*v(k, j, i)
              ke(k, j, i) = ke(k, j, i) + u2 + u22
              u2 = u22
           enddo
        enddo
     enddo
  endif
end subroutine kin
