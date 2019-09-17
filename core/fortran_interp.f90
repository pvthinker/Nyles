! f2py -m fortran_interp -c fortran_interp.f90 --opt='-O3 -fPIC'
!----------------------------------------
  subroutine adv_upwind(u,v,du,dv,nh,m,n)

    implicit none

    integer,intent(in):: n,m,nh
    real*8,dimension(m,n) :: u,v,du,dv

!f2py intent(inplace)::u,v,du,dv
    
    integer::i,j

    real*8:: idx2, idy2
    real*8:: omega, omega_i, omega_j, vu, uv

    idx2 = 1.
    idy2 = 1.
    
    ! main and single double do-loop
    do j=nh,m-nh
       do i=nh,n-nh
          ! vorticity
          omega = (v(j, i+1)-v(j, i))-(u(j+1, i)-u(j, i))
          ! v at u-point ! v is the contravariant component
          vu = 0.25*(v(j, i)+v(j, i+1)+v(j-1, i)+v(j-1, i+1))
          ! u at v-point ! u is the contravariant component
          uv = 0.25*(u(j, i)+u(j, i-1)+u(j+1, i)+u(j+1, i-1))
          omega_j = omega ! upwinded interpolation at u point along j
          omega_i = omega ! upwinded interpolation at v point along i

          du(j, i) = du(j, i) + vu * omega_j * idy2
          dv(j, i) = dv(j, i) - uv * omega_i * idx2
       enddo
    enddo
    
  end subroutine adv_upwind

