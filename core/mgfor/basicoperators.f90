module basicoperators
  use types
  implicit none
contains
  subroutine frestrict_centers2d(xf, xc, coef, nxf,nyf,nzf,nxc,nyc,nzc,nh)
    implicit none
    integer::nh
    integer:: nxf,nyf,nzf
    integer:: nxc,nyc,nzc
    REAL,dimension(1-nh:nxf+nh, 1-nh:nyf+nh, 1:nzf):: xf
    REAL,dimension(1-nh:nxc+nh, 1-nh:nyc+nh, 1:nzc):: xc, coef

    intent(in):: xf,coef
    intent(out):: xc

    integer:: i,j,ic,jc,k

    do k = 1, nzc
       do jc = 1, nyc
          j = 1+(jc-1)*2
          do ic = 1, nxc
             i = 1+(ic-1)*2
             xc(ic,jc,k) = coef(ic,jc,k)*(&
                  xf(i,j  ,k)+xf(i+1,j  ,k)&
                  +xf(i,j+1,k)+xf(i+1,j+1,k)&
                  )
          enddo
       enddo
    enddo
  end subroutine frestrict_centers2d

  subroutine frestrict_centers3d(xf, xc, coef, nxf,nyf,nzf,nxc,nyc,nzc,nh)
    implicit none
    integer::nh
    integer:: nxf,nyf,nzf
    integer:: nxc,nyc,nzc
    REAL,dimension(1-nh:nxf+nh, 1-nh:nyf+nh, 1:nzf):: xf
    REAL,dimension(1-nh:nxc+nh, 1-nh:nyc+nh, 1:nzc):: xc, coef

    intent(in):: xf,coef
    intent(out):: xc

    integer:: i,j,k,ic,jc,kc

    do kc = 1+nh, nzc-nh
       k = 1+nh+(kc-1-nh)*2
       do jc = 1, nyc
          j = 1+(jc-1)*2
          do ic = 1, nxc
             i = 1+(ic-1)*2
             xc(ic,jc,kc) = coef(ic,jc,kc)*(&
                  xf(i,j  ,k)+xf(i+1,j  ,k)&
                  +xf(i,j+1,k)+xf(i+1,j+1,k)&
                  +xf(i,j  ,k+1)+xf(i+1,j  ,k+1)&
                  +xf(i,j+1,k+1)+xf(i+1,j+1,k+1)&
                  )
          enddo
       enddo
    enddo
  end subroutine frestrict_centers3d

  subroutine frestrict_vertices2d(xf, xc, coef, nxf,nyf,nzf,nxc,nyc,nzc,nh)
    implicit none
    integer::nh
    integer:: nxf,nyf,nzf
    integer:: nxc,nyc,nzc
    REAL,dimension(1-nh:nxf+nh, 1-nh:nyf+nh, 1:nzf):: xf
    REAL,dimension(1-nh:nxc+nh, 1-nh:nyc+nh, 1:nzc):: xc, coef

    intent(in):: xf,coef
    intent(out):: xc

    integer:: i,j,ic,jc,k

    do k = 1, nzc
       do jc = 1, nyc
          j = 1+(jc-1)*2
          do ic = 1, nxc
             i = 1+(ic-1)*2
             xc(ic,jc,k) = coef(ic,jc,k)*(&
                  +  xf(i-1,j-1,k)+2*xf(i,j-1,k)+  xf(i+1,j-1,k)&
                  +2*xf(i-1,j  ,k)+4*xf(i,j  ,k)+2*xf(i+1,j  ,k)&
                  +  xf(i-1,j+1,k)+2*xf(i,j+1,k)+  xf(i+1,j+1,k)&
                  )
          enddo
       enddo
    enddo

  end subroutine frestrict_vertices2d

  subroutine frestrict_vertices3d(xf, xc, coef, nxf,nyf,nzf,nxc,nyc,nzc,nh)
    implicit none
    integer::nh
    integer:: nxf,nyf,nzf
    integer:: nxc,nyc,nzc
    REAL,dimension(1-nh:nxf+nh, 1-nh:nyf+nh, 1:nzf):: xf
    REAL,dimension(1-nh:nxc+nh, 1-nh:nyc+nh, 1:nzc):: xc, coef

    intent(in):: xf,coef
    intent(out):: xc

    integer:: i,j,k,ic,jc,kc

    do kc = 1+nh, nzc-nh
       k = 1+nh+(kc-1-nh)*2
       do jc = 1, nyc
          j = 1+(jc-1)*2
          do ic = 1, nxc
             i = 1+(ic-1)*2
             xc(ic,jc,kc) = coef(ic,jc,kc)*(&
                  +  xf(i-1,j-1,k-1)+2*xf(i,j-1,k-1)+  xf(i+1,j-1,k-1)&
                  +2*xf(i-1,j  ,k-1)+4*xf(i,j  ,k-1)+2*xf(i+1,j  ,k-1)&
                  +  xf(i-1,j+1,k-1)+2*xf(i,j+1,k-1)+  xf(i+1,j+1,k-1)&
                  !
                  +2*xf(i-1,j-1,k)+4*xf(i,j-1,k)+2*xf(i+1,j-1,k)&
                  +4*xf(i-1,j  ,k)+8*xf(i,j  ,k)+4*xf(i+1,j  ,k)&
                  +2*xf(i-1,j+1,k)+4*xf(i,j+1,k)+2*xf(i+1,j+1,k)&
                  !
                  +  xf(i-1,j-1,k+1)+2*xf(i,j-1,k+1)+  xf(i+1,j-1,k+1)&
                  +2*xf(i-1,j  ,k+1)+4*xf(i,j  ,k+1)+2*xf(i+1,j  ,k+1)&
                  +  xf(i-1,j+1,k+1)+2*xf(i,j+1,k+1)+  xf(i+1,j+1,k+1)&
                  )
          enddo
       enddo
    enddo

  end subroutine frestrict_vertices3d

  subroutine fprolongation_centers2d(xf, xc, coef, nxf,nyf,nzf,nxc,nyc,nzc,nh)
    implicit none
    integer::nh
    integer:: nxf,nyf,nzf
    integer:: nxc,nyc,nzc

    REAL,dimension(1-nh:nxf+nh, 1-nh:nyf+nh, 1:nzf):: xf, coef
    REAL,dimension(1-nh:nxc+nh, 1-nh:nyc+nh, 1:nzc):: xc
    !REAL,dimension(:,:,:):: xf, coef
    !REAL,dimension(:,:,:):: xc


    intent(in):: xc, coef
    intent(out):: xf

    integer:: i,j,ic,jc,k

    do k = 1, nzf
       do j = 1, nyf-1, 2
          jc = (j+1)/2
          do i = 1, nxf-1, 2
             ic = (i+1)/2
             xf(i,  j  ,k) = xf(i,  j  ,k)+coef(i  ,j  ,k)*(&
                  9*xc(ic,jc  ,k)+3*xc(ic-1,jc  ,k)&
                  +3*xc(ic,jc-1,k)+  xc(ic-1,jc-1,k))

             xf(i+1,j  ,k) = xf(i+1,j  ,k)+coef(i+1,j  ,k)*(&
                  9*xc(ic,jc  ,k)+3*xc(ic+1,jc  ,k)&
                  +3*xc(ic,jc-1,k)+  xc(ic+1,jc-1,k))

             xf(i  ,j+1,k) = xf(i  ,j+1,k)+coef(i  ,j+1,k)*(&
                  9*xc(ic,jc  ,k)+3*xc(ic-1,jc  ,k)&
                  +3*xc(ic,jc+1,k)+  xc(ic-1,jc+1,k))

             xf(i+1,j+1,k) = xf(i+1,j+1,k)+coef(i+1,j+1,k)*(&
                  9*xc(ic,jc  ,k)+3*xc(ic+1,jc  ,k)&
                  +3*xc(ic,jc+1,k)+  xc(ic+1,jc+1,k))

          enddo
       enddo
    enddo

  end subroutine fprolongation_centers2d

  subroutine fprolongation_centers3d(xf, xc, coef, nxf,nyf,nzf,nxc,nyc,nzc,nh)
    implicit none
    integer::nh
    integer:: nxf,nyf,nzf
    integer:: nxc,nyc,nzc

    REAL,dimension(1-nh:nxf+nh, 1-nh:nyf+nh, 1:nzf):: xf, coef
    REAL,dimension(1-nh:nxc+nh, 1-nh:nyc+nh, 1:nzc):: xc
    !REAL,dimension(:,:,:):: xf, coef
    !REAL,dimension(:,:,:):: xc

    intent(in):: xc, coef
    intent(out):: xf

    integer:: i,j,k,ic,jc,kc

    real:: a,b,c

    do k = 1+nh, nzf-nh, 2
       kc = nh+(k+1-nh)/2
       do j = 1, nyf-1, 2
          jc = (j+1)/2
          do i = 1, nxf-1, 2
             ic = (i+1)/2

             a = 9*xc(ic,jc,kc-1)+3*xc(ic-1,jc,kc-1)+3*xc(ic,jc-1,kc-1)+xc(ic-1,jc-1,kc-1)
             b = 9*xc(ic,jc,kc)+3*xc(ic-1,jc,kc)+3*xc(ic,jc-1,kc)+xc(ic-1,jc-1,kc)
             c = 9*xc(ic,jc,kc+1)+3*xc(ic-1,jc,kc+1)+3*xc(ic,jc-1,kc+1)+xc(ic-1,jc-1,kc+1)

             xf(i,  j  ,k  ) = xf(i ,j,k) +coef(i,j,k)*(3*b+a)
             xf(i,  j  ,k+1) = xf(i,j,k+1)+coef(i,j,k+1)*(3*b+c)

             a = 9*xc(ic,jc  ,kc-1)+3*xc(ic+1,jc  ,kc-1)+3*xc(ic,jc-1,kc-1)+  xc(ic+1,jc-1,kc-1)
             b = 9*xc(ic,jc  ,kc)+3*xc(ic+1,jc  ,kc)+3*xc(ic,jc-1,kc)+  xc(ic+1,jc-1,kc)
             c = 9*xc(ic,jc  ,kc+1)+3*xc(ic+1,jc  ,kc+1)+3*xc(ic,jc-1,kc+1)+  xc(ic+1,jc-1,kc+1)

             xf(i+1,j  ,k)   = xf(i+1,j  ,k)  +coef(i+1,j  ,k)*(3*b+a)
             xf(i+1,j  ,k+1) = xf(i+1,j  ,k+1)+coef(i+1,j  ,k+1)*(3*b+c)


             a = 9*xc(ic,jc  ,kc-1)+3*xc(ic-1,jc  ,kc-1)+3*xc(ic,jc+1,kc-1)+  xc(ic-1,jc+1,kc-1)
             b = 9*xc(ic,jc  ,kc)+3*xc(ic-1,jc  ,kc)+3*xc(ic,jc+1,kc)+  xc(ic-1,jc+1,kc)
             c = 9*xc(ic,jc  ,kc+1)+3*xc(ic-1,jc  ,kc+1)+3*xc(ic,jc+1,kc+1)+  xc(ic-1,jc+1,kc+1)

             xf(i  ,j+1,k)   = xf(i  ,j+1,k)  +coef(i  ,j+1,k)*(3*b+a)
             xf(i  ,j+1,k+1) = xf(i  ,j+1,k+1)+coef(i  ,j+1,k+1)*(3*b+c)

             a = 9*xc(ic,jc  ,kc-1)+3*xc(ic+1,jc  ,kc-1)+3*xc(ic,jc+1,kc-1)+  xc(ic+1,jc+1,kc-1)
             b = 9*xc(ic,jc  ,kc)+3*xc(ic+1,jc  ,kc)+3*xc(ic,jc+1,kc)+  xc(ic+1,jc+1,kc)
             c = 9*xc(ic,jc  ,kc+1)+3*xc(ic+1,jc  ,kc+1)+3*xc(ic,jc+1,kc+1)+  xc(ic+1,jc+1,kc+1)

             xf(i+1,j+1,k)   = xf(i+1,j+1,k)  +coef(i+1,j+1,k)*(3*b+a)
             xf(i+1,j+1,k+1) = xf(i+1,j+1,k+1)+coef(i+1,j+1,k+1)*(3*b+c)

          enddo
       enddo
    enddo

  end subroutine fprolongation_centers3d

  subroutine fprolongation_vertices2d(xf, xc, coef, nxf,nyf,nzf,nxc,nyc,nzc,nh)
    implicit none

    integer::nh
    integer:: nxf,nyf,nzf
    integer:: nxc,nyc,nzc

    REAL,dimension(1-nh:nxf+nh, 1-nh:nyf+nh, 1:nzf):: xf, coef
    REAL,dimension(1-nh:nxc+nh, 1-nh:nyc+nh, 1:nzc):: xc
    !REAL,dimension(:,:,:):: xf, coef
    !REAL,dimension(:,:,:):: xc


    intent(in):: xc, coef
    intent(out):: xf

    integer:: i,j,ic,jc,k

    do k = 1, nzf
       do j = 1, nyf-1, 2
          jc = (j+1)/2
          do i = 1, nxf-1, 2
             ic = (i+1)/2
             xf(i,  j  ,k) = xf(i,  j  ,k)&
                  +coef(i  ,j  ,k)* xc(ic,jc  ,k)
             xf(i+1,j  ,k) = xf(i+1,j  ,k)&
                  +coef(i+1,j  ,k)*(xc(ic,jc  ,k)+xc(ic+1,jc,k))
          enddo
          do i = 1, nxf-1, 2
             ic = (i+1)/2
             xf(i  ,j+1,k) = xf(i  ,j+1,k)&
                  +coef(i  ,j+1,k)*(xc(ic,jc  ,k)+xc(ic  ,jc+1,k))
             xf(i+1,j+1,k) = xf(i+1,j+1,k)&
                  +coef(i+1,j+1,k)*(xc(ic,jc  ,k)+xc(ic+1,jc  ,k)&
                  +xc(ic,jc+1,k)+xc(ic+1,jc+1,k))
          enddo
       enddo
    enddo

  end subroutine fprolongation_vertices2d



  subroutine fresidual2d(x,b,r,msk,diag,nx,ny,nz,nh)
    integer::nh
    integer:: nx,ny,nz
    REAL,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz):: b,x,r,diag
    integer,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz):: msk
    intent(out):: r
    intent(in)::x,b,msk,diag

    integer:: i,j,k

    do k=1,nz
       do j=1,ny
          do i=1,nx
             r(i,j,k) = msk(i,j,k)*(b(i,j,k)&
                  +diag(i,j,k)*x(i,j,k)&
                  -(x(i-1,j,k)+x(i+1,j,k)+x(i,j-1,k)+x(i,j+1,k)))

          enddo
       enddo
    enddo

  end subroutine fresidual2d


  subroutine fresidual3d(x,b,r,msk,diag,nx,ny,nz,nh)
    integer::nh
    integer:: nx,ny,nz
    REAL,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz):: b,x,r,diag
    integer,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz):: msk
    intent(out):: r
    intent(in)::x,b,msk,diag

    integer:: i,j,k

    do k=1+nh,nz-nh
       do j=1,ny
          do i=1,nx
             r(i,j,k) = msk(i,j,k)*(b(i,j,k)&
                  +diag(i,j,k)*x(i,j,k)&
                  -(x(i-1,j,k)+x(i+1,j,k)&
                  +x(i,j-1,k)+x(i,j+1,k)&
                  +x(i,j,k-1)+x(i,j,k+1)))

          enddo
       enddo
    enddo

  end subroutine fresidual3d

  subroutine fsmoother2d(x,b,y,idiag,omega,nx,ny,nz,nh)
    integer::nh
    integer:: nx,ny,nz
    REAL,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::x,y,b,idiag
    !integer,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::msk

    real::omega

    intent(inout)::x,y
    intent(in)::b,idiag!,msk

    integer:: i,j,k
    real::cff1

    cff1 = 1.-omega

    do k = 1, nz
       do j = 0, ny+1
          do i = 0, nx+1
             y(i,j,k) = (cff1*x(i,j,k)+omega*(&
                  (x(i-1,j,k)+x(i+1,j,k)+x(i,j-1,k)+x(i,j+1,k))&
                  -b(i,j,k))*idiag(i,j,k))
          enddo
       enddo
    enddo

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             x(i,j,k) = (cff1*y(i,j,k)+omega*(&
                  (y(i-1,j,k)+y(i+1,j,k)+y(i,j-1,k)+y(i,j+1,k))&
                  -b(i,j,k))*idiag(i,j,k))
          enddo
       enddo
    enddo
  end subroutine fsmoother2d


  subroutine fsmoother3d(x,b,y,idiag,omega,nx,ny,nz,nh)
    integer::nh
    integer:: nx,ny,nz
    REAL,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::x,y,b,idiag
    !integer,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::msk

    real::omega

    intent(inout)::x,y
    intent(in)::b,idiag!,msk

    integer:: i,j,k
    real::cff1

    cff1 = 1.-omega

    do k = nh, nz+1-nh
       do j = 0, ny+1
          do i = 0, nx+1
             y(i,j,k) = (cff1*x(i,j,k)+omega*(&
                  (x(i-1,j,k)+x(i+1,j,k)+x(i,j-1,k)+x(i,j+1,k)&
                  +x(i,j,k-1)+x(i,j,k+1))&
                  -b(i,j,k))*idiag(i,j,k))
          enddo
       enddo
    enddo

    do k = nh+1, nz-nh
       do j = 1, ny
          do i = 1, nx
             x(i,j,k) = (cff1*y(i,j,k)+omega*(&
                  (y(i-1,j,k)+y(i+1,j,k)+y(i,j-1,k)+y(i,j+1,k)&
                  +y(i,j,k-1)+y(i,j,k+1))&
                  -b(i,j,k))*idiag(i,j,k))
          enddo
       enddo
    enddo
  end subroutine fsmoother3d

  function fnorm2d(msk, x, nx,ny,nz,nh) result(sum)
    integer:: i,j,k
    integer::nh
    integer:: nx,ny,nz
    REAL,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::x
    integer,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::msk
    REAL:: sum

    sum = zero

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             sum = sum+msk(i,j,k)*(x(i,j,k)*x(i,j,k))
          enddo
       enddo
    enddo

  end function fnorm2d

  function fnorm3d(msk, x, nx,ny,nz,nh) result(sum)
    integer:: i,j,k
    integer::nh
    integer:: nx,ny,nz
    REAL,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::x
    integer,dimension(1-nh:nx+nh, 1-nh:ny+nh, 1:nz)::msk
    REAL:: sum

    sum = zero

    do k = 1+nh, nz-nh
       do j = 1, ny
          do i = 1, nx
             sum = sum+msk(i,j,k)*(x(i,j,k)*x(i,j,k))
          enddo
       enddo
    enddo

  end function fnorm3d

end module basicoperators
