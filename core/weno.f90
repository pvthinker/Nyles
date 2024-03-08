function weno3(qm,q0,qp) result(q)
  real*8:: qm,q0,qp
  real*8:: q

  real*8:: eps,beta1,beta2,tau
  real*8:: qi1,qi2,w1,w2

  eps = 1e-14

  qi1 = (-qm + 3*q0)*0.5
  qi2 = (q0 + qp)*0.5

  beta1 = (q0-qm)**2
  beta2 = (qp-q0)**2
  tau = abs(beta2-beta1)

  w1 = (1. + tau / (beta1 + eps))
  w2 = (1. + tau / (beta2 + eps))*2

  q = (w1*qi1 + w2*qi2) / (w1 + w2)

end function weno3


function weno5(qmm,qm,q0,qp,qpp) result(q)
  real*8:: qmm,qm,q0,qp,qpp
  real*8:: q

  real*8:: eps,beta1,beta2,beta3,tau
  real*8:: k1, k2
  real*8:: qi1,qi2,qi3,w1,w2,w3

  eps = 1e-16

  qi1 = 1./3.*qmm - 7./6.*qm + 11./6.*q0
  qi2 = -1./6.*qm + 5./6.*q0 + 1./3.*qp
  qi3 = 1./3.*q0 + 5./6.*qp - 1./6.*qpp

  k1 = 13./12.
  k2 = .25
  beta1 = k1 * (qmm-2*qm+q0)**2 + k2 * (qmm-4*qm+3*q0)**2
  beta2 = k1 * (qm-2*q0+qp)**2  + k2 * (qm-qp)**2
  beta3 = k1 * (q0-2*qp+qpp)**2 + k2 * (3*q0-4*qp+qpp)**2

  tau5 = abs(beta1 - beta3)
  !g1, g2, g3 = 0.1, 0.6, 0.3

  w1 =     (1. + tau5 / (beta1 + eps))
  w2 = 6 * (1. + tau5 / (beta2 + eps))
  w3 = 3 * (1. + tau5 / (beta3 + eps))

  q = (w1*qi1 + w2*qi2 + w3*qi3) / (w1 + w2 + w3)

end function weno5


subroutine flux1d_test(u, q, flux, n)
  implicit none

  integer, intent(in):: n

  real*8, dimension(n), intent(in) :: q, u
  real*8, dimension(n), intent(out) :: flux
  real*8:: weno3, weno5
  integer:: i

  if(u(1).gt.0) then
     flux(1) = u(1)*q(1)
  else
     flux(1) = u(1)*q(2)!weno3(q(3),q(2),q(1))
  endif

  if(u(2).gt.0) then
     flux(2) = u(2)*q(2)!weno3(q(1),q(2),q(3))
  else
     flux(2) = u(2)*q(3)!weno5(q(5),q(4),q(3),q(2),q(1))
  endif

  do i = 3, n-3
     if(u(i).gt.0) then
        flux(i) = u(i)*q(i)!weno5(q(i-2),q(i-1),q(i),q(i+1),q(i+2))
     else
        flux(i) = u(i)*q(i+1)!weno5(q(i+3),q(i+2),q(i+1),q(i),q(i-1))
     endif
  enddo

  i = n-2
  if(u(i).gt.0) then
     flux(i) = u(i)*q(i)!weno5(q(i-2),q(i-1),q(i),q(i+1),q(i+2))
  else
     flux(i) = u(i)*q(i+1)!weno3(q(i+2),q(i+1),q(i))
  endif

  i = n-1
  if(u(i).gt.0) then
     flux(i) = u(i)*q(i)!weno3(q(i-1),q(i),q(i+1))
  else
     flux(i) = u(i)*q(i+1)
  endif

  flux(n) = 0.


end subroutine flux1d_test

subroutine flux1d(u, q, flux, n)
  implicit none

  integer, intent(in):: n

  real*8, dimension(n), intent(in) :: q, u
  real*8, dimension(n), intent(out) :: flux
  real*8:: weno3, weno5
  integer:: i

  if(u(1).gt.0) then
     flux(1) = u(1)*q(1)
  else
     flux(1) = u(1)*weno3(q(3),q(2),q(1))
  endif

  if(u(2).gt.0) then
     flux(2) = u(2)*weno3(q(1),q(2),q(3))
  else
     flux(2) = u(2)*weno5(q(5),q(4),q(3),q(2),q(1))
  endif

  do i = 3, n-3
     if(u(i).gt.0) then
        flux(i) = u(i)*weno5(q(i-2),q(i-1),q(i),q(i+1),q(i+2))
     else
        flux(i) = u(i)*weno5(q(i+3),q(i+2),q(i+1),q(i),q(i-1))
     endif
  enddo

  i = n-2
  if(u(i).gt.0) then
     flux(i) = u(i)*weno5(q(i-2),q(i-1),q(i),q(i+1),q(i+2))
  else
     flux(i) = u(i)*weno3(q(i+2),q(i+1),q(i))
  endif

  i = n-1
  if(u(i).gt.0) then
     flux(i) = u(i)*weno3(q(i-1),q(i),q(i+1))
  else
     flux(i) = u(i)*q(i+1)
  endif

  flux(n) = 0.


end subroutine flux1d
