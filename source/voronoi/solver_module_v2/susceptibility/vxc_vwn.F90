  subroutine vxc_vwn(den,mag,vxc,bxc,exc)
! xc potential from VWN parametrization
! densities with 4pi and without r^2
! check definition of rs

  implicit none

! Charge and magnetization densities
  real(kind=r8b), intent(in)  :: den, mag
! xc-potentials, charge and magnetic
  real(kind=r8b), intent(out) :: vxc, bxc, exc
! Parameters of the VWN parametrization (taken from the vosko.f routine)
  real(kind=r8b), parameter :: tol = 1.d-10
  real(kind=r8b), parameter :: fourpi = 16.d0*atan(1.d0)
  real(kind=r8b), parameter :: ap = 0.0621814d0, xp0 = -0.10498d0, bp = 3.72744d0
  real(kind=r8b), parameter :: cp = 12.9352d0, qp = 6.1519908d0
  real(kind=r8b), parameter :: cp1= 1.2117833d0, cp2= 1.1435257d0, cp3 = -0.031167608d0
  real(kind=r8b), parameter :: af = 0.0310907d0, xf0=-0.32500d0, bf = 7.06042d0
  real(kind=r8b), parameter :: cf = 18.0578d0, qf = 4.7309269d0
  real(kind=r8b), parameter :: cf1 = 2.9847935d0, cf2= 2.7100059d0, cf3= -0.1446006d0
! --------------------------------------------------------------------------------------
  real(kind=r8b) :: rs, zeta
  real(kind=r8b) :: fzeta, dfzeta, cbrt1, cbrt2, zeta3, zeta4, x, x2, x3, x4, xpx, xfx
  real(kind=r8b) :: beta, dbeta, atanp, atanf, ecp, ecf, ec
  real(kind=r8b) :: tp1, tf1, ucp, ucf, uc0, uc10, uc20, duc, duc1, duc2, uc1, uc2
  

! cut the densities for low values (vacuum)  
  if (den < tol) then
    vxc = 0.d0
    bxc = 0.d0
    return
  end if

! density and relative spin polarization
!  rs = (3.d0/(fourpi*den))**(1.d0/3.d0)
  rs = (3.d0/den)**(1.d0/3.d0)
!  rs = (3.d0/den)**(1.d0/3.d0)
  x = sqrt(rs); x2 = x*x; x3 = x*x2; x4 = x2*x2
  zeta = mag/den
  zeta3 = zeta**3
  zeta4 = zeta**4 - 1.d0

! exchange enhancement
  fzeta = ((1.d0 + zeta)**(4.d0/3.d0) + (1.d0 - zeta)**(4.d0/3.d0) - 2.d0)/(2.d0**(4.d0/3.d0) - 2.d0)
  cbrt1 = (1.d0 + zeta)**(1.d0/3.d0)
  cbrt2 = (1.d0 - zeta)**(1.d0/3.d0)
  dfzeta = (4.d0/3.d0)*(cbrt1 - cbrt2)/(2.d0**(4.d0/3.d0) - 2.d0)

! build up the correlation energy
  xpx = x2 + bp*x + cp
  xfx = x2 + bf*x + cf
  beta = 1.d0/(2.74208d0 + 3.182d0*x + 0.09873d0*x2 + 0.18268d0*x3)
  dbeta = -(0.27402d0*x + 0.09873d0 + 1.591d0/x)*beta**2
  atanp = atan(qp/(2.d0*x + bp))
  atanf = atan(qf/(2.d0*x + bf))
  ecp = ap*(log(x2/xpx) + cp1*atanp - cp3*(log((x - xp0)**2/xpx) + cp2*atanp))
  ecf = af*(log(x2/xfx) + cf1*atanf - cf3*(log((x - xf0)**2/xfx) + cf2*atanf))
  ec = ecp + fzeta*(ecf - ecp)*(1.d0 + zeta4*beta)
  exc = ec - 0.9163306d0/rs - 0.2381735d0/rs*fzeta

! build up the correlation potential
  tp1 = (x2 + bp*x)/xpx
  tf1 = (x2 + bf*x)/xfx
  ucp = ecp - ap/3.d0*(1.d0 - tp1 - cp3*(x/(x - xp0) - tp1 - xp0*x/xpx))
  ucf = ecf - af/3.d0*(1.d0 - tf1 - cf3*(x/(x - xf0) - tf1 - xf0*x/xfx))
  uc0 = ucp + (ucf - ucp)*fzeta
  uc10 = uc0 - (ecf - ecp)*(zeta - 1.d0)*dfzeta
  uc20 = uc0 - (ecf - ecp)*(zeta + 1.d0)*dfzeta
  duc = (ucf - ucp)*beta*zeta4*fzeta + (ecf - ecp)*(-rs/3.d0)*dbeta*zeta4*fzeta
  duc1 = duc - (ecf - ecp)*beta*(zeta - 1.d0)*(4.d0*zeta3*fzeta + zeta4*dfzeta)
  duc2 = duc - (ecf - ecp)*beta*(zeta + 1.d0)*(4.d0*zeta3*fzeta + zeta4*dfzeta)
  uc1 = uc10 + duc1
  uc2 = uc20 + duc2

! xc potentials
  vxc = 0.5d0*(uc1 + uc2) - 0.5d0*(cbrt1 + cbrt2)*1.221774d0/rs
  bxc = 0.5d0*(uc1 - uc2) - 0.5d0*(cbrt1 - cbrt2)*1.221774d0/rs

! All done!
  end subroutine vxc_vwn

