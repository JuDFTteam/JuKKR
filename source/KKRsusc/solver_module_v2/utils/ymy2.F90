  subroutine ymy2(ctheta,phi,nl,nlm,ylm)
! this subroutine calculates real spherical harmonics with the
! normalization : <y|y> =1
!
! pilfered and abused in 2012 by MdSD

  implicit none

! spherical angles
  real(kind=r8b),    intent(in)  :: ctheta, phi
! dimensions
  integer(kind=i4b), intent(in)  :: nl, nlm
! real spherical harmonics
  real(kind=r8b),    intent(out) :: ylm(nlm)
! -----------------------------------------------------------------
! parameters
  real(kind=r8b), parameter :: invr4pi = 0.25d0/sqrt(atan(1.d0))
  real(kind=r8b), parameter :: rtwo = sqrt(2.d0)
  real(kind=r8b), parameter :: snull = 1.d-10
! -----------------------------------------------------------------
  real(kind=r8b)    :: xy, xyz, cth, sth, cph, sph, sthm
  integer(kind=i4b) :: i, l, m
  real(kind=r8b)    :: p(0:nl,0:nl), c(0:nl), s(0:nl)
  real(kind=r8b)    :: dfac, root, root1, root2


! I can't be bothered doing the extra if's
  if (nl < 2) stop 'check ymy2'


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            calculate sin and cos of theta and phi
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cth = ctheta
  sth = sqrt(1.d0 - ctheta*ctheta)
! pathological thetas
  if (abs(ctheta) < snull) then
    cth =  0.d0; sth = 1.d0
  else if (abs(ctheta - 1.d0) < snull) then
    cth =  1.d0; sth = 0.d0
  else if (abs(ctheta + 1.d0) < snull) then
    cth = -1.d0; sth = 0.d0
  end if
  cph = cos(phi)
  sph = sin(phi)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     generate normalized associated legendre functions for m >= 0
!                see NR section 6.7 (only 2007)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! loop over m values
  dfac = 1; root = 1.d0; sthm = 1.d0; p(0,0) = 1.d0!invr4pi
  do m=0,nl
    if (m > 0) then
      sthm = sthm*sth
      dfac = dfac*sqrt((2*m-1)**2/(2.d0*m*(2*m-1)))
      root = sqrt((2.d0*m+1.d0))!*invr4pi
      p(m,m) = dfac*root*sthm  ! there should be a minus sign here for Condon-Shortley
    end if
    if (m < nl) p(m+1,m) = sqrt(2*m+3.d0)*cth*p(m,m)
!   upward recursion in l
    do l=m+2,nl
      root1 = sqrt((4*l*l-1.d0)/(l*l-m*m))
      root2 = sqrt(((l-1)*(l-1)-m*m)/(4*(l-1)*(l-1)-1.d0))
      p(l,m) = root1*(cth*p(l-1,m) - root2*p(l-2,m))
    end do
  end do
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            recursions for sin(m phi) and cos(m phi)
!                     see NR section 5.4
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  s(0) = 0.d0; s(1) = sph
  c(0) = 1.d0; c(1) = cph
  do m=2,nl
    s(m) = sin(m*phi)
    c(m) = cos(m*phi)
  end do
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         multiply legendre functions with cosines and sines
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  i = 0
  do l=0,nl
    i = i + l + 1
    ylm(i) = p(l,0)  ! m = 0
    do m=1,l
      ylm(i+m) = rtwo*p(l,m)*c(m)  ! m > 0 cosine type
      ylm(i-m) = rtwo*p(l,m)*s(m)  ! m < 0 sine type
    end do
    i = i + l
  end do
! All done!
  end subroutine ymy2

