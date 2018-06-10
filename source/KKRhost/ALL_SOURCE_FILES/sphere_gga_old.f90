subroutine sphere_gga(lmax, yr, wtyr, rij, ijd, lmmaxd, thet, ylm, dylmt1, &
  dylmt2, dylmf1, dylmf2, dylmtf)
!-----------------------------------------------------------------------
!     generate an angular mesh and spherical harmonics at those
!     mesh points. For an angular integration the weights are ge-
!     rated .

!     R. Zeller      Feb. 1996
!     Small change for GGA implementation
!     Nikos          Dec. 1996
!-----------------------------------------------------------------------
  implicit none
!.. Scalar Arguments ..
  integer :: ijd, lmax, lmmaxd
!..
!.. Local Scalars ..
  double precision :: dx1, dx2, dx3, f0, pi, r, r1, r2, r3
  integer :: ij, lm1
!..
!.. External Subroutines ..
  external :: cylm02, ymy
!..
!.. Array Arguments ..
  double precision :: dylmf1(ijd, lmmaxd), dylmf2(ijd, lmmaxd), &
    dylmt1(ijd, lmmaxd), dylmt2(ijd, lmmaxd), dylmtf(ijd, lmmaxd), &
    rij(ijd, 3), thet(ijd), wtyr(ijd, *), ylm(ijd, lmmaxd), yr(ijd, *)
!..
!.. Local Arrays ..

  double precision :: cosfi(ijd), cosx(ijd), fai(ijd), nd(3, 3), sinfi(ijd), &
    wght, y(1000)
!..
!.. Intrinsic Functions ..
  intrinsic :: abs, acos, atan, cos, sin, sqrt
!     ..
  pi = 4.d0*atan(1.d0)
  write (6, *) 'SPHERE for GGA: read LEBEDEV mesh'
  if (ijd>1000) stop 'SPHERE'


  do ij = 1, ijd
    call lebedev(ij, r1, r2, r3, wght)

!      make a small rotation

    f0 = 0.08d0
    nd(1, 1) = cos(f0)
    nd(1, 2) = 0d0
    nd(1, 3) = sin(f0)
    nd(2, 1) = 0d0
    nd(2, 2) = 1d0
    nd(2, 3) = 0d0
    nd(3, 1) = -sin(f0)
    nd(3, 2) = 0d0
    nd(3, 3) = cos(f0)

    dx1 = nd(1, 1)*r1 + nd(2, 1)*r2 + nd(3, 1)*r3
    dx2 = nd(1, 2)*r1 + nd(2, 2)*r2 + nd(3, 2)*r3
    dx3 = nd(1, 3)*r1 + nd(2, 3)*r2 + nd(3, 3)*r3

    r1 = dx1
    r2 = dx2
    r3 = dx3

    rij(ij, 1) = r1
    rij(ij, 2) = r2
    rij(ij, 3) = r3

    call ymy(r1, r2, r3, r, y, lmax)
    do lm1 = 1, (lmax+1)**2
      yr(ij, lm1) = y(lm1)
    end do

!---> multiply the spherical harmonics with the weights

    do lm1 = 1, (lmax+1)**2
      wtyr(ij, lm1) = yr(ij, lm1)*wght*pi*4.d0
    end do

!---> produce what is needed for GGA

    cosx(ij) = r3
    if (abs(r3)/=1.d0) then
      cosfi(ij) = r1/sqrt(1.d0-r3*r3)
      sinfi(ij) = r2/sqrt(1.d0-r3*r3)
      if (abs(cosfi(ij))>1.d0) cosfi(ij) = cosfi(ij)/abs(cosfi(ij))
      if (abs(sinfi(ij))>1.d0) sinfi(ij) = sinfi(ij)/abs(sinfi(ij))
      fai(ij) = acos(cosfi(ij))
    else if (sinfi(ij)==0.d0) then
      fai(ij) = pi/2.d0
    else
      cosfi(ij) = r1
      sinfi(ij) = r2
      if (abs(cosfi(ij))>1.d0) cosfi(ij) = cosfi(ij)/abs(cosfi(ij))
    end if
    fai(ij) = acos(cosfi(ij))
  end do

  call cylm02(lmax, cosx, fai, 2*lmax+1, lmmaxd, thet, ylm, dylmt1, dylmt2, &
    dylmf1, dylmf2, dylmtf)

end subroutine
