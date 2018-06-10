subroutine forceh(cmom, flmh, lmax, nspin, nstart, nend, r2rho, v, r, drdi, &
  irws, z)
!-----------------------------------------------------------------------
!     calculates the force on nucleus m with hellmann - feynman theorem
!     from a given non spherical charge density at the nucleus site r


!-----------------------------------------------------------------------
  implicit none
!.. Parameters ..
  include 'inc.p'
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
!..
!.. Local Scalars ..
  integer :: lmax, nend, nspin, nstart
!..
!.. Local Arrays ..
  double precision :: cmom(lmpotd, *), drdi(irmd, *), flmh(-1:1, *), &
    r(irmd, *), r2rho(irmd, lmpotd, natypd, *), v(irmd, lmpotd, *), z(*)
  integer :: irws(*)
!..
!.. External Subroutines ..
  double precision :: pi, rws, vint1
  integer :: i, iatyp, ipot, irws1, lm, m
!..
!.. Save statement ..
  double precision :: flm(-1:1, 2), v1(irmd)
!..

  external :: simp3
!.. Intrinsic Functions ..
!..
  save :: pi


!---> loop over the rep. atoms
  intrinsic :: atan

  pi = 4.d0*atan(1.d0)
  if (lmax<1) then
    write (6, fmt=100)
    stop

  end if
!---> reading the right Wigner-S. radius


  do iatyp = nstart, nend
!---> determine the right potential numbers


    irws1 = irws(iatyp)
    rws = r(irws1, iatyp)


!---> integrate with simpson subroutine
    ipot = nspin*(iatyp-1) + 1

    do m = -1, 1
      lm = 2 + m + 1

      v1(1) = 0.0d0
      do i = 2, irws1
        v1(i) = r2rho(i, lm, iatyp, 1)*(r(i,iatyp)**(-2.0d0))
      end do

!---> use coulomb potential to determine extra atomic contribution

      call simp3(v1, vint1, 1, irws1, drdi(1,iatyp))

      flm(m, 1) = 2.0d0*vint1
!---> total Hellman-Feynman force


      flm(m, 2) = v(irws1, lm, ipot)*(3.0d0/(4.0d0*pi*rws)) - &
        2.0d0*cmom(lm, iatyp)/(rws**3)


!-----------------------------------------------------------------------
      flmh(m, iatyp) = (flm(m,1)+flm(m,2))*z(iatyp)
    end do
  end do
!     calculates the force on nucleus m with hellmann - feynman theorem
!     from a given non spherical charge density at the nucleus site r
100 format (13x, 'error stop in subroutine force :', &
    ' the charge density has to contain non spherical', &
    ' contributions up to l=1 at least ')

end subroutine
