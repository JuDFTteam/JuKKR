module mod_forceh

contains

subroutine forceh(cmom, flmh, lmax, nspin, nstart, nend, r2rho, v, r, drdi, &
  irws, z)
  use :: mod_datatypes, only: dp
  ! -----------------------------------------------------------------------
  ! calculates the force on nucleus m with hellmann - feynman theorem
  ! from a given non spherical charge density at the nucleus site r


  ! -----------------------------------------------------------------------
  use global_variables
  implicit none
  ! ..
  ! .. Local Scalars ..
  integer :: lmax, nend, nspin, nstart
  ! ..
  ! .. Local Arrays ..
  real (kind=dp) :: cmom(lmpotd, *), drdi(irmd, *), flmh(-1:1, *), r(irmd, *), &
    r2rho(irmd, lmpotd, natypd, *), v(irmd, lmpotd, *), z(*)
  integer :: irws(*)
  ! ..
  ! .. External Subroutines ..
  real (kind=dp) :: pi, rws, vint1
  integer :: i, iatyp, ipot, irws1, lm, m
  ! ..
  ! .. Save statement ..
  real (kind=dp) :: flm(-1:1, 2), v1(irmd)
  ! ..

  external :: simp3
  ! .. Intrinsic Functions ..
  ! ..
  save :: pi


  ! ---> loop over the rep. atoms
  intrinsic :: atan

  pi = 4.e0_dp*atan(1.e0_dp)
  if (lmax<1) then
    write (6, fmt=100)
    stop

  end if
  ! ---> reading the right Wigner-S. radius


  do iatyp = nstart, nend
    ! ---> determine the right potential numbers


    irws1 = irws(iatyp)
    rws = r(irws1, iatyp)


    ! ---> integrate with simpson subroutine
    ipot = nspin*(iatyp-1) + 1

    do m = -1, 1
      lm = 2 + m + 1

      v1(1) = 0.0e0_dp
      do i = 2, irws1
        v1(i) = r2rho(i, lm, iatyp, 1)*(r(i,iatyp)**(-2.0e0_dp))
      end do

      ! ---> use coulomb potential to determine extra atomic contribution

      call simp3(v1, vint1, 1, irws1, drdi(1,iatyp))

      flm(m, 1) = 2.0e0_dp*vint1
      ! ---> total Hellman-Feynman force


      flm(m, 2) = v(irws1, lm, ipot)*(3.0e0_dp/(4.0e0_dp*pi*rws)) - &
        2.0e0_dp*cmom(lm, iatyp)/(rws**3)


      ! -----------------------------------------------------------------------
      flmh(m, iatyp) = (flm(m,1)+flm(m,2))*z(iatyp)
    end do
  end do
  ! calculates the force on nucleus m with hellmann - feynman theorem
  ! from a given non spherical charge density at the nucleus site r
100 format (13x, 'error stop in subroutine force :', &
    ' the charge density has to contain non spherical', &
    ' contributions up to l=1 at least ')

end subroutine forceh

end module mod_forceh
