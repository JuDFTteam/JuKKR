!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_forceh

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates fore on nucleaus with Hellmann-Feynmann theorem
  !> Author: 
  !> Category: KKRhost, physical-observables
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates the force on nucleus m with Hellmann - Feynman theorem
  !> from a given non spherical charge density at the nucleus site r
  !-------------------------------------------------------------------------------
  subroutine forceh(cmom, flmh, lmax, nspin, nstart, nend, r2rho, v, r, drdi, irws, z)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: lmpotd, irmd, natypd
    use :: mod_simp3, only: simp3
    use :: mod_constants, only: pi
    implicit none

    integer :: lmax, nend, nspin, nstart
    real (kind=dp) :: cmom(lmpotd, *), drdi(irmd, *), flmh(-1:1, *), r(irmd, *), r2rho(irmd, lmpotd, natypd, *), v(irmd, lmpotd, *), z(*)
    integer :: irws(*)
    real (kind=dp) :: rws, vint1
    integer :: i, iatyp, ipot, irws1, lm, m
    real (kind=dp) :: flm(-1:1, 2), v1(irmd)


    if (lmax<1) then
      write (6, fmt=100)
      stop
    end if

    ! ---> loop over the rep. atoms
    do iatyp = nstart, nend

      ! ---> reading the right Wigner-S. radius

      irws1 = irws(iatyp)
      rws = r(irws1, iatyp)

      ! ---> determine the right potential numbers

      ipot = nspin*(iatyp-1) + 1

      do m = -1, 1
        lm = 2 + m + 1

        v1(1) = 0.0e0_dp
        do i = 2, irws1
          v1(i) = r2rho(i, lm, iatyp, 1)*(r(i,iatyp)**(-2.0e0_dp))
        end do

        ! ---> integrate with simpson subroutine

        call simp3(v1, vint1, 1, irws1, drdi(1,iatyp))

        flm(m, 1) = 2.0e0_dp*vint1

        ! ---> use coulomb potential to determine extra atomic contribution

        flm(m, 2) = v(irws1, lm, ipot)*(3.0e0_dp/(4.0e0_dp*pi*rws)) - 2.0e0_dp*cmom(lm, iatyp)/(rws**3)

        ! ---> total Hellman-Feynman force

        flmh(m, iatyp) = (flm(m,1)+flm(m,2))*z(iatyp)

      end do
    end do

100 format (13x, 'error stop in subroutine force :', ' the charge density has to contain non spherical', ' contributions up to l=1 at least ')

  end subroutine forceh

end module mod_forceh
