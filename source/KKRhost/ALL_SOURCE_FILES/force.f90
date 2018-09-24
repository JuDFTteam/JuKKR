module mod_force

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates force on nucleus with core contribution (Coulomb contribution)
  !> Author: 
  !> Category: KKRhost, physical-observables
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates the force on nucleus m
  !> from a given non spherical charge density at the nucleus site r
  !> with core correction (coulomb contribution)
  !-------------------------------------------------------------------------------
  subroutine force(flm, flmc, lmax, nspin, nstart, nend, rhoc, v, r, drdi, irws)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: irmd, lmpotd, natypd
    use :: mod_simp3, only: simp3
    use :: mod_constants, only: pi
    implicit none

    integer :: lmax, nend, nspin, nstart
    ! ..
    real (kind=dp) :: drdi(irmd, *), flm(-1:1, *), flmc(-1:1, *), r(irmd, *), rhoc(irmd, *), v(irmd, lmpotd, *)
    integer :: irws(*)
    ! ..
    real (kind=dp) :: dv, rws, vint1
    real (kind=dp), parameter :: fac = sqrt((4.0e0_dp*pi)/3.0e0_dp)
    integer :: i, iatyp, ipot, irws1, ispin, lm, m
    ! ..
    real (kind=dp) :: flmh(-1:1, natypd), v1(irmd)

    intrinsic :: atan, sqrt


    if (lmax<1) then
      write (6, fmt=100)
      stop

    end if

    ! ---> loop over rep. atoms

    do iatyp = nstart, nend

      irws1 = irws(iatyp)
      rws = r(irws1, iatyp)

      do m = -1, 1
        lm = 2 + m + 1

        ! ---> initialize v1

        do i = 1, irws1
          v1(i) = 0.0e0_dp
        end do

        do ispin = 1, nspin

          ! ---> determine the right potential numbers

          ipot = nspin*(iatyp-1) + ispin

          ! ---> determine the derivative of the potential using a 5-point formular

          dv = (-3.0e0_dp*v(1,lm,ipot)-10.0e0_dp*v(2,lm,ipot)+18.0e0_dp*v(3,lm,ipot)-6.0e0_dp*v(4,lm,ipot)+v(5,lm,ipot))/(12.0e0_dp*drdi(2,iatyp))

          v1(2) = rhoc(2, ipot)*(2.0e0_dp*v(2,lm,ipot)/r(2,iatyp)+dv)/(4.0e0_dp*pi) + v1(2)

          do i = 3, irws1 - 2

            dv = (v(i-2,lm,ipot)-v(i+2,lm,ipot)+8.0e0_dp*(v(i+1,lm,ipot)-v(i-1,lm,ipot)))/(12.0e0_dp*drdi(i,iatyp))

            v1(i) = rhoc(i, ipot)*(2.0e0_dp*v(i,lm,ipot)/r(i,iatyp)+dv)/(4.0e0_dp*pi) + v1(i)
          end do

          dv = (-v(irws1-4,lm,ipot)+6.0e0_dp*v(irws1-3,lm,ipot)-18.0e0_dp*v(irws1-2,lm,ipot)+10.0e0_dp*v(irws1-1,lm,ipot)+3.0e0_dp*v(irws1,lm,ipot))/(12.0e0_dp*drdi(irws1-1,iatyp))
          v1(irws1-1) = rhoc(irws1-1, ipot)*(2.0e0_dp*v(irws1-1,lm,ipot)/r(irws1-1,iatyp)+dv)/(4.0e0_dp*pi) + v1(irws1-1)

          dv = (3.0e0_dp*v(irws1-4,lm,ipot)-16.0e0_dp*v(irws1-3,lm,ipot)+36.0e0_dp*v(irws1-2,lm,ipot)-48.0e0_dp*v(irws1-1,lm,ipot)+25.0e0_dp*v(irws1,lm,ipot))/ &
            (12.0e0_dp*drdi(irws1,iatyp))

          v1(irws1) = rhoc(irws1, ipot)*(2.0e0_dp*v(irws1,lm,ipot)/r(irws1,iatyp)+dv)/(4.0e0_dp*pi) + v1(irws1)
        end do

        ! ---> integrate with simpson subroutine

        call simp3(v1, vint1, 1, irws1, drdi(1,iatyp))

        flmh(m, iatyp) = fac*flm(m, iatyp)
        flmc(m, iatyp) = -fac*vint1
        flm(m, iatyp) = flmh(m, iatyp) + flmc(m, iatyp)

      end do

    end do

100 format (13x, 'error stop in subroutine force :', ' the charge density has to contain non spherical', ' contributions up to l=1 at least ')

  end subroutine force

end module mod_force
