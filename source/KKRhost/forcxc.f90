!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_forcxc

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates forrce on nucleus with core correction (xc-contribution)
  !> Author: 
  !> Category: KKRhost, physical-observables, xc-potential
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates the force on nucleus m
  !> from a given non spherical charge density at the nucleus site r
  !> with core correction(exchange contribution)
  !>
  !> @warning
  !> BEWARE!!! RM commented away!!! -->Dipole Tensor is useless
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine forcxc(flm, flmc, lmax, nspin, nstart, nend, rhoc, v, r, alat, drdi, irws, natref)
    !>>>>>BEWARE!!! RM commented away!!! -->Dipole Tensor is useless
    ! SUBROUTINE FORCXC(FLM,FLMC,LMAX,NSPIN,NSTART,NEND,RHOC,V,R,ALAT,
    ! +                  RM,NSHELL,DRDI,IRWS,NATREF)
    use :: mod_types, only: t_inc
    use :: mod_datatypes, only: dp
    use :: global_variables, only: irmd, lmpotd, natypd
    use :: mod_simp3, only: simp3
    use :: mod_constants, only: pi
    implicit none

    real (kind=dp), parameter :: fac = sqrt((4.0_dp*pi)/3.0_dp)
    real (kind=dp) :: alat
    integer :: lmax, natref, nend, nspin, nstart
    real (kind=dp) :: drdi(irmd, *), flm(-1:1, *), flmc(-1:1, *), r(irmd, *), rhoc(irmd, *), v(irmd, lmpotd, *)
    integer :: irws(*)
    real (kind=dp) :: dv, vint1
    integer :: i, iatyp, iper, ipot, irws1, ispin, lm, m
    real (kind=dp) :: flmh(-1:1, natypd), flmxc(-1:1, natypd), v1(irmd)
    intrinsic :: dsqrt

    if (lmax<1) then
      write (6, fmt=100)
      stop
    end if

    if (t_inc%i_write>0) write (1337, fmt=130)
    if (t_inc%i_write>0) write (1337, fmt=120)
    if (t_inc%i_write>0) write (1337, fmt=130)

    do iatyp = nstart, nend

      iper = iatyp - natref
      if (t_inc%i_write>0) write (1337, fmt=140) iper

      irws1 = irws(iatyp)

      do m = -1, 1
        lm = 2 + m + 1

        do i = 1, irws1
          v1(i) = 0.0_dp
        end do

        do ispin = 1, nspin
          ! ---> determine the right potential numbers
          ipot = nspin*(iatyp-1) + ispin
          dv = (-3.0_dp*v(1,lm,ipot)-10.0_dp*v(2,lm,ipot)+18.0_dp*v(3,lm,ipot)-6.0_dp*v(4,lm,ipot)+v(5,lm,ipot))/(12.0_dp*drdi(2,iatyp))
          v1(2) = rhoc(2, ipot)*(2.0_dp*v(2,lm,ipot)/r(2,iatyp)+dv)/(4.0_dp*pi) + v1(2)

          do i = 3, irws1 - 2
            dv = (v(i-2,lm,ipot)-v(i+2,lm,ipot)+8.0_dp*(v(i+1,lm,ipot)-v(i-1,lm,ipot)))/(12.0_dp*drdi(i,iatyp))
            v1(i) = rhoc(i, ipot)*(2.0_dp*v(i,lm,ipot)/r(i,iatyp)+dv)/(4.0_dp*pi) + v1(i)
          end do

          dv = (-v(irws1-4,lm,ipot)+6.0_dp*v(irws1-3,lm,ipot)-18.0_dp*v(irws1-2,lm,ipot)+10.0_dp*v(irws1-1,lm,ipot)+3.0_dp*v(irws1,lm,ipot))/(12.0_dp*drdi(irws1-1,iatyp))
          v1(irws1-1) = rhoc(irws1-1, ipot)*(2.0_dp*v(irws1-1,lm,ipot)/r(irws1-1,iatyp)+dv)/(4.0_dp*pi) + v1(irws1-1)

          dv = (3.0_dp*v(irws1-4,lm,ipot)-16.0_dp*v(irws1-3,lm,ipot)+36.0_dp*v(irws1-2,lm,ipot)-48.0_dp*v(irws1-1,lm,ipot)+25.0_dp*v(irws1,lm,ipot))/(12.0_dp*drdi(irws1,iatyp))
          v1(irws1) = rhoc(irws1, ipot)*(2.0_dp*v(irws1,lm,ipot)/r(irws1,iatyp)+dv)/(4.0_dp*pi) + v1(irws1)
        end do

        ! ---> integrate with simpson subroutine
        call simp3(v1, vint1, 1, irws1, drdi(1,iatyp))

        flmh(m, iatyp) = flm(m, iatyp) - flmc(m, iatyp)
        flmxc(m, iatyp) = -fac*vint1 - flmc(m, iatyp)
        flm(m, iatyp) = flm(m, iatyp) + flmxc(m, iatyp)

      end do

      if (t_inc%i_write>0) then
        write (1337, fmt=150) flmh(1, iatyp), flmc(1, iatyp), flmxc(1, iatyp), flm(1, iatyp)
        write (1337, fmt=160) flmh(-1, iatyp), flmc(-1, iatyp), flmxc(-1, iatyp), flm(-1, iatyp)
        write (1337, fmt=170) flmh(0, iatyp), flmc(0, iatyp), flmxc(0, iatyp), flm(0, iatyp)
      end if

    end do

    if (t_inc%i_write>0) then
      write (1337, fmt=130)
      write (1337, fmt=110)
      write (1337, fmt=130)
    end if

100 format (13x, 'error stop in subroutine force :', ' the charge density has to contain non spherical', ' contributions up to l=1 at least ')
110 format (1x, 81('-'))
120 format (1x, 33('-'), ' force on the nucleus ', 33('-'), /, 34x, ' in units Ry/(a(Bohr) ')
130 format (1x, '>')
140 format (3x, i5, 'th shell')
150 format (7x, 'fhx=', e13.6, 2x, 'fcx=', e13.6, 2x, 'fxcx=', e13.6, 2x, 'fx=', e13.6, ' Ry/(a(Bohr))')
160 format (7x, 'fhy=', e13.6, 2x, 'fcy=', e13.6, 2x, 'fxcy=', e13.6, 2x, 'fy=', e13.6, ' Ry/(a(Bohr))')
170 format (7x, 'fhz=', e13.6, 2x, 'fcz=', e13.6, 2x, 'fxcz=', e13.6, 2x, 'fz=', e13.6, ' Ry/(a(Bohr))')

  end subroutine forcxc

end module mod_forcxc
