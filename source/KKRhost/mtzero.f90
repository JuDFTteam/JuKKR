!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Determine muffin tin zero and shift potential to muffin tin zero
!> Author: 
!> Determine muffin tin zero and shift potential to muffin tin zero. For spin 
!> polarized calculations muffin tin zero is related to the average of the 2 spins.
!------------------------------------------------------------------------------------
module mod_mtzero

contains

  !-------------------------------------------------------------------------------
  !> Summary: Determine muffin tin zero and shift potential to muffin tin zero
  !> Author: 
  !> Category: potential, KKRhost
  !> Deprecated: False 
  !> Determine muffin tin zero and shift potential to muffin tin zero. For spin 
  !> polarized calculations muffin tin zero is related to the average of the 2 spins.
  !-------------------------------------------------------------------------------
  subroutine mtzero(lmpot,natyp,conc,nspin,v,vbc,z,r,drdi,imt,ircut,ipan,ntcell,    &
use :: mod_runoptions, only: use_decimation --manopt-- 
    lmsp,ifunm,thetas,irws,eshift,ishift,nshell,lsurf)

    use :: mod_datatypes, only: dp
    use :: global_variables
    use :: mod_simp3
    use :: mod_simpk
    use :: mod_constants, only : pi
    implicit none
    ! ..
    ! .. Input variables
    integer, intent(in) :: lmpot !! (LPOT+1)**2 
    integer, intent(in) :: natyp !! Number of kinds of atoms in unit cell
    integer, intent(in) :: nspin !! Counter for spin directions
    integer, intent(in) :: ishift
    real (kind=dp), intent(in) :: eshift
    logical, intent(in) :: lsurf !! If True a matching with semi-inifinite surfaces must be performed
    integer, dimension(*), intent(in) :: imt    !! R point at MT radius
    integer, dimension(*), intent(in) :: ipan   !! Number of panels in non-MT-region
    integer, dimension(*), intent(in) :: irws   !! R point at WS radius
    integer, dimension(*), intent(in) :: ntcell !! Index for WS cell
    integer, dimension(0:nsheld), intent(in) :: nshell !! Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
    integer, dimension(natypd,*), intent(in)  :: lmsp  !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension(0:ipand,*), intent(in) :: ircut !! R points of panel borders
    integer, dimension(natypd,*), intent(in)  :: ifunm
    real (kind=dp), dimension(*), intent(in) :: z  !! Nuclear charge
    real (kind=dp), dimension(natypd), intent(in) :: conc !! Concentration of a given atom
    real (kind=dp), dimension(irmd,*), intent(in) :: r !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension(irmd,*), intent(in) :: drdi !! Derivative dr/di
    real (kind=dp), dimension(irid,nfund,*), intent(in) :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    ! .. In/Out variables
    real (kind=dp), dimension(*), intent(inout) :: vbc !! Potential constants
    real (kind=dp), dimension(irmd,lmpotd,*), intent(inout) :: v !! output potential (nonspherical VONS)
    ! .. Local variables
    real (kind=dp) :: fpi, rfpi, vav0, vol0, zzor
    integer :: icell, ifun, ih, imt1, ipan1, ipot, ir, irc1, irh, is, lm
    real (kind=dp), dimension(2) :: vav1, vol1
    real (kind=dp), dimension(irmd) :: v1, v2

    logical, external :: test, opt

    intrinsic :: sqrt

    fpi = 4.0_dp*pi
    rfpi = sqrt(fpi)
    ! ---  >     muffin tin or atomic sphere calculation
    vav0 = 0.0e0_dp
    vol0 = 0.0e0_dp
    vav1(1) = 0.e0_dp
    vav1(2) = 0.e0_dp
    vol1(1) = 0.e0_dp
    vol1(2) = 0.e0_dp
    do ih = 1, natyp

      do ir = 1, irmd
        v1(ir) = 0.0e0_dp
        v2(ir) = 0.0e0_dp
      end do
      do is = 1, nspin
        ipot = nspin*(ih-1) + is
        ipan1 = ipan(ih)
        imt1 = imt(ih)

        if (ipan1==1) then

          ! (IPAN1.EQ.1)

          irc1 = irws(ih)
          do ir = imt1, irc1
            v2(ir) = fpi*r(ir, ih)**2
            zzor = 2.0e0_dp*z(ih)/r(ir, ih)
            v1(ir) = (v(ir,1,ipot)/rfpi-zzor)*v2(ir)
          end do
          ! ---  >     full potential calculation
          call simp3(v1, vav1(is), imt1, irc1, drdi(1,ih))
          call simp3(v2, vol1(is), imt1, irc1, drdi(1,ih))

        else

          irc1 = ircut(ipan1, ih)
          icell = ntcell(ih)
          imt1 = imt(ih)
          do ir = imt1 + 1, irc1
            v2(ir) = r(ir, ih)**2*thetas(ir-imt1, 1, icell)*rfpi
            zzor = 2.0e0_dp*z(ih)/r(ir, ih)
            v1(ir) = (v(ir,1,ipot)/rfpi-zzor)*v2(ir)
          end do
          do lm = 2, lmpot
            if (lmsp(icell,lm)>0) then
              ifun = ifunm(icell, lm)

              do ir = imt1 + 1, irc1
                irh = ir - imt1
                v1(ir) = v1(ir) + r(ir, ih)**2*v(ir, lm, ipot)*thetas(irh, ifun, icell)
              end do
              ! (IPAN1.EQ.1)
            end if

          end do
          ! SPIN LOOP
          call simpk(v1, vav1(is), ipan1, ircut(0,ih), drdi(1,ih))
          call simpk(v2, vol1(is), ipan1, ircut(0,ih), drdi(1,ih))

        end if
        ! 19.5.99   Nikos
      end do                       ! This way it is compatible with old kkr and tb-kkr
      if (nspin==1) then
        vav1(2) = vav1(1)
        vol1(2) = vol1(1)
      end if

      ! added 10.11.99 to fix vbc

      if (lsurf .and. (ih==1)) write (1337, *) 'Vacancies are ignored for VBC'
      if (lsurf .and. (z(ih)<1.e0_dp)) cycle
      ! ---  > shift potential to muffin tin zero
      vav0 = vav0 + conc(ih)*nshell(ih)*(vav1(1)+vav1(2))/2.e0_dp
      vol0 = vol0 + conc(ih)*nshell(ih)*(vol1(1)+vol1(2))/2.e0_dp
    end do
    if (.not. (use_decimation)) then
      vbc(1) = 0.0e0_dp
      if (abs(vav0)>1e-10_dp) vbc(1) = -vav0/vol0
      if (ishift>0) vbc(1) = vbc(1) + eshift
    end if

    write (1337, fmt=100) vol0, vav0, vbc(1)
    vbc(2) = vbc(1)

    ! ************************************************************************
    do is = 1, nspin
      do ih = 1, natyp
        ipot = nspin*(ih-1) + is
        do ir = 1, ircut(ipan(ih), ih)
          v(ir, 1, ipot) = v(ir, 1, ipot) + rfpi*vbc(is)
        end do
        ! ************************************************************************
      end do
    end do

    return
    ! determine muffin tin zero and shift potential to muffin tin zero
100 format ('  VOL INT.', f16.9, '  VAV INT.', f16.9, '  VMT ZERO', f16.9)
110 format ('  ATOM ', i4, ' VMT ZERO :', f16.9)
  end subroutine mtzero

end module mod_mtzero
