!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Performs a straight mixing of the potential 
!> Author: 
!> Performs a straight mixing between the potential from a previous iteration
!> and the one from the current iteration to speed up convergence. 
!------------------------------------------------------------------------------------
module mod_mixstr

contains

  !-------------------------------------------------------------------------------  
  !> Summary: Performs a straight mixing of the potential 
  !> Author: 
  !> Category: potential, solver, KKRhost
  !> Deprecated: False 
  !> Performs a straight mixing between the potential from a previous iteration 
  !> and the one from the current iteration to speed up convergence. 
  !-------------------------------------------------------------------------------  
  subroutine mixstr(rmsavq,rmsavm,ins,lpot,lmpot,natref,nshell,nstart,nend,conc,    &
    nspin,itc,rfpi,fpi,ipf,mixing,fcm,irc,irmin,r,drdi,vons,visp,vins,vspsmo,vspsme,&
    lsmear)
  
    use :: global_variables
    use :: mod_datatypes, only: dp
    
    implicit none

    ! .. Input variables
    integer, intent(in) :: ins    !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent(in) :: ipf    !! Not real used, IPFE should be 0
    integer, intent(in) :: itc    !! Current iteration number
    integer, intent(in) :: lpot   !! Maximum l component in potential expansion
    integer, intent(in) :: nend   !! Final atom in the loop
    integer, intent(in) :: lmpot  !! (LPOT+1)**2
    integer, intent(in) :: nspin  !! Counter for spin directions
    integer, intent(in) :: natref
    integer, intent(in) :: nstart !! First atom in the loop
    integer, intent(in) :: lsmear
    real (kind=dp), intent(in) :: fcm  !! Factor for increased linear mixing of magnetic part of potential compared to non-magnetic part.
    real (kind=dp), intent(in) :: fpi  !! 4 \(\pi\)
    real (kind=dp), intent(in) :: rfpi !! \(\sqrt{4\pi}\)
    real (kind=dp), intent(in) :: mixing !! Magnitude of the mixing parameter
    integer, dimension(natypd), intent(in) :: irc  !! R point for potential cutting
    integer, dimension(natypd), intent(in) :: irmin  !! Max R for spherical treatment
    integer, dimension(0:nsheld), intent(in) :: nshell !! Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
    real (kind=dp), dimension(natypd), intent(in) :: conc !! Concentration of a given atom
    real (kind=dp), dimension(irmd, natypd), intent(in) :: r !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension(irmd, natypd), intent(in) :: drdi !! Derivative dr/di
    real (kind=dp), dimension(irmd, *), intent(in)      :: visp !! Spherical part of the potential
    real (kind=dp), dimension(irmind:irmd, lmpotd, *), intent(in) :: vins !! Non-spherical part of the potential
    ! .. Output variables 
    real (kind=dp), intent(out) :: rmsavm !! RMS value of the magnetization
    real (kind=dp), intent(out) :: rmsavq !! RMS value of the charge
    ! .. In/Out variables
    real (kind=dp), dimension(irmd, nspotd), intent(inout) :: vspsmo
    real (kind=dp), dimension(irmd, nspotd), intent(inout) :: vspsme
    real (kind=dp), dimension(irmd, lmpotd, *), intent(inout) :: vons !! output potential (nonspherical VONS)
    ! .. Local variables
    real (kind=dp) :: fac, rmserm, rmserq, vmn, vnm, vnp, voldm, voldp, vpn
    real (kind=dp) :: natom
    integer :: i, ih, ihp1, irc1, irmin1, j, lm, np
    ! ---> final construction of the potentials
    ! attention : the spherical averaged potential is the lm=1
    intrinsic :: mod, real, sqrt
    ! component of vons times sqrt(4 pi).
    rmsavq = 0.0e0_dp
    rmsavm = 0.0e0_dp

    ! first mixing scheme : straight mixing
    ! ---> determination of the root mean sqare error

    natom = 0.0e0_dp

    do np = nstart, nend

      i = np - natref
      natom = natom + real(nshell(i), kind=dp)*conc(i)
      ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      if (nspin==2) then
        ih = 2*np - 1
        ihp1 = ih + 1
      else
        ih = np
        ihp1 = ih
      end if

      irc1 = irc(np)
      rmserq = 0.0e0_dp
      rmserm = 0.0e0_dp
      fac = 0.5e0_dp/rfpi
      ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      do j = 1, irc1
        vnp = fac*(vons(j,1,ih)+vons(j,1,ihp1))
        vnm = fac*(vons(j,1,ih)-vons(j,1,ihp1))
        voldp = 0.5e0_dp*(visp(j,ih)+visp(j,ihp1))
        voldm = 0.5e0_dp*(visp(j,ih)-visp(j,ihp1))
        rmserq = rmserq + 2.0e0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r(j, np)*drdi(j, np)*(vnp-voldp)*(vnp-voldp)
        rmserm = rmserm + 2.0e0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r(j, np)*drdi(j, np)*(vnm-voldm)*(vnm-voldm)
        vpn = voldp + mixing*(vnp-voldp)
        vmn = voldm + fcm*mixing*(vnm-voldm)
        vons(j, 1, ihp1) = vpn - vmn
        vons(j, 1, ih) = vpn + vmn
      end do


      if (lsmear>=3) then
        do j = 1, irc1
          vnp = 0.5e0_dp*(vspsmo(j,ih)+vspsmo(j,ihp1))
          vnm = 0.5e0_dp*(vspsmo(j,ih)-vspsmo(j,ihp1))
          voldp = 0.5e0_dp*(vspsme(j,ih)+vspsme(j,ihp1))
          voldm = 0.5e0_dp*(vspsme(j,ih)-vspsme(j,ihp1))
          vpn = voldp + mixing*(vnp-voldp)
          vmn = voldm + fcm*mixing*(vnm-voldm)
          vspsmo(j, ihp1) = vpn - vmn
          vspsmo(j, ih) = vpn + vmn
        end do
      end if

      if ((lsmear==1) .or. (lsmear==2)) then
        do j = 1, irc1
          vspsme(j, ihp1) = vspsmo(j, ihp1)
          vspsme(j, ih) = vspsmo(j, ih)
        end do
      end if


      rmserq = rmserq/(r(irc1,np)**3)
      rmserm = rmserm/(r(irc1,np)**3)
      rmsavq = rmsavq + rmserq*nshell(i)*conc(i)
      rmsavm = rmsavm + rmserm*nshell(i)*conc(i)

      if (nspin==2) then
        write (ipf, fmt=100) i, sqrt(rmserq), sqrt(rmserm)
      else
        write (ipf, fmt=120) i, sqrt(rmserq)
      end if

      if (ins/=0 .and. lpot>0) then

        rmserq = 0.0e0_dp
        rmserm = 0.0e0_dp
        irmin1 = irmin(np)
        do lm = 2, lmpot
          do j = irmin1, irc1
            vnp = 0.5e0_dp*(vons(j,lm,ih)+vons(j,lm,ihp1))
            vnm = 0.5e0_dp*(vons(j,lm,ih)-vons(j,lm,ihp1))
            voldp = 0.5e0_dp*(vins(j,lm,ih)+vins(j,lm,ihp1))
            voldm = 0.5e0_dp*(vins(j,lm,ih)-vins(j,lm,ihp1))
            rmserq = rmserq + 2.0e0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r(j, np)*drdi(j, np)*(vnp-voldp)*(vnp-voldp)
            rmserm = rmserm + 2.0e0_dp*real(1+mod(j,2), kind=dp)*r(j, np)*r(j, np)*drdi(j, np)*(vnm-voldm)*(vnm-voldm)
            vpn = voldp + mixing*(vnp-voldp)
            vmn = voldm + fcm*mixing*(vnm-voldm)
            vons(j, lm, ihp1) = vpn - vmn
            vons(j, lm, ih) = vpn + vmn
          end do
        end do
        rmserq = rmserq/(r(irc1,np)**3)/fpi
        rmserm = rmserm/(r(irc1,np)**3)/fpi
        rmsavq = rmsavq + rmserq*nshell(i)*conc(i)
        rmsavm = rmsavm + rmserm*nshell(i)*conc(i)

        if (nspin==2) then
          write (ipf, fmt=110) i, sqrt(rmserq), sqrt(rmserm)
        else
          write (ipf, fmt=130) i, sqrt(rmserq)
        end if

      end if

    end do
    ! 13.10.95 ***************************************************************
    ! ************************************************************************
    rmsavq = sqrt(rmsavq/natom)
    rmsavm = sqrt(rmsavm/natom)
    ! .. Parameters ..
    write (1337, '(79("-"),/)')
    if (nspin==2) then
      write (ipf, fmt=140) itc, rmsavq, rmsavm
      write (6, fmt=140) itc, rmsavq, rmsavm
    else
      write (ipf, fmt=150) itc, rmsavq
      write (6, fmt=150) itc, rmsavq
    end if
    write (1337, '(79("-"))')
    ! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
100 format (5x, ' rms-error for atom', i3, 1x, ':', 'v+ + v- = ', 1p, d11.4, 2x, ',', 2x, 'v+ - v- = ', 1p, d11.4)
110 format (5x, ' rms-error non spherical contribution for atom ', i3, 1x, ':', 'v+ + v- = ', 1p, d11.4, 02x, ',', 2x, 'v+ - v- = ', 1p, d11.4)
120 format (5x, ' rms-error for atom', i3, 1x, ':', 'v+ + v- = ', 1p, d11.4)
130 format (5x, ' rms-error non spherical contribution for atom ', i3, 1x, ':', 'v+ + v- = ', 1p, d11.4)
140 format ('      ITERATION', i4, ' average rms-error : v+ + v- = ', 1p, d11.4, /, 39x, ' v+ - v- = ', 1p, d11.4)
150 format ('      ITERATION', i4, ' average rms-error : v+ + v- = ', 1p, d11.4)
  end subroutine mixstr

end module mod_mixstr
