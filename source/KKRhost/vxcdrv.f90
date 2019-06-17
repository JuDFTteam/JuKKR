!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Driver for the exchange-correlation potential and energy calculation
!> Author: 
!> Driver for the exchange-correlation potential and energy calculation. It wraps
!> all the different exchange-correlation potentials and make sure to call the 
!> appropriate subroutines depending on the type of exchange correlation potential
!> indicated in the `inputcard`
!------------------------------------------------------------------------------------
module mod_vxcdrv

contains

  !-------------------------------------------------------------------------------
  !> Summary: Driver for the exchange-correlation potential and energy calculation
  !> Author:
  !> Category: xc-potential, KKRhost 
  !> Deprecated: False 
  !> Driver for the exchange-correlation potential and energy calculation. It wraps
  !> all the different exchange-correlation potentials and make sure to call the 
  !> appropriate subroutines depending on the type of exchange correlation potential
  !> indicated in the `inputcard`
  !-------------------------------------------------------------------------------
  subroutine vxcdrv(exc,kte,kxc,lpot,nspin,nstart,nend,rho2ns,vons,r,drdi,a,irws,   &
    ircut,ipan,ntcell,kshape,gsh,ilm,imaxsh,ifunm,thetas,lmsp)

    use :: global_variables, only: lmpotd, ngshd, ipand, natypd, irmd, irid, nfund, lmxspd, krel, lpotd
    use :: mod_datatypes, only: dp
    use :: mod_sphere_nogga, only: sphere_nogga 
    use :: mod_sphere_gga, only: sphere_gga
    use :: mod_vxcgga, only: vxcgga
    use :: mod_vxclm, only: vxclm
    use :: mod_wunfiles, only: t_params
    implicit none
    ! Parameters ..
    integer :: ijd 
    parameter (ijd=434)

    ! .. Input variables 
    integer, intent(in) :: kte    !! Calculation of the total energy On/Off (1/0)
    integer, intent(in) :: kxc    !! Type of xc-potential 0=vBH 1=MJW 2=VWN 3=PW91
    integer, intent(in) :: nend
    integer, intent(in) :: lpot   !! Maximum l component in potential expansion
    integer, intent(in) :: nspin  !! Counter for spin directions
    integer, intent(in) :: nstart
    integer, intent(in) :: kshape !! Exact treatment of WS cell
    integer, dimension(*), intent(in) :: irws !! R point at WS radius
    integer, dimension(*), intent(in) :: ipan !! Number of panels in non-MT-region
    integer, dimension(*), intent(in) :: ntcell !! Index for WS cell
    integer, dimension(0:lmpotd), intent(in)    :: imaxsh
    integer, dimension(ngshd, 3), intent(in)    :: ilm
    integer, dimension(0:ipand, *), intent(in)  :: ircut !! R points of panel borders
    integer, dimension(natypd, *), intent(in)   :: lmsp  !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension(natypd, *), intent(in)   :: ifunm
    real (kind=dp), dimension(natypd), intent(in)   :: a
    real (kind=dp), dimension(*), intent(in)        :: gsh
    real (kind=dp), dimension(irmd, *), intent(in)  :: r    !! Set of real space vectors (in a.u.)
    real (kind=dp), dimension(irmd, *), intent(in)  :: drdi !! Derivative dr/di
    real (kind=dp), dimension(irid, nfund, *), intent(in)   :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    real (kind=dp), dimension(irmd, lmpotd, natypd, *), intent(in) :: rho2ns !! radial density
    ! .. In/Out variables
    real (kind=dp), dimension(0:lpotd, *), intent(inout) :: exc !! xc-energy
    real (kind=dp), dimension(irmd, lmpotd, *), intent(inout)  :: vons   !! output potential (nonspherical VONS)

    ! Local Arrays ..
    integer, dimension(lmxspd) :: ifunmiat, lmspiat
    real (kind=dp), dimension(ijd)              :: thet
    real (kind=dp), dimension(ijd, 3)           :: rij
    real (kind=dp), dimension(irmd, lmpotd, 2)  :: rho2iat
    real (kind=dp), dimension(ijd,lmpotd) :: dylmf1,dylmf2,dylmt1,dylmt2,dylmtf,wtyr,ylm,yr 
    integer :: kxc_tmp

    ! Local Scalars ..
    integer :: iatyp, icell, ipot, lmx1

    if (kxc<3) then
      call sphere_nogga(lpot,yr,wtyr,rij,ijd)
    else
      call sphere_gga(lpot,yr,wtyr,rij,ijd,lmpotd,thet,ylm,dylmt1,dylmt2,dylmf1,    &
        dylmf2,dylmtf)
    end if
    do iatyp = nstart, nend
      icell = ntcell(iatyp)
      ipot = nspin*(iatyp-1) + 1
      kxc_tmp = kxc
      if (kxc>=3 .and. abs(t_params%zat(iatyp))<1.0e-8_dp) then
        kxc_tmp = 2 ! use LDA also for empty sites
        write(1337, '(A)') 'WARNING: use LDA (VWN) for empty cells instead of GGA!'
      end if
      do lmx1 = 1, lmxspd
        ifunmiat(lmx1) = ifunm(icell, lmx1)
        lmspiat(lmx1) = lmsp(icell, lmx1)
      end do
      call dcopy(irmd*lmpotd, rho2ns(1,1,iatyp,1), 1, rho2iat(1,1,1), 1)
      if (nspin==2 .or. krel==1) then
        call dcopy(irmd*lmpotd, rho2ns(1,1,iatyp,2), 1, rho2iat(1,1,2), 1)
      end if
      if (kxc_tmp<3) then 
        call vxclm(exc,kte,kxc_tmp,lpot,nspin,iatyp,rho2iat,vons(1,1,ipot),r(1,iatyp),  &
          drdi(1,iatyp),irws(iatyp),ircut(0,iatyp),ipan(iatyp),kshape,gsh,ilm,      &
          imaxsh,ifunmiat,thetas(1,1,icell),yr,wtyr,ijd,lmspiat)
      else
        !----------------------------------------------------------------------------
        ! GGA EX-COR POTENTIAL
        !----------------------------------------------------------------------------
        call vxcgga(exc,kte,kxc_tmp,lpot,nspin,iatyp,rho2iat,vons(1,1,ipot),r(1,iatyp), &
          drdi(1,iatyp),a(iatyp),irws(iatyp),ircut(0,iatyp),ipan(iatyp),kshape,gsh, &
          ilm,imaxsh,ifunmiat,thetas(1,1,icell),wtyr,ijd,lmspiat,thet,ylm,dylmt1,   &
          dylmt2,dylmf1,dylmf2,dylmtf)
      end if
    end do
  end subroutine vxcdrv

end module mod_vxcdrv
