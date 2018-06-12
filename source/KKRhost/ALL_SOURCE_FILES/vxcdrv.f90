! -------------------------------------------------------------------------------
! SUBROUTINE: VXCDRV
! > @brief Wrapper for the calculation of the exchange-correlation energy, for
! > the different treatments of the xc-potential.
! > @note
! > - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to
! Fortran90
! -------------------------------------------------------------------------------
subroutine vxcdrv(exc, kte, kxc, lpot, nspin, nstart, nend, rho2ns, vons, r, &
  drdi, a, irws, ircut, ipan, ntcell, kshape, gsh, ilm_map, imaxsh, ifunm, &
  thetas, lmsp, lmpot, natyp)

  use :: global_variables
  use :: mod_datatypes, only: dp

  implicit none

  ! .. Input variables
  integer, intent (in) :: kte      ! < Calculation of the total energy On/Off
                                   ! (1/0)
  integer, intent (in) :: kxc      ! < Type of xc-potential 0=vBH 1=MJW 2=VWN
                                   ! 3=PW91
  integer, intent (in) :: nend
  integer, intent (in) :: lpot     ! < Maximum l component in potential
                                   ! expansion
  integer, intent (in) :: natyp    ! < Number of kinds of atoms in unit cell
  integer, intent (in) :: lmpot    ! < (LPOT+1)**2
  integer, intent (in) :: nspin    ! < Counter for spin directions
  integer, intent (in) :: nstart
  integer, intent (in) :: kshape   ! < Exact treatment of WS cell
  ! .. Array Arguments ..
  real (kind=dp), dimension (natyp), intent (in) :: a ! < Constants for
                                                      ! exponential R mesh
  real (kind=dp), dimension (ngshd), intent (in) :: gsh
  real (kind=dp), dimension (irmd, natyp), intent (in) :: r ! < Radial mesh (
                                                            ! in units a Bohr)
  real (kind=dp), dimension (irmd, natyp), intent (in) :: drdi ! < Derivative
                                                               ! dr/di
  real (kind=dp), dimension (irid, nfund, ncelld), intent (in) :: thetas ! <
                                                                         ! shape
                                                                         ! function
                                                                         ! THETA=0
                                                                         ! outer
                                                                         ! space
                                                                         ! THETA
                                                                         ! =1
                                                                         ! inside
                                                                         ! WS
                                                                         ! cell
                                                                         ! in
                                                                         ! spherical
                                                                         ! harmonics
                                                                         ! expansion
  real (kind=dp), dimension (irmd, lmpot, natyp, 2), intent (in) :: rho2ns
  ! < radial density
  integer, dimension (natyp), intent (in) :: irws ! < R point at WS radius
  integer, dimension (natyp), intent (in) :: ipan ! < Number of panels in
                                                  ! non-MT-region
  integer, dimension (natyp), intent (in) :: ntcell ! < index for WS cell
  integer, dimension (0:lmpot), intent (in) :: imaxsh
  integer, dimension (ngshd, 3), intent (in) :: ilm_map
  integer, dimension (natyp, lmxspd), intent (in) :: lmsp ! < 0,1 :
                                                          ! non/-vanishing
                                                          ! lm=(l,m) component
                                                          ! of non-spherical
                                                          ! potential
  integer, dimension (natyp, lmxspd), intent (in) :: ifunm
  integer, dimension (0:ipand, natyp), intent (in) :: ircut ! < r points of
                                                            ! panel borders
  ! .. In/Out variables
  real (kind=dp), dimension (0:lpot, natyp), intent (inout) :: exc ! <
                                                                   ! exchange
                                                                   ! correlation
                                                                   ! energy
  ! .. Output variables
  real (kind=dp), dimension (irmd, lmpot, npotd), intent (out) :: vons ! <
                                                                       ! output
                                                                       ! potential
                                                                       ! (nonspherical
                                                                       ! VONS)
  ! ..
  ! .. External Subroutines
  external :: dcopy, sphere_gga, sphere_nogga, vxcgga, vxclm
  ! .. Local Scalars
  integer :: iatyp, icell, ipot, lmx1
  integer :: ijd
  ! .. Parameters
  parameter (ijd=434)
  ! .. Local Arrays ..
  integer, dimension (lmxspd) :: lmspiat
  integer, dimension (lmxspd) :: ifunmiat
  real (kind=dp), dimension (ijd) :: thet
  real (kind=dp), dimension (ijd, lmpot) :: yr
  real (kind=dp), dimension (ijd, 3) :: rij
  real (kind=dp), dimension (ijd, lmpot) :: ylm
  real (kind=dp), dimension (ijd, lmpot) :: wtyr
  real (kind=dp), dimension (ijd, lmpot) :: dylmf1
  real (kind=dp), dimension (ijd, lmpot) :: dylmf2
  real (kind=dp), dimension (ijd, lmpot) :: dylmt1
  real (kind=dp), dimension (ijd, lmpot) :: dylmt2
  real (kind=dp), dimension (ijd, lmpot) :: dylmtf
  real (kind=dp), dimension (irmd, lmpot, 2) :: rho2iat
  ! ..
  if (kxc<3) then
    call sphere_nogga(lpot, yr, wtyr, rij, ijd)
  else
    call sphere_gga(lpot, yr, wtyr, rij, ijd, lmpot, thet, ylm, dylmt1, &
      dylmt2, dylmf1, dylmf2, dylmtf)
  end if
  do iatyp = nstart, nend
    icell = ntcell(iatyp)
    ipot = nspin*(iatyp-1) + 1
    do lmx1 = 1, lmxspd
      ifunmiat(lmx1) = ifunm(icell, lmx1)
      lmspiat(lmx1) = lmsp(icell, lmx1)
    end do
    call dcopy(irmd*lmpot, rho2ns(1,1,iatyp,1), 1, rho2iat(1,1,1), 1)
    if (nspin==2 .or. krel==1) then
      call dcopy(irmd*lmpot, rho2ns(1,1,iatyp,2), 1, rho2iat(1,1,2), 1)
    end if
    if (kxc<3) then
      call vxclm(exc, kte, kxc, lpot, nspin, iatyp, rho2iat, vons(1,1,ipot), &
        r(1,iatyp), drdi(1,iatyp), irws(iatyp), ircut(0,iatyp), ipan(iatyp), &
        kshape, gsh, ilm_map, imaxsh, ifunmiat, thetas(1,1,icell), yr, wtyr, &
        ijd, lmspiat, lmpot, lmmaxd, lpot, natyp)
    else
      ! ----------------------------------------------------------------------
      ! GGA EX-COR POTENTIAL
      ! ----------------------------------------------------------------------
      call vxcgga(exc, kte, kxc, lpot, nspin, iatyp, rho2iat, vons(1,1,ipot), &
        r(1,iatyp), drdi(1,iatyp), a(iatyp), irws(iatyp), ircut(0,iatyp), &
        ipan(iatyp), kshape, gsh, ilm_map, imaxsh, ifunmiat, &
        thetas(1,1,icell), wtyr, ijd, lmspiat, thet, ylm, dylmt1, dylmt2, &
        dylmf1, dylmf2, dylmtf, lmpot, lmmaxd, lpot, natyp)
    end if
  end do
end subroutine vxcdrv
