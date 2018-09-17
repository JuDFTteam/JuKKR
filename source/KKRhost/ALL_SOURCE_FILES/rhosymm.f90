module mod_rhosymm

contains

  ! -------------------------------------------------------------------------------
  ! SUBROUTINE: RHOSYMM
  ! > @brief Symmetrize the charge densities and magnetic moments of
  ! > atoms which are magnetic 'antisymmetric'
  ! > (dependencies in IXIPOL(*))

  ! > @author P. Zahn
  ! > @date Aug. 1996
  ! > @note -Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to
  ! Fortran90
  ! -------------------------------------------------------------------------------
  subroutine rhosymm(lmpot, nspin, nstart, nend, rho2ns, ixipol, irws, ircut, ipan, kshape, natyp, irm)

    use :: global_variables
    use :: mod_datatypes, only: dp

    implicit none
    ! .. Input variables

    integer, intent (in) :: irm    ! < Maximum number of radial points
    integer, intent (in) :: nend
    integer, intent (in) :: natyp  ! < Number of kinds of atoms in unit cell
    integer, intent (in) :: lmpot  ! < (LPOT+1)**2
    integer, intent (in) :: nspin  ! < Counter for spin directions
    integer, intent (in) :: kshape ! < Exact treatment of WS cell
    integer, intent (in) :: nstart
    integer, dimension (*), intent (in) :: ipan ! < Number of panels in
    ! non-MT-region
    integer, dimension (*), intent (in) :: irws ! < R point at WS radius
    integer, dimension (*), intent (in) :: ixipol ! < Constraint of spin pol.
    integer, dimension (0:ipand, *), intent (in) :: ircut ! < R points of panel
    ! borders
    ! .. In/Out variables
    real (kind=dp), dimension (irm, lmpot, natyp, *), intent (inout) :: rho2ns
    ! < radial density
    ! .. Local variables
    integer :: i, iatyp, iatyp1, irc, irc1, lm
    real (kind=dp) :: fac
    ! .. Intrinsic Functions
    intrinsic :: abs
    ! ----------------------------------------------------------------------------

    do iatyp = nstart, nend

      iatyp1 = abs(ixipol(iatyp))

      fac = 1.e0_dp
      if (ixipol(iatyp)<0) fac = -1.e0_dp

      if (iatyp1>=iatyp) then
        write (1337, *) 'Symmetrize atom ', iatyp, ' with ', iatyp1, '.'
        if (kshape/=0) then
          irc = ircut(ipan(iatyp), iatyp)
          irc1 = ircut(ipan(iatyp1), iatyp1)
        else
          irc = irws(iatyp)
          irc1 = irws(iatyp1)
        end if

        if (irc/=irc1) then
          write (6, *) 'Error in RHOSYMM : ***********************'
          write (6, *) 'Radial mesh of atoms ', iatyp, ' and ', iatyp1, ' are not equal.'
        end if

        do lm = 1, lmpot
          do i = 1, irc1
            rho2ns(i, lm, iatyp, 1) = (rho2ns(i,lm,iatyp,1)+rho2ns(i,lm,iatyp1,1))/2.e0_dp
            rho2ns(i, lm, iatyp1, 1) = rho2ns(i, lm, iatyp, 1)
            if (nspin>1) then
              rho2ns(i, lm, iatyp, 2) = (rho2ns(i,lm,iatyp,2)+fac*rho2ns(i,lm,iatyp1,2))/2.e0_dp
              rho2ns(i, lm, iatyp1, 2) = fac*rho2ns(i, lm, iatyp, 2)
            end if
          end do                   ! I =1,IRC1
        end do                     ! LM =1,LMPOT
      end if                       ! (IATYP1.GT.IATYP)
    end do                         ! IATYP=NSTART,NEND

    return

  end subroutine rhosymm

end module mod_rhosymm
