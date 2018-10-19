!------------------------------------------------------------------------------------
!> Summary: Calculates angles of the **local frame** with respect to the **global frame**
!> Author: 
!> After correcting `MVEVI` with the Fermi energy value `MVEVIEF` (outside 
!> this routine) it calculates the new angles of the **local frame** quantization 
!> axis with respect to the **global frame** 
!------------------------------------------------------------------------------------
module mod_mdirnewang
  use :: mod_datatypes, only: dp
  use :: constants, only : pi
  private :: dp

contains

  !-------------------------------------------------------------------------------  
  !> Summary: Calculates angles of the **local frame** with respect to the **global frame**
  !> Author:
  !> Category: physical-observables, KKRhost 
  !> Deprecated: False 
  !> After correcting `MVEVI` with the Fermi energy value `MVEVIEF` (outside 
  !> this routine) it calculates the new angles of the **local frame** quantization 
  !> axis with respect to the **global frame**
  !-------------------------------------------------------------------------------  
  !> @note This routine has been build up from the last part of the original 
  !> Munich `CALCMVEC()` routine.
  !> @endnote                                                                       
  !------------------------------------------------------------------------------- 
  subroutine mdirnewang(it,nmvec,mvevi,mvphi,mvtet,mvgam,natypd,lmaxd,nmvecmax)

    implicit none

    ! Parameter definitions
    integer :: lmaxdloc
    parameter (lmaxdloc=8)

    ! Scalar Arguments
    integer, intent(in) :: it     !! Current atom type
    integer, intent(in) :: nmvec 
    integer, intent(in) :: natypd !! Number of kinds of atoms in unit cell
    integer, intent(in) :: lmaxd  !! Maximum l component in wave function expansion
    integer, intent(in) :: nmvecmax

    ! Array Arguments
    complex (kind=dp), dimension(natypd, 3, nmvecmax), intent(in) :: mvevi
    real (kind=dp), dimension(natypd,nmvecmax), intent(out) :: mvphi
    real (kind=dp), dimension(natypd,nmvecmax), intent(out) :: mvtet
    real (kind=dp), dimension(natypd,nmvecmax), intent(out) :: mvgam

    ! Local Scalars
    real (kind=dp) :: mv, mvx, mvxy, mvy, mvz
    integer :: i, imv, icall

    ! Local Arrays
    real (kind=dp), dimension(3,nmvecmax) :: mvglo

    ! Intrinsic Functions
    intrinsic :: abs, atan

    ! Data Statements
    data icall/0/

    ! Save Statements
    save :: icall

    icall = icall + 1
    ! =======================================================================
    if (icall==1) then

      if (lmaxd>lmaxdloc) then
        write (6, *)
        write (6, *) ' Please increase parameter LMAXDLOC to ', lmaxd
        write (6, *) ' in the < MVECGLOBAL > routine.'
        stop ' < TBKKR2 > '
      end if

    end if
    ! =======================================================================

    do imv = 1, nmvec

      do i = 1, 3
        mvglo(i, imv) = aimag(mvevi(it,i,imv))
      end do

      mvx = mvglo(1, imv)
      mvy = mvglo(2, imv)
      mvz = mvglo(3, imv)

      mv = sqrt(mvx**2+mvy**2+mvz**2)
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (mv<1e-8_dp) then
        mvphi(it, imv) = 0e0_dp
        mvtet(it, imv) = 0e0_dp
        mvgam(it, imv) = 0e0_dp
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      else
        ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        mvxy = sqrt(mvx**2+mvy**2)
        ! ----------------------------------------------------------------------
        if (abs(mvxy)<1e-8_dp) then
          mvphi(it, imv) = 0e0_dp
          ! ----------------------------------------------------------------------
        else
          ! ----------------------------------------------------------------------
          if (mvy>=0e0_dp) then
            mvphi(it, imv) = acos(mvx/mvxy)
          else if (mvx<0e0_dp) then
            mvphi(it, imv) = pi + acos(-mvx/mvxy)
          else
            mvphi(it, imv) = 2*pi - acos(mvx/mvxy)
          end if
          mvphi(it, imv) = mvphi(it, imv)*180e0_dp/pi
          if (abs(mvphi(it,imv)-360.0e0_dp)<1e-8_dp) mvphi(it, imv) = 0e0_dp
        end if
        ! ----------------------------------------------------------------------
        if (mvphi(it,imv)>=345.e0_dp) mvphi(it, imv) = 360.e0_dp - mvphi(it, imv)
        mvtet(it, imv) = acos(mvz/mv)*180e0_dp/pi
        mvgam(it, imv) = 0e0_dp
      end if
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end do
    ! =======================================================================
  end subroutine mdirnewang

end module mod_mdirnewang
