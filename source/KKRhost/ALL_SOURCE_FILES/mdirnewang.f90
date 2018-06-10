subroutine mdirnewang(it, nmvec, mvevi, mvphi, mvtet, mvgam, natypd, lmaxd, &
  nmvecmax)
!   ********************************************************************
!   *                                                                  *
!   *  this routine has been build up from the last part of the        *
!   *  original Munich CALCMVEC routine.                               *
!   *  After correcting MVEVI with the Fermi energy value MVEVIEF      *
!   *  (outside this routine) it calculates the new angles of the      *
!   *  LOCAL FRAME quantisation axis with respect to the GLOBAL FRAME  *
!   *                                                                  *
!   ********************************************************************
  implicit none

!Parameter definitions
  integer :: lmaxdloc
  parameter (lmaxdloc=8)

!Scalar Arguments
  integer :: it, nmvec, natypd, lmaxd, nmvecmax

!Array Arguments
  double complex :: mvevi(natypd, 3, nmvecmax)
  double precision :: mvphi(natypd, nmvecmax), mvtet(natypd, nmvecmax), &
    mvgam(natypd, nmvecmax)

!Local Scalars
  double precision :: mv, mvx, mvxy, mvy, mvz, pi
  integer :: i, imv, icall

!Local Arrays
  double precision :: mvglo(3, nmvecmax)

!Intrinsic Functions
  intrinsic :: abs, atan

!Data Statements
  data icall/0/

!Save Statements
  save :: icall, pi

  icall = icall + 1
!=======================================================================
  if (icall==1) then

    if (lmaxd>lmaxdloc) then
      write (6, *)
      write (6, *) ' Please increase parameter LMAXDLOC to ', lmaxd
      write (6, *) ' in the < MVECGLOBAL > routine.'
      stop ' < TBKKR2 > '
    end if

    pi = 4.d0*atan(1.d0)

  end if
!=======================================================================

  do imv = 1, nmvec

    do i = 1, 3
      mvglo(i, imv) = dimag(mvevi(it,i,imv))
    end do

    mvx = mvglo(1, imv)
    mvy = mvglo(2, imv)
    mvz = mvglo(3, imv)

    mv = sqrt(mvx**2+mvy**2+mvz**2)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (mv<1d-8) then
      mvphi(it, imv) = 0d0
      mvtet(it, imv) = 0d0
      mvgam(it, imv) = 0d0
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      mvxy = sqrt(mvx**2+mvy**2)
! ----------------------------------------------------------------------
      if (abs(mvxy)<1d-8) then
        mvphi(it, imv) = 0d0
! ----------------------------------------------------------------------
      else
! ----------------------------------------------------------------------
        if (mvy>=0d0) then
          mvphi(it, imv) = acos(mvx/mvxy)
        else if (mvx<0d0) then
          mvphi(it, imv) = pi + acos(-mvx/mvxy)
        else
          mvphi(it, imv) = 2*pi - acos(mvx/mvxy)
        end if
        mvphi(it, imv) = mvphi(it, imv)*180d0/pi
        if (abs(mvphi(it,imv)-360.0d0)<1d-8) mvphi(it, imv) = 0d0
      end if
! ----------------------------------------------------------------------
      if (mvphi(it,imv)>=345.d0) mvphi(it, imv) = 360.d0 - mvphi(it, imv)
      mvtet(it, imv) = acos(mvz/mv)*180d0/pi
      mvgam(it, imv) = 0d0
    end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do
!=======================================================================
end subroutine
