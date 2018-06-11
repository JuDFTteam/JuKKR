    Subroutine mdirnewang(it, nmvec, mvevi, mvphi, mvtet, mvgam, natypd, &
      lmaxd, nmvecmax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *  this routine has been build up from the last part of the        *
!   *  original Munich CALCMVEC routine.                               *
!   *  After correcting MVEVI with the Fermi energy value MVEVIEF      *
!   *  (outside this routine) it calculates the new angles of the      *
!   *  LOCAL FRAME quantisation axis with respect to the GLOBAL FRAME  *
!   *                                                                  *
!   ********************************************************************
      Implicit None

!Parameter definitions
      Integer :: lmaxdloc
      Parameter (lmaxdloc=8)

!Scalar Arguments
      Integer :: it, nmvec, natypd, lmaxd, nmvecmax

!Array Arguments
      Complex (Kind=dp) :: mvevi(natypd, 3, nmvecmax)
      Real (Kind=dp) :: mvphi(natypd, nmvecmax), mvtet(natypd, nmvecmax), &
        mvgam(natypd, nmvecmax)

!Local Scalars
      Real (Kind=dp) :: mv, mvx, mvxy, mvy, mvz, pi
      Integer :: i, imv, icall

!Local Arrays
      Real (Kind=dp) :: mvglo(3, nmvecmax)

!Intrinsic Functions
      Intrinsic :: abs, atan

!Data Statements
      Data icall/0/

!Save Statements
      Save :: icall, pi

      icall = icall + 1
!=======================================================================
      If (icall==1) Then

        If (lmaxd>lmaxdloc) Then
          Write (6, *)
          Write (6, *) ' Please increase parameter LMAXDLOC to ', lmaxd
          Write (6, *) ' in the < MVECGLOBAL > routine.'
          Stop ' < TBKKR2 > '
        End If

        pi = 4.E0_dp*atan(1.E0_dp)

      End If
!=======================================================================

      Do imv = 1, nmvec

        Do i = 1, 3
          mvglo(i, imv) = aimag(mvevi(it,i,imv))
        End Do

        mvx = mvglo(1, imv)
        mvy = mvglo(2, imv)
        mvz = mvglo(3, imv)

        mv = sqrt(mvx**2+mvy**2+mvz**2)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        If (mv<1E-8_dp) Then
          mvphi(it, imv) = 0E0_dp
          mvtet(it, imv) = 0E0_dp
          mvgam(it, imv) = 0E0_dp
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        Else
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          mvxy = sqrt(mvx**2+mvy**2)
! ----------------------------------------------------------------------
          If (abs(mvxy)<1E-8_dp) Then
            mvphi(it, imv) = 0E0_dp
! ----------------------------------------------------------------------
          Else
! ----------------------------------------------------------------------
            If (mvy>=0E0_dp) Then
              mvphi(it, imv) = acos(mvx/mvxy)
            Else If (mvx<0E0_dp) Then
              mvphi(it, imv) = pi + acos(-mvx/mvxy)
            Else
              mvphi(it, imv) = 2*pi - acos(mvx/mvxy)
            End If
            mvphi(it, imv) = mvphi(it, imv)*180E0_dp/pi
            If (abs(mvphi(it,imv)-360.0E0_dp)<1E-8_dp) mvphi(it, imv) = 0E0_dp
          End If
! ----------------------------------------------------------------------
          If (mvphi(it,imv)>=345.E0_dp) mvphi(it, imv) = 360.E0_dp - &
            mvphi(it, imv)
          mvtet(it, imv) = acos(mvz/mv)*180E0_dp/pi
          mvgam(it, imv) = 0E0_dp
        End If
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      End Do
!=======================================================================
    End Subroutine
