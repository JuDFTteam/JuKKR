    Subroutine rhocore(nsra, ispin, nspin, i1, drdi, r, visp, a, b, zat, &
      ircut, rhoc, ecore, ncore, lcore, cscl, vtrel, btrel, rmrel, drdirel, &
      r2drdirel, zrel, jwsrel, irshift, ecorerel, nkcore, kapcore)
      Use mod_datatypes, Only: dp
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
      Implicit None
!.. Parameters ..
      Include 'inc.p'
! ===================================================================
!  RELATIVISTIC TREATMENT OF CORE ELECTRONS   July/2002
      Real (Kind=dp) :: a, b, zat
      Integer :: jwsrel, zrel, irshift
      Integer :: i1, ispin, ncore, nspin, nsra
!  SEE ROUTINE <DRVCORE> FOR A SHORT DESCRIPTION OF THE VARIABLES
!
      Real (Kind=dp) :: drdi(irmd), ecore(20*(krel+1)), r(irmd), &
        rhoc(irmd, 2), visp(irmd)
      Integer :: ircut(0:ipand), lcore(20*(krel+1))
! ===================================================================
!..
!.. Local Scalars ..
!..
      Real (Kind=dp) :: ecorerel(krel*20+(1-krel), 2)
      Integer :: nkcore(20), kapcore(20*2)
      Real (Kind=dp) :: cscl(krel*lmaxd+1)
      Real (Kind=dp) :: vtrel(irmd*krel+(1-krel))
      Real (Kind=dp) :: btrel(irmd*krel+(1-krel))
      Real (Kind=dp) :: drdirel(irmd*krel+(1-krel)), &
        r2drdirel(irmd*krel+(1-krel)), rmrel(irmd*krel+(1-krel))
!.. External Subroutines ..
! --------------------------------------------------------------
!     ipr=0 : do not write state dependent information
      Real (Kind=dp) :: qc, qc1, rmax
      Integer :: ipr, nr
      Save :: qc
!     ipr=1 : write something
!     ipr=2 : write all (for debugging)
      External :: corel, drvcore
! --------------------------------------------------------------

!=======================================================================
! non/scalar-relativistic OR relativistic

      ipr = 0

      If (ispin==1) qc = 0.0E0_dp
      nr = ircut(1)
      rmax = r(nr)


!=======================================================================
      If (krel==0) Then
!=======================================================================

        Call corel(nsra, ipr, i1, rhoc(1,ispin), visp, ecore, lcore, ncore, &
          drdi, zat, qc1, a, b, ispin, nspin, nr, rmax, irmd)
! non/scalar-relativistic OR relativistic
        If (ipr/=0) Write (1337, Fmt=100) i1
        qc = qc + qc1
        If (ispin==nspin) Write (1337, Fmt=110) zat, qc
!=======================================================================
      Else
! *********************************************************************
        Call drvcore(ipr, i1, lcore, ncore, cscl, vtrel, btrel, rmrel, a, b, &
          drdirel, r2drdirel, zrel, jwsrel, irshift, rhoc, ecorerel, nkcore, &
          kapcore, ecore, lmaxd, irmd)
      End If
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
      Return
100   Format (1X, 5('*'), ' core-relaxation for ', I3, 'th cell', &
        ' was done ', 5('*'))
110   Format (4X, 'nuclear charge  ', F10.6, 9X, 'core charge =   ', F10.6)
    End Subroutine
