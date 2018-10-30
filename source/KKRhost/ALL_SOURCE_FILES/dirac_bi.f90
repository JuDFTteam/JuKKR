!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------
!> Summary: EMPTY
!> Deprecated: True ! This needs to be set to True for deprecated subroutines
!>
!> @warning everythong is commented out! @ endwarning
!-------------------------------------------------------------------------------
module mod_dirac_bi


  !-------------------------------------------------------------------------------
  !> Summary: EMPTY
  !> Author: 
  !> Category: KKRhost, dirac
  !> Deprecated: True ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  ! SUBROUTINE DIRABMBI(GETIRRSOL,C,IT,E,L,MJ,KAP1,KAP2,PIS,CG1,CG2,
  ! &                    CG4,CG5,CG8,AMEBI1,AMEBI2,V,B,AT,Z,NUCLEUS,R,
  ! &                    DRDI,DOVR,NMESH,PR,QR,PI,QI,DP,DQ,AP,AQ,NTMAX,
  ! &                    NLAMAX,NKMMAX,NRMAX)
  ! C   ********************************************************************
  ! C   *                                                                  *
  ! C   *   ROUTINE TO SOLVE THE SPIN-POLARISED RADIAL DIRAC EQUATIONS     *
  ! C   *                                                                  *
  ! C   *            including a vector potential A(L,M)                   *
  ! C   *                                                                  *
  ! C   *   the outward integration is started by a power expansion        *
  ! C   *   and continued by ADAMS-BASHFORTH-MOULTON - pred./corr.-method  *
  ! C   *   NABM = 4(5) selects the 4(5)-point formula                     *
  ! C   *                                                                  *
  ! C   *   the inward integration is started analytically                 *
  ! C   *                                                                  *
  ! C   *   returns the wave functions up to the mesh point NMESH          *
  ! C   *   PR,QR and PI,QI  with   P=r*g and Q=r*c*f                      *
  ! C   *   and    R/I standing for regular/irregular solution             *
  ! C   *                                                                  *
  ! C   *  26/01/95  HE                                                    *
  ! C   ********************************************************************
  ! IMPLICIT COMPLEX*16(A-H,O-Z)
  ! C
  ! C
  ! C PARAMETER definitions
  ! C
  ! INTEGER NLMAXLOC,NTMAXLOC,NRMAXLOC,MPSMAX,NPEMAX,NABM
  ! PARAMETER (NLMAXLOC=4,NTMAXLOC=2,NRMAXLOC=750,MPSMAX=40,NPEMAX=4,
  ! &           NABM=4)
  ! COMPLEX*16 CZ
  ! PARAMETER (CZ=(0.0D0,0.0D0))
  ! REAL*8 TOL
  ! PARAMETER (TOL=1.0D-6)
  ! INTEGER ITMAX
  ! PARAMETER (ITMAX=50)
  ! C
  ! C Dummy arguments
  ! C
  ! REAL*8 C,CG1,CG2,CG4,CG5,CG8,MJ
  ! COMPLEX*16 E,PIS
  ! LOGICAL GETIRRSOL
  ! INTEGER IT,KAP1,KAP2,L,NKMMAX,NLAMAX,NMESH,NRMAX,NTMAX,NUCLEUS,Z
  ! REAL*8 AMEBI1(NKMMAX,NKMMAX,NLAMAX,-1:+1),
  ! &       AMEBI2(NKMMAX,NKMMAX,NLAMAX,-1:+1),AP(2,2,NRMAXLOC),
  ! &       AQ(2,2,NRMAXLOC),AT(NRMAXLOC,NLAMAX,-1:+1,NTMAXLOC),
  ! &       B(NRMAX),DOVR(NRMAX),DRDI(NRMAX),R(NRMAX),V(NRMAX)
  ! COMPLEX*16 DP(2,2,NRMAX),DQ(2,2,NRMAX),PI(2,2,NRMAX),PR(2,2,NRMAX)
  ! &           ,QI(2,2,NRMAX),QR(2,2,NRMAX)
  ! C
  ! ! C Local variables
  ! ! C
  ! !       REAL*8 ACORR0(0:NABM-1),APRED0(NABM),
  ! !      &       ASTEP
  ! ! C
  ! ! C     DATA APRED0 / 1901.0D0, -2774.0D0, 2616.0D0, -1274.0D0, 251.0D0 /
  ! ! C     DATA ACORR0 /  251.0D0,  +646.0D0, -264.0D0,  +106.0D0, -19.0D0 /
  ! ! C     DATA ASTEP  /  720.0D0 /
  ! !       DATA APRED0/55.0D0, - 59.0D0, + 37.0D0, - 9.0D0/
  ! !       DATA ACORR0/9.0D0, + 19.0D0, - 5.0D0, + 1.0D0/
  ! !       DATA ASTEP/24.0D0/
  ! C
  ! C#######################################################################

  ! stop ' < DIRABMBI > : Not implemented. Set SOLVER=BS in inputcard'

  ! END

end module mod_dirac_bi
