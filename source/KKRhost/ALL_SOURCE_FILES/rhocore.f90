SUBROUTINE rhocore(nsra,ispin,nspin,i1,drdi,r,visp,a,b,zat,ircut,  &
    rhoc,ecore,ncore,lcore,cscl,vtrel,btrel,rmrel,  &
    drdirel,r2drdirel,zrel,jwsrel,irshift,ecorerel, nkcore,kapcore)
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
      IMPLICIT NONE
!.. Parameters ..
      INCLUDE 'inc.p'
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION A,B,ZAT
      INTEGER JWSREL,ZREL,IRSHIFT
      INTEGER I1,ISPIN,NCORE,NSPIN,NSRA
!..
!.. Array Arguments ..
DOUBLE PRECISION DRDI(IRMD),ECORE(20*(KREL+1)), &
                 R(IRMD),RHOC(IRMD,2), &
                 VISP(IRMD)
INTEGER IRCUT(0:IPAND), &
        LCORE(20*(KREL+1))
! ===================================================================
!  RELATIVISTIC TREATMENT OF CORE ELECTRONS   July/2002
!  SEE ROUTINE <DRVCORE> FOR A SHORT DESCRIPTION OF THE VARIABLES
!
DOUBLE PRECISION ECOREREL(KREL*20+(1-KREL),2)
INTEGER NKCORE(20),KAPCORE(20*2)
DOUBLE PRECISION CSCL(KREL*LMAXD+1)
DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL))
DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL))
DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL)), &
                 R2DRDIREL(IRMD*KREL+(1-KREL)), &
                 RMREL(IRMD*KREL+(1-KREL))
! ===================================================================
!..
!.. Local Scalars ..
      DOUBLE PRECISION QC,QC1,RMAX
      INTEGER IPR,NR
      SAVE QC
!..
!.. External Subroutines ..
      EXTERNAL COREL,DRVCORE
! --------------------------------------------------------------
!     ipr=0 : do not write state dependent information
!     ipr=1 : write something
!     ipr=2 : write all (for debugging)
! --------------------------------------------------------------
ipr = 0

IF ( ispin == 1 ) qc = 0.0D0
nr = ircut(1)
rmax = r(nr)
!=======================================================================
! non/scalar-relativistic OR relativistic

IF ( krel == 0 ) THEN
  
  
  CALL corel(nsra,ipr,i1,rhoc(1,ispin),visp,ecore,lcore,ncore,  &
      drdi,zat,qc1,a,b,ispin,nspin,nr,rmax,irmd)
  
  IF (ipr /= 0) WRITE (1337,FMT=99001) i1
  qc = qc + qc1
  IF (ispin == nspin) WRITE (1337,FMT=99002) zat,qc
!=======================================================================
ELSE
!=======================================================================
  CALL drvcore(ipr,i1,lcore,ncore,cscl,vtrel,btrel,rmrel,a,b,  &
      drdirel,r2drdirel,zrel,jwsrel,irshift,rhoc,  &
      ecorerel,nkcore,kapcore,ecore,lmaxd,irmd)
END IF

! non/scalar-relativistic OR relativistic
!=======================================================================
RETURN
99001 FORMAT (1X,5('*'),' core-relaxation for ',i3,'th cell',  &
    ' was done ',5('*'))
99002 FORMAT (4X,'nuclear charge  ',f10.6,9X,'core charge =   ',f10.6)
END SUBROUTINE rhocore
