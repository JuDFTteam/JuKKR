C*==rhocore.f    processed by SPAG 6.05Rc at 11:29 on 10 May 2004
      SUBROUTINE RHOCORE(NSRA,ISPIN,NSPIN,I1,DRDI,R,VISP,A,B,ZAT,IRCUT,
     &                   RHOC,ECORE,NCORE,LCORE,CSCL,VTREL,BTREL,RMREL,
     &                   DRDIREL,R2DRDIREL,ZREL,JWSREL,IRSHIFT,ECOREREL,
     &                   NKCORE,KAPCORE)
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *                                                                   *
C *********************************************************************
C
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'inc.p'
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,ZAT
      INTEGER JWSREL,ZREL,IRSHIFT
      INTEGER I1,ISPIN,NCORE,NSPIN,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD),ECORE(20*(KREL+1)),
     +                 R(IRMD),RHOC(IRMD,2),
     +                 VISP(IRMD)
      INTEGER IRCUT(0:IPAND),
     +        LCORE(20*(KREL+1))
C ===================================================================
C  RELATIVISTIC TREATMENT OF CORE ELECTRONS   July/2002
C  SEE ROUTINE <DRVCORE> FOR A SHORT DESCRIPTION OF THE VARIABLES
C
      DOUBLE PRECISION ECOREREL(KREL*20+(1-KREL),2)
      INTEGER NKCORE(20),KAPCORE(20*2)
      DOUBLE PRECISION CSCL(KREL*LMAXD+1)
      DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL))
      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL))
      DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL)),
     &                 R2DRDIREL(IRMD*KREL+(1-KREL)),
     &                 RMREL(IRMD*KREL+(1-KREL))
C ===================================================================
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION QC,QC1,RMAX
      INTEGER IPR,NR
      SAVE QC
C     ..
C     .. External Subroutines ..
      EXTERNAL COREL,DRVCORE
C --------------------------------------------------------------
C     ipr=0 : do not write state dependent information
C     ipr=1 : write something
C     ipr=2 : write all (for debugging)
C --------------------------------------------------------------
      IPR = 0
C     
      IF ( ISPIN.EQ.1 ) QC = 0.0D0
      NR = IRCUT(1)
      RMAX = R(NR)
C=======================================================================
C non/scalar-relativistic OR relativistic
C
      IF ( KREL.EQ.0 ) THEN
C     

         CALL COREL(NSRA,IPR,I1,RHOC(1,ISPIN),VISP,ECORE,LCORE,NCORE,
     +        DRDI,ZAT,QC1,A,B,ISPIN,NSPIN,NR,RMAX,IRMD)
C     
         IF (IPR.NE.0) WRITE (1337,FMT=99001) I1
         QC = QC + QC1
         IF (ISPIN.EQ.NSPIN) WRITE (1337,FMT=99002) ZAT,QC
C=======================================================================
      ELSE
C=======================================================================
         CALL DRVCORE(IPR,I1,LCORE,NCORE,CSCL,VTREL,BTREL,RMREL,A,B,
     &                DRDIREL,R2DRDIREL,ZREL,JWSREL,IRSHIFT,RHOC,
     &                ECOREREL,NKCORE,KAPCORE,ECORE,LMAXD,IRMD)
      END IF
C
C non/scalar-relativistic OR relativistic
C=======================================================================
      RETURN
99001 FORMAT (1x,5('*'),' core-relaxation for ',i3,'th cell',
     &        ' was done ',5('*'))
99002 FORMAT (4X,'nuclear charge  ',F10.6,9X,'core charge =   ',F10.6)
      END
