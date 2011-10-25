      SUBROUTINE RHOCORE(EBOT,NSRA,ISPIN,NSPIN,I1,DRDI,R,VISP,A,B,ZAT,
     &              IRCUT,RHOC,QC,ECORE,NCORE,LCORE,
C                   new input parameters after inc.p removal
     &              irmd, ipand)
C
      IMPLICIT NONE

      INTEGER irmd
      INTEGER ipand
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,B,ZAT,EBOT
      INTEGER I1,ISPIN,NCORE,NSPIN,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD),ECORE(20),
     +                 R(IRMD),RHOC(IRMD,2),
     +                 VISP(IRMD)
      INTEGER IRCUT(0:IPAND),
     +        LCORE(20)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION QC,QC1,RMAX
      INTEGER IPR,NR
C     ..
C     .. External Subroutines ..
      EXTERNAL COREL
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
C     
         CALL COREL(NSRA,IPR,I1,RHOC(1,ISPIN),VISP,ECORE,LCORE,NCORE,
     +        DRDI,ZAT,QC1,A,B,ISPIN,NSPIN,NR,RMAX,IRMD,EBOT)
C     
         IF (IPR.NE.0) WRITE (6,FMT=99001) I1
         QC = QC + QC1

      RETURN
99001 FORMAT (1x,5('*'),' core-relaxation for ',i3,'th cell',
     &        ' was done ',5('*'))
      END
