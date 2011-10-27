      SUBROUTINE VXCDRV(EXC,KTE,KXC,LPOT,NSPIN,IATYP,RHO2NS,VONS,
     +                  R,DRDI,A,IRWS,IRCUT,IPAN,ICELL,GSH,ILM,
     +                  IMAXSH,IFUNM,THETAS,LMSP,
C                       new input parameters after inc.p removal
     &                  naez, irmd, irid, nfund, ngshd, ipand)
      IMPLICIT NONE

      INTEGER naez
      INTEGER irmd
      INTEGER irid
      INTEGER nfund
      INTEGER ngshd
      INTEGER ipand

C     .. Parameters ..
      INTEGER IJD
C     INTEGER LMPOTD,LMXSPD
C     PARAMETER (LMPOTD= (LPOTD+1)**2,LMXSPD= (2*LPOTD+1)**2)

C     FIX: hardcoded
      PARAMETER (IJD = 434)
C     ..
C     .. Scalar Arguments ..
      INTEGER KTE,KXC,LPOT,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
C     DOUBLE PRECISION A(NAEZD),DRDI(IRMD,*),EXC(0:LPOTD),GSH(*),
C    +                 R(IRMD,*),RHO2NS(IRMD,LMPOTD,2),
C    +                 THETAS(IRID,NFUND,*),VONS(IRMD,LMPOTD,2)
C     INTEGER IFUNM(*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
C    +        IRCUT(0:IPAND,*),IRWS(*),LMSP(*)

      DOUBLE PRECISION A(NAEZ)
      DOUBLE PRECISION DRDI(IRMD,*)
      DOUBLE PRECISION EXC(0:LPOT)
      DOUBLE PRECISION GSH(*)
      DOUBLE PRECISION R(IRMD,*)
      DOUBLE PRECISION RHO2NS(IRMD,(LPOT+1)**2,2)
      DOUBLE PRECISION THETAS(IRID,NFUND,*)
      DOUBLE PRECISION VONS(IRMD,(LPOT+1)**2,2)
      INTEGER IFUNM(*)
      INTEGER ILM(NGSHD,3)
      INTEGER IMAXSH(0:(LPOT+1)**2)
      INTEGER IPAN(*)
      INTEGER IRCUT(0:IPAND,*)
      INTEGER IRWS(*)
      INTEGER LMSP(*)
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOPY,SPHERE_GGA,SPHERE_NOGGA,VXCGGA,VXCLM
C     ..
C     .. Local Arrays .. Fortran 90 automatic arrays
C     DOUBLE PRECISION DYLMF1(IJD,LMPOTD),DYLMF2(IJD,LMPOTD),
C    +                 DYLMT1(IJD,LMPOTD),DYLMT2(IJD,LMPOTD),
C    +                 DYLMTF(IJD,LMPOTD),
C    +                 RIJ(IJD,3),THET(IJD),WTYR(IJD,LMPOTD),
C    +                 YLM(IJD,LMPOTD),YR(IJD,LMPOTD)
C     INTEGER IFUNMIAT(LMXSPD)

      DOUBLE PRECISION DYLMF1(IJD,(LPOT+1)**2)
      DOUBLE PRECISION DYLMF2(IJD,(LPOT+1)**2)
      DOUBLE PRECISION DYLMT1(IJD,(LPOT+1)**2)
      DOUBLE PRECISION DYLMT2(IJD,(LPOT+1)**2)
      DOUBLE PRECISION DYLMTF(IJD,(LPOT+1)**2)
      DOUBLE PRECISION RIJ(IJD,3)
      DOUBLE PRECISION THET(IJD)
      DOUBLE PRECISION WTYR(IJD,(LPOT+1)**2)
      DOUBLE PRECISION YLM(IJD,(LPOT+1)**2)
      DOUBLE PRECISION YR(IJD,(LPOT+1)**2)
      INTEGER IFUNMIAT((2*LPOT+1)**2)
C     ..
C     .. Local Scalars ..
      INTEGER IATYP,ICELL,LMX1

      INTEGER LMPOTD

      LMPOTD = (LPOT+1)**2
C     ..
      IF (KXC.LT.3) THEN
        CALL SPHERE_NOGGA(LPOT,YR,WTYR,RIJ,IJD)
      ELSE
        CALL SPHERE_GGA(LPOT,YR,WTYR,RIJ,IJD,LMPOTD,THET,YLM,DYLMT1,
     +                  DYLMT2,DYLMF1,DYLMF2,DYLMTF)
      END IF

        IF (KXC.LT.3) THEN

          CALL VXCLM(EXC,KTE,KXC,LPOT,NSPIN,IATYP,RHO2NS,
     +               VONS,R(1,IATYP),DRDI(1,IATYP),
     +               IRWS(IATYP),IRCUT(0,IATYP),IPAN(IATYP),
     +               GSH,ILM,IMAXSH,IFUNM,THETAS(1,1,ICELL),
     +               YR,WTYR,IJD,LMSP,
     &               irmd, irid, nfund, ngshd, ipand)
        ELSE
c
c GGA EX-COR POTENTIAL
c
          CALL VXCGGA(EXC,KTE,KXC,LPOT,NSPIN,IATYP,RHO2NS,
     +                VONS,R(1,IATYP),DRDI(1,IATYP),A(IATYP),
     +                IRWS(IATYP),IRCUT(0,IATYP),IPAN(IATYP),
     +                GSH,ILM,IMAXSH,IFUNM,THETAS(1,1,ICELL),
     +                YR,WTYR,IJD,LMSP,THET,YLM,DYLMT1,DYLMT2,
     +                DYLMF1,DYLMF2,DYLMTF,
     &                irmd, irid, nfund, ngshd, ipand)
        END IF
      END
