      SUBROUTINE VXCDRV(EXC,KTE,KXC,LPOT,NSPIN,IATYP,RHO2NS,VONS,
     +                  R,DRDI,A,IRWS,IRCUT,IPAN,ICELL,GSH,ILM,
     +                  IMAXSH,IFUNM,THETAS,LMSP)
      IMPLICIT NONE
      INCLUDE 'inc.p'
C     .. Parameters ..
      INTEGER IJD,LMPOTD,LMXSPD
      PARAMETER (LMPOTD= (LPOTD+1)**2,LMXSPD= (2*LPOTD+1)**2)
      PARAMETER (IJD = 434)
C     ..
C     .. Scalar Arguments ..
      INTEGER KTE,KXC,LPOT,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NAEZD),DRDI(IRMD,*),EXC(0:LPOTD),GSH(*),
     +                 R(IRMD,*),RHO2NS(IRMD,LMPOTD,2),
     +                 THETAS(IRID,NFUND,*),VONS(IRMD,LMPOTD,2)
      INTEGER IFUNM(*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
     +        IRCUT(0:IPAND,*),IRWS(*),LMSP(*)
C     ..
C     .. External Subroutines ..
      EXTERNAL DCOPY,SPHERE_GGA,SPHERE_NOGGA,VXCGGA,VXCLM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DYLMF1(IJD,LMPOTD),DYLMF2(IJD,LMPOTD),
     +                 DYLMT1(IJD,LMPOTD),DYLMT2(IJD,LMPOTD),
     +                 DYLMTF(IJD,LMPOTD),
     +                 RIJ(IJD,3),THET(IJD),WTYR(IJD,LMPOTD),
     +                 YLM(IJD,LMPOTD),YR(IJD,LMPOTD)
      INTEGER IFUNMIAT(LMXSPD)
C     ..
C     .. Local Scalars ..
      INTEGER IATYP,ICELL,LMX1
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
     +               YR,WTYR,IJD,LMSP)
        ELSE
c
c GGA EX-COR POTENTIAL
c
          CALL VXCGGA(EXC,KTE,KXC,LPOT,NSPIN,IATYP,RHO2NS,
     +                VONS,R(1,IATYP),DRDI(1,IATYP),A(IATYP),
     +                IRWS(IATYP),IRCUT(0,IATYP),IPAN(IATYP),
     +                GSH,ILM,IMAXSH,IFUNM,THETAS(1,1,ICELL),
     +                YR,WTYR,IJD,LMSP,THET,YLM,DYLMT1,DYLMT2,
     +                DYLMF1,DYLMF2,DYLMTF)
        END IF
      END
