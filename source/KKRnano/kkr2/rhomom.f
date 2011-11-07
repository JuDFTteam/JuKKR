c 13.10.95 ***************************************************************
      SUBROUTINE RHOMOM(CMOM,CMINST,LMAX,IATYP,RHO2NS,R,
     +                  DRDI,IRCUT,IPAN,ICELL,ILM,IFUNM,
     +                  IMAXSH,GSH,THETAS,LMSP,
C                       new input parameters after inc.p removal
     &                  irmd, irid, nfund, ipand, ngshd)
c ************************************************************************
c     calculate charge moments of given charge densities 
c
c                             rcut
c              cmom(lm,i) =    s dr' r'** l rho2ns(r',lm,i,1)
c                              0
c-----------------------------------------------------------------------

      INTEGER irmd
      INTEGER irid
      INTEGER nfund
      INTEGER ipand
      INTEGER ngshd

C     INTEGER LMPOTD
C     PARAMETER (LMPOTD= (LPOTD+1)**2) ! = (2*LMAX+1)**2
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAX
C     ..
C     .. Array Arguments ..
C     DOUBLE PRECISION CMINST(LMPOTD),CMOM(LMPOTD),DRDI(IRMD,*),
C    +                 GSH(*),R(IRMD,*),RHO2NS(IRMD,LMPOTD),
C    +                 THETAS(IRID,NFUND,*)
C     INTEGER IFUNM(*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
C    +        IRCUT(0:IPAND,*),IRWS(*),LMSP(*)
      DOUBLE PRECISION CMINST((2*LMAX+1)**2)
      DOUBLE PRECISION CMOM((2*LMAX+1)**2)
      DOUBLE PRECISION DRDI(IRMD,*)
      DOUBLE PRECISION GSH(*)
      DOUBLE PRECISION R(IRMD,*)
      DOUBLE PRECISION RHO2NS(IRMD,(2*LMAX+1)**2)
      DOUBLE PRECISION THETAS(IRID,NFUND,*)
      INTEGER IFUNM(*)
      INTEGER ILM(NGSHD,3)
      INTEGER IMAXSH(0:(2*LMAX+1)**2)
      INTEGER IPAN(*)
      INTEGER IRCUT(0:IPAND,*)
      INTEGER LMSP(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FAC,PI,RL
      INTEGER I,IATYP,ICELL,IEND,IFUN,IRC1,IRS1,ISTART,J,L,LM,LM2,
     +        LM3,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION V1(IRMD),VINT1(IRMD)
      INTEGER IRCUTM(0:IPAND)
C     ..
C     .. External Subroutines ..
      EXTERNAL SINWK,SOUTK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,REAL
C     ..
      PI = 4.D0*ATAN(1.D0)
c
        IRS1 = IRCUT(1,IATYP)
        IRC1 = IRCUT(IPAN(IATYP),IATYP)
        DO 10 I = 0,IPAN(IATYP)
          IRCUTM(I) = IRCUT(I,IATYP)
   10   CONTINUE

        DO 100 L = 0,LMAX
          FAC = 8.0D0*PI/REAL(2*L+1)
          DO 90 M = -L,L
            LM = L*L + L + M + 1
c
c---> set up of the integrands v1 and v2
c
            V1(1) = 0.0D0
            DO 20 I = 2,IRS1
              RL = R(I,IATYP)**L
              V1(I) = RHO2NS(I,LM)*RL*DRDI(I,IATYP)
   20       CONTINUE
c
c---> convolute charge density of interstial with shape function
c
            DO 30 I = IRS1 + 1,IRC1
              V1(I) = 0.0D0
   30       CONTINUE
            ISTART = IMAXSH(LM-1) + 1
            IEND = IMAXSH(LM)
            DO 50 J = ISTART,IEND
              LM2 = ILM(J,2)
              LM3 = ILM(J,3)
              IF (LMSP(LM3).GT.0) THEN
                IFUN = IFUNM(LM3)
                DO 40 I = IRS1 + 1,IRC1
                  V1(I) = V1(I) + GSH(J)*RHO2NS(I,LM2)*
     +                    THETAS(I-IRS1,IFUN,ICELL)
   40           CONTINUE
              END IF
   50       CONTINUE

            DO 60 I = IRS1 + 1,IRC1
              RL = R(I,IATYP)**L
              V1(I) = V1(I)*RL*DRDI(I,IATYP)
   60       CONTINUE
c
c---> now integrate
c
            CALL SOUTK(V1,VINT1,IPAN(IATYP),IRCUTM)

            CMOM(LM) = VINT1(IRS1)
            CMINST(LM) = VINT1(IRC1) - VINT1(IRS1)

   90     CONTINUE

  100   CONTINUE

      RETURN

      END
