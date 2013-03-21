c 13.10.95 ***************************************************************
      SUBROUTINE VINTRAS_NEW(LPOT,NSPIN,RHO2NS,VONS,R,
     +                   DRDI,IRCUT,IPAN,ILM,IFUNM,
     +                   IMAXSH,GSH,THETAS,LMSP,
     &                   irmd, irid, nfund, ngshd, ipand)
c ************************************************************************
c     calculate the electron-intracell-potentials and the charge-
c     moments of given charge densities . ( for each spin-direc-
c     tion the potential is the same in the polarized case . )
c     initialize the potential v with the electron-intracell-potentials
c     the intracell-potential is expanded into spherical harmonics .
c     the lm-term of the intracell-potential of the representive atom i
c     is given by
c                    8pi        r      r'** l
c      v(r,lm,i) =  ----- *  (  s dr' --------   rho2ns(r',lm,i,1)
c                   2*l+1       0     r **(l+1)
c
c                                 rcut    r ** l
c                               +  s dr' ---------   rho2ns(r',lm,i,1) )
c                                  r     r' **(l+1)
c
c             (see notes by b.drittler and u.klemradt)
c
c     attention : rho2ns(...,1) is the real charge density times r**2
c                 developed into spherical harmonics . (see deck rholm)
c
c                               b.drittler   may 1987
c-----------------------------------------------------------------------
      IMPLICIT NONE
C     .. Parameters ..

      INTEGER irmd
      INTEGER irid
      INTEGER nfund
      INTEGER ngshd
      INTEGER ipand

C
C     INTEGER LMPOTD
C     PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT,NSPIN
C     ..
C     .. Array Arguments ..

C     DOUBLE PRECISION DRDI(IRMD,*),
C    +                 GSH(*),R(IRMD,*),RHO2NS(IRMD,LMPOTD),
C    +                 THETAS(IRID,NFUND,*),VONS(IRMD,LMPOTD,2)
C     INTEGER IFUNM(*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
C    +        IRCUT(0:IPAND,*),IRWS(*),LMSP(*)

      DOUBLE PRECISION DRDI(IRMD)
      DOUBLE PRECISION GSH(*)
      DOUBLE PRECISION R(IRMD)
C     DOUBLE PRECISION RHO2NS(IRMD,LMPOTD)
      DOUBLE PRECISION RHO2NS(IRMD,(LPOT+1)**2)

      DOUBLE PRECISION THETAS(IRID,NFUND)
C     DOUBLE PRECISION VONS(IRMD,LMPOTD,2)
      DOUBLE PRECISION VONS(IRMD,(LPOT+1)**2,2)

      INTEGER IFUNM(*)
      INTEGER ILM(NGSHD,3)
C     INTEGER IMAXSH(0:LMPOTD)
      INTEGER IMAXSH(0:(LPOT+1)**2)
      INTEGER IPAN
      INTEGER IRCUT(0:IPAND)
      INTEGER LMSP(*)

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FAC,PI,RL
      INTEGER I,IEND,IFUN,IPOT,IRC1,IRS1,ISTART,J,L,LM,LM2,
     +        LM3,M
C     ..
C     .. Local Arrays ..
C     Fortran 90 automatic arrays
      DOUBLE PRECISION V1(IRMD),V2(IRMD),VINT1(IRMD),VINT2(IRMD)
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
      IRS1 = IRCUT(1)
      IRC1 = IRCUT(IPAN)
      DO 10 I = 0,IPAN
        IRCUTM(I) = IRCUT(I)
   10 CONTINUE

c---> determine the right potential numbers
      IPOT = NSPIN

      DO 100 L = 0,LPOT
        FAC = 8.0D0*PI/REAL(2*L+1)
        DO 90 M = -L,L
          LM = L*L + L + M + 1
c
c---> set up of the integrands v1 and v2
c
          V1(1) = 0.0D0
          V2(1) = 0.0D0
          DO 20 I = 2,IRS1
            RL = R(I)**L
            V1(I) = RHO2NS(I,LM)*RL*DRDI(I)
            V2(I) = RHO2NS(I,LM)/R(I)/RL*DRDI(I)
   20     CONTINUE
c
c---> convolute charge density of interstial with shape function
c
          DO 30 I = IRS1 + 1,IRC1
            V1(I) = 0.0D0
   30     CONTINUE
          ISTART = IMAXSH(LM-1) + 1
          IEND = IMAXSH(LM)
          DO 50 J = ISTART,IEND
            LM2 = ILM(J,2)
            LM3 = ILM(J,3)
            IF (LMSP(LM3).GT.0) THEN
              IFUN = IFUNM(LM3)
              DO 40 I = IRS1 + 1,IRC1
                V1(I) = V1(I) + GSH(J)*RHO2NS(I,LM2)*
     +                  THETAS(I-IRS1,IFUN)
   40         CONTINUE
            END IF
   50     CONTINUE

          DO 60 I = IRS1 + 1,IRC1
            RL = R(I)**L
            V2(I) = V1(I)/R(I)/RL*DRDI(I)
            V1(I) = V1(I)*RL*DRDI(I)
   60     CONTINUE
c
c---> now integrate v1 and v2
c
          CALL SOUTK(V1,VINT1,IPAN,IRCUTM)
          CALL SINWK(V2,VINT2,IPAN,IRCUTM)
c
c---> gather all parts
c
          IF (LM.EQ.1) THEN
            VONS(1,LM,IPOT) = FAC*VINT2(1)

          ELSE

            VONS(1,LM,IPOT) = 0.0D0
          END IF

          DO 70 I = 2,IRC1
            RL = R(I)**L
            VONS(I,LM,IPOT) = FAC* (VINT1(I)/R(I)/RL+VINT2(I)*RL)
   70     CONTINUE
c
          IF (NSPIN.EQ.2) THEN
            DO 80 I = 1,IRC1
              VONS(I,LM,IPOT-1) = VONS(I,LM,IPOT)
   80       CONTINUE
          END IF

   90   CONTINUE

  100 CONTINUE

      RETURN

      END
