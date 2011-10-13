      SUBROUTINE CYLM02(LMAX,COSX,FAI,LPOT2P,LMMAXD,THET,YLM,DYLMT1,
     +                  DYLMT2,DYLMF1,DYLMF2,DYLMTF)
c.....------------------------------------------------------------------
c     preparation of cylm0(=ylm(ip,i)), cylmt1(=dylm/dtheta),
c     cylmt2(=d2ylm/dt2),
c     cylmf1, cylmf2 are for fai.
c     cylmtf=d2ylm/dfdt
c     i=1,2,....,(lmax+1)**2
c.....------------------------------------------------------------------
c
C     .. Parameters ..
      INTEGER IJD
      PARAMETER (IJD=434)
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAX,LMMAXD,LPOT2P
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION COSX(IJD),DYLMF1(IJD,LMMAXD),DYLMF2(IJD,LMMAXD),
     +                 DYLMT1(IJD,LMMAXD),DYLMT2(IJD,LMMAXD),
     +                 DYLMTF(IJD,LMMAXD),FAI(IJD),THET(IJD),
     +                 YLM(IJD,LMMAXD)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX CI,EM1F,EM2F,EP1F,EP2F
      DOUBLE PRECISION AAA,CCC,DI,FI,ONE,PI,SSS
      INTEGER I,IP,L,LLMAX,LM,LM1,LM1M,LM2,LMM,LMM1,LMM1M,LMM2,M,MM
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CYLM0(LMMAXD),CYLMF1(LMMAXD),CYLMF2(LMMAXD),
     +               CYLMT1(LMMAXD),CYLMT2(LMMAXD),CYLMTF(LMMAXD)
      DOUBLE PRECISION BB1(LMMAXD),YL(LPOT2P)
C     ..
C     .. External Subroutines ..
      EXTERNAL SPHER,TRAREA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ACOS,ATAN,CMPLX,CONJG,COS,DBLE,SIN,SQRT
C     ..
      CI = CMPLX(0.d0,1.d0)
      ONE = 1.d0
      PI = 4.d0*ATAN(ONE)
      LLMAX = (LMAX+1)**2

      DO 120 IP = 1,IJD

        THET(IP) = ACOS(COSX(IP))
        FI = FAI(IP)
        DI = 2*FAI(IP)
        EP1F = CMPLX(COS(FI),SIN(FI))
        EM1F = CONJG(EP1F)
        EP2F = CMPLX(COS(DI),SIN(DI))
        EM2F = CONJG(EP2F)

        DO 50 L = 0,LMAX
c
          CALL SPHER(YL,L,COSX(IP))
          DO 10 M = -L,L
            MM = L + M + 1
            I = (L+1)**2 - L + M
            AAA = M*FAI(IP)
            CCC = COS(AAA)
            SSS = SIN(AAA)
            CYLM0(I) = YL(MM)*CMPLX(CCC,SSS)
   10     CONTINUE

          DO 20 M = -L,L
            I = (L+1)**2 - L + M
            CYLMT1(I) = 0.D0
            CYLMT2(I) = 0.D0
            CYLMTF(I) = 0.D0
   20     CONTINUE

          DO 30 M = -L,L
            I = (L+1)**2 - L + M

            LMM1M = L - M - 1
            LMM = L - M
            LMM1 = L - M + 1
            LMM2 = L - M + 2
            LM1M = L + M - 1
            LM = L + M
            LM1 = L + M + 1
            LM2 = L + M + 2

            CYLMT2(I) = CYLMT2(I) - (LMM*LM1+LMM1*LM)/4.D0*CYLM0(I)

            IF (M+2.LE.L) CYLMT2(I) = CYLMT2(I) +
     +                                SQRT(DBLE(LMM1M*LMM*LM1*LM2))/4*
     +                                CYLM0(I+2)*EM2F

            IF (M+1.LE.L) CYLMT1(I) = CYLMT1(I) +
     +                                SQRT(DBLE(LMM*LM1))/2*CYLM0(I+1)*
     +                                EM1F

            IF (M-1.GE.-L) CYLMT1(I) = CYLMT1(I) -
     +                                 SQRT(DBLE(LM*LMM1))/2*CYLM0(I-1)*
     +                                 EP1F

            IF (M-2.GE.-L) CYLMT2(I) = CYLMT2(I) +
     +                                 SQRT(DBLE(LMM1*LMM2*LM1M*LM))/4*
     +                                 CYLM0(I-2)*EP2F

   30     CONTINUE

          DO 40 M = -L,L
            I = (L+1)**2 - L + M
            CYLMF1(I) = CI*M*CYLM0(I)
            CYLMF2(I) = -M*M*CYLM0(I)
            CYLMTF(I) = CI*M*CYLMT1(I)
   40     CONTINUE

   50   CONTINUE
c
c        calculate real spherical harmonics differenciated
c
c
c        write(6,9005) (cylm0(i),i=1,5)
c9005 format(1x,' cylm0',4f10.5)
        CALL TRAREA(CYLM0,BB1,LMAX)

        DO 60 M = 1,LLMAX
          YLM(IP,M) = BB1(M)
   60   CONTINUE
c
c        write(6,9006) (ylm(ip,i),i=1,5)
c9006 format(1x,' ylm',10f10.5)
c
c
        CALL TRAREA(CYLMT1,BB1,LMAX)
        DO 70 M = 1,LLMAX
          DYLMT1(IP,M) = BB1(M)
   70   CONTINUE
C
        CALL TRAREA(CYLMT2,BB1,LMAX)
        DO 80 M = 1,LLMAX
          DYLMT2(IP,M) = BB1(M)
   80   CONTINUE
c
        CALL TRAREA(CYLMF1,BB1,LMAX)
        DO 90 M = 1,LLMAX
          DYLMF1(IP,M) = BB1(M)
   90   CONTINUE
C
        CALL TRAREA(CYLMF2,BB1,LMAX)
        DO 100 M = 1,LLMAX
          DYLMF2(IP,M) = BB1(M)
  100   CONTINUE
C
        CALL TRAREA(CYLMTF,BB1,LMAX)
        DO 110 M = 1,LLMAX
          DYLMTF(IP,M) = BB1(M)
  110   CONTINUE
C
  120 CONTINUE
      RETURN
      END
