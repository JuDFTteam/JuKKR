      SUBROUTINE CINTHFF(AG,AF,BG,BF,RMEHF,NKA,NKB,JTOP,FX,R,DRDI,NRMAX)
C   ********************************************************************
C   *                                                                  *
C   *  routine to calculate the radial hyperfine matrixelement         *
C   *  by extrapolating the lower integration boundary to r -> 0       *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER JTOP,NKA,NKB,NRMAX
      DOUBLE COMPLEX AF(NRMAX,NKA),AG(NRMAX,NKA),BF(NRMAX,NKB),
     &               BG(NRMAX,NKB),FX(JTOP),RMEHF(2,2)
      DOUBLE PRECISION DRDI(NRMAX),R(NRMAX)
C
C Local variables
C
      DOUBLE PRECISION AI,AR,BI,BR,DELTA,HIBAR,R1BAR,RELDIF,
     &                 RLIM1,RLIM2,S,SI,SR,SX,SXI,SXR,SXX,YI,YR
      
      DOUBLE COMPLEX HI,HIF,Z(NRMAX)
      INTEGER I,IMAX,IMIN,KA,KB,N
C
      DO KB = 1,NKB
         DO KA = 1,NKA
C
            DO I = 1,JTOP
               FX(I) = (AG(I,KA)*BF(I,KB)+AF(I,KA)*BG(I,KB))*DRDI(I)
            END DO
C
            CALL CINT4PTS(FX,JTOP,Z)
C
            RLIM1 = 0.5D-5
            RLIM2 = 0.5D-4
C            RLIM1 = 1D-4
C            RLIM2 = 5D-4
            IF ( R(1).GT.RLIM1 ) THEN
               RLIM1 = R(1)
               RLIM2 = R(20)
            END IF
C
            IMIN = 0
            IMAX = 0
            DO I = 1,JTOP
               IF ( R(I).LE.RLIM1 ) IMIN = I
               IF ( R(I).GE.RLIM2 ) THEN
                  IMAX = I
                  GOTO 20
               END IF
            END DO
 20         CONTINUE
            IF ( IMIN.EQ.0 ) STOP '<CINTHFF> IMIN = 0'
            IF ( IMIN.GT.JTOP ) STOP '<CINTHFF> IMIN > JTOP'
            IF ( IMAX.EQ.0 ) STOP '<CINTHFF> IMAX = 0'
            IF ( IMAX.GT.JTOP ) STOP '<CINTHFF> IMAX > JTOP'
C
            N = 0
            S = 0.0D0
            SX = 0.0D0
            SXX = 0.0D0
            SR = 0.0D0
            SI = 0.0D0
            SXR = 0.0D0
            SXI = 0.0D0
C
            DO I = IMIN,IMAX
               YR = DREAL(Z(JTOP)-Z(I))
               YI = DIMAG(Z(JTOP)-Z(I))
               N = N + 1
               S = S + 1.0D0
               SX = SX + R(I)
               SXX = SXX + R(I)**2
               SR = SR + YR
               SI = SI + YI
               SXR = SXR + YR*R(I)
               SXI = SXI + YI*R(I)
            END DO
C
            DELTA = S*SXX - SX*SX
            AR = (SXX*SR-SX*SXR)/DELTA
            AI = (SXX*SI-SX*SXI)/DELTA
            BR = (S*SXR-SX*SR)/DELTA
            BI = (S*SXI-SX*SI)/DELTA
C
            HIBAR = 0.0D0
            R1BAR = 0.0D0
            DO I = IMIN,IMAX
               HIBAR = HIBAR + DBLE(Z(JTOP)-Z(I))
               R1BAR = R1BAR + R(I)
            END DO
            HIBAR = HIBAR/DBLE(IMAX-IMIN+1)
            R1BAR = R1BAR/DBLE(IMAX-IMIN+1)
C
            SXX = 0.0D0
            SXR = 0.0D0
            SXI = 0.0D0
            RELDIF = 0.0D0
            DO I = IMIN,IMAX
               HI = Z(JTOP) - Z(I)
               HIF = DCMPLX(AR,AI) + DCMPLX(BR,BI)*R(I)
               IF ( ABS(HIF).NE.0.0D0 )
     &              RELDIF = MAX(RELDIF,ABS(1.0D0-HI/HIF))
            END DO
C
            RMEHF(KA,KB) = DCMPLX(AR,AI)
C
         END DO
      END DO
C
      END
