C ************************************************************************
      SUBROUTINE VLLMAT(IRMIN,NRMAXD,IRC,LMMAX,LMMAXSO,VNSPLL0,VINS,
     +                  CLEB,ICLEB,IEND,NSPIN,Z,RNEW,USE_SRATRICK)
C ************************************************************************
C     .. Parameters ..
      IMPLICIT NONE
      include 'inc.p'
      INTEGER LMPOTD,USE_SRATRICK
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION Z
      INTEGER IEND,IRC,IRMIN,LMMAX,ISP,LMMAXSO,NSPIN,NRMAXD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CLEB(*),VINS(IRMIN:IRC,LMPOTD,NSPIND),
     +                 RNEW(IRMIN:NRMAXD),
     +                 VNSPLL(LMMAX,LMMAX,IRMIN:IRC,NSPIND)
      DOUBLE COMPLEX   VNSPLL0(LMMAXSO,LMMAXSO,IRMIN:IRC)
      INTEGER ICLEB(NCLEB,4)
C     ..
C     .. Local Scalars ..
      INTEGER I,IR,J,LM1,LM2,LM3
C     ..
      DO ISP=1,NSPIN
        DO 30 LM1 = 1,LMMAX
          DO 20 LM2 = 1,LM1
            DO 10 IR = IRMIN,IRC
              VNSPLL(LM1,LM2,IR,ISP) = 0.0D0
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE

        DO 50 J = 1,IEND
          LM1 = ICLEB(J,1)
          LM2 = ICLEB(J,2)
          LM3 = ICLEB(J,3)
          DO 40 I = IRMIN,IRC
            VNSPLL(LM1,LM2,I,ISP) = VNSPLL(LM1,LM2,I,ISP) + 
     +                                 CLEB(J)*VINS(I,LM3,ISP)
   40     CONTINUE
   50   CONTINUE
c
c---> use symmetry of the gaunt coef.
c
        DO 80 LM1 = 1,LMMAX
          DO 70 LM2 = 1,LM1 - 1
            DO 60 I = IRMIN,IRC
              VNSPLL(LM2,LM1,I,ISP) = VNSPLL(LM1,LM2,I,ISP)
   60       CONTINUE
   70     CONTINUE
   80   CONTINUE
        
        IF (USE_SRATRICK.EQ.0) THEN
        DO LM1=1,LMMAX
         DO I=IRMIN,IRC
          VNSPLL(LM1,LM1,I,ISP)=VNSPLL(LM1,LM1,I,ISP)+
     +       VINS(I,1,ISP)-2d0*Z/RNEW(I)
         ENDDO
        ENDDO
        ENDIF

      END DO !NSPIN
      
c set vnspll as twice as large

       VNSPLL0(1:LMMAX,1:LMMAX,IRMIN:IRC)=
     +    VNSPLL(1:LMMAX,1:LMMAX,IRMIN:IRC,1)

       VNSPLL0(LMMAX+1:LMMAXSO,LMMAX+1:LMMAXSO,IRMIN:IRC)=
     +    VNSPLL(1:LMMAX,1:LMMAX,IRMIN:IRC,NSPIN)
      END