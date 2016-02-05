C ************************************************************************
      SUBROUTINE VLLMAT(IRMIN,IRC,LMMAX,LMMAXSO,VNSPLL0,VINS,CLEB,
     +                  ICLEB,IEND,NSPIN,Z,RNEW,USE_SRATRICK,I1)
C ************************************************************************
C     .. Parameters ..
      IMPLICIT NONE
      include 'inc.p'
      INTEGER LMPOTD,USE_SRATRICK
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION Z
      INTEGER IEND,IRC,IRMIN,LMMAX,ISP,LMMAXSO,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CLEB(*),VINS(IRMIN:IRC,LMPOTD,NSPIND),
     +                 RNEW(IRMIN:IRC),
     +                 VNSPLL(LMMAX,LMMAX,IRMIN:IRC,NSPIND)
      DOUBLE COMPLEX   VNSPLL0(LMMAXSO,LMMAXSO,IRMIN:IRC)
      INTEGER ICLEB(NCLEB,4)
C     ..
C     .. Local Scalars ..
      INTEGER I,IR,J,LM1,LM2,LM3,I1
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

cc shift d orbital for Ag
c       IF (I1.GE.4.AND.I1.LE.13) THEN
c        DO LM1=1,LMMAX
c         DO LM2=1,LMMAX
c          IF (LM1.EQ.LM2) THEN
c           IF (LM2.EQ.5.OR.LM2.EQ.6.OR.LM2.EQ.7
c     +                 .OR.LM2.EQ.8.OR.LM2.EQ.9) THEN
c            DO I=IRMIN,IRC
c             VNSPLL(LM2,LM1,I,ISP)=VNSPLL(LM2,LM1,I,ISP)+0.18374517
c            ENDDO
c           ENDIF
c          ENDIF
c         ENDDO
c        ENDDO
c       ENDIF

      END DO !NSPIN
      
c set vnspll as twice as large

       VNSPLL0(1:LMMAX,1:LMMAX,IRMIN:IRC)=
     +    VNSPLL(1:LMMAX,1:LMMAX,IRMIN:IRC,1)

       VNSPLL0(LMMAX+1:LMMAXSO,LMMAX+1:LMMAXSO,IRMIN:IRC)=
     +    VNSPLL(1:LMMAX,1:LMMAX,IRMIN:IRC,NSPIN)
      END
