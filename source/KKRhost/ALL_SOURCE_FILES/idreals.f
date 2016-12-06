      SUBROUTINE IDREALS(DARRY,NARRY,IPRINT)
      IMPLICIT NONE
C
C PARAMETER definitions
C
      INTEGER NSQR,NMUL,DIVMAX
      PARAMETER (NSQR=7,NMUL=5,DIVMAX=15)
      DOUBLE PRECISION TOL
      PARAMETER (TOL=1D-6)
C
C Dummy arguments
C
      INTEGER IPRINT,NARRY
      DOUBLE PRECISION DARRY(NARRY)
C
C Local variables
C
      DOUBLE PRECISION DABS,DBLE,DSQRT,DSIGN
      INTEGER DIV,I1,I2,IDONE(NARRY),IMUL(NMUL),ISQR(NSQR)
      DOUBLE PRECISION DSQ,X,XN
      INTEGER IABS,IDNINT
C
      DATA ISQR/2,3,5,6,7,8,10/
      DATA IMUL/3,7,11,13,17/
C
C --> mark all numbers as unchecked
C
      DO I1 = 1,NARRY
         IDONE(I1) = 0
      END DO
C
C --> check darry**2/i integer?, i=1,divmax
C
      DO DIV = 1,DIVMAX
         DSQ = DBLE(DIV)
         DO I2 = 1,NARRY
            IF ( IDONE(I2).EQ.0 ) THEN
               X = DARRY(I2)*DARRY(I2)*DSQ
               XN = DNINT(X)
               IF ( DABS(X-XN)/DSQ.LT.TOL .AND. XN.NE.0.D0 ) THEN
                  IF (IPRINT.GT.4) WRITE (1337,99000) 
     &                 DABS(DARRY(I2)),NINT(X),DIV
                  DARRY(I2) = DSIGN(1D0,DARRY(I2))*DSQRT(XN/DSQ)
                  IDONE(I2) = 1
               END IF
            END IF
         END DO
      END DO
C
C --> check darry/sqrt(n) =?=  i/j 
C        n=2,3,5,6,7,8,10      i=1,divmax j=i*n
C
      DO I1 = 1,NSQR
         DO DIV = 1,DIVMAX
            DSQ = DSQRT(DBLE(DIV*DIV*ISQR(I1)))
            DO I2 = 1,NARRY
               IF ( IDONE(I2).EQ.0 ) THEN
                  X = DARRY(I2)*DSQ
                  XN = DNINT(X)
                  IF ( DABS(X-XN)/DSQ.LT.TOL .AND. XN.NE.0.D0 ) THEN
                     IF (IPRINT.GT.4) WRITE (1337,99001) 
     &                    DABS(DARRY(I2)),ISQR(I1),
     &                    IABS(IDNINT(XN)),IABS(ISQR(I1)*DIV)
                     DARRY(I2) = XN/DSQ
                     IDONE(I2) = 1
                  END IF
               END IF
            END DO
         END DO
      END DO
C
C --> check darry = j/i * n ?
C        n=3,7,11,13,17
C                  
      DO I1 = 1,NMUL
         DO DIV = 1,DIVMAX
            DSQ = DBLE(DIV*IMUL(I1))
            DO I2 = 1,NARRY
               IF ( IDONE(I2).EQ.0 ) THEN
                  X = DARRY(I2)*DSQ
                  XN = DNINT(X)
                  IF ( DABS(X-XN)/DSQ.LT.TOL .AND. XN.NE.0.D0 ) THEN
                     IF (IPRINT.GT.4) WRITE(1337,99002) 
     &                    DABS(DARRY(I2)),IMUL(I1),IABS(IDNINT(XN)),DIV
                     DARRY(I2) = XN/DSQ
                     IDONE(I2) = 1
                  END IF
               END IF
            END DO
         END DO
      END DO
      RETURN
C
99000 FORMAT (8X,'< IDREALS > : identify ',f12.8,
     &     ' as dsqrt(',i3,'/',i3,')')
99001 FORMAT (8X,'< IDREALS > : identify ',f12.8,
     &     ' as dsqrt(',i2,')*',i3,'/',i3)
99002 FORMAT (8X,'< IDREALS > : identify ',f12.8,' as 1/',i2,
     &     ' * ',i2,'/',i1)
C
      END
