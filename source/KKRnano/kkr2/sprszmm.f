      SUBROUTINE SPRSZMM(
     >                   IAT,GLLH,NUMN0,INDN0,X,DONE,OMEGA,DELTA,
     <                   AX)
C
      IMPLICIT NONE
C     .. Parameters ..
      include 'inc.p'
      include 'inc.cls'
C
      INTEGER           LMMAXD
      PARAMETER        (LMMAXD= (LMAXD+1)**2)
      INTEGER           NDIM,NAEZ
      PARAMETER        (NAEZ=NAEZD,NDIM=NAEZD*LMMAXD)
      INTEGER           NGTBD
      PARAMETER        (NGTBD = NACLSD*LMMAXD)
C
      DOUBLE COMPLEX    CONE,CZERO
      PARAMETER        (CONE=(1.0D0,0.0D0),CZERO=(0.0D0,0.0D0))
C     ..
C     ... Global Scalars ..
      INTEGER          IAT
      DOUBLE COMPLEX   OMEGA,DELTA  ! scalars in Matrix-Matrix-Mult.
C     ... Global Arrays ..
      DOUBLE COMPLEX   X(NDIM,LMMAXD),AX(NDIM,LMMAXD),
     +                  GLLH(LMMAXD,NGTBD,NAEZD)
      INTEGER           NUMN0(NAEZD),
     +                  INDN0(NAEZD,NACLSD)
      LOGICAL           DONE(LMMAXD)
C .. 
C     .. External ..
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX   SPRSX(NGTBD,LMMAXD)
C ..
C     .. Local Scalars ..
      INTEGER          I1,I2,I3,LM2,IL1B,I2H,I3H
C#ifndef FPP_OMP
      INTEGER          MYTHRD
!$    INTEGER          OMP_GET_THREAD_NUM
C#endif
C ..
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C#ifndef FPP_OMP
!$    CALL OMP_SET_NUM_THREADS(NTHRDS)
!$OMP PARALLEL PRIVATE (I1,LM2,I2,I3H,I2H,
!$OMP&                  IL1B,SPRSX,MYTHRD)
      MYTHRD = 0
!$    MYTHRD = OMP_GET_THREAD_NUM()
C#endif
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      DO I1=1,NAEZ
C#ifndef FPP_OMP
      IF (MOD(I1,NTHRDS).EQ.MYTHRD) THEN
C#endif
C
        DO LM2=1,LMMAXD
        IF (.NOT.DONE(LM2)) THEN
          DO I2=1,NUMN0(IAT)
            I3=INDN0(I1,I2)
            I3H = (I3-1)*LMMAXD + 1
            I2H = (I2-1)*LMMAXD + 1
            CALL ZCOPY(LMMAXD,X(I3H,LM2),1,
     +                 SPRSX(I2H,LM2),1)
          ENDDO
        ENDIF
        ENDDO
C
        IL1B=LMMAXD*(I1-1)
C
        CALL ZGEMM('N','N',LMMAXD,LMMAXD,NUMN0(IAT)*LMMAXD,
     +             OMEGA,GLLH(1,1,I1),LMMAXD,
     +             SPRSX,NGTBD,
     +             DELTA,AX(IL1B+1,1),NDIM)
C
C#ifndef FPP_OMP
      ENDIF
C#endif
      ENDDO
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C add OMP-parallelization here
C end
C#ifndef FPP_OMP
!$OMP END PARALLEL
C#endif
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      END
C
C
