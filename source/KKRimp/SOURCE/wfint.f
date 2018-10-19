!------------------------------------------------------------------------------------
!> Summary: Determines the integrands `CDER`, `DDER` or `ADER`, `BDER`
!> Author: R. Zeller
!> Determines the integrands `CDER`, `DDER` or `ADER`, `BDER` in the integral equations 
!> for the non-spherical wavefunctions from the non-spherical contributions of the 
!> potential `vinsPLL`.
!------------------------------------------------------------------------------------
      MODULE mod_wfint
      CONTAINS
  !-------------------------------------------------------------------------------
  !> Summary: Determines the integrands `CDER`, `DDER` or `ADER`, `BDER`
  !> Author: R. Zeller
  !> Category: numerical-tools, single-site, KKRimp
  !> Deprecated: False 
  !> determines the integrands `CDER`, `DDER` or `ADER`, `BDER` in the integral equations 
  !> for the non-spherical wavefunctions from the non-spherical contributions of the 
  !> potential `vinsPLL`. 
  !-------------------------------------------------------------------------------      
      SUBROUTINE WFINT(QNS,CDER,DDER,QZEKDR,PZEKDR,VNSPLL,NSRA,IRMIND,
     +                   IRMD,LMMAXD)
c     Implicit None
C     .. Scalar Arguments ..
      INTEGER IRMD,IRMIND,LMMAXD,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               DDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2)
      DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. Local Scalars ..
      INTEGER IR,LM1,LM2
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION QNSI(LMMAXD,LMMAXD),QNSR(LMMAXD,LMMAXD),
     +                 VTQNSI(LMMAXD,LMMAXD),VTQNSR(LMMAXD,LMMAXD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCMPLX,DIMAG,DBLE
C     ..
c
      DO 90 IR = IRMIND,IRMD
        DO 20 LM2 = 1,LMMAXD
          DO 10 LM1 = 1,LMMAXD
            QNSR(LM1,LM2) = DBLE(QNS(LM1,LM2,IR,1))
            QNSI(LM1,LM2) = DIMAG(QNS(LM1,LM2,IR,1))
   10     CONTINUE
   20   CONTINUE
        CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +             LMMAXD,QNSR,LMMAXD,0.D0,VTQNSR,LMMAXD)
        CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +             LMMAXD,QNSI,LMMAXD,0.D0,VTQNSI,LMMAXD)
        DO 40 LM1 = 1,LMMAXD
          DO 30 LM2 = 1,LMMAXD
            CDER(LM1,LM2,IR) = QZEKDR(LM1,IR,1)*
     +                         DCMPLX(VTQNSR(LM1,LM2),VTQNSI(LM1,LM2))
            DDER(LM1,LM2,IR) = PZEKDR(LM1,IR,1)*
     +                         DCMPLX(VTQNSR(LM1,LM2),VTQNSI(LM1,LM2))
   30     CONTINUE
   40   CONTINUE
        IF (NSRA.EQ.2) THEN
          DO 60 LM2 = 1,LMMAXD
            DO 50 LM1 = 1,LMMAXD
              QNSR(LM1,LM2) = DBLE(QNS(LM1,LM2,IR,2))
              QNSI(LM1,LM2) = DIMAG(QNS(LM1,LM2,IR,2))
   50       CONTINUE
   60     CONTINUE
          CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +               LMMAXD,QNSR,LMMAXD,0.D0,VTQNSR,LMMAXD)
          CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +               LMMAXD,QNSI,LMMAXD,0.D0,VTQNSI,LMMAXD)
          DO 80 LM2 = 1,LMMAXD
            DO 70 LM1 = 1,LMMAXD
              CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) +
     +                           QZEKDR(LM1,IR,2)*DCMPLX(VTQNSR(LM1,
     +                           LM2),VTQNSI(LM1,LM2))
              DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) +
     +                           PZEKDR(LM1,IR,2)*DCMPLX(VTQNSR(LM1,
     +                           LM2),VTQNSI(LM1,LM2))
   70       CONTINUE
   80     CONTINUE
        END IF

   90 CONTINUE
      END SUBROUTINE
      END MODULE mod_wfint