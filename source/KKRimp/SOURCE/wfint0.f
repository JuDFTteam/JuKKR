!-------------------------------------------------------------------------------
!> Summary: Determines the integrands `CDER`, `DDER` or `ADER`, `BDER`
!> Author: R. Zeller
!> Determines the integrands `CDER`, `DDER` or `ADER`, `BDER` in the integral equations 
!> for the non-spherical wavefunctions from the non-spherical contributions of the 
!> potential `vinsPLL`. (This subroutine is used in zeroth order Born approximation,
!> otherwise subroutine WFINT must be used)
!-------------------------------------------------------------------------------
      MODULE mod_wfint0
      CONTAINS
  !-------------------------------------------------------------------------------
  !> Summary: Determines the integrands `CDER`, `DDER` or `ADER`, `BDER`
  !> Author: R. Zeller
  !> Category: numerical-tools, single-site, KKRimp
  !> Deprecated: False 
  !> Determines the integrands `CDER`, `DDER` or `ADER`, `BDER` in the integral equations 
  !> for the non-spherical wavefunctions from the non-spherical contributions of the 
  !> potential `vinsPLL`. (This subroutine is used in zeroth order Born approximation,
  !> otherwise subroutine WFINT must be used)
  !-------------------------------------------------------------------------------      
      SUBROUTINE WFINT0(CDER,DDER,QZLM,QZEKDR,PZEKDR,VNSPLL,NSRA,
     +                    IRMIND,IRMD,LMMAXD)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c      determines the integrands CDER, DDER or ADER, BDER in the
c        integral equations for the non-spherical wavefunctions from
c        the non-spherical contributions of the potential vinsPLL.
c        (This subroutine is used in zeroth order Born approximation,
c         otherwise subroutine WFINT must be used)
c      R. Zeller      Aug. 1994
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IRMD,IRMIND,LMMAXD,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               DDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZLM(LMMAXD,IRMIND:IRMD,2)
      DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX V1
      INTEGER IR,LM1,LM2
C     ..
c
      DO 50 IR = IRMIND,IRMD
        DO 20 LM2 = 1,LMMAXD
          DO 10 LM1 = 1,LMMAXD
            V1 = VNSPLL(LM1,LM2,IR)*QZLM(LM2,IR,1)
            CDER(LM1,LM2,IR) = QZEKDR(LM1,IR,1)*V1
            DDER(LM1,LM2,IR) = PZEKDR(LM1,IR,1)*V1
   10     CONTINUE
   20   CONTINUE
        IF (NSRA.EQ.2) THEN
          DO 40 LM2 = 1,LMMAXD
            DO 30 LM1 = 1,LMMAXD
              V1 = VNSPLL(LM1,LM2,IR)*QZLM(LM2,IR,2)
              CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) + QZEKDR(LM1,IR,2)*V1
              DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) + PZEKDR(LM1,IR,2)*V1
   30       CONTINUE
   40     CONTINUE
        END IF

   50 CONTINUE
      END SUBROUTINE
      END MODULE mod_wfint0