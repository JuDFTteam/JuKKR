SUBROUTINE wfint0(cder,dder,qzlm,qzekdr,pzekdr,vnspll,nsra,  &
        irmind,irmd,lmmaxd,irmin,irmax)                 ! Added IRMIN,IRMAX 1.7.2014
!-----------------------------------------------------------------------
!      determines the integrands CDER, DDER or ADER, BDER in the
!        integral equations for the non-spherical wavefunctions from
!        the non-spherical contributions of the potential vinsPLL.
!        (This subroutine is used in zeroth order Born approximation,
!         otherwise subroutine WFINT must be used)
!      R. Zeller      Aug. 1994
!-----------------------------------------------------------------------
implicit none
!.. Scalar Arguments ..
      INTEGER IRMD,IRMIND,LMMAXD,NSRA,IRMIN,IRMAX
!..
!.. Array Arguments ..
DOUBLE COMPLEX CDER(LMMAXD,LMMAXD,IRMIND:IRMD), &
               DDER(LMMAXD,LMMAXD,IRMIND:IRMD), &
               PZEKDR(LMMAXD,IRMIND:IRMD,2), &
               QZEKDR(LMMAXD,IRMIND:IRMD,2), &
               QZLM(LMMAXD,IRMIND:IRMD,2)
DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
!..
!.. Local Scalars ..
      DOUBLE COMPLEX V1
      INTEGER IR,LM1,LM2

DO  ir = irmin,irmax
  DO  lm2 = 1,lmmaxd
    DO  lm1 = 1,lmmaxd
      v1 = vnspll(lm1,lm2,ir)*qzlm(lm2,ir,1)
      cder(lm1,lm2,ir) = qzekdr(lm1,ir,1)*v1
      dder(lm1,lm2,ir) = pzekdr(lm1,ir,1)*v1
    END DO
  END DO
  IF (nsra == 2) THEN
    DO  lm2 = 1,lmmaxd
      DO  lm1 = 1,lmmaxd
        v1 = vnspll(lm1,lm2,ir)*qzlm(lm2,ir,2)
        cder(lm1,lm2,ir) = cder(lm1,lm2,ir) + qzekdr(lm1,ir,2)*v1
        dder(lm1,lm2,ir) = dder(lm1,lm2,ir) + pzekdr(lm1,ir,2)*v1
      END DO
    END DO
  END IF
  
END DO
END SUBROUTINE wfint0
