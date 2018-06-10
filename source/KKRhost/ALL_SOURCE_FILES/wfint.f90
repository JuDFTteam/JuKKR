SUBROUTINE wfint(qns,cder,dder,qzekdr,pzekdr,vnspll,nsra,irmind,  &
        irmd,lmmaxd,irmin,irmax)                          ! Added IRMIN,IRMAX 1.7.2014
!-----------------------------------------------------------------------
!      determines the integrands CDER, DDER or ADER, BDER in the
!        integral equations for the non-spherical wavefunctions from
!        the non-spherical contributions of the potential vinsPLL.

!      R. Zeller      Aug. 1994
!-----------------------------------------------------------------------
IMPLICIT NONE
!.. Scalar Arguments ..
      INTEGER IRMD,IRMIND,LMMAXD,NSRA,IRMIN,IRMAX
!..
!.. Array Arguments ..
      DOUBLE COMPLEX CDER(LMMAXD,LMMAXD,IRMIND:IRMD), &
                     DDER(LMMAXD,LMMAXD,IRMIND:IRMD), &
                     PZEKDR(LMMAXD,IRMIND:IRMD,2), &
                     QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2), &
                     QZEKDR(LMMAXD,IRMIND:IRMD,2)
      DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
!..
!.. Local Scalars ..
      INTEGER IR,LM1,LM2
!..
!.. External Subroutines ..
      EXTERNAL DGEMM
!..
!.. Local Arrays ..
      DOUBLE PRECISION QNSI(LMMAXD,LMMAXD),QNSR(LMMAXD,LMMAXD), &
                       VTQNSI(LMMAXD,LMMAXD),VTQNSR(LMMAXD,LMMAXD)
!..
!.. Intrinsic Functions ..
      INTRINSIC DCMPLX,DIMAG,DBLE
!     ..

DO  ir = irmin,irmax
  DO  lm2 = 1,lmmaxd
    DO  lm1 = 1,lmmaxd
      qnsr(lm1,lm2) = DBLE(qns(lm1,lm2,ir,1))
      qnsi(lm1,lm2) = DIMAG(qns(lm1,lm2,ir,1))
    END DO
  END DO
  CALL dgemm('N','N',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),  &
      lmmaxd,qnsr,lmmaxd,0.d0,vtqnsr,lmmaxd)
  CALL dgemm('N','N',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),  &
      lmmaxd,qnsi,lmmaxd,0.d0,vtqnsi,lmmaxd)
  DO  lm1 = 1,lmmaxd
    DO  lm2 = 1,lmmaxd
      cder(lm1,lm2,ir) = qzekdr(lm1,ir,1)*  &
          DCMPLX(vtqnsr(lm1,lm2),vtqnsi(lm1,lm2))
      dder(lm1,lm2,ir) = pzekdr(lm1,ir,1)*  &
          DCMPLX(vtqnsr(lm1,lm2),vtqnsi(lm1,lm2))
    END DO
  END DO
  IF (nsra == 2) THEN
    DO  lm2 = 1,lmmaxd
      DO  lm1 = 1,lmmaxd
        qnsr(lm1,lm2) = DBLE(qns(lm1,lm2,ir,2))
        qnsi(lm1,lm2) = DIMAG(qns(lm1,lm2,ir,2))
      END DO
    END DO
    CALL dgemm('N','N',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),  &
        lmmaxd,qnsr,lmmaxd,0.d0,vtqnsr,lmmaxd)
    CALL dgemm('N','N',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),  &
        lmmaxd,qnsi,lmmaxd,0.d0,vtqnsi,lmmaxd)
    DO  lm2 = 1,lmmaxd
      DO  lm1 = 1,lmmaxd
        cder(lm1,lm2,ir) = cder(lm1,lm2,ir) +  &
            qzekdr(lm1,ir,2)*DCMPLX(vtqnsr(lm1, lm2),vtqnsi(lm1,lm2))
        dder(lm1,lm2,ir) = dder(lm1,lm2,ir) +  &
            pzekdr(lm1,ir,2)*DCMPLX(vtqnsr(lm1, lm2),vtqnsi(lm1,lm2))
      END DO
    END DO
  END IF
  
END DO
END SUBROUTINE wfint
