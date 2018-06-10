SUBROUTINE wftsca(drdi,efac,pz,qz,fz,sz,nsra,pzlm,qzlm,pzekdr,  &
        qzekdr,ek,loflm,irmind,irmd,irmin,irmax, &     ! Added IRMIN,IRMAX 1.7.2014  &
        lmaxd,lmmaxd)
!-----------------------------------------------------------------------
!                 R. Zeller      Oct. 1993
!-----------------------------------------------------------------------
implicit none
!.. Parameters ..
DOUBLE COMPLEX CONE
PARAMETER (CONE= (1.D0,0.D0))
!..
!.. Scalar Arguments ..
DOUBLE COMPLEX EK
INTEGER IRMD,IRMIND,LMAXD,LMMAXD,NSRA,IRMIN,IRMAX
!..
!.. Array Arguments ..
DOUBLE COMPLEX EFAC(LMMAXD),FZ(IRMD,0:LMAXD),PZ(IRMD,0:LMAXD), &
               PZEKDR(LMMAXD,IRMIND:IRMD,2), &
               PZLM(LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD), &
               QZEKDR(LMMAXD,IRMIND:IRMD,2), &
               QZLM(LMMAXD,IRMIND:IRMD,2),SZ(IRMD,0:LMAXD)
DOUBLE PRECISION DRDI(*)
INTEGER LOFLM(*)
!..
!.. Local Scalars ..
DOUBLE COMPLEX EFAC1,V1
INTEGER IR,J,L,L1,LM,LM1,M
!..
!.. Intrinsic Functions ..
INTRINSIC DBLE
!..


!---> set up array efac : efac(lm) = sqrt(e)**l/(2l - 1)!!

efac(1) = cone
v1 = cone
DO  l = 1,lmaxd
  v1 = v1*ek/DBLE(2*l-1)
  DO  m = -l,l
    lm = l* (l+1) + m + 1
    efac(lm) = v1
  END DO
END DO


!---> get wfts of same magnitude by scaling with efac

DO  lm1 = 1,lmmaxd
  l1 = loflm(lm1)
  efac1 = efac(lm1)
  DO  ir = irmin,irmax
    pzlm(lm1,ir,1) = pz(ir,l1)/efac1
    qzlm(lm1,ir,1) = qz(ir,l1)*efac1
  END DO
  IF (nsra == 2) THEN
    DO  ir = irmin,irmax
      pzlm(lm1,ir,nsra) = fz(ir,l1)/efac1
      qzlm(lm1,ir,nsra) = sz(ir,l1)*efac1
    END DO
  END IF
  
  DO  j = 1,nsra
    DO  ir = irmin,irmax
      pzekdr(lm1,ir,j) = pzlm(lm1,ir,j)*ek*drdi(ir)
      qzekdr(lm1,ir,j) = qzlm(lm1,ir,j)*ek*drdi(ir)
    END DO
  END DO
END DO


END SUBROUTINE wftsca
