SUBROUTINE pnstmat(drdi,ek,icst,pz,qz,fz,sz,pns,tmatll,vins,irmin,  &
        ipan,ircut,nsra,cleb,icleb,iend,loflm,tmat,lkonv, & ! Added IRMIN 1.7.2014  &
        idoldau,lopt,lmlo,lmhi,wldau,wldauav,cutoff,  &
        alpha0)  ! LLY
IMPLICIT NONE
!     .. Parameters ..
INCLUDE 'inc.p'
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *  LDA+U implementation     Mar. 2002-Dec.2004                      *
! *                           ph.mavropoulos, h. ebert, v. popescu    *
! *                                                                   *
! *********************************************************************
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER MMAXD
      PARAMETER (MMAXD=2*LMAXD+1)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+1) * (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO=(0.D0,0.D0))
!..
!.. Scalar Arguments ..
      DOUBLE COMPLEX EK
      INTEGER ICST,IDOLDAU,IEND,IPAN,LKONV,LOPT,NSRA,LMLO,LMHI,IRMIN
      DOUBLE PRECISION WLDAUAV
!..
!.. Array Arguments ..
DOUBLE COMPLEX FZ(IRMD,0:LMAXD),PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2), &
               PZ(IRMD,0:LMAXD),QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD), &
               TMAT(0:LMAXD),TMATLL(LMMAXD,LMMAXD), &
               ALPHA0(LMMAXD,LMMAXD)  ! LLY
DOUBLE PRECISION CLEB(NCLEB,2),DRDI(IRMD),VINS(IRMIND:IRMD,LMPOTD)
DOUBLE PRECISION WLDAU(MMAXD,MMAXD),CUTOFF(IRMD)
INTEGER ICLEB(NCLEB,4),IRCUT(0:IPAND),LOFLM(*)
!..
!.. Local Scalars ..
      INTEGER I,IR,LM1,LM2,LMMKONV,M1,M2,IRMAX
!..
!.. Local Arrays ..
DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CMAT(LMMAXD,LMMAXD,IRMIND:IRMD), &
               DMAT(LMMAXD,LMMAXD,IRMIND:IRMD),EFAC(LMMAXD), &
               PZEKDR(LMMAXD,IRMIND:IRMD,2), &
               PZLM(LMMAXD,IRMIND:IRMD,2), &
               QZEKDR(LMMAXD,IRMIND:IRMD,2), &
               QZLM(LMMAXD,IRMIND:IRMD,2)
DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
!..
!.. External Subroutines ..
      EXTERNAL REGNS,VLLNS,WFTSCA,ZGEMM

irmax = ircut(ipan)

CALL vllns(vnspll,vins,cleb,icleb,iend, irmd,ncleb,lmpotd,irmind,lmmaxd)
IF (lkonv /= lmaxd) THEN
  lmmkonv = (lkonv+1)* (lkonv+1)
  DO lm1 = 1,lmmaxd
    DO lm2 = lmmkonv + 1,lmmaxd
      DO i = irmind,irmd
        vnspll(lm2,lm1,i) = 0.0D0
        vnspll(lm1,lm2,i) = 0.0D0
      END DO
    END DO
  END DO
ELSE
  lmmkonv = lmmaxd
END IF
! ======================================================================
! LDA+U
! Add WLDAU to non-spherical porential VINS in case of LDA+U
! Use the average wldau (=wldauav) and calculate the deviation
! of wldau from this. Use the deviation in the Born series
! for the non-spherical wavefunction, while the average is
! used for the spherical wavefunction.

IF ( idoldau == 1.AND.lopt >= 0 ) THEN
  DO ir = irmind,irmd
    
! -> First add wldau to all elements ...
    
    DO lm2 = lmlo,lmhi
      m2 = lm2 - lmlo + 1
      DO lm1 = lmlo,lmhi
        m1 = lm1 - lmlo + 1
        vnspll(lm1,lm2,ir) =  vnspll(lm1,lm2,ir) + wldau(m1,m2) * cutoff(ir)
      END DO
      
! ... and then subtract average from diag. elements
      
      vnspll(lm2,lm2,ir) =  vnspll(lm2,lm2,ir) - wldauav * cutoff(ir)
    END DO
  END DO
END IF

! LDA+U
! ======================================================================
pzlm(:,irmind:irmd,:) = czero
qzlm(:,irmind:irmd,:) = czero
pzekdr(:,irmind:irmd,:) = czero
qzekdr(:,irmind:irmd,:) = czero
cmat(:,:,irmind:irmd) = czero
dmat(:,:,irmind:irmd) = czero

!---> get wfts of same magnitude by scaling with efac

CALL wftsca(drdi,efac,pz,qz,fz,sz,nsra,pzlm,qzlm,pzekdr,qzekdr,  &
    ek,loflm,irmind,irmd,irmin,irmax,lmaxd,lmmaxd)      ! Added IRMIN,IRMAX 1.7.2014

!---> determine the regular non sph. wavefunction

CALL regns(ar,tmatll,efac,pns,vnspll,icst,ipan,ircut,pzlm,qzlm,  &
    pzekdr,qzekdr,ek,pns(1,1,irmind,1),cmat,  &
    pns(1,1,irmind,2),dmat,nsra,irmind,irmd,irmin,irmax, & ! Added IRMIN,IRMAX 1.7.2014  &
    ipand,lmmaxd)


DO lm1 = 1,lmmkonv
  tmatll(lm1,lm1) = tmatll(lm1,lm1) + tmat(loflm(lm1))
END DO


DO lm2 = 1,lmmkonv
  DO lm1 = 1,lmmkonv
    ar(lm1,lm2) = alpha0(lm1,lm1) * ar(lm1,lm2) ! LLY non-spher. contribution to alpha matrix
  END DO                                          ! LLY Drittler PhD eq. 3.106
END DO
alpha0(1:lmmaxd,1:lmmaxd) = ar(1:lmmaxd,1:lmmaxd) ! on output.


RETURN
END SUBROUTINE pnstmat
