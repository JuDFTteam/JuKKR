      SUBROUTINE RHONS(DEN,DF,DRDI,GMAT,EK,RHO2NS,IPAN,IRCUT,THETAS,
     +                   IFUNM,LMSP,NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB,
     +                   ICLEB,JEND,IEND,EKL,
     &                   lmax, irmd, irnsd, irid, ipand, nfund, ncleb)
c-----------------------------------------------------------------------
c
c     the charge density is developed in spherical harmonics :
c
c             rho(r) =   { rho(lm,r) * y(r,lm) }       (summed over lm)
c
c          rho(lm,r) =   { do rho(r) * y(r,lm)         (integrated over
c                                                           unit sphere)
c     in the case of spin-polarization :
c       the spin density is developed in spherical harmonics :
c
c            sden(r) =   { sden(lm,r) * y(r,lm) }      (summed over lm)
c
c         sden(lm,r) =   { do sden(r) * y(r,lm)        (integrated over
c                                                           unit sphere)
c     n(r,e) is developed in
c
c        n(r,e) = { y(r,l'm') * n(l'm',lm,r,e) * y(r,lm) }
c
c     therefore a faltung of n(l'm',lm,r,e) with the gaunt coeffients
c     has to be used to calculate the lm-contribution of the charge
c     density .
c
c
c     calculate the valence density of states , in the spin-polarized
c      case spin dependent .
c     recognize that the density of states is always complex also in
c      the case of "real-energy-integation" (ief>0) since in that case
c      the energy integration is done parallel to the real energy axis
c      but not on the real energy axis .
c     in the last energy-spin loop the l-contribution of the valence
c      charge is calculated .
c
c                               b.drittler   aug. 1988
c
c     modified for the use of shape functions
c
c     attention : irmin + 3 has to be less then imt
c                 if shape functions are used
c
c                               b.drittler   july 1989
c-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER lmax
      INTEGER irmd
      INTEGER ncleb
      INTEGER irnsd
      INTEGER irid
      INTEGER ipand
      INTEGER nfund

C     INTEGER IRMIND
C     PARAMETER (IRMIND=IRMD-IRNSD)
C     INTEGER LMPOTD,LMMAXD
C     PARAMETER (LMPOTD= (LPOTD+1)**2) ! LMPOTD= (2*LMAX+1)**2
C     PARAMETER (LMMAXD= (LMAXD+1)**2)
C     INTEGER LMAXD1
C     PARAMETER (LMAXD1= LMAXD+1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,EK
      INTEGER IEND,IPAN,NSRA
C     ..
C     .. Array Arguments ..
C     DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),DEN(0:LMAXD1),
C    +               EKL(0:LMAXD),FZ(IRMD,0:LMAXD),GMAT(LMMAXD,LMMAXD),
C    +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD),
C    +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD),
C    +               SZ(IRMD,0:LMAXD)
C     DOUBLE PRECISION CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD),
C    +                 THETAS(IRID,NFUND)
C     INTEGER ICLEB(NCLEB,3),IFUNM(*),IRCUT(0:IPAND),
C    +        JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(*)

      DOUBLE COMPLEX AR((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE COMPLEX CR((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE COMPLEX DEN(0:LMAX+1)
      DOUBLE COMPLEX EKL(0:LMAX)
      DOUBLE COMPLEX FZ(IRMD,0:LMAX)
      DOUBLE COMPLEX GMAT((LMAX+1)**2,(LMAX+1)**2)
      DOUBLE COMPLEX PNS((LMAX+1)**2,(LMAX+1)**2,IRMD-IRNSD:IRMD,2)
      DOUBLE COMPLEX PZ(IRMD,0:LMAX)
      DOUBLE COMPLEX QNS((LMAX+1)**2,(LMAX+1)**2,IRMD-IRNSD:IRMD,2)
      DOUBLE COMPLEX QZ(IRMD,0:LMAX)
      DOUBLE COMPLEX SZ(IRMD,0:LMAX)

C     DOUBLE PRECISION CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD),
C    +                 THETAS(IRID,NFUND)
C     INTEGER ICLEB(NCLEB,3),IFUNM(*),IRCUT(0:IPAND),
C    +        JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(*)

      DOUBLE PRECISION CLEB(*)
      DOUBLE PRECISION DRDI(IRMD)
      DOUBLE PRECISION RHO2NS(IRMD,(2*LMAX+1)**2)
      DOUBLE PRECISION THETAS(IRID,NFUND)

      INTEGER ICLEB(NCLEB,3)
      INTEGER IFUNM(*)
      INTEGER IRCUT(0:IPAND)
      INTEGER JEND((2*LMAX+1)**2,0:LMAX,0:LMAX)
      INTEGER LMSP(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX DENNS,V1
      INTEGER IMT1,L,LM,M
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CDEN(IRMD,0:LMAX)
      DOUBLE COMPLEX CDENNS(IRMD)
      DOUBLE COMPLEX EFAC((LMAX+1)**2)
C     ..
C     .. External Subroutines ..
      EXTERNAL CSIMPK,RHOIN,RHOOUT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE

      INTEGER IRMIND
      IRMIND=IRMD-IRNSD
c
c---> set up efac(lm) = sqrt(e))**l/(2l - 1)!!
c
      EFAC(1) = 1.0D0
      V1 = 1.0D0
      DO 20 L = 1,LMAX
        V1 = V1*EK/DBLE(2*L-1)
        DO 10 M = -L,L
          LM = L* (L+1) + M + 1
          EFAC(LM) = V1
   10   CONTINUE
   20 CONTINUE
c
      IMT1 = IRCUT(1)
c

      CALL RHOOUT(CDEN,DF,GMAT,EK,PNS,QNS,RHO2NS,THETAS,IFUNM,IPAN,
     +            IMT1,LMSP,CDENNS,NSRA,CLEB,ICLEB,IEND,
     &            lmax, irmd, irnsd, irid, ipand, nfund, ncleb)

c
      CALL RHOIN(AR,CDEN,CR,DF,GMAT,EK,RHO2NS,IRMIND,NSRA,EFAC,PZ,FZ,
     +             QZ,SZ,CLEB,ICLEB,JEND,IEND,EKL,
     &             lmax, irmd, ncleb)

c
c---> calculate complex density of states
c
      DO 30 L = 0,LMAX
c
c---> call integration subroutine
c
        CALL CSIMPK(CDEN(1,L),DEN(L),IPAN,IRCUT,DRDI)
   30 CONTINUE

      IF (IPAN.GT.1) THEN
        CALL CSIMPK(CDENNS,DENNS,IPAN,IRCUT,DRDI)
        DEN(LMAX+1) = DENNS
      END IF

      RETURN
      END
