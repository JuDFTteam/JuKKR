      SUBROUTINE RHONS(DEN,DF,DRDI,GMAT,EK,RHO2NS,IPAN,IRCUT,THETAS,
     +                   IFUNM,LMSP,NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB,
     +                   ICLEB,JEND,IEND,EKL)
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
C     .. Parameters ..
      INCLUDE 'inc.p'
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER LMPOTD,LMMAXD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER LMAXD1
      PARAMETER (LMAXD1= LMAXD+1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,EK
      INTEGER IEND,IPAN,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),DEN(0:LMAXD1),
     +               EKL(0:LMAXD),FZ(IRMD,0:LMAXD),GMAT(LMMAXD,LMMAXD),
     +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD),
     +               SZ(IRMD,0:LMAXD)
      DOUBLE PRECISION CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD),
     +                 THETAS(IRID,NFUND)
      INTEGER ICLEB(NCLEB,3),IFUNM(*),IRCUT(0:IPAND),
     +        JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX DENNS,V1
      INTEGER IMT1,L,LM,M
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CDEN(IRMD,0:LMAXD),CDENNS(IRMD),EFAC(LMMAXD)
C     ..
C     .. External Subroutines ..
      EXTERNAL CSIMPK,RHOIN,RHOOUT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
c
c---> set up efac(lm) = sqrt(e))**l/(2l - 1)!!
c
      EFAC(1) = 1.0D0
      V1 = 1.0D0
      DO 20 L = 1,LMAXD
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
     +              IMT1,LMSP,CDENNS,NSRA,CLEB,ICLEB,IEND)

c
      CALL RHOIN(AR,CDEN,CR,DF,GMAT,EK,RHO2NS,IRMIND,NSRA,EFAC,PZ,FZ,
     +             QZ,SZ,CLEB,ICLEB,JEND,IEND,EKL)

c
c---> calculate complex density of states
c
      DO 30 L = 0,LMAXD
c
c---> call integration subroutine
c
        CALL CSIMPK(CDEN(1,L),DEN(L),IPAN,IRCUT,DRDI)
   30 CONTINUE

      IF (IPAN.GT.1) THEN
        CALL CSIMPK(CDENNS,DENNS,IPAN,IRCUT,DRDI)
        DEN(LMAXD1) = DENNS
      END IF

      RETURN
      END
