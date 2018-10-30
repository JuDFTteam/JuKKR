      MODULE MOD_RHONS
      CONTAINS
!-------------------------------------------------------------------------
!> Summary: Driver for valence charge density for non-spherical potential
!> Category: physical-observables, KKRimp
!>
!>     the charge density is developed in spherical harmonics :
!>
!>             rho(r) =   { rho(lm,r) * y(r,lm) }       (summed over lm)
!>
!>          rho(lm,r) =   { do rho(r) * y(r,lm)         (integrated over
!>                                                           unit sphere)
!>     in the case of spin-polarization :
!>       the spin density is developed in spherical harmonics :
!>
!>            sden(r) =   { sden(lm,r) * y(r,lm) }      (summed over lm)
!>
!>         sden(lm,r) =   { do sden(r) * y(r,lm)        (integrated over
!>                                                           unit sphere)
!>     n(r,e) is developed in
!>
!>        n(r,e) = { y(r,l'm') * n(l'm',lm,r,e) * y(r,lm) }
!>
!>     therefore a faltung of n(l'm',lm,r,e) with the gaunt coeffients
!>     has to be used to calculate the lm-contribution of the charge
!>     density .
!>
!>
!>     calculate the valence density of states , in the spin-polarized
!>      case spin dependent .
!>     recognize that the density of states is always complex also in
!>      the case of "real-energy-integation" (ief>0) since in that case
!>      the energy integration is done parallel to the real energy axis
!>      but not on the real energy axis .
!>     in the last energy-spin loop the l-contribution of the valence
!>      charge is calculated .
!>
!>                               b.drittler   aug. 1988
!>
!>     modified for the use of shape functions
!>
!>     attention : irmin + 3 has to be less then imt
!>                 if shape functions are used
!>
!>                               b.drittler   july 1989
!>-----------------------------------------------------------------------
      SUBROUTINE RHONS(DEN,DENLM,DF,DRDI,GMAT,EK,RHO2NS,IPAN,IRCUT,
     +                 THETAS,
     +                   IFUNM,LMSP,NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB,
     +                   ICLEB,JEND,IEND,EKL,
     +                   IRID,NFUND,IRMIND,
     +                   IRMD,NCLEB,LMAXD,LMMAXD,LMPOTD)
c-----------------------------------------------------------------------
C     .. Parameters ..
!       INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *                                                                   *
C *********************************************************************
C
      USE MOD_CSIMPK
      USE MOD_RHOIN
      USE MOD_RHOOUT
      IMPLICIT NONE
      INTEGER IRMIND,NCLEB
      INTEGER IRMD,IRID,NFUND
!       PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER LMAXD,LMPOTD,LMMAXD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
!       parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
!       INTEGER LMAXD1
!       PARAMETER (LMAXD1= LMAXD+1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,EK
      INTEGER IEND,IPAN,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),DEN(0:LMAXD+1),
     +               EKL(0:LMAXD),FZ(IRMD,0:LMAXD),GMAT(LMMAXD,LMMAXD),
     +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD),
     +               SZ(IRMD,0:LMAXD)
     +              ,DENLM(LMMAXD),DENLM2(LMMAXD),ENERG  ! lm-dos
      DOUBLE PRECISION CLEB(NCLEB),DRDI(IRMD),RHO2NS(IRMD,LMPOTD),
     +                 THETAS(IRID,NFUND)
      INTEGER ICLEB(NCLEB,4),IFUNM(*),IRCUT(0:IPAN),
     +        JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX DENNS,V1
      INTEGER IMT1,L,LM,M
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CDEN(IRMD,0:LMAXD),CDENNS(IRMD),EFAC(LMMAXD)
     +              ,CDENLM(IRMD,LMMAXD)  ! lm-dos
C     ..
C     .. External Subroutines ..
!       EXTERNAL RHOIN,RHOOUT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
      DOUBLE PRECISION PI
      PI = 4.D0*DATAN(1.D0)
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
     +              IMT1,LMSP,CDENNS,NSRA,CLEB,ICLEB,IEND
     +             ,CDENLM, ! lm-dos
     +              NCLEB,LMAXD,LMMAXD,LMPOTD,IRMD,IRMIND,IRID,
     +              NFUND)
c
      CALL RHOIN(AR,CDEN,CR,DF,GMAT,EK,RHO2NS,IRMIND,NSRA,EFAC,PZ,FZ,
     +             QZ,SZ,CLEB,ICLEB,JEND,IEND,EKL
     +            ,CDENLM,  ! lm-dos
     +             NCLEB,LMAXD,LMMAXD,LMPOTD,IRMD)
c
c---> calculate complex density of states
c
      DO 30 L = 0,LMAXD
c
c---> call integration subroutine
c
        CALL CSIMPK(CDEN(1,L),DEN(L),IPAN,IRCUT,DRDI)
   30 CONTINUE

      DO 40 L = 1,LMMAXD  ! lm-dos
        CALL CSIMPK(CDENLM(1,L),DENLM2(L),IPAN,IRCUT,DRDI)  ! lm-dos
        DENLM(L)=DENLM2(L)
   40 CONTINUE  ! lm-dos

c Energy depends on EK and NSRA:                            ! lm-dos
c     IF (NSRA.EQ.1) EK = SQRT(E)                           ! lm-dos
c     IF (NSRA.EQ.2) EK = SQRT(E+E*E/ (CVLIGHT*CVLIGHT))    ! lm-dos
c     CVLIGHT=274.0720442D0                                 ! lm-dos 
c Therefore the following is a good approximation           ! lm-dos
c for energies of a few Ryd:                                ! lm-dos
      ENERG = EK**2                                         ! lm-dos

!       WRITE(30,9000) DREAL(ENERG),(-DIMAG(DENLM(LM))/PI,LM=1,LMMAXD)
!  9000 FORMAT(30E12.4)


      IF (IPAN.GT.1) THEN
        CALL CSIMPK(CDENNS,DENNS,IPAN,IRCUT,DRDI)
        DEN(LMAXD+1) = DENNS
      END IF

      RETURN
      END SUBROUTINE RHONS
      END MODULE MOD_RHONS