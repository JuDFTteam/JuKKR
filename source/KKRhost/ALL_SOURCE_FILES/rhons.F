      SUBROUTINE RHONS(DEN,DF,DRDI,GMAT,EK,RHO2NS,IPAN,IRCUT,IRMIN,    ! Added IRMIN 1.7.2014
     +                   THETAS,IFUNM,LMSP,NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,
     +                   SZ,CLEB,ICLEB,JEND,IEND,EKL,DENLM,GFLLE_PART)
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
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER LMPOTD,LMMAXD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
      INTEGER LMAXD1
      PARAMETER (LMAXD1= LMAXD+1)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,EK
      INTEGER IEND,IPAN,NSRA,IRMIN
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),DEN(0:LMAXD1),
     +               EKL(0:LMAXD),FZ(IRMD,0:LMAXD),GMAT(LMMAXD,LMMAXD),
     +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD),
     +               SZ(IRMD,0:LMAXD)
     +              ,DENLM(LMMAXD),ENERG  ! lm-dos
      DOUBLE PRECISION CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD),
     +                 THETAS(IRID,NFUND)
     +                ,CVLIGHT              ! lm-dos
      INTEGER ICLEB(NCLEB,4),IFUNM(*),IRCUT(0:IPAND),
     +        JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX DENNS,V1
      INTEGER IMT1,L,LM,M,IRMAX
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CDEN(IRMD,0:LMAXD),CDENNS(IRMD),EFAC(LMMAXD)
     +              ,CDENLM(IRMD,LMMAXD),CWR(IRMD,LMMAXD,LMMAXD)  ! lm-dos
     +              ,GFLLE_PART(LMMAXD,LMMAXD)
C     ..
C     .. Arrays in Common ..
      CHARACTER*8 OPTC(32)            ! qdos
C     ..
C     .. Common blocks ..
      COMMON /OPTC/OPTC                         ! qdos
C     .. External Functions ..
      LOGICAL OPT                          ! qdos
      EXTERNAL OPT                         ! qdos
C     ..
C     .. External Subroutines ..
      EXTERNAL CSIMPK,RHOIN,RHOOUT
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
      IRMAX = IRCUT(IPAN)
c
      CALL RHOOUT(CDEN,DF,GMAT,EK,PNS,QNS,RHO2NS,THETAS,IFUNM,IPAN,
     +              IMT1,IRMIN,IRMAX,LMSP,CDENNS,NSRA,CLEB,ICLEB,IEND      ! Added IRMIN,IRMAX 1.7.2014
     +             ,CDENLM,CWR)  ! lm-dos
c
      CALL RHOIN(AR,CDEN,CR,DF,GMAT,EK,RHO2NS,IRMIN,NSRA,EFAC,PZ,FZ,   ! Changed from IRMIND to IRMIN 1.7.2014
     +             QZ,SZ,CLEB,ICLEB,JEND,IEND,EKL
     +            ,CDENLM)  ! lm-dos  ! Attention, cwr does not go into rhoin, does lmlm-dos work properly?

c
c---> calculate complex density of states
c
      DO 30 L = 0,LMAXD
c
c---> call integration subroutine
c
        CALL CSIMPK(CDEN(1,L),DEN(L),IPAN,IRCUT,DRDI)
   30 CONTINUE

      DO 40 LM1 = 1,LMMAXD  ! lm-dos
        CALL CSIMPK(CDENLM(1,LM1),DENLM(LM1),IPAN,IRCUT,DRDI)  ! lm-dos
        IF (OPT('lmlm-dos').OR.OPT('qdos    ').OR.           ! lmlm-dos
     +      OPT('LDA+U   ')) THEN                            ! LDAU
          DO 41 LM2 = 1,LMMAXD                               ! lmlm-dos
            CALL CSIMPK(CWR(1,LM1,LM2),GFLLE_PART(LM1,LM2),  ! lmlm-dos
     +                  IPAN,IRCUT,DRDI)                     ! lmlm-dos
   41     CONTINUE                                           ! lmlm-dos
        ENDIF                                                ! lmlm-dos
   40 CONTINUE  ! lm-dos

c Energy depends on EK and NSRA:                            ! lm-dos
c     IF (NSRA.EQ.1) EK = SQRT(E)                           ! lm-dos
c     IF (NSRA.EQ.2) EK = SQRT(E+E*E/ (CVLIGHT*CVLIGHT))    ! lm-dos
c     CVLIGHT=274.0720442D0                                 ! lm-dos 
c Therefore the following is a good approximation           ! lm-dos
c for energies of a few Ryd:                                ! lm-dos
#ifndef CPP_MPI
      IF (.NOT.OPT('qdos    ')) THEN
      ENERG = EK**2                                         ! lm-dos
      WRITE(30,9000) DREAL(ENERG),(-DIMAG(DENLM(LM))/PI,LM=1,LMMAXD)
 9000 FORMAT(30E12.4)
      ENDIF  ! not qdos option
#endif


      IF (IPAN.GT.1) THEN
        CALL CSIMPK(CDENNS,DENNS,IPAN,IRCUT,DRDI)
        DEN(LMAXD1) = DENNS
      END IF

      RETURN
      END
