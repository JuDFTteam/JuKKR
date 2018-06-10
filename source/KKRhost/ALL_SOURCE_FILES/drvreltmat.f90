SUBROUTINE drvreltmat(eryd,tmatll,vt,bt,r,drdi,r2drdi,zat,jws,  &
    solver,soctl,ctl,lmmaxd,lmaxd,irmd)
!   ********************************************************************
!   *                                                                  *
!   * driving routine to call relativistic < SSITE > routine           *
!   * only to calculate the single-site t matrix                       *
!   * v.popescu, munich, may 2004                                      *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

! PARAMETER definitions
INTEGER NRMAX
PARAMETER ( NRMAX=900 )
INTEGER NLAMAX,NQMAX,NTMAX,NMMAX
PARAMETER (NLAMAX=1,NQMAX=1,NTMAX=1,NMMAX=1)
INTEGER NLMAX,NKMMAX,NMUEMAX,NKMPMAX,NKMAX,LINMAX
PARAMETER ( NLMAX = 5 ) ! this should be >= LMAXD + 1
PARAMETER ( NKMMAX = 2*NLMAX**2, NKMAX = 2*NLMAX-1 )
PARAMETER ( NKMPMAX = NKMMAX+2*NLMAX, NMUEMAX = 2*NLMAX)
PARAMETER ( LINMAX = 2*NLMAX*(2*NLMAX-1) )

! Dummy arguments
INTEGER LMAXD,LMMAXD,IRMD
INTEGER ZAT(NTMAX),JWS(NMMAX)
DOUBLE COMPLEX TMATLL(LMMAXD,LMMAXD)
DOUBLE PRECISION SOCTL(NLMAX)
DOUBLE PRECISION CTL(NLMAX)
DOUBLE PRECISION VT(NRMAX),BT(NRMAX)
DOUBLE PRECISION R(NRMAX,NMMAX),R2DRDI(NRMAX,NMMAX)
DOUBLE PRECISION DRDI(NRMAX,NMMAX)

! Local variables
REAL*8 AMEOPO(NKMMAX,NKMMAX,NLAMAX,3),AT(NRMAX,NLAMAX,3,NTMAX)
COMPLEX*16 BZJ(LINMAX,NTMAX),BZZ(LINMAX,NTMAX), &
           DZJ(LINMAX,NTMAX), &
           DZZ(LINMAX,NTMAX),ERYD, &
           MSST(NKMMAX,NKMMAX,NTMAX), &
           OZJ(LINMAX,NTMAX), &
           OZZ(LINMAX,NTMAX),P, &
           QZJ(LINMAX,NTMAX),QZZ(LINMAX,NTMAX), &
           SZJ(LINMAX,NTMAX),SZZ(LINMAX,NTMAX)
COMPLEX*16 TSST(NKMMAX,NKMMAX,NTMAX),TSSTLIN(LINMAX,NTMAX), &
           TZJ(LINMAX,NTMAX),TZZ(LINMAX,NTMAX)
COMPLEX*16 OZZS(LINMAX,NTMAX,2),OZJS(LINMAX,NTMAX,2)
LOGICAL CALCINT,GETIRRSOL
REAL*8 CGC(NKMPMAX,2)
INTEGER I,IHYPER,IKM1LIN(LINMAX),IKM2LIN(LINMAX),IL, &
        IMT(NTMAX),IMUE,IPRINT,IQ,IQAT(NQMAX,NTMAX), &
        IT,IWRIRRWF,IWRREGWF,J,LOPT(NTMAX), &
        MMAX,NKM,NKMQ(NQMAX),NL, &
        NLINQ(NQMAX),NLQ(NQMAX),NT,NUCLEUS
INTEGER NFILCBWF
INTEGER NSOLLM(NLMAX,NMUEMAX),LTAB(NMUEMAX), &
        KAPTAB(NMUEMAX),NMUETAB(NMUEMAX)
CHARACTER*10 SOLVER
INTEGER ICALL

DATA ICALL / 0 /

SAVE ICALL,IKM1LIN,IKM2LIN,LOPT,NLQ,NKMQ, &
     IQAT,IMT, &
     NKM,IHYPER,IPRINT,IT,NT,NUCLEUS, &
     IWRREGWF,IWRIRRWF,CALCINT,GETIRRSOL,NFILCBWF

icall = icall + 1

!=======================================================================
!       initialise relativistic and dummy variables and SAVE them
!=======================================================================
IF ( icall == 1 ) THEN
  
  IF ( lmaxd > nlmax-1) THEN
    WRITE(6,*) ' LMAXD = ',lmaxd, ' > NLMAX-1 = ',nlmax - 1
    STOP  ' Increase NLMAX in < DRVRELTMAT > '
  END IF
  
  IF ( irmd > nrmax ) THEN
    WRITE (6,*) ' IRMD = ',irmd,' > NRMAX = ',nrmax
    WRITE (6,*) ' Increase NRMAX in < sprkkr_rmesh.dim > '
    STOP ' and in < DRVRELTMAT > '
  END IF
  
  nl = lmaxd+1           ! no need to save, used only here once
  iprint = 0
  
  DO i = 1,nmuemax
    ltab(i) = i/2
    IF( 2*ltab(i) == i ) THEN
      kaptab(i) =  ltab(i)
    ELSE
      kaptab(i) = -ltab(i) - 1
    END IF
    nmuetab(i) = 2*ABS(kaptab(i))
  END DO
  
  DO il = 1,nlmax
    mmax = 2*il
    DO imue = 1,mmax
      IF ( (imue == 1) .OR. (imue == mmax) ) THEN
        nsollm(il,imue) = 1
      ELSE
        nsollm(il,imue) = 2
      END IF
    END DO
  END DO
  
  CALL ikmlin(iprint,nsollm,ikm1lin,ikm2lin,nlmax,nmuemax, linmax,nlmax)
  
  CALL calccgc(ltab,kaptab,nmuetab,cgc,nkmax,nmuemax,nkmpmax)
  
  DO it = 1,ntmax
    imt(it) = 1
    lopt(it) = -1       ! this should change for Brooks' OP
  END DO
  
  DO iq = 1,nqmax
    nlq(iq) = nl
    nkmq(iq) = lmmaxd
    nlinq(iq) = 2*nlq(iq)*(2*nlq(iq)-1)
    iqat(iq,1) = 1
  END DO
  
  nkm = lmmaxd
  ihyper = 0
  nt = 1
  it = 1
  nucleus = 0
  
  iwrregwf = 0
  iwrirrwf = 0
  calcint = .false.
  getirrsol = .false.
  nfilcbwf = 0
END IF                    ! ICALL.EQ.1
!=======================================================================

CALL ssite(iwrregwf,iwrirrwf,nfilcbwf,calcint,getirrsol,soctl,ctl,  &
    eryd,p,ihyper,iprint,ikm1lin,ikm2lin,nlq,nkmq,nlinq,nt,  &
    nkm,iqat,tsst,msst,tsstlin,dzz,dzj,szz,szj,ozz,ozj,bzz,  &
    bzj,qzz,qzj,tzz,tzj,vt,bt,at,zat,nucleus,r,drdi,r2drdi,  &
    jws,imt,ameopo,lopt,solver,cgc,ozzs,ozjs,nlmax,nqmax,  &
    linmax,nrmax,nmmax,ntmax,nkmmax,nkmpmax,nlamax)

DO j = 1,nkm
  CALL zcopy(nkm,tsst(1,j,it),1,tmatll(1,j),1)
END DO

RETURN
END SUBROUTINE drvreltmat
