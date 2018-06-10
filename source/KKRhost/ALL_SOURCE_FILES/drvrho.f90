SUBROUTINE drvrho_qdos(ldorhoef,rho2ns,r2nef,den,dmuorb,rhotborb,  &
    iecurr,eryd,we,ielast, gmatll,vt,bt,r,drdi,r2drdi,zat,  &
    jws,ishift,solver,soctl,ctl,qmtet,qmphi,  &
    itermvdir,mvevil,mvevilef,lmmaxd,lmaxd,irmd, lmpotd,iemxd,nmvecmax,  &
    i1,nqdos)                       ! qdos ruess
!   ********************************************************************
!   *                                                                  *
!   * driving routine to call relativistic routines                    *
!   *          < SSITE >, < SCFCHRDNS >, < CALCMVEC >                  *
!   * to calculate the charge and spin density in the REL mode         *
!   * v.popescu, munich, may 2004                                      *
!   *                                                                  *
!   ********************************************************************
use mod_types, only: t_tgmat
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
COMPLEX*16 CONE,CZERO
PARAMETER ( CONE=(1.0D0,0.0D0), CZERO = (0.0D0,0.0D0))
DOUBLE PRECISION DZERO
PARAMETER ( DZERO=0.0D0 )

! Dummy arguments
INTEGER LMAXD,LMMAXD,IRMD,IELAST
INTEGER ZAT(NTMAX),JWS(NMMAX),ISHIFT
INTEGER LMPOTD,IEMXD,I1
LOGICAL LDORHOEF
COMPLEX*16 WE,ERYD
DOUBLE PRECISION RHO2NS(IRMD,LMPOTD,2),R2NEF(IRMD,LMPOTD,2)
!  DOUBLE PRECISION VT(NRMAX,NTMAX),BT(NRMAX,NTMAX)
DOUBLE PRECISION VT(NRMAX),BT(NRMAX)
DOUBLE PRECISION R(NRMAX,NMMAX),R2DRDI(NRMAX,NMMAX)
DOUBLE PRECISION DRDI(NRMAX,NMMAX),SOCTL(NTMAX,NLMAX)
DOUBLE PRECISION CTL(NTMAX,NLMAX)
DOUBLE COMPLEX GMATLL(LMMAXD,LMMAXD,IEMXD), &
     DEN(0:LMAXD+1,2*IELAST)
! l-resolved orbital polarisation 
COMPLEX*16 DMUORB(0:LMAXD,3)
! orbital density
REAL*8 RHOTBORB(IRMD)

! Local variables
REAL*8 AMEOPO(NKMMAX,NKMMAX,NLAMAX,3),AT(NRMAX,NLAMAX,3,NTMAX), &
       BCOR(NTMAX),BCORS(NTMAX),CONC(NTMAX), &
       DOS(NTMAX),DOSI(NTMAX), &
       EFERMI,HFF(NTMAX), &
       HFFI(NTMAX),MUEORB,MUESPN,NVALTOT,OMT(NTMAX),OMTI(NTMAX), &
       QEL(NTMAX), &
       RHOORB(NRMAX,NTMAX),RHOCHR(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX)
REAL*8 SHFTEF,SMT(NTMAX),SMTI(NTMAX),PI,SQPI,TOTDOS
COMPLEX*16 BZJ(LINMAX,NTMAX),BZZ(LINMAX,NTMAX), &
           DOSINT(NLMAX,NTMAX),DOSL0(NLMAX,NTMAX), &
           DOSM(NMUEMAX),DZJ(LINMAX,NTMAX), &
           DZZ(LINMAX,NTMAX),EBAND,EBANDT(NTMAX), &
           HFFINT(NLMAX,NTMAX),HFFL0(NLMAX,NTMAX),HFFM(NMUEMAX), &
           MSST(NKMMAX,NKMMAX,NTMAX),OMTINT(NLMAX,NTMAX), &
           OMTL0(NLMAX,NTMAX),OMTM(NMUEMAX),OZJ(LINMAX,NTMAX), &
           OZZ(LINMAX,NTMAX),P, &
           QZJ(LINMAX,NTMAX),QZZ(LINMAX,NTMAX), &
           SMTINT(NLMAX,NTMAX),SMTL0(NLMAX,NTMAX),SMTM(NMUEMAX), &
           SZJ(LINMAX,NTMAX),SZZ(LINMAX,NTMAX)
COMPLEX*16 TAUT(NKMMAX,NKMMAX,NTMAX),OMTLS0(NLMAX,NTMAX,2)
COMPLEX*16 OZZS(LINMAX,NTMAX,2),OZJS(LINMAX,NTMAX,2)
COMPLEX*16 TAUTLIN(LINMAX,NTMAX), &
           TSST(NKMMAX,NKMMAX,NTMAX),TSSTLIN(LINMAX,NTMAX), &
           TZJ(LINMAX,NTMAX),TZZ(LINMAX,NTMAX)
LOGICAL CALCINT,GETIRRSOL
REAL*8 CGC(NKMPMAX,2)
REAL*8 GDIA(NKMMAX),GMDIA(NKMMAX),GOFF(NKMMAX),GMOFF(NKMMAX)
REAL*8 FDIA(NKMMAX),FMDIA(NKMMAX),FOFF(NKMMAX),FMOFF(NKMMAX)
INTEGER I,IECURR,IHYPER,IKM1LIN(LINMAX),IKM2LIN(LINMAX),IL, &
        IMT(NTMAX),IMUE,IP,IPRINT,IQ,IQAT(NQMAX,NTMAX),IREL, &
        IT,IWRIRRWF,IWRREGWF,J,LIN,LOPT(NTMAX), &
        MMAX,NAT(NTMAX),NETAB,NKM,NKMQ(NQMAX),NL, &
        NLINQ(NQMAX),NLQ(NQMAX),NT,NUCLEUS
INTEGER NFILCBWF,IOL
INTEGER NSOLLM(NLMAX,NMUEMAX),LTAB(NMUEMAX),LBTAB(NMUEMAX), &
     KAPTAB(NMUEMAX),NMUETAB(NMUEMAX)
CHARACTER*10 SOLVER
CHARACTER*4 TXTT(NTMAX)
DOUBLE COMPLEX W1(LMMAXD,LMMAXD)
INTEGER ICALL,IEC
!qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos
COMPLEX*16 GMAT0(LMMAXD,LMMAXD)                    !qdos ruess
INTEGER NQDOS,IREC,IPOINT                          !qdos ruess 
!qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos
INTRINSIC ATAN,SQRT

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ITERMDIR

LOGICAL ITERMVDIR,SPLITSS
INTEGER NMVECMAXD,NMVECMAX
PARAMETER (NMVECMAXD=4)
REAL*8 AMEMVEC(NKMMAX,NKMMAX,3,NMVECMAXD),FACT(0:100)
INTEGER IMKMTAB(NKMMAX),IKMLLIM1(NKMMAX),IKMLLIM2(NKMMAX)
CHARACTER*1 TXTL(0:NLMAX)
INTEGER IGRID(2),IEPATH,NEPATH

REAL*8 QMTET,QMPHI        ! ARG. LIST
REAL*8 QMPHILOC(NQMAX),QMTETLOC(NQMAX) ! DUMMY

COMPLEX*16 BMVEVDL0(NLMAX,NTMAX,3,NMVECMAX), &
     BMVEVIL1(NLMAX,NTMAX,3,NMVECMAX), &
     MVEVDL0(NLMAX,NTMAX,3,NMVECMAX), &
     MVEVIL1(NLMAX,NTMAX,3,NMVECMAX)
COMPLEX*16 MVEVIL(0:LMAXD,3,NMVECMAX) ! OUTPUT
COMPLEX*16 MVEVILEF(0:LMAXD,3,NMVECMAX) ! OUTPUT
!.. dummy arrays
COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMVECMAX), &
           MEZZ(NKMMAX,NKMMAX,NTMAX,NMVECMAX)

! ITERMDIR
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!..
!.. External Subroutines ..
EXTERNAL AMEMAGVEC,CALCCGC,CALCGF,CALCMVEC,CINIT,IKMLIN,RINIT, &
         SCFCHRDNS,SSITE,ZCOPY,ZGEMM

DATA ICALL / 0 /

SAVE ICALL,IKM1LIN,IKM2LIN,GDIA,GMDIA,GOFF,LOPT,NLQ,NKMQ, &
     IQAT,IREL,BCOR,BCORS,QEL,NAT,CONC,TXTT,IMT,SHFTEF, &
     NVALTOT,NKM,IHYPER,IPRINT,IT,IQ,NL,NT,NUCLEUS,CGC, &
     IWRREGWF,IWRIRRWF,CALCINT,GETIRRSOL,NFILCBWF,PI,SQPI, &
     AMEMVEC,IMKMTAB,IKMLLIM1,IKMLLIM2,FACT,SPLITSS, &
     TXTL,IGRID,IEPATH,NEPATH

icall = icall + 1

!=======================================================================
!       initialise relativistic and dummy variables and SAVE them
!=======================================================================
IF ( icall == 1 ) THEN
  
  IF ( lmaxd > nlmax-1) THEN
    WRITE(6,*) ' LMAXD = ',lmaxd, ' > NLMAX-1 = ',nlmax - 1
    STOP  ' Increase NLMAX in < DRVRHO > '
  END IF
  
  IF ( irmd > nrmax ) THEN
    WRITE(6,*) ' IRMD = ',irmd, ' > NRMAX = ',nrmax
    WRITE(6,*) ' Increase NRMAX in < sprkkr_rmesh.dim > '
    STOP ' In < DRVRHO > '
  END IF
  
  IF ( nmvecmax > nmvecmaxd ) THEN
    WRITE (6,*) ' NMVECMAX = ',nmvecmax,' > NMVECMAXD ', nmvecmaxd
    WRITE (6,*) ' Increase NVECMAXD in < DRVRHO > ',  &
        'or reduce NMVECMAX in < main1c > '
    STOP ' In < DRVRHO > '
  END IF
  
  iprint = 0
  nl = lmaxd + 1
  
  DO i = 1,nmuemax
    ltab(i) = i/2
    IF( 2*ltab(i) == i ) THEN
      lbtab(i)  = ltab(i) - 1
      kaptab(i) = ltab(i)
    ELSE
      lbtab(i)  =  ltab(i) + 1
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
  
  CALL calcgf(nkmax,cgc,gdia,gmdia,goff,gmoff,fdia,fmdia,  &
      foff,fmoff,ltab,lbtab,kaptab,nmuetab, nmuemax,nkmmax,nkmpmax)
  
  DO it = 1,ntmax
    bcor(it) = 0D0
    bcors(it) = 0D0
    qel(it) = 0D0
    nat(it) = 1
    conc(it) = 1D0
    txtt(it) = '    '
    imt(it) = 1
    lopt(it) = -1       ! this should change for Brooks' OP
  END DO
  
  DO iq = 1,nqmax
    nlq(iq) = nl
    nkmq(iq) = lmmaxd
    nlinq(iq) = 2*nlq(iq)*(2*nlq(iq)-1)
    iqat(iq,1) = 1
  END DO
  
  irel = 3
  shftef = 0D0
  efermi = 0D0
  nvaltot = 0
  pi = 4D0*ATAN(1D0)
  sqpi = SQRT(pi)
  
  nkm = lmmaxd
  ihyper = 0
  it = 1
  nt = 1
  iq = 1
  nucleus = 0
  
  iwrregwf = 1
  iwrirrwf = 1
  calcint = .true.
  getirrsol = .true.
  nfilcbwf = 87
!     Length in Bytes
  iol = 8*4 + 3 + (16*4*nrmax)
  OPEN (nfilcbwf,STATUS='SCRATCH',FORM='UNFORMATTED',  &
      ACCESS='DIRECT',RECL=iol)
  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ITERMDIR
  
  IF ( itermvdir ) THEN
    splitss = .false.
    fact(0) = 1.0D0
    DO i=1,100
      fact(i) = fact(i-1)*DBLE(i)
    END DO
    
    CALL amemagvec(irel,iprint+1,nkm,amemvec,ikmllim1,ikmllim2,  &
        imkmtab,cgc,nlmax,nkmmax,nkmpmax,nmvecmax)
    
    DO i = 0,nlmax
      txtl(i) = ' '
    END DO
    igrid(1) = 5
    iepath = 1
    nepath = 1
  END IF
  
!      ITERMDIR
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
END IF                    ! ICALL.EQ.1
!=======================================================================

CALL ssite(iwrregwf,iwrirrwf,nfilcbwf,calcint,getirrsol,soctl,ctl,  &
    eryd,p,ihyper,iprint,ikm1lin,ikm2lin,nlq,nkmq,nlinq,nt,  &
    nkm,iqat,tsst,msst,tsstlin,dzz,dzj,szz,szj,ozz,ozj,bzz,  &
    bzj,qzz,qzj,tzz,tzj,vt,bt,at,zat,nucleus,r,drdi,r2drdi,  &
    jws,imt,ameopo,lopt,solver,cgc,ozzs,ozjs,nlmax,nqmax,  &
    linmax,nrmax,nmmax,ntmax,nkmmax,nkmpmax,nlamax)

!-----------------------------------------------------------------------
!     get charge density
!-----------------------------------------------------------------------

netab = iecurr + 1
iec = iecurr

! Loop over all qdos points specified in qvec.dat
DO  ipoint = 1,nqdos                                        ! qdos ruess
!                                                                    ! qdos ruess
! Read in Green function; remember that for the rel. case, nspin = 1 ! qdos ruess
! (without qdos, IPOINT=NQDOS=1)                                     ! qdos ruess
  irec = ipoint + nqdos * (iecurr-1) +  nqdos * ielast * (i1-1)  ! qdos ruess
  IF (t_tgmat%gmat_to_file) THEN
    READ(69,REC=irec) gmat0                                        ! qdos ruess
  ELSE
    gmat0(:,:) = t_tgmat%gmat(:,:,irec)
  END IF
  gmatll(:,:,iecurr) = gmat0(:,:)                                ! qdos ruess
!                                                                    ! qdos ruess
  
!-------- GET TAU MATRIX ------------------------
!         TAUT = t G t + t
  
! ---> taut = t
  
  DO j = 1,nkm
    CALL zcopy(nkm,tsst(1,j,it),1,taut(1,j,it),1)
  END DO
  
! ---> w1 = G * t
  
  CALL zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,gmatll(1,1,iecurr),  &
      lmmaxd,tsst(1,1,it),nkmmax,czero,w1,lmmaxd)
  
! ---> taut = t * G * t + t = t * w1 + taut
  
  CALL zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,tsst(1,1,it),nkmmax,  &
      w1,lmmaxd,cone,taut(1,1,it),nkmmax)
  
! ---> store taut in linear array tautlin
  
  DO lin = 1,nlinq(iq)
    tautlin(lin,it) = taut(ikm1lin(lin),ikm2lin(lin),it)
  END DO
  
  CALL rinit(nrmax,rhochr(1,it))
  CALL rinit(nrmax,rhospn(1,it))
  CALL rinit(nrmax,rhoorb(1,it))
  CALL cinit(nlmax,omtl0(1,it))
  DO lin = 1,2
    CALL cinit(nlmax,omtls0(1,it,lin))
  END DO
  
  CALL scfchrdns(nfilcbwf,r2drdi,jws,imt,shftef,totdos,muespn,  &
      mueorb,irel,iprint,nt,nl,nkm,eryd,we,efermi,iec,  &
      netab,dos,smt,omt,hff,dosi,smti,omti,hffi,dosm,  &
      dosl0,dosint,smtm,smtl0,smtint,omtm,omtl0,omtint,  &
      hffm,hffl0,hffint,bcor,bcors,dzz,dzj,szz,szj,ozz,  &
      ozj,bzz,bzj,ozzs,ozjs,omtls0,tautlin,nvaltot,txtt,  &
      conc,nat,rhochr,rhospn,rhoorb,qel,gdia,gmdia,goff,  &
      ntmax,nlmax,nmuemax,linmax,nrmax,nmmax,nkmmax, eband,ebandt)
  
  DO i = 1,ishift
    rho2ns(i,1,1) = dzero
    rho2ns(i,1,2) = dzero
    rhotborb(i) = dzero
  END DO
  
  DO i = 1,jws(it)
    ip = i + ishift
    rho2ns(ip,1,1) = rho2ns(ip,1,1)  &
        - 0.5D0 * sqpi * rhochr(i,it) * (r(i,1)**2)
    rho2ns(ip,1,2) = rho2ns(ip,1,2)  &
        - 0.5D0 * sqpi * rhospn(i,it) * (r(i,1)**2)
    rhotborb(ip) = rhotborb(ip) - 0.5D0 * sqpi * rhoorb(i,it) * (r(i,1)**2)
  END DO
  
  DO il = 1,nl
    den(il-1,iecurr+ielast) = -0.5D0 * (dosl0(il,it)+smtl0(il,it)) * pi
    
    den(il-1,iecurr) = -0.5D0 * (dosl0(il,it)-smtl0(il,it)) * pi
    
    DO i = 1,2
      dmuorb(il-1,i) = -omtls0(il,it,i) * pi
    END DO
    
    dmuorb(il-1,3) = -omtl0(il,it) * pi
  END DO
  den(nl,iecurr+ielast) = czero
  den(nl,iecurr) = czero
  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ITERMDIR
  
  IF ( itermvdir ) THEN
    
    qmphiloc(iq) = qmphi
    qmtetloc(iq) = qmtet
    
    CALL cinit(nlmax*ntmax*3*nmvecmax,mvevdl0)
    CALL cinit(nlmax*ntmax*3*nmvecmax,bmvevdl0)
    CALL cinit(nlmax*ntmax*3*nmvecmax,mvevil1)
    CALL cinit(nlmax*ntmax*3*nmvecmax,bmvevil1)
    
    CALL calcmvec(nfilcbwf,splitss,iepath,nepath,irel,  &
        iprint,nt,nl,mezz,mezj,taut,tsst,iqat,nkmq,nkm,  &
        iec,netab,igrid(iepath),we,mvevdl0, mvevil1,bmvevdl0,bmvevil1,  &
        r2drdi,jws,imt,amemvec, ikmllim1,ikmllim2,imkmtab,ntmax,nlmax,nmuemax,  &
        nqmax,nkmmax,nmmax,nmvecmax,nrmax)
    
    DO i=1,nmvecmax
      DO j=1,3
        DO il=1,nl
          mvevil(il-1,j,i) = mvevil(il-1,j,i) - mvevdl0(il,it,j,i) * pi * we
        END DO
      END DO
    END DO
  END IF
  
!      ITERMDIR
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  
  enddo ! IPOINT = 1,NQDOS

IF ( (iecurr /= ielast) .OR. (.NOT.ldorhoef) ) RETURN

! ======================================================================
!     get the charge at the Fermi energy (IELAST)
!     call SCFCHRDNS with the energy weight CONE --> not overwrite WE

CALL rinit(nrmax,rhochr(1,it))
CALL rinit(nrmax,rhospn(1,it))

CALL scfchrdns(nfilcbwf,r2drdi,jws,imt,shftef,totdos,muespn,  &
    mueorb,irel,iprint,nt,nl,nkm,eryd,cone,efermi,iec,  &
    netab,dos,smt,omt,hff,dosi,smti,omti,hffi,dosm,  &
    dosl0,dosint,smtm,smtl0,smtint,omtm,omtl0,omtint,  &
    hffm,hffl0,hffint,bcor,bcors,dzz,dzj,szz,szj,ozz,  &
    ozj,bzz,bzj,ozzs,ozjs,omtls0,tautlin,nvaltot,txtt,  &
    conc,nat,rhochr,rhospn,rhoorb,qel,gdia,gmdia,goff,  &
    ntmax,nlmax,nmuemax,linmax,nrmax,nmmax,nkmmax, eband,ebandt)

DO i = 1,ishift
  r2nef(i,1,1) = dzero
  r2nef(i,1,2) = dzero
END DO

DO i = 1,jws(it)
  ip = i + ishift
  r2nef(ip,1,1) = - 0.5D0 * sqpi * rhochr(i,it) * (r(i,1)**2)
  r2nef(ip,1,2) = - 0.5D0 * sqpi * rhospn(i,it) * (r(i,1)**2)
END DO

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ITERMDIR

IF( itermvdir ) THEN
  DO i=1,nmvecmax
    DO j=1,3
      DO il=1,nl
        mvevilef(il-1,j,i) = -mvevdl0(il,it,j,i) * pi
      END DO
    END DO
  END DO
END IF

!      ITERMDIR
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


! ======================================================================

END SUBROUTINE drvrho_qdos
