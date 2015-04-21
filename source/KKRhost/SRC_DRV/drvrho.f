C*==drvrho.f    processed by SPAG 6.05Rc at 11:40 on 10 May 2004
      SUBROUTINE DRVRHO_QDOS(LDORHOEF,RHO2NS,R2NEF,DEN,DMUORB,RHOTBORB,
     &                  IECURR,ERYD,WE,IELAST,
     &                  GMATLL,VT,BT,R,DRDI,R2DRDI,ZAT,
     &                  JWS,ISHIFT,SOLVER,SOCTL,CTL,QMTET,QMPHI,
     &                  ITERMVDIR,MVEVIL,MVEVILEF,LMMAXD,LMAXD,IRMD,
     &                  LMPOTD,IEMXD,NMVECMAX,
     &                  I1,QVEC,NQDOS)                       ! qdos ruess 
C   ********************************************************************
C   *                                                                  *
C   * driving routine to call relativistic routines                    *
C   *          < SSITE >, < SCFCHRDNS >, < CALCMVEC >                  *
C   * to calculate the charge and spin density in the REL mode         *
C   * v.popescu, munich, may 2004                                      *
C   *                                                                  *
C   ********************************************************************
      use mod_types, only: t_tgmat
      IMPLICIT NONE
C
C PARAMETER definitions
C
      INTEGER NRMAX
      PARAMETER ( NRMAX=750 )
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
C
C Dummy arguments
C
      INTEGER LMAXD,LMMAXD,IRMD,IELAST
      INTEGER ZAT(NTMAX),JWS(NMMAX),ISHIFT
      INTEGER LMPOTD,IEMXD,I1
      LOGICAL LDORHOEF
      COMPLEX*16 WE,ERYD
      DOUBLE PRECISION RHO2NS(IRMD,LMPOTD,2),R2NEF(IRMD,LMPOTD,2)
      DOUBLE PRECISION VT(NRMAX,NTMAX),BT(NRMAX,NTMAX)
      DOUBLE PRECISION R(NRMAX,NMMAX),R2DRDI(NRMAX,NMMAX)
      DOUBLE PRECISION DRDI(NRMAX,NMMAX),SOCTL(NTMAX,NLMAX)
      DOUBLE PRECISION CTL(NTMAX,NLMAX)
      DOUBLE COMPLEX GMATLL(LMMAXD,LMMAXD,IEMXD),
     &     DEN(0:LMAXD+1,2*IEMXD)
C l-resolved orbital polarisation 
      COMPLEX*16 DMUORB(0:LMAXD,3)
C orbital density
      REAL*8 RHOTBORB(IRMD)
C
C Local variables
C
      REAL*8 AMEOPC(NKMMAX,NKMMAX,NLAMAX,3),
     &       AMEOPO(NKMMAX,NKMMAX,NLAMAX,3),AT(NRMAX,NLAMAX,3,NTMAX),
     &       BCOR(NTMAX),BCORS(NTMAX),CONC(NTMAX),
     &       DOS(NTMAX),DOSI(NTMAX),
     &       EFERMI,HFF(NTMAX),
     &       HFFI(NTMAX),MUEORB,MUESPN,NVALTOT,OMT(NTMAX),OMTI(NTMAX),
     &       QEL(NTMAX),
     &       RHOORB(NRMAX,NTMAX),RHOCHR(NRMAX,NTMAX),RHOSPN(NRMAX,NTMAX)
      REAL*8 SHFTEF,SMT(NTMAX),SMTI(NTMAX),PI,SQPI,TOTDOS
      COMPLEX*16 BZJ(LINMAX,NTMAX),BZZ(LINMAX,NTMAX),
     &           DOSINT(NLMAX,NTMAX),DOSL0(NLMAX,NTMAX),
     &           DOSM(NMUEMAX),DZJ(LINMAX,NTMAX),
     &           DZZ(LINMAX,NTMAX),EBAND,EBANDT(NTMAX),
     &           HFFINT(NLMAX,NTMAX),HFFL0(NLMAX,NTMAX),HFFM(NMUEMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),OMTINT(NLMAX,NTMAX),
     &           OMTL0(NLMAX,NTMAX),OMTM(NMUEMAX),OZJ(LINMAX,NTMAX),
     &           OZZ(LINMAX,NTMAX),P,
     &           QZJ(LINMAX,NTMAX),QZZ(LINMAX,NTMAX),
     &           SMTINT(NLMAX,NTMAX),SMTL0(NLMAX,NTMAX),SMTM(NMUEMAX),
     &           SZJ(LINMAX,NTMAX),SZZ(LINMAX,NTMAX)
      COMPLEX*16 TAUT(NKMMAX,NKMMAX,NTMAX),OMTLS0(NLMAX,NTMAX,2)
      COMPLEX*16 OZZS(LINMAX,NTMAX,2),OZJS(LINMAX,NTMAX,2)
      COMPLEX*16 TAUTLIN(LINMAX,NTMAX),
     &           TSST(NKMMAX,NKMMAX,NTMAX),TSSTLIN(LINMAX,NTMAX),
     &           TZJ(LINMAX,NTMAX),TZZ(LINMAX,NTMAX)
      LOGICAL CALCINT,GETIRRSOL
      REAL*8 CGC(NKMPMAX,2)
      REAL*8 GDIA(NKMMAX),GMDIA(NKMMAX),GOFF(NKMMAX),GMOFF(NKMMAX)
      REAL*8 FDIA(NKMMAX),FMDIA(NKMMAX),FOFF(NKMMAX),FMOFF(NKMMAX)
      INTEGER I,IECURR,IHYPER,IKM1LIN(LINMAX),IKM2LIN(LINMAX),IL,
     &        IMT(NTMAX),IMUE,IP,IPRINT,IQ,IQAT(NQMAX,NTMAX),IREL,
     &        IT,IWRIRRWF,IWRREGWF,J,LIN,LOPT(NTMAX),
     &        MMAX,NAT(NTMAX),NETAB,NKM,NKMQ(NQMAX),NL,
     &        NLINQ(NQMAX),NLQ(NQMAX),NT,NUCLEUS
      INTEGER NFILCBWF,IOL
      INTEGER NSOLLM(NLMAX,NMUEMAX),LTAB(NMUEMAX),LBTAB(NMUEMAX),
     &     KAPTAB(NMUEMAX),NMUETAB(NMUEMAX)
      CHARACTER*10 SOLVER
      CHARACTER*4 TXTT(NTMAX)
      DOUBLE COMPLEX W1(LMMAXD,LMMAXD)
      INTEGER ICALL,IEC
C     qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos
      COMPLEX*16 GMAT0(LMMAXD,LMMAXD)                    !qdos ruess
      INTEGER NQDOS,IREC,IPOINT                          !qdos ruess 
      DOUBLE PRECISION QVEC(3,NQDOS)                     !qdos ruess 
      COMPLEX*16 DENTOT1,DENTOT2          !dummy arrays  !qdos ruess 
C     .. Arrays in Common ..                             !qdos ruess 
      CHARACTER*8 OPTC(32)                                !qdos ruess 
C     .. Common blocks ..                                !qdos ruess 
      COMMON /OPTC/OPTC                                  !qdos ruess 
C     qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos
      INTRINSIC ATAN,SQRT

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ITERMDIR
C
      LOGICAL ITERMVDIR,SPLITSS
      INTEGER NMVECMAXD,NMVECMAX
      PARAMETER (NMVECMAXD=4)
      REAL*8 AMEMVEC(NKMMAX,NKMMAX,3,NMVECMAXD),FACT(0:100)
      INTEGER IMKMTAB(NKMMAX),IKMLLIM1(NKMMAX),IKMLLIM2(NKMMAX)
      CHARACTER*1 TXTL(0:NLMAX)
      INTEGER IGRID(2),IEPATH,NEPATH
C
      REAL*8 MVGAM(NTMAX,NMVECMAX),MVPHI(NTMAX,NMVECMAX),
     &     MVTET(NTMAX,NMVECMAX) ! DUMMY 
C     
      REAL*8 QMTET,QMPHI        ! ARG. LIST
      REAL*8 QMPHILOC(NQMAX),QMTETLOC(NQMAX) ! DUMMY
C     
      COMPLEX*16 BMVEVDL0(NLMAX,NTMAX,3,NMVECMAX),
     &     BMVEVIL1(NLMAX,NTMAX,3,NMVECMAX),
     &     MVEVDL0(NLMAX,NTMAX,3,NMVECMAX),
     &     MVEVIL1(NLMAX,NTMAX,3,NMVECMAX)
      COMPLEX*16 MVEVIL(0:LMAXD,3,NMVECMAX) ! OUTPUT
      COMPLEX*16 MVEVILEF(0:LMAXD,3,NMVECMAX) ! OUTPUT
C     .. dummy arrays
      COMPLEX*16 MEZJ(NKMMAX,NKMMAX,NTMAX,NMVECMAX),
     &           MEZZ(NKMMAX,NKMMAX,NTMAX,NMVECMAX)
C
C      ITERMDIR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ..
C     .. External Subroutines ..
      EXTERNAL AMEMAGVEC,CALCCGC,CALCGF,CALCMVEC,CINIT,IKMLIN,RINIT,
     +         SCFCHRDNS,SSITE,ZCOPY,ZGEMM
C
      DATA ICALL / 0 /
C
      SAVE ICALL,IKM1LIN,IKM2LIN,GDIA,GMDIA,GOFF,LOPT,NLQ,NKMQ,
     &     IQAT,IREL,BCOR,BCORS,QEL,NAT,CONC,TXTT,IMT,SHFTEF,
     &     NVALTOT,NKM,IHYPER,IPRINT,IT,IQ,NL,NT,NUCLEUS,CGC,
     &     IWRREGWF,IWRIRRWF,CALCINT,GETIRRSOL,NFILCBWF,PI,SQPI,
     &     AMEMVEC,IMKMTAB,IKMLLIM1,IKMLLIM2,FACT,SPLITSS,
     &     TXTL,IGRID,IEPATH,NEPATH
C
      ICALL = ICALL + 1
C
C=======================================================================
C       initialise relativistic and dummy variables and SAVE them
C=======================================================================
      IF ( ICALL.EQ.1 ) THEN
C
         IF ( LMAXD.GT.NLMAX-1) THEN
            WRITE(6,*) ' LMAXD = ',LMAXD, ' > NLMAX-1 = ',NLMAX - 1
            STOP  ' Increase NLMAX in < DRVRHO > '
         END IF
C
         IF ( IRMD.GT.NRMAX ) THEN
            WRITE(6,*) ' IRMD = ',IRMD, ' > NRMAX = ',NRMAX
            WRITE(6,*) ' Increase NRMAX in < sprkkr_rmesh.dim > '
            STOP ' In < DRVRHO > '
         END IF
C
         IF ( NMVECMAX.GT.NMVECMAXD ) THEN
            WRITE (6,*) ' NMVECMAX = ',NMVECMAX,' > NMVECMAXD ',
     &                  NMVECMAXD
            WRITE (6,*) ' Increase NVECMAXD in < DRVRHO > ',
     &                  'or reduce NMVECMAX in < main1c > '
            STOP ' In < DRVRHO > '
         END IF
C
         IPRINT = 0
         NL = LMAXD + 1
C
         DO I = 1,NMUEMAX
            LTAB(I) = I/2
            IF( 2*LTAB(I).EQ.I ) THEN
               LBTAB(I)  = LTAB(I) - 1
               KAPTAB(I) = LTAB(I)
            ELSE
               LBTAB(I)  =  LTAB(I) + 1
               KAPTAB(I) = -LTAB(I) - 1
            END IF
            NMUETAB(I) = 2*ABS(KAPTAB(I))
         END DO
C          
         DO IL = 1,NLMAX
            MMAX = 2*IL
            DO IMUE = 1,MMAX
               IF ( (IMUE.EQ.1) .OR. (IMUE.EQ.MMAX) ) THEN
                  NSOLLM(IL,IMUE) = 1
               ELSE
                  NSOLLM(IL,IMUE) = 2
               END IF
            END DO
         END DO
C     
         CALL IKMLIN(IPRINT,NSOLLM,IKM1LIN,IKM2LIN,NLMAX,NMUEMAX,
     &               LINMAX,NLMAX)
C
         CALL CALCCGC(LTAB,KAPTAB,NMUETAB,CGC,NKMAX,NMUEMAX,NKMPMAX)
C
         CALL CALCGF(NKMAX,CGC,GDIA,GMDIA,GOFF,GMOFF,FDIA,FMDIA,
     &               FOFF,FMOFF,LTAB,LBTAB,KAPTAB,NMUETAB,
     &               NMUEMAX,NKMMAX,NKMPMAX)
C     
         DO IT = 1,NTMAX
            BCOR(IT) = 0D0
            BCORS(IT) = 0D0
            QEL(IT) = 0D0
            NAT(IT) = 1
            CONC(IT) = 1D0
            TXTT(IT) = '    '
            IMT(IT) = 1
            LOPT(IT) = -1       ! this should change for Brooks' OP
         END DO
C          
         DO IQ = 1,NQMAX
            NLQ(IQ) = NL
            NKMQ(IQ) = LMMAXD
            NLINQ(IQ) = 2*NLQ(IQ)*(2*NLQ(IQ)-1)
            IQAT(IQ,1) = 1
         END DO
C     
         IREL = 3
         SHFTEF = 0D0
         EFERMI = 0D0
         NVALTOT = 0
         PI = 4D0*ATAN(1D0)
         SQPI = SQRT(PI)
C          
         NKM = LMMAXD
         IHYPER = 0
         IT = 1
         NT = 1
         IQ = 1
         NUCLEUS = 0
C     
         IWRREGWF = 1
         IWRIRRWF = 1
         CALCINT = .TRUE.
         GETIRRSOL = .TRUE.
         NFILCBWF = 87
C     Length in Bytes
         IOL = 8*4 + 3 + (16*4*NRMAX)
         OPEN (NFILCBWF,STATUS='SCRATCH',FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=IOL)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ITERMDIR
C
         IF ( ITERMVDIR ) THEN
            SPLITSS = .FALSE.
            FACT(0) = 1.0D0
            DO I=1,100
               FACT(I) = FACT(I-1)*DBLE(I)
            END DO
C     
            CALL AMEMAGVEC(IREL,IPRINT+1,NKM,AMEMVEC,IKMLLIM1,IKMLLIM2,
     &                     IMKMTAB,CGC,NLMAX,NKMMAX,NKMPMAX,NMVECMAX)
C
            DO I = 0,NLMAX
               TXTL(I) = ' '
            END DO
            IGRID(1) = 5
            IEPATH = 1 
            NEPATH = 1 
         END IF
C
C      ITERMDIR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END IF                    ! ICALL.EQ.1
C=======================================================================
C     
      CALL SSITE(IWRREGWF,IWRIRRWF,NFILCBWF,CALCINT,GETIRRSOL,SOCTL,CTL,
     &     ERYD,P,IHYPER,IPRINT,IKM1LIN,IKM2LIN,NLQ,NKMQ,NLINQ,NT,
     &     NKM,IQAT,TSST,MSST,TSSTLIN,DZZ,DZJ,SZZ,SZJ,OZZ,OZJ,BZZ,
     &     BZJ,QZZ,QZJ,TZZ,TZJ,VT,BT,AT,ZAT,NUCLEUS,R,DRDI,R2DRDI,
     &     JWS,IMT,AMEOPC,AMEOPO,LOPT,SOLVER,CGC,OZZS,OZJS,NLMAX,NQMAX,
     &     LINMAX,NRMAX,NMMAX,NTMAX,NKMMAX,NKMPMAX,NLAMAX)
C     
C-----------------------------------------------------------------------
C     get charge density
C-----------------------------------------------------------------------
C     
      NETAB = IECURR + 1
      IEC = IECURR
C
C Loop over all qdos points specified in qvec.dat
      DO 200 IPOINT = 1,NQDOS                                        ! qdos ruess 
C                                                                    ! qdos ruess 
C Read in Green function; remember that for the rel. case, nspin = 1 ! qdos ruess 
c (without qdos, IPOINT=NQDOS=1)                                     ! qdos ruess 
      IREC = IPOINT + NQDOS * (IECURR-1) +  NQDOS * IELAST * (I1-1)  ! qdos ruess 
      if (t_tgmat%gmat_to_file) then
         READ(69,REC=IREC) GMAT0                                        ! qdos ruess 
      else
         GMAT0(:,:) = t_tgmat%gmat(:,:,irec)
      end if
      GMATLL(:,:,IECURR) = GMAT0(:,:)                                ! qdos ruess 
C                                                                    ! qdos ruess 
C
C-------- GET TAU MATRIX ------------------------
C         TAUT = t G t + t
C
C ---> taut = t
C
      DO J = 1,NKM
         CALL ZCOPY(NKM,TSST(1,J,IT),1,TAUT(1,J,IT),1)
      END DO
C
C ---> w1 = G * t
C
      CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,GMATLL(1,1,IECURR),
     &           LMMAXD,TSST(1,1,IT),NKMMAX,CZERO,W1,LMMAXD)
C
C ---> taut = t * G * t + t = t * w1 + taut
C
      CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,TSST(1,1,IT),NKMMAX,
     &           W1,LMMAXD,CONE,TAUT(1,1,IT),NKMMAX)
C
C ---> store taut in linear array tautlin
C     
      DO LIN = 1,NLINQ(IQ)
         TAUTLIN(LIN,IT) = TAUT(IKM1LIN(LIN),IKM2LIN(LIN),IT)
      END DO
C     
      CALL RINIT(NRMAX,RHOCHR(1,IT))
      CALL RINIT(NRMAX,RHOSPN(1,IT))
      CALL RINIT(NRMAX,RHOORB(1,IT))
      CALL CINIT(NLMAX,OMTL0(1,IT))
      DO LIN = 1,2
         CALL CINIT(NLMAX,OMTLS0(1,IT,LIN))
      END DO
C
      CALL SCFCHRDNS(NFILCBWF,R2DRDI,JWS,IMT,SHFTEF,TOTDOS,MUESPN,
     &               MUEORB,IREL,IPRINT,NT,NL,NKM,ERYD,WE,EFERMI,IEC,
     &               NETAB,DOS,SMT,OMT,HFF,DOSI,SMTI,OMTI,HFFI,DOSM,
     &               DOSL0,DOSINT,SMTM,SMTL0,SMTINT,OMTM,OMTL0,OMTINT,
     &               HFFM,HFFL0,HFFINT,BCOR,BCORS,DZZ,DZJ,SZZ,SZJ,OZZ,
     &               OZJ,BZZ,BZJ,OZZS,OZJS,OMTLS0,TAUTLIN,NVALTOT,TXTT,
     &               CONC,NAT,RHOCHR,RHOSPN,RHOORB,QEL,GDIA,GMDIA,GOFF,
     &               NTMAX,NLMAX,NMUEMAX,LINMAX,NRMAX,NMMAX,NKMMAX,
     &               EBAND,EBANDT)
C     
      DO I = 1,ISHIFT
         RHO2NS(I,1,1) = DZERO
         RHO2NS(I,1,2) = DZERO
         RHOTBORB(I) = DZERO
      END DO
C
      DO I = 1,JWS(IT)
         IP = I + ISHIFT
         RHO2NS(IP,1,1) = RHO2NS(IP,1,1) 
     &        - 0.5D0 * SQPI * RHOCHR(I,IT) * (R(I,1)**2)
         RHO2NS(IP,1,2) = RHO2NS(IP,1,2) 
     &        - 0.5D0 * SQPI * RHOSPN(I,IT) * (R(I,1)**2)
         RHOTBORB(IP) = RHOTBORB(IP)
     &        - 0.5D0 * SQPI * RHOORB(I,IT) * (R(I,1)**2)
      END DO
C
      DO IL = 1,NL
         DEN(IL-1,IECURR+IEMXD) = 
     &        -0.5D0 * (DOSL0(IL,IT)+SMTL0(IL,IT)) * PI
C
         DEN(IL-1,IECURR) = 
     &        -0.5D0 * (DOSL0(IL,IT)-SMTL0(IL,IT)) * PI
C     
         DO I = 1,2
            DMUORB(IL-1,I) = -OMTLS0(IL,IT,I) * PI
         END DO
C
         DMUORB(IL-1,3) = -OMTL0(IL,IT) * PI
      END DO
      DEN(NL,IECURR+IEMXD) = CZERO
      DEN(NL,IECURR) = CZERO
c
c Write out qdos
      DENTOT1 = DCMPLX(0.D0,0.D0)                                         ! qdos ruess 
      DENTOT2 = DCMPLX(0.D0,0.D0)                                         ! qdos ruess 
      DO IL = 0,NL                                                        ! qdos ruess 
         DENTOT1 = DENTOT1 + DEN(IL,IECURR)                               ! qdos ruess 
         DENTOT2 = DENTOT2 + DEN(IL,IECURR+IEMXD)                         ! qdos ruess 
      ENDDO                                                               ! qdos ruess 
      WRITE(31,9000) ERYD,QVEC(1,IPOINT),QVEC(2,IPOINT),QVEC(3,IPOINT),   ! qdos ruess 
     &     -DIMAG(DENTOT1)/PI,(-DIMAG(DEN(IL,IECURR))/PI,IL=0,LMAXD+1),   ! qdos ruess 
     &     -DREAL(DENTOT1)/PI,(-DREAL(DEN(IL,IECURR))/PI,IL=0,LMAXD+1)    ! qdos ruess 
      WRITE(32,9000) ERYD,QVEC(1,IPOINT),QVEC(2,IPOINT),QVEC(3,IPOINT),   ! qdos ruess 
     &-DIMAG(DENTOT2)/PI,(-DIMAG(DEN(IL,IECURR+IEMXD))/PI,IL=0,LMAXD+1),  ! qdos ruess 
     &-DREAL(DENTOT2)/PI,(-DREAL(DEN(IL,IECURR+IEMXD))/PI,IL=0,LMAXD+1)   ! qdos ruess 
 9000 FORMAT(5F10.6,40E16.8)
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ITERMDIR
C
      IF ( ITERMVDIR ) THEN 
C
         QMPHILOC(IQ) = QMPHI
         QMTETLOC(IQ) = QMTET
C
         CALL CINIT(NLMAX*NTMAX*3*NMVECMAX,MVEVDL0)
         CALL CINIT(NLMAX*NTMAX*3*NMVECMAX,BMVEVDL0)
         CALL CINIT(NLMAX*NTMAX*3*NMVECMAX,MVEVIL1)
         CALL CINIT(NLMAX*NTMAX*3*NMVECMAX,BMVEVIL1)
C
         CALL CALCMVEC(TXTL,NFILCBWF,SHFTEF,SPLITSS,IEPATH,NEPATH,IREL,
     &                 IPRINT,NT,NL,MEZZ,MEZJ,TAUT,TSST,IQAT,NKMQ,NKM,
     &                 IEC,NETAB,IGRID(IEPATH),WE,TXTT,FACT,MVEVDL0,
     &                 MVEVIL1,BMVEVDL0,BMVEVIL1,MVPHI,MVTET,MVGAM,
     &                 QMTETLOC,QMPHILOC,R2DRDI,JWS,IMT,AMEMVEC,
     &                 IKMLLIM1,IKMLLIM2,IMKMTAB,NTMAX,NLMAX,NMUEMAX,
     &                 NQMAX,NKMMAX,NMMAX,NMVECMAX,NRMAX)
C     
         DO I=1,NMVECMAX
            DO J=1,3
               DO IL=1,NL
                  MVEVIL(IL-1,J,I) = MVEVIL(IL-1,J,I) 
     &                               - MVEVDL0(IL,IT,J,I) * PI * WE
               END DO
            END DO
         END DO
      END IF
C
C      ITERMDIR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
 200  ENDDO ! IPOINT = 1,NQDOS
C
      IF ( (IECURR.NE.IELAST) .OR. (.NOT.LDORHOEF) ) RETURN
C
C ======================================================================
C     get the charge at the Fermi energy (IELAST)
C     call SCFCHRDNS with the energy weight CONE --> not overwrite WE
C
      CALL RINIT(NRMAX,RHOCHR(1,IT))
      CALL RINIT(NRMAX,RHOSPN(1,IT))
C
      CALL SCFCHRDNS(NFILCBWF,R2DRDI,JWS,IMT,SHFTEF,TOTDOS,MUESPN,
     &               MUEORB,IREL,IPRINT,NT,NL,NKM,ERYD,CONE,EFERMI,IEC,
     &               NETAB,DOS,SMT,OMT,HFF,DOSI,SMTI,OMTI,HFFI,DOSM,
     &               DOSL0,DOSINT,SMTM,SMTL0,SMTINT,OMTM,OMTL0,OMTINT,
     &               HFFM,HFFL0,HFFINT,BCOR,BCORS,DZZ,DZJ,SZZ,SZJ,OZZ,
     &               OZJ,BZZ,BZJ,OZZS,OZJS,OMTLS0,TAUTLIN,NVALTOT,TXTT,
     &               CONC,NAT,RHOCHR,RHOSPN,RHOORB,QEL,GDIA,GMDIA,GOFF,
     &               NTMAX,NLMAX,NMUEMAX,LINMAX,NRMAX,NMMAX,NKMMAX,
     &               EBAND,EBANDT)
C
      DO I = 1,ISHIFT
         R2NEF(I,1,1) = DZERO
         R2NEF(I,1,2) = DZERO
      END DO
C
      DO I = 1,JWS(IT)
         IP = I + ISHIFT
         R2NEF(IP,1,1) = 
     &        - 0.5D0 * SQPI * RHOCHR(I,IT) * (R(I,1)**2)
         R2NEF(IP,1,2) = 
     &        - 0.5D0 * SQPI * RHOSPN(I,IT) * (R(I,1)**2)
      END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ITERMDIR
C
      IF( ITERMVDIR ) THEN 
         DO I=1,NMVECMAX
            DO J=1,3
               DO IL=1,NL
                  MVEVILEF(IL-1,J,I) = -MVEVDL0(IL,IT,J,I) * PI
               END DO
            END DO
         END DO
      END IF
C     
C      ITERMDIR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C ======================================================================
C
      END
