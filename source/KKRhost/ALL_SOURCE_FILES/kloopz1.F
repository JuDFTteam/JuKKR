c*==kloopz1.f    processed by SPAG 6.05Rc at 16:33 on 10 Oct 2001
c 13.10.95 *************************************************************
      SUBROUTINE KLOOPZ1_QDOS(ERYD,GMATLL,INS,ALAT,IE,IGF,
     &                   NSHELL,NAEZ,NOFKS,VOLBZ,BZKP,VOLCUB,CLS,NACLS,
     &                   NACLSMAX,RR,RBASIS,EZOA,ATOM,RCLS,ICC,GINP,   
     &                   IDECI,LEFTTINVLL,RIGHTTINVLL,VACFLAG,NLBASIS,
     &                   NRBASIS,FACTL,NATOMIMP,NSYMAT,DSYMLL,RATOM,
     &                   RROT,NSH1,NSH2,IJTABSYM,IJTABSH,ICHECK,INVMOD,
     &                   REFPOT,TREFLL,TSST,MSST,CFCTOR,          
     &                   CFCTORINV,CREL,RC,RREL,SRREL,IRREL,NRREL,DROTQ,
     &                   SYMUNITARY,KMROT,NATYP,NCPA,ICPA,ITCPAMAX,
     &                   CPATOL,NOQ,IQAT,ITOQ,CONC,IPRINT,ICPAFLAG,
     &                   ISPIN,NSPIN,
     &                   TQDOS,IQDOSRUN,  !qdos ruess
     &                   DTREFLL,DTMATLL,DGINP,LLY_GRTR,TRACET,LLY)! LLY Lloyd
c **********************************************************************
      use mod_types, only: t_inc
      use mod_mympi, only: myrank, master
      IMPLICIT NONE
C     .. Parameters ..
      include 'inc.p'
C
C **********************************************************************
C * For KREL = 1 (relativistic mode)                                   *
C *                                                                    *
C *  NPOTD = 2 * NATYPD                                                *
C *  LMMAXD = 2 * (LMAXD+1)^2                                          *
C *  NSPIND = 1                                                        *
C *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green      *
C *          function, set up in the spin-independent non-relativstic  *
C *          (l,m_l)-representation                                    *
C *                                                                    *
C **********************************************************************
C
      INTEGER NEMBD1
      PARAMETER (NEMBD1=NEMBD+1)
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+KORBIT+1) * (LMAXD+1)**2)
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAXD+1)**2)
      INTEGER LINMAX
      PARAMETER (LINMAX=1)
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER ( CONE  = ( 1.0D0,0.0D0), CZERO = ( 0.0D0,0.0D0 ) )
      DOUBLE PRECISION TOLMSSQ
      PARAMETER ( TOLMSSQ=1.0D-6 )
C     ..
C     .. Scalar Arguments ..
      INTEGER ICC,ICPAFLAG,IDECI,IE,IGF,INS,INVMOD,IPRINT,ITCPAMAX,
     +        KMROT,NAEZ,NATOMIMP,NATYP,NCPA,NLBASIS,NRBASIS,NOFKS,
     +        NSYMAT,ISPIN,NSPIN,NACLSMAX,
     +        IQDOSRUN            !qdos ruess: counts qdos run
     &       ,LLY ! LLY <> 0 => use Lloyd formula
      DOUBLE PRECISION ALAT,VOLBZ,CPATOL
      DOUBLE COMPLEX ERYD,CFCTOR,CFCTORINV
      DOUBLE COMPLEX LLY_GRTR ! Trace Eq.5.38 PhD Thiess (k-integrated)! LLY Lloyd 
      DOUBLE COMPLEX TRACET   ! Tr[ (t-tref)^-1 d(t-tref)/dE ]
C     ..
C     .. Array Arguments ..
      INTEGER ATOM(NACLSD,*),CLS(*),EZOA(NACLSD,*),NACLS(*),
     +        NSHELL(0:NSHELD),NSH1(*),NSH2(*)
      INTEGER ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD),
     +        IJTABSYM(*),IJTABSH(*),REFPOT(*)
      INTEGER IRREL(2,2,LMMAXD),NRREL(2,LMMAXD)
      INTEGER ICPA(NAEZD),NOQ(NAEZD),IQAT(NATYPD),
     +        ITOQ(NATYPD,NAEZD)
      DOUBLE PRECISION BZKP(3,KPOIBZ),CONC(NATYPD),VOLCUB(KPOIBZ)
      DOUBLE PRECISION RATOM(3,*),RBASIS(3,*),RCLS(3,NACLSD,*),
     +                 RR(3,0:NRD),RROT(48,3,*)
      DOUBLE COMPLEX GMATLL(LMMAXD,LMMAXD,*),
     +               GINP(LMGF0D*NACLSMAX,LMGF0D,*),
     &               DGINP(LMGF0D*NACLSMAX,LMGF0D,*) ! LLY Energy deriv. of GINP (=Gref)
      DOUBLE COMPLEX LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPIN),
     &               RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPIN)
      DOUBLE COMPLEX TSST(LMMAXD,LMMAXD,NATYPD),
     +               MSST(LMMAXD,LMMAXD,NATYPD),
     +               TREFLL(LMMAXD,LMMAXD,NREFD),
     &               DTREFLL(LMMAXD,LMMAXD,NREFD), ! LLY dtref/dE
     &               DTMATLL(LMMAXD,LMMAXD,NAEZD)  ! LLY  dt/dE (should be av.-tmatrix in CPA)
      DOUBLE COMPLEX FACTL(LMMAXD,LMMAXD)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,*),DROTQ(LMMAXD,LMMAXD,NAEZD)
      DOUBLE COMPLEX CREL(LMMAXD,LMMAXD),RC(LMMAXD,LMMAXD),
     &               RREL(LMMAXD,LMMAXD),SRREL(2,2,LMMAXD)
      DOUBLE COMPLEX TQDOS(LMMAXD,LMMAXD,NAEZD)  ! qdos : Read-in inverse t-matrix
      LOGICAL VACFLAG(2),SYMUNITARY(*)
C     ..
C     .. Local Scalars ..
      INTEGER IH,LM1,LM2,NS,NSDIA,ICALL,IREC
      INTEGER IQ,JQ,IT,I,J,IQTAU,ICPASTART,ITCPA,IU,NSMAX
      DOUBLE PRECISION CPAERRL,CPAERR,CPACORR,CPACHNG
      DOUBLE COMPLEX EZ,CNSYMAT,TAUVBZ
      LOGICAL LDIA
      CHARACTER*4 STR4
      CHARACTER*10 STR10
C     ..
C     .. Local Arrays ..
      INTEGER IKM1LIN(LINMAX),IKM2LIN(LINMAX)
      INTEGER NKMQ(NAEZD),NLINQ(NAEZD)
      INTEGER ISUMT(NSYMAXD,NATYPD),ISUMQ(NSYMAXD,NAEZD),ISUMG(NSYMAXD)
      DOUBLE COMPLEX GLL(LMMAXD,LMMAXD)
      DOUBLE COMPLEX TAUTLIN(LINMAX,NATYPD)                
C     ..
C     .. effective (site-dependent) Delta_t^(-1) matrix ..
      DOUBLE COMPLEX MSSQ(LMMAXD,LMMAXD,NAEZD)
      DOUBLE COMPLEX XC(LMMAXD,LMMAXD),W1(LMMAXD,LMMAXD),
     &               W2(LMMAXD,LMMAXD)
C     .. 
CF90-------------------------------------------------------------
      DOUBLE COMPLEX GS(:,:,:,:),TAUDELQ(:,:,:),TAUDELT(:,:,:)
      DOUBLE COMPLEX DMSSQ(:,:,:)
      ALLOCATABLE GS,TAUDELQ,TAUDELT,DMSSQ
CF90-------------------------------------------------------------
#ifdef CPP_MPI
      INTEGER :: irec0
#endif
C
CF90CF90-------------------------------------------------------------
CF90      DOUBLE COMPLEX GS(:,:,:,:),TAUDELQ(:,:,:),TAUDELT(:,:,:)
CF90      DOUBLE COMPLEX DMSSQ(:,:,:)
CF90      ALLOCATABLE GS,TAUDELQ,TAUDELT,DMSSQ
CF90CF90-------------------------------------------------------------
CF77CF77-------------------------------------------------------------
CF77      DOUBLE COMPLEX GS(LMMAXD,LMMAXD,NSYMAXD,NSHELD)
CF77      DOUBLE COMPLEX TAUDELQ(LMMAXD,LMMAXD,NSHELD)
CF77      DOUBLE COMPLEX TAUDELT(LMMAXD,LMMAXD,NSHELD)
CF77      DOUBLE COMPLEX DMSSQ(LMMAXD,LMMAXD,NAEZD)
CF77CF77-------------------------------------------------------------
C     ..
C     .. External Functions/Subroutines ..
      LOGICAL TEST,OPT
      EXTERNAL TEST,OPT,KKRMAT01,ROTGLL,ROTATE,SYMETRMAT,
     &         ZCOPY,ZGEMM,ZSCAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
C     .. Data statement ..
      DATA ICALL / 0 /
C     ..
C     .. Save statement ..
      SAVE ICALL,ISUMT,ISUMQ,ISUMG,CNSYMAT
C     ..

C **********************************************************************
      IF ( TEST('flow    ') ) WRITE (6,*) 
     &       '>>> KLOOPZ1: invert delta_t and do Fourier transformation'
      ICALL = ICALL + 1
C
C     Reinitialise the ICALL counter for second run of kloopz1   !qdos ruess
      IF ( (OPT('qdos    ').AND.(IQDOSRUN.EQ.1)) ) ICALL = 1     !qdos ruess
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C --> the arrays ISUM are used in the symmetrisation routine SYMETRMAT 
C     symmetrising single-site : same matrix for each symmetry
C     symmetrising G matrix    : pick G(ISYM) for symmetry ISYM
C
      IF ( ICALL.EQ.1 ) THEN
         CNSYMAT = CONE/DBLE(NSYMAT)
         DO IU = 1,NSYMAXD
            DO IT = 1,NATYPD
               ISUMT(IU,IT) = IT
            END DO
            DO IQ = 1,NAEZD
               ISUMQ(IU,IQ) = IQ
            END DO
            ISUMG(IU) = IU
         END DO
      END IF
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C ======================================================================
C
C     TSST in the LOCAL frame is used to set up
C          MSST = (TSST-TREF)^(-1) in the LOCAL frame to be used
C                                  in <CPAMILLSX> and < PROJTAU > below
C          
C          MSSQ = the inverse of the effective (on-site) Delta_t matrix
C                 in the GLOBAL frame; 
C                 the Average T-matrix Approximation (ATA) (ICPASTART=1)
C                 is used 
C


      ICPASTART = 1
      CALL MSSINIT(NCPA,ICPASTART,TSST,MSST,MSSQ,TREFLL,DROTQ,
     &             REFPOT,IQAT,ITOQ,NOQ,CONC,
     &             KMROT,NATYP,NAEZ,LMMAXD)  ! nrefd was taken out of calling list 1.2.2012



C    VIRTUAL ATOMS:
C    Be careful! in case of  OPT('VIRATOMS')==1 MSSQ is the Tmatrix 
C    not the inverse T-matrix!!


C
C ----------------------------------------------------------------------
C    Output now:
C    MSSQ is the (Delta_t)^(-1) matrix in the GLOBAL frame
C         and refers to a site IQ occupied by NOQ(IQ) atoms having
C         the occupancies CONC(ITOQ(1..NOQ(IQ))
C
C    TSST is the t-matrix in the LOCAL frame 
C         for each of the NATYP atoms
C    MSST is the Delta_t^(-1) matrix in the LOCAL frame 
C         for each of the NATYP atoms
C ----------------------------------------------------------------------
C
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         BEGIN CPA - LOOP  (if required)
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  ikm1lin,ikm2lin,nlinq --> dummy settings for <PROJTAU>
C
      IKM1LIN(1) = 1
      IKM2LIN(1) = 1
      DO IQ = 1,NAEZ
        NKMQ(IQ) = LMMAXD
        NLINQ(IQ) = 1
      END DO
C
      CPAERRL = 1.0D+6
      ITCPA = 0
      EZ = ERYD
C
      IF ( NCPA.NE.0 ) THEN
CF90-------------------------------------------------------------
         ALLOCATE (DMSSQ(LMMAXD,LMMAXD,NAEZD),STAT=LM1)
         IF ( LM1.NE.0 ) STOP '      ERROR: <kloopz1> allocate DMSSQ'
         DMSSQ(:,:,:) = (0.d0,0.d0)
CF90-------------------------------------------------------------
         CALL CINIT(LMMAXD*LMMAXD*NAEZD,DMSSQ)
      END IF
C
C
C --->  < GLL2K >  incorporated now
C
C
      TAUVBZ = 1.D0/VOLBZ
C
C=======================================================================
C
C --> convert inverted delta_t-matrices to p.u. 
C
      DO I = 1,NAEZ

         IF ( .not. OPT('VIRATOMS') ) THEN
            CALL ZSCAL(LMMAXD*LMMAXD,CFCTOR,MSSQ(1,1,I),1)
         ELSE
            CALL ZSCAL(LMMAXD*LMMAXD,CFCTORINV,MSSQ(1,1,I),1)
         END IF   !( .not. OPT('VIRATOMS') ) 
C
        IF (.NOT.OPT('NEWSOSOL') ) THEN
         IF (KMROT.EQ.0) THEN
            DO LM2 = 1,LMMAXD
               DO LM1 = 1,LM2
                  MSSQ(LM1,LM2,I) = 0.5D0 * 
     &                  ( MSSQ(LM1,LM2,I) + MSSQ(LM2,LM1,I) )
                  MSSQ(LM2,LM1,I) = MSSQ(LM1,LM2,I) 
               END DO
            END DO
         END IF
       ENDIF
      END DO
C
      DO I = 1,NATYP
         CALL ZSCAL(LMMAXD*LMMAXD,CFCTOR,MSST(1,1,I),1)
      END DO
C ======================================================================
C
C --> symmetrise the delta_t^(-1) matrices 
C
C ----------------------------------------------------------------------
      IF (.NOT.OPT('NEWSOSOL')) THEN
      IF ( ( KREL.EQ.1 ).OR.( INS.NE.0 ) ) THEN
         DO IQ = 1,NAEZ
            CALL SYMETRMAT(NSYMAT,CNSYMAT,DSYMLL,SYMUNITARY,MSSQ,
     &                     ISUMQ(1,IQ),MSSQ(1,1,IQ),LMMAXD)
            IF ( KMROT.EQ.0 ) THEN
               DO I = 1,NOQ(IQ)
                  IT = ITOQ(I,IQ)
                  CALL SYMETRMAT(NSYMAT,CNSYMAT,DSYMLL,SYMUNITARY,MSST,
     &                           ISUMT(1,IT),MSST(1,1,IT),LMMAXD)
               END DO
            END IF
         END DO
      END IF
      ENDIF
C
      EZ = EZ*CFCTOR*CFCTOR
      NSDIA = MAX(NAEZ,NATYP)
C
CF90-------------------------------------------------------------
      ALLOCATE (TAUDELQ(LMMAXD,LMMAXD,NSHELL(0)),STAT=LM1)
      IF ( LM1.NE.0 ) STOP '      ERROR: <kloopz1> allocate TAUDELQ'
      TAUDELQ(:,:,:) = (0.d0,0.d0)
CF90-------------------------------------------------------------
C
 100  CONTINUE
      ITCPA = ITCPA + 1
C
CF90-------------------------------------------------------------
      ALLOCATE (GS(LMMAXD,LMMAXD,NSYMAXD,NSHELL(0)),STAT=LM1)
      IF ( LM1.NE.0 ) STOP '      ERROR: <kloopz1> allocate GS'
      GS(:,:,:,:) = (0.d0,0.d0)
CF90-------------------------------------------------------------
C
C     copy read-in cpa t-matrix but only after fort.37 was created in first run !qdos ruess
      IF (OPT('readcpa ').OR.(OPT('qdos    ').AND.                              !qdos ruess
     &                               (IQDOSRUN.EQ.1))) THEN                     !qdos ruess
         MSSQ(:,:,:) = TQDOS(:,:,:) ! lmmaxd,lmmaxd,naezd                       !qdos ruess
      ENDIF
C
      CALL KKRMAT01(BZKP,NOFKS,GS,VOLCUB,MSSQ,RROT,NSHELL(0),NSDIA,
     &              ALAT,NSYMAT,NAEZ,CLS,NACLS,NACLSMAX,RR,EZOA,ATOM,
     &              NSH1,NSH2,GINP,RBASIS,RCLS,LEFTTINVLL(1,1,1,ISPIN),
     &              RIGHTTINVLL(1,1,1,ISPIN),VACFLAG,
     &              NLBASIS,NRBASIS,FACTL,ICHECK,INVMOD,IDECI,
     &              SRREL,IRREL,NRREL,
     &          DTREFLL,DTMATLL,DGINP,REFPOT,LLY_GRTR,TRACET,CFCTOR,LLY) ! LLY
C
      NSMAX = NSHELL(0) 
        
C
C=========================================================== NS=1,NSMAX
C
C For qdos calculation, do not symmetrize the GF matrix, but call routine 
C symetrmat anyway because of factor tauvbz:
      DO NS = 1,NSMAX
C
C --> symmetrise GS, get GLL as sum over all symmetry wedges of GS
C     (isumg(i) = i, see above)
C
         CALL SYMETRMAT(NSYMAT,TAUVBZ,DSYMLL,SYMUNITARY,
     &                  GS(1,1,1,NS),ISUMG,GLL,LMMAXD)
C
         IF (NS.LE.NATYP) THEN
            IQTAU = IQAT(NS)
         ELSE 
            IQTAU = NS
         END IF
C
         TAUDELQ(1:LMMAXD,1:LMMAXD,IQTAU) = -GLL(1:LMMAXD,1:LMMAXD)
      END DO
C=========================================================== NS=1,NSMAX
CF90-------------------------------------------------------------
      DEALLOCATE (GS,STAT=LM1)
      IF ( LM1.NE.0 ) STOP '      ERROR: <kloopz1> deallocate GS'
CF90-------------------------------------------------------------

C=== CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF ( NCPA.GT.0 ) THEN
C
C  --> do one CPA iteration, the output is a new MSSQ = (Delta_t)^(-1)
C      in the GLOBAL frame
C
         ICPAFLAG = 0
         IF (OPT('readcpa ').OR.(OPT('qdos    ')
     &                          .AND.IQDOSRUN.GT.0)) THEN   ! copy read-in cpa t-matrix in second run
            CALL CPAMILLSX(ITCPA,CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,
     &                  NAEZ,NKMQ,NOQ,ITOQ,CONC,MSSQ,MSST,TAUDELQ,DMSSQ,
     &                  KMROT,DROTQ,NATYPD,NAEZD,LMMAXD)
            MSSQ(:,:,:) = TQDOS(:,:,:) ! lmmaxd,lmmaxd,naezd   ! qdos 
            ITCPAMAX = 0
            CPAERR = 0.D0

         ELSE
            CALL CPAMILLSX(ITCPA,CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,
     &                  NAEZ,NKMQ,NOQ,ITOQ,CONC,MSSQ,MSST,TAUDELQ,DMSSQ,
     &                  KMROT,DROTQ,NATYPD,NAEZD,LMMAXD)
C
C -->  symmetrise m-CPA
C
            IF (.NOT.OPT('NEWSOSOL')) THEN
               DO IQ = 1,NAEZ
                  IF ( ICPA(IQ).NE.0 ) 
     &                 CALL SYMETRMAT(NSYMAT,CNSYMAT,DSYMLL,
     &                 SYMUNITARY,MSSQ,ISUMQ(1,IQ),MSSQ(1,1,IQ),LMMAXD)
               END DO
            ENDIF
C
            IF ( IPRINT.GE.1 ) WRITE (6,99004) CPAERR,CPACORR,CPACHNG
C
            IF ( CPAERR.LE.CPATOL ) THEN
               IF ( IPRINT.GT.0 ) WRITE (6,99001) ITCPA,CPAERR,CPACORR,
     &                                            CPACHNG
            ELSE IF ( ITCPA.GT.ITCPAMAX ) THEN
               WRITE (6,99002) ITCPA,CPAERR,CPACORR,CPACHNG
               ICPAFLAG = 1
            ELSE IF ( CPAERR.GT.20*CPAERRL ) THEN
               WRITE (6,99003) ITCPA
               WRITE (6,99004) CPAERR,CPACORR,CPACHNG
               ICPAFLAG = 2
            ELSE
C
C  --> go to the next CPA iteration if not converged
C
               CPAERRL = CPAERR
               GOTO 100
            END IF              !   ( CPAERR.LE.CPATOL )
         END IF                 !  (OPT('readcpa ').OR.OPT('qdos    ')) 
CF90-------------------------------------------------------------
         DEALLOCATE (DMSSQ,STAT=LM1)
         IF ( LM1.NE.0 ) STOP '      ERROR: <kloopz1> deallocate DMSSQ'
CF90-------------------------------------------------------------
C
      END IF
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                         END CPA - LOOP
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C=======================================================================
C  The inverse of the Delta_t(CPA)-matrix is calculated, now write
C  this on a file in case of decimation output. 
C  attention: new format Dec. 2004
C
C  In first qdos run write out t matrix which is read in for the calculation for every k point
      IF ( OPT('deci-out').OR.(IQDOSRUN.EQ.0) ) THEN  ! qdos ruess
         DO IH = 1,NAEZ
#ifdef CPP_MPI
            irec0 = LMMAXD**2*(ih-1)+LMMAXD**2*naez*(ie-1)+
     &              LMMAXD**2*t_inc%IELAST*naez*(ispin-1)
            IF (OPT('deci-out'))  WRITE (37,99005) IE,ERYD,IH
#else
            WRITE (37,99005) IE,ERYD,IH
#endif
            DO LM1 = 1,LMMAXD
               DO LM2 = 1,LMMAXD
                  IF ( LM1.EQ.LM2 ) THEN
#ifdef CPP_MPI
                     IF (.not.OPT('deci-out')) then
                     irec = irec0 + LM2 + LMMAXD*(LM1-1)
                     WRITE(37,rec=irec) MSSQ(LM1,LM2,IH)*CFCTORINV
                     else
                     WRITE(37,99006) LM1,LM2,MSSQ(LM1,LM2,IH)*CFCTORINV
                     end if
#else
                     WRITE(37,99006) LM1,LM2,MSSQ(LM1,LM2,IH)*CFCTORINV
#endif
                  ELSE
#ifdef CPP_MPI
                     irec = irec0 + LM2 + LMMAXD*(LM1-1)
                     IF ( CDABS(MSSQ(LM1,LM2,IH)/MSSQ(LM1,LM1,IH))
     &                    .GT.TOLMSSQ ) THEN
                         if (.not.OPT('deci-out')) then
                           WRITE(37,rec=irec) MSSQ(LM1,LM2,IH)*CFCTORINV
                         else
                           WRITE(37,99006) 
     &                                LM1,LM2,MSSQ(LM1,LM2,IH)*CFCTORINV
                         end if
                     ELSE
                         if (.not.OPT('deci-out')) then
                            WRITE(37,rec=irec) (0.d0, 0.d0) 
                         end if
                     END IF
#else
                     IF ( CDABS(MSSQ(LM1,LM2,IH)/MSSQ(LM1,LM1,IH))
     &                    .GT.TOLMSSQ ) WRITE(37,99006) 
     &                    LM1,LM2,MSSQ(LM1,LM2,IH)*CFCTORINV
#endif
                  END IF
               END DO
            END DO
         END DO
      END IF
C=======================================================================
C
C --> calculate the component-projected site-diagonal 
C     TAU-matrices TAUDELT(IT).
C     - there are NSMAX = NAEZ/NATYP (NCPA=0/1) site-diagonal
C       elements TAUDELQ(1..NSMAX) in the GLOBAL (crystal) frame of 
C       reference
C     - they are NSMAX site-diagonal elements projected on atomic types 
C       in the array TAUDELT(1..NSMAX) in the LOCAL frame of reference
C
CF90-------------------------------------------------------------
      ALLOCATE (TAUDELT(LMMAXD,LMMAXD,NSHELL(0)),STAT=LM1)
      IF ( LM1.NE.0 ) STOP '      ERROR: <kloopz1> allocate TAUDELT'
      TAUDELT(:,:,:) = (0.d0,0.d0)
CF90-------------------------------------------------------------
      IF ( NCPA.GT.0) THEN
         NSMAX = NATYP
         CALL PROJTAU(ICPAFLAG,CPACHNG,KMROT,.FALSE.,.FALSE.,9,CZERO,
     &        NATYP,NAEZ,NKMQ,MSST,MSSQ,NLINQ,IQAT,CONC,TAUDELQ,TAUDELT,
     &        TAUTLIN,IKM1LIN,IKM2LIN,DROTQ,NSHELD,NAEZD,LMMAXD,LINMAX)
      ELSE
         NSMAX = NAEZ
         DO NS=1,NSMAX
            CALL ZCOPY(LMMAXD*LMMAXD,TAUDELQ(1,1,NS),1,
     &           TAUDELT(1,1,NS),1)
            IF ( KMROT.NE.0 ) THEN
               DO J = 1,LMMAXD
                  CALL ZCOPY(LMMAXD,TAUDELT(1,J,NS),1,W1(1,J),1)
               END DO
               CALL ROTATE(W1,'G->L',TAUDELT(1,1,NS),LMMAXD,
     &              DROTQ(1,1,NS),LMMAXD)
            END IF
         END DO
      END IF
C
C -> off-diagonal elements (if present) are in the same array 
C    TAUDELQ in the range (NSMAX+1,..,NSHELL(0)). 
C    The G_ij's are left UNPROJECTED and in the GLOBAL frame of 
C    reference (remember this for further use!), i.e. in case of CPA, 
C    G_ij(CPA) is stored. The component-projected G-elements can
C    be obtained through the projection matrices (see below)
C
      DO NS = NSMAX+1, NSHELL(0)
         CALL ZCOPY(LMMAXD*LMMAXD,TAUDELQ(1,1,NS),1,TAUDELT(1,1,NS),1)
      END DO
C ----------------------------------------------------------------------
C
C -> switch from TAU to G multiplying with MSST on the site-diagonal,
C    system-matrix (1..NSMAX) and with MSSQ on the rest of the elements
C
C ----------------------------------------------------------------------
C
      IF ( TEST('Gmat    ') ) WRITE (6,'(/,4X,70(1H-),/,4X,A,I4)')
     &                   'system G_ii matrix for i = 1,',NSMAX
      DO NS = 1,NSHELL(0)
C
           
         LDIA = (DABS(RATOM(1,NS)**2+RATOM(2,NS)**2+RATOM(3,NS)**2)
     &           .LT.1D-6)
C
C --->  GLL = -TAU
C
         CALL ZCOPY(LMMAXD*LMMAXD,TAUDELT(1,1,NS),1,GLL,1)
         CALL ZSCAL(LMMAXD*LMMAXD,-CONE,GLL,1)
C
C --->  XC = (Delta_t(I))^(-1) * GLL
C
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         IF ( NS.LE.NSMAX ) THEN 
C 
C deal with the site-diagonal GFUN of the system - use MSST(I)
C on both sides of TAU multiplication
C
            DO J = 1,LMMAXD
               CALL ZCOPY(LMMAXD,MSST(1,J,NS),1,W1(1,J),1)
            END DO
            CALL ZCOPY(LMMAXD*LMMAXD,W1,1,W2,1)
         ELSE
C 
C deal with the Gij of the cluster - use MSSQ(IQ) and MSSQ(JQ)
C
            IQ = NSH1(NS)
            JQ = NSH2(NS)
            DO J = 1,LMMAXD
               CALL ZCOPY(LMMAXD,MSSQ(1,J,IQ),1,W1(1,J),1)
            END DO
            DO J = 1,LMMAXD
               CALL ZCOPY(LMMAXD,MSSQ(1,J,JQ),1,W2(1,J),1)
            END DO
         END IF
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::: NS.LT.NSMAX
C
      IF ( .not. OPT('VIRATOMS') ) THEN
       IF ( .not. TEST('testgmat') ) THEN


         CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,W1,
     &              LMMAXD,GLL,LMMAXD,CZERO,XC,LMMAXD)
C
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        IF ( LDIA ) THEN
C
C --->  GLL = - (Delta_t(I))^(-1) 
C                          - (Delta_t(I))^(-1) * GLL * (Delta_t(I))^(-1)
C
           CALL ZCOPY(LMMAXD*LMMAXD,W2,1,GLL,1)
           CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,-CONE,XC,LMMAXD,
     &                W2,LMMAXD,-CONE,GLL,LMMAXD)
        ELSE
C
C --->    GLL =  - (Delta_t(I))^(-1)  * GLL * (Delta_t(J))^(-1)
C
           CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,-CONE,XC,LMMAXD,
     &                W2,LMMAXD,CZERO,GLL,LMMAXD)
        END IF

       END IF !( .not. TEST('testgmat') ) THEN
      END IF !( .not. OPT('VIRATOMS') ) THEN


C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: LDIA
C
C --->  GMATLL = GLL/RFCTOR
C
      GMATLL(1:LMMAXD,1:LMMAXD,NS) = GLL(1:LMMAXD,1:LMMAXD)*CFCTORINV
C
        IF ( ( NS.LE.NSMAX ).AND.( TEST('Gmat    ') ) ) THEN
           WRITE(STR4,'(I4)') NS
           STR10 = '   i ='//STR4(1:4)
           CALL CMATSTR(STR10,10,GMATLL(1,1,NS),LMMAXD,LMMAXD,
     &                  2*KREL+1,2*KREL+1,0,1d-8,6)
        END IF
      END DO
C ----------------------------------------------------------------------
      IF ( TEST('Gmat    ') ) WRITE (6,'(/,4X,70(1H-))')
C ----------------------------------------------------------------------
C --> it calculates the rest of the G n n' matrix from the
C     knowledge of the representative pairs (shells) using the
C     real space symmetries (added 23.2.2000)
C

      IF ( ICC.GT.0 ) then
          IREC = 1+IE+t_inc%IELAST*(ISPIN-1) ! added for mpi run
          CALL ROTGLL(GMATLL,NATOMIMP,IJTABSYM,IJTABSH,
     &                DSYMLL,SYMUNITARY,IGF,RC,CREL,RREL,
     &                KREL,LMMAXD,IREC)
      END IF


!       IF ( OPT('VIRATOMS') ) THEN
!         write(*,*) 'VIRTUAL ATOM OPTION : stop calculation '
!         stop
!       END IF

C ----------------------------------------------------------------------
C
C ----------------------------------------------------------------------
C --> in the case of NCPA.NE.0 and NSHELL(0).GT.NATYP the projection 
C     matrices DMAT and DTIL which are used to get
C                     ij            ij    _ 
C                    G     =  D  * G    * D
C                     ab       a    CPA    b
C   - with a/b the atom of type a/b sitting on site i/j - are calculated
C   and stored for later use.  the allocated work space for 
C   TSST (DMAT) and MSST (DTIL) is used.
C   for an atom having occupancy 1, DMAT/DTIL = unit matrix
C 
      IF (( NCPA.NE.0 ).AND. ( NSHELL(0).GT.NSMAX )) THEN
         DO IT = 1,NATYP
            IQ = IQAT(IT)
            IH = REFPOT(IQ)
C
            IF ( KMROT.NE.0 ) THEN
               CALL ROTATE(TSST(1,1,IT),'L->G',
     &              W1,LMMAXD,DROTQ(1,1,IQ),LMMAXD)
            ELSE
               CALL ZCOPY(LMMAXD*LMMAXD,TSST(1,1,IT),1,W1,1)
            END IF
C
            DO J = 1,LMMAXD
               CALL ZAXPY(LMMAXD,-CONE,TREFLL(1,J,IH),1,W1(1,J),1)
            END DO
C
            CALL GIJDMAT(TAUDELQ(1,1,IQ),W1,MSSQ(1,1,IQ),
     &                   TSST(1,1,IT),MSST(1,1,IT),CFCTORINV,IPRINT,
     &                   IE,IT,KREL,LMMAXD)
         END DO
      END IF
C ----------------------------------------------------------------------
      IF ( TEST('flow    ') ) WRITE (6,*) '<<< KLOOPZ1'

CF90-------------------------------------------------------------
      DEALLOCATE (TAUDELQ,TAUDELT,STAT=LM1)
      IF ( LM1.NE.0 ) STOP '      ERROR: <kloopz1> deallocate TAUs'
CF90-------------------------------------------------------------
C
      RETURN

99001 FORMAT (' CPA converged after',I3,' iterations   ','ERR:',F9.6,
     &        ' CORR:',F9.6,' CHNG:',F9.6)
99002 FORMAT (10X,10('!'),' CPA-cycle  NOT  converged after ',I3,
     &        ' iterations ',10('!'),/,14X,'ERROR ',F12.8,
     &        ' CORRECTION ',F15.8,' CHANGE ',F15.8)
99003 FORMAT (' CPA: ERROR increased by more than ','20*TOL for ITCPA=',
     &        i4,' >>>> iteration stopped ',5X,10('!'))
99004 FORMAT (' CPA: ERROR ',F12.8,'    CORRECTION ',F15.8,
     &        '    CHANGE ',F15.8)
C
99005 FORMAT (/,80('*')/,'ENERGY ',I5,2D16.8,' SITE ',I3)
99006 FORMAT (2I5,1P,2D22.14)
      END
