      PROGRAM MAIN1B
      IMPLICIT NONE
CMPI  INCLUDE 'mpif.h'
      INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
C *          function, set up in the spin-independent non-relativstic *
C *          (l,m_l)-representation                                   *
C *                                                                   *
C *********************************************************************
C     ..
C     .. Parameters ..
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAXD+1)**2)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+KORBIT+1)* (LMAXD+1)**2)
      INTEGER NSPINDD
      PARAMETER (NSPINDD=NSPIND-KORBIT)
C     ..
C parameter nembd1 avoids zero sized arrays.(2.1.01 R.Zeller)
      INTEGER NEMBD1
      PARAMETER (NEMBD1=NEMBD+1)
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      INTEGER NOFGIJD
      PARAMETER (NOFGIJD = NATOMIMPD*NATOMIMPD+1)
      INTEGER MAXMSHD
      PARAMETER (MAXMSHD=30)
      INTEGER LRECGRF,LRECTMT,LRECGREEN
      PARAMETER (LRECGRF=WLENGTH*4*NACLSD*LMGF0D*LMGF0D*NCLSD) ! 4 words = 16 bytes / complex number (in ifort 4; in gfort 16)
c        word/byte distiction moved to subroutine opendafile to be the same for all unformatted files
      PARAMETER (LRECTMT=WLENGTH*4*LMMAXD*LMMAXD)
      PARAMETER (LRECGREEN=WLENGTH*2*NATOMIMPD*LMMAXD*NATOMIMPD*LMMAXD) 
      INTEGER LRECTRA ! LLY Lloyd
      PARAMETER (LRECTRA=WLENGTH*4) ! LLY Lloyd
C     ..
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO=(0.0D0,0.0D0))
      DOUBLE COMPLEX CONE
      PARAMETER (CONE=(1.0D0,0.0D0))
      DOUBLE PRECISION CVLIGHT
      PARAMETER (CVLIGHT = 274.0720442D0 )
C     ..
C     .. Local Scalars ..
      INTEGER NLBASIS,NRBASIS                     
      INTEGER IELAST,INS,LMAX,NATYP,NREF,NSPIN,NSRA
      INTEGER KMROT,INVMOD,ICC,IGF,IPRINT,IDECI
      INTEGER MAXMESH,NAEZ,I,ISPIN,I1,IE,IREC,L,LM1,LM2,L1
      INTEGER NATOMIMP,NSYMAT,NCPA,NCPAFAIL,ICPAFLAG,ITCPAMAX
      INTEGER NMESH,ID,NOFGIJ,NQCALC
      INTEGER ITMPDIR,ILTMP
      INTEGER IQ,NQDOS ! qdos ruess:number of qdos points
      INTEGER IX,ISITE ! qdos ruess
      INTEGER IQDOSRUN ! qdos ruess: counter to organise qdos run
      INTEGER NCLS,NACLSMAX,IC,LRECGRF1
      DOUBLE COMPLEX TREAD ! qdos ruess
C     RFCTOR=A/(2*PI) conversion factor to p.u.
      DOUBLE PRECISION CPATOL,ALAT,RFCTOR,THETA,PHI
      DOUBLE PRECISION PI
      DOUBLE COMPLEX DF,ERYD,EK,CFCTOR,CFCTORINV
      LOGICAL OPT,TEST,LCPAIJ,LREAD
      CHARACTER*80 TMPDIR
      CHARACTER*80 TEXT !qdos ruess
C     .. 
C     .. Local Arrays
      INTEGER NSHELL(0:NSHELD),REFPOT(NAEZD+NEMBD),ATOMIMP(NATOMIMPD)
      INTEGER KAOEZ(NATYPD,NAEZD+NEMBD)
      DOUBLE PRECISION QVEC(:,:)       ! qdos ruess, q-vectors for qdos
      ALLOCATABLE QVEC                 ! qdos ruess
      DOUBLE COMPLEX TQDOS(LMMAXD,LMMAXD,NAEZD)  ! qdos ruess
C     .. TB-Cluster arrays
      INTEGER ATOM(NACLSD,NAEZD+NEMBD),CLS(NAEZD+NEMBD),NACLS(NCLSD),
     +        EZOA(NACLSD,NAEZD+NEMBD)       
C     ..
      DOUBLE PRECISION RMTREF(NREFD),VREF(NREFD),RBASIS(3,NAEZD+NEMBD),
     +                 RCLS(3,NACLSD,NCLSD),RR(3,0:NRD)
C     ..
C     .. GMATLL = diagonal elements of the G matrix (system)
C     .. GINP   = cluster GF (ref. syst.)
      DOUBLE COMPLEX EZ(IEMXD),WEZ(IEMXD)
      DOUBLE COMPLEX
     +     GMATLL(LMMAXD,LMMAXD,NSHELD),GMAT0(LMMAXD,LMMAXD)
      DOUBLE COMPLEX,ALLOCATABLE:: 
     &     GINP(:,:,:), ! cluster GF (ref syst.)              ! GINP(NACLSD*LMGF0D,LMGF0D,NCLSD)
     &     DGINP(:,:,:) ! LLY Lloyd Energy derivative of GINP ! DGINP(NACLSD*LMGF0D,LMGF0D,NCLSD)
      DOUBLE COMPLEX TMAT(LMMAXD,LMMAXD),TSST(LMMAXD,LMMAXD,NATYPD),
     +               TREFLL(LMMAXD,LMMAXD,NREFD),
     &               DTREFLL(LMMAXD,LMMAXD,NREFD),  ! LLY Lloyd dtref/dE
     &               DTMATLL(LMMAXD,LMMAXD,NAEZD),  ! LLY Lloyd  dt/dE
     &               LLY_GRTR(IEMXD,NSPIND),        ! LLY Lloyd  Trace[ M^-1 dM/dE ], Eq.5.38 PhD Thiess
     &               LLY_G0TR(IEMXD),               ! LLY Lloyd  Trace[ X ], Eq.5.27 PhD Thiess
     &               TRALPHA(IEMXD,NSPIND),TRALPHA1,                   ! LLY Lloyd  
     &               TRACET(IEMXD,NSPIND), ! Tr[ (t-tref)^-1 d(t-tref)/dE ]  ! LLY Lloyd 
     &               ALPHAREF(0:LMAXD,NREFD),DALPHAREF(0:LMAXD,NREFD), ! LLY Lloyd Alpha matrix and deriv.
     &               TRALPHAREF(IEMXD),                                ! LLY Lloyd  
     &               CDOS_LLY(IEMXD,NSPIND),CDOSREF_LLY(IEMXD)         ! LLY Lloyd
C     .. effective (site-dependent) Delta_t^(-1) matrix ..
      DOUBLE COMPLEX WN2(LMGF0D,LMGF0D)                  ! LLY
      DOUBLE COMPLEX MSST(LMMAXD,LMMAXD,NATYPD)
      DOUBLE COMPLEX LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD),
     &               RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD)
      DOUBLE COMPLEX W1(LMMAXD,LMMAXD),WN1(LMGF0D,LMGF0D)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD)
      DOUBLE COMPLEX FACTL(LMMAXD,LMMAXD)
      DOUBLE PRECISION RATOM(3,NSHELD),RROT(48,3,NSHELD)
      INTEGER NSH1(NSHELD),NSH2(NSHELD)
      INTEGER IJTABCALC(NOFGIJD),ISH(NSHELD,NOFGIJD),JSH(NSHELD,NOFGIJD)
      INTEGER IJTABSYM(NOFGIJD),IJTABSH(NOFGIJD),IQCALC(NAEZD)
      DOUBLE COMPLEX CREL(LMMAXD,LMMAXD),RC(LMMAXD,LMMAXD),
     &               RREL(LMMAXD,LMMAXD),SRREL(2,2,LMMAXD)
      INTEGER IRREL(2,2,LMMAXD),NRREL(2,LMMAXD)
      INTEGER ICPA(NAEZD),IECPAFAIL(IEMXD)
      INTEGER NOQ(NAEZD),IQAT(NATYPD),ITOQ(NATYPD,NAEZD)
      DOUBLE PRECISION CONC(NATYPD)
      COMPLEX*8 GIMP(LMMAXD*LMMAXD)
      INTEGER ILM
      INTEGER LLY ! LLY <> 0 : apply Lloyd's formula
C     ..
C     .. relativistic mode 
C     .. rotation matrices, unitary/antiunitary symmetry flag
      DOUBLE COMPLEX DROTQ(LMMAXD,LMMAXD,NAEZD)
      LOGICAL SYMUNITARY(NSYMAXD)
C     ..
      LOGICAL VACFLAG(2)
      INTEGER KMESH(IEMXD),NOFKS(MAXMSHD)
      DOUBLE PRECISION BZKP(3,KPOIBZ,MAXMSHD)
      DOUBLE PRECISION VOLCUB(KPOIBZ,MAXMSHD),VOLBZ(MAXMSHD)
      DOUBLE COMPLEX TR1(IEMXD,NATYPD)
      INTEGER ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD)
      CHARACTER*35 INVALG(0:2)
C     ..
C     .. Arrays in Common ..
      CHARACTER*8 OPTC(32),TESTC(32)
C     ..
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     for conductivity calculation
C     INTEGER NCPAIRD
C     PARAMETER(NCPAIRD=10)
C     INTEGER IATCONDL(NCPAIRD),IATCONDR(NCPAIRD),NCONDPAIR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN
C     ..
C     .. External Subroutines ..
      EXTERNAL CALCTREF13,CHANGEREP,KLOOPZ1_QDOS,OPT,CINIT
C     ..
C     .. Common blocks ..
      COMMON /OPTC/OPTC
      COMMON /TESTC/TESTC
C     ..
C     .. Data statements
      DATA INVALG /'FULL MATRIX                        ',
     &             'BANDED MATRIX (slab)               ',
     &             'BANDED + CORNERS MATRIX (supercell)' /
      DATA IPRINT / 1 /
C     ..
C     .. MPI variables ..
CMPI  INTEGER MYRANK,NROFNODES,IERR
CMPI  COMMON /MPI/MYRANK,NROFNODES
CT3E  INTEGER ICLKTCK,ICLOCKS,IE,MEND,MSTART
CT3E  DOUBLE PRECISION TIME1,TIME2,TIME3,TIME4
CMPI  INTEGER MAPBLOCK
CMPI  EXTERNAL MPI_ALLREDUCE,MPI_COMM_RANK,MPI_COMM_SIZE,
CMPI +         MPI_FINALIZE,MPI_INIT
C     ..
C

Consistency check
      IF ( (KREL.LT.0) .OR. (KREL.GT.1) )
     &     STOP ' set KREL=0/1 (non/fully) relativistic mode in inc.p'
      IF ( (KREL.EQ.1) .AND. (NSPIND.EQ.2) ) 
     &   STOP ' set NSPIND = 1 for KREL = 1 in inc.p'

      PI = 4.D0*DATAN(1.D0)
C    
CT3E  CALL SYSTEM_CLOCK(MSTART)
CMPI  CALL MPI_INIT(IERR)
CMPI  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)
CMPI  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NROFNODES,IERR)
CMPI  IF ( MYRANK.NE.0 ) OPEN (6,STATUS='scratch',FORM='formatted')
C
C ======================================================================
C =             read in variables from unformatted files               =
C ======================================================================
C
C -------------------------------------------------------------- input1b
C
      OPEN (67,FILE='input1b.unformatted',FORM='unformatted')
      READ (67) NSRA,INS,NATYP,NAEZ,NSPIN,LMAX,NREF,ICC,IGF,
     &          NLBASIS,NRBASIS,NCPA,ICPA,ITCPAMAX,CPATOL
C ......................................................................
Consistency check 
C
      IF ( (KREL.EQ.1) .AND. (INS.NE.0) ) THEN
         WRITE(6,*)
     &        ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         STOP ' set INS = 0 in the input'
      END IF
C     
      IF ( NSRA.LE.2 ) THEN
         IF ( KREL.EQ.1 ) STOP
     &        ' KVREL <= 1 in input, but relativistic program used'
      ELSE
         IF ( KREL.EQ.0 ) STOP
     &        ' KVREL > 1 in input, but non-relativistic program used'
      END IF
      IDECI = 0
C ......................................................................
      READ (67) ALAT,RBASIS,REFPOT,RMTREF,VREF,RCLS,RR,
     &          ATOM,CLS,NCLS,EZOA,NACLS,NSHELL,KMROT,KAOEZ,IQAT,NOQ,
     &          CONC,KMESH,MAXMESH,TESTC,OPTC,LLY
      READ (67) NSYMAT,NATOMIMP,NOFGIJ,NQCALC,RATOM,RROT,NSH1,NSH2,DROTQ
      READ (67) IJTABCALC
      READ(67) ((ISH(I,L),L=1,NOFGIJ),I=1,NSHELL(0))
      READ(67) ((JSH(I,L),L=1,NOFGIJ),I=1,NSHELL(0))
      READ (67) IJTABSYM,IJTABSH,IQCALC
      READ (67) DSYMLL,INVMOD,ICHECK,ATOMIMP,SYMUNITARY
      READ (67) TMPDIR,ITMPDIR,ILTMP
      IF ( KREL.EQ.1 ) READ (67) RC,CREL,RREL,SRREL,NRREL,IRREL
      IF ( OPT('DECIMATE') ) THEN
         READ (67) LEFTTINVLL,RIGHTTINVLL,VACFLAG
         IDECI = 1
      END IF
      CLOSE (67)
C ------------------------------------------------------------- k-points
C
      OPEN (52,FILE='kpoints',FORM='formatted')
      REWIND (52)
      DO L = 1,MAXMESH
         READ (52,FMT='(I8,f15.10)') NOFKS(L),VOLBZ(L)
         DO I=1,NOFKS(L)
            READ (52,FMT=*) (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L)
         END DO
      END DO
      CLOSE (52)
C ---------------------------------------------------------- energy_mesh
C
      OPEN (67,FILE='energy_mesh',FORM='unformatted')
      READ (67) IELAST,EZ,WEZ
      CLOSE (67)
      IF (TEST('gmatasci')) 
     &           OPEN(298347,FILE='gmat.ascii',FORM='formatted')
C -------------------------------------------------------------- itermdir
C
      IF (OPT('ITERMDIR')) THEN
         OPEN(67,FILE='itermdir.unformatted',FORM='unformatted')
         READ (67)
         READ (67) DROTQ
         CLOSE(67)
         IF ( KMROT.EQ.0 ) KMROT = 1
      END IF
C ======================================================================
C =                     End read in variables                          =
C ======================================================================
C
C     If qdos option is used set IQDOSRUN so that in a first run the 
C     (t(E)-t_ref(E))^-1 matrix (fort.37) and the gref matrix can be 
C     written out for one k point, in a second run these matrices are 
C     read in to continue the calculation with the k points specified by 
C     the user in the qvec.dat file
      IF ( OPT('qdos    ') ) THEN                       ! qdos ruess
         IQDOSRUN=0                                     ! qdos ruess
      ELSE                                              ! qdos ruess
         IQDOSRUN=-1                                    ! qdos ruess
      ENDIF                                             ! qdos ruess
C Jump back here to continue with second run if qdos option is selected
  210 CONTINUE                                          ! qdos ruess
C     Reset GMATLL for calculation in second run
      IF ( IQDOSRUN.EQ.1 ) THEN                         ! qdos ruess
         DO I1 = 1,NSHELL(0)                            ! qdos ruess
            GMATLL(1:LMMAXD,1:LMMAXD,I1) = CZERO        ! qdos ruess
         ENDDO                                          ! qdos ruess
      ENDIF                                             ! qdos ruess
      OPEN(37)                                          ! qdos ruess
C                                                       ! qdos ruess
C
CT3E  CALL SYSTEM_CLOCK(MEND)
CT3E  TIME1 = (MEND-MSTART)
CT3E  TIME2 = 0.0D0
CT3E  CALL SYSTEM_CLOCK(MSTART)
C
      DO I=1,NAEZ
         DO L=1,NOQ(I)
            ITOQ(L,I) = KAOEZ(L,I)
         END DO
      END DO
      RFCTOR = ALAT/(8.D0*ATAN(1.0D0))           ! = ALAT/(2*PI)
      CFCTOR = CONE*RFCTOR
      CFCTORINV = CONE/RFCTOR

      CALL SETFACTL(FACTL,LMAX,KREL,LMMAXD)
C     
CMPI  IF( MYRANK.EQ.0 ) THEN
      WRITE(6,2100) INVALG(INVMOD)
CMPI  END IF
C


      NACLSMAX = 1
      DO IC = 1,NCLS
         IF (NACLS(IC).GT.NACLSMAX) NACLSMAX = NACLS(IC)
      ENDDO
      LRECGRF1 = WLENGTH*4*NACLSMAX*LMGF0D*LMGF0D*NCLS 

      IF (.NOT.ALLOCATED(GINP)) 
     &     ALLOCATE(  GINP(NACLSMAX*LMGF0D,LMGF0D,NCLS) )
      IF (.NOT.ALLOCATED(DGINP)) 
     &     ALLOCATE( DGINP(NACLSMAX*LMGF0D,LMGF0D,NCLS) )


      CALL OPENDAFILE(68,'gref',4,LRECGRF1,TMPDIR,ITMPDIR,ILTMP)
      CALL OPENDAFILE(69,'tmat',4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)
      CALL OPENDAFILE(70,'gmat',4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)
      IF (LLY.NE.0) THEN                                               ! LLY Lloyd
         CALL OPENDAFILE(681,'dgrefde',7,LRECGRF1,TMPDIR,ITMPDIR,ILTMP) ! LLY Lloyd: derivative of Gref
         OPEN(682,FILE='lly_g0tr_ie.ascii',FORM='FORMATTED')           ! LLY Lloyd: trace eq.5.27 PhD Thiess
         CALL OPENDAFILE(691,'dtmatde',7,LRECTMT,TMPDIR,ITMPDIR,ILTMP) ! LLY Lloyd: derivative of t-matrix
         CALL OPENDAFILE(692,'tralpha',7,LRECTRA,TMPDIR,ITMPDIR,ILTMP) ! LLY Lloyd: Tr[alpha^{-1} dalpha/dE]
      ENDIF                                                            ! LLY Lloyd

C
      LCPAIJ = .FALSE.
      IF ( ( NCPA.NE.0 ).AND.( NSHELL(0).GT.NATYP ) ) LCPAIJ = .TRUE.
C
      IF ( LCPAIJ ) 
     +     OPEN (71,ACCESS='direct',RECL=2*LRECTMT,
     +           FILE='dmatproj.unformatted',FORM='unformatted')
C
      IF ( IGF.NE.0 ) THEN

          IF ( ( OPT('GPLAIN  ') ) ) THEN
            OPEN (8888,FILE='kkrflex_green.dat')
          END IF

          IF ( ( OPT('KKRFLEX ') ) ) THEN
              OPEN (888,ACCESS='direct',
     &                   RECL=WLENGTH*2*NATOMIMP*LMMAXD*NATOMIMP*LMMAXD,
     &                   FILE='kkrflex_green',FORM='unformatted')
          END IF

          OPEN (88,ACCESS='direct',RECL=LRECGREEN,
     &             FILE='green',FORM='unformatted')
          IREC=1

          IF ( ( OPT('KKRFLEX ') ) ) THEN
            WRITE(888,REC=IREC) IELAST,NSPIN,NATOMIMP,
     &                NATOMIMP,LMMAXD,KORBIT,
     &             (EZ(IE),IE=1,IELAST),(WEZ(IE),IE=1,IELAST)
            IF ( ( OPT('GPLAIN  ') ) ) THEN
              WRITE(8888,'(5I,50000F)') IELAST,NSPIN,NATOMIMP,NATOMIMP,
     &               (LMAX+1)**2,
     &               (EZ(IE),IE=1,IELAST),(WEZ(IE),IE=1,IELAST)
            END IF
          END IF
          IF ( (.not. OPT('KKRFLEX ') ) ) THEN
            WRITE(88,REC=IREC) IELAST,NSPIN,
     &         (EZ(IE),IE=1,IELAST),(WEZ(IE),IE=1,IELAST),
     &         NATOMIMPD*LMMAXD
          END IF
      END IF

C Value of NQDOS changes to a read-in value if option qdos is applied, otherwise:
      NQDOS = 1                                         ! qdos ruess
      IF (OPT('qdos    ').AND.(IQDOSRUN.EQ.1)) THEN     ! qdos ruess
C        Read BZ path for qdos calculation:
         OPEN(67,FILE='qvec.dat')                       ! qdos ruess
         READ(67,*) NQDOS                               ! qdos ruess
         DEALLOCATE(QVEC)                               ! qdos ruess: deallocate in first run allocated array to change it
         ALLOCATE(QVEC(3,NQDOS))                        ! qdos ruess
         DO IQ = 1,NQDOS                                ! qdos ruess
            READ(67,*) (QVEC(IX,IQ),IX=1,3)             ! qdos ruess
         ENDDO                                          ! qdos ruess
         CLOSE(67)                                      ! qdos ruess
C Prepare k-mesh information to be appropriate for qdos calculation.
C The idea is that subr. KLOOPZ1 is called for only one point at a time,
C with weight equal to the full BZ; in this way we avoid changing the 
C calling list or the contents of kloopz1.
         KMESH(1:IELAST) = 1                            ! qdos ruess
         NOFKS(1) = 1                                   ! qdos ruess
         VOLCUB(1,1) = VOLBZ(1)                         ! qdos ruess
         NSYMAT = 1
      ELSEIF (OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) THEN ! qdos ruess
C Call the k loop just once with one k point to write out the fort.37 file
         ALLOCATE(QVEC(3,NQDOS))                        ! qdos ruess
         QVEC(1:3,1) = 0.D0                             ! qdos ruess
         KMESH(1:IELAST) = 1                            ! qdos ruess
         NOFKS(1) = 1                                   ! qdos ruess
         VOLCUB(1,1) = VOLBZ(1)                         ! qdos ruess
      END IF                                            ! qdos ruess


! Initialize trace for Lloyd formula                    ! LLY Lloyd
      LLY_GRTR(:,:) = CZERO ! 1:IEMXD,1:NSPIND          ! LLY Lloyd

      IF (.NOT.OPT('NEWSOSOL')) THEN
C
C    ----------------------------------------------------------------
C    |          BEGIN do loop over spins and energies               |
C    ----------------------------------------------------------------
C
      DO 370 ISPIN = 1,NSPIN
C
         NCPAFAIL = 0
C
         DO 360 IE = 1,IELAST
CMPI        IF( MYRANK.EQ.MAPBLOCK(IE,1,IELAST,1,0,NROFNODES-1) ) THEN
C
            READ (68,REC=IE) GINP
            IF (LLY.NE.0) READ (681,REC=IE) DGINP   ! LLY Lloyd

            ERYD = EZ(IE)
            NMESH = KMESH(IE)
            WRITE (6,'(A,I3,A,2(1X,F10.6),A,I3)') 
     &             ' ************ IE = ',IE,' ENERGY =',EZ(IE),
     &                ' KMESH = ', NMESH
C ********************************************************** I1 = 1,NREF
C -> calculate t(ll') of the reference system (on clusters)
C
            IF ( KREL.EQ.0 ) THEN
               DO I1 = 1,NREF
                  CALL CALCTREF13(ERYD,VREF(I1),RMTREF(I1),LMAX,LM1,
     &                    TREFLL(1,1,I1),DTREFLL(1,1,I1),                 ! LLY Lloyd
     &                    ALPHAREF(0,I1),DALPHAREF(0,I1),LMAXD+1,LMMAXD)  ! LLY Lloyd
               END DO
            ELSE
               DO I1 = 1,NREF
                  CALL CALCTREF13(ERYD,VREF(I1),RMTREF(I1),LMAX,LM1,
     &                          WN1,WN2,               ! LLY Lloyd
     &                          ALPHAREF(0,I1),DALPHAREF(0,I1),    ! LLY Lloyd
     &                          LMAXD+1,LMGF0D)
C------------------------------------------------------------
C add second spin-block for relativistic calculation and transform
C from NREL to REL representation
C------------------------------------------------------------
                  CALL CINIT(LMMAXD*LMMAXD,W1)
                  IF (LMMAXD.NE.LM1*2) STOP 'LMMAXD <> LM1*2 '
                  DO I=1,LM1
                     W1(I,I) = WN1(I,I)
                     W1(LM1+I,LM1+I) = WN1(I,I)
                  END DO
                  CALL CHANGEREP(W1,'RLM>REL',TREFLL(1,1,I1),LMMAXD,
     &                           LMMAXD,RC,CREL,RREL,'TREFLL',0)
               END DO
            END IF
C
C **********************************************************************
            TRALPHA(IE,ISPIN) = CZERO                                    ! LLY
            TRALPHAREF(IE) = CZERO                                       ! LLY
            DO I1 = 1,NATYP

               IREC = IE + IELAST* (ISPIN-1) + IELAST*NSPIN* (I1-1) 
               READ (69,REC=IREC) TMAT
               TSST(1:LMMAXD,1:LMMAXD,I1)=TMAT(1:LMMAXD,1:LMMAXD)

               IF (LLY.NE.0) THEN                                         ! LLY
                  READ (691,REC=IREC) TMAT                                ! LLY dt/dE
                  DTMATLL(1:LMMAXD,1:LMMAXD,I1) =TMAT(1:LMMAXD,1:LMMAXD)  ! LLY
                  READ (692,REC=IREC) TRALPHA1                            ! LLY
                  TRALPHA(IE,ISPIN) = TRALPHA(IE,ISPIN) + TRALPHA1        ! LLY Tr[ alpha^{-1} dalpha/dE]

                  IF (ISPIN.EQ.1) THEN  ! Ref. system is spin-independent ! LLY
                     TRALPHA1 = CZERO                                     ! LLY
                     DO L1 = 0,LMAX                                       ! LLY
                        TRALPHA1 = TRALPHA1 + (2*L1 + 1) *                ! LLY
     &                                        DALPHAREF(L1,REFPOT(I1)) /  ! LLY
     &                                         ALPHAREF(L1,REFPOT(I1))    ! LLY
                     ENDDO
                     TRALPHAREF(IE) = TRALPHAREF(IE) + TRALPHA1           ! LLY Tr[ alpharef^{-1} dalpharef/dE
                  ENDIF
              ENDIF                                                       ! LLY

            END DO
            
            IF (LLY.NE.0.AND.ISPIN.EQ.1)                            ! LLY
     &                      READ(682,FMT='(2E24.16)') LLY_G0TR(IE)  ! LLY
C ------------------------------------------------------------------------
C
C --> setting up of Delta_t moved to < KLOOPZ1 >
C
            IF (OPT('readcpa ').OR.
     &        (OPT('qdos    ').AND.(IQDOSRUN.EQ.1))) THEN     ! qdos ruess: read in cpa t-matrix
               DO ISITE = 1,NAEZ                              ! qdos ruess
                  TQDOS(:,:,ISITE) = CZERO                    ! qdos ruess
                  READ(37,*) TEXT                             ! qdos ruess
                  READ(37,*) TEXT                             ! qdos ruess
 9921             CONTINUE                                    ! qdos ruess
                  READ(37,99013) LM1,LM2,TREAD                ! qdos ruess
99013             FORMAT (2I5,1P,2D22.14)                     ! qdos ruess
                  IF ( (LM1+LM2).NE.0 ) THEN                  ! qdos ruess
                     TQDOS(LM1,LM2,ISITE) = TREAD / CFCTORINV ! qdos ruess
                     IF ( (LM1+LM2).LT.2*LMMAXD ) GOTO 9921   ! qdos ruess
                  END IF                                      ! qdos ruess
               ENDDO                                          ! qdos ruess
            END IF                                            ! qdos ruess
C  Loop over all QDOS points and change volume for KLOOPZ run accordingly
            DO 200 IQ = 1,NQDOS                               ! qdos ruess
            IF (OPT('qdos    ')) BZKP(:,1,1) = QVEC(:,IQ)     ! qdos ruess: Set q-point x,y,z
C
            CALL KLOOPZ1_QDOS(ERYD,GMATLL,INS,ALAT,IE,IGF,NSHELL,NAEZ,
     &                 NOFKS(NMESH),VOLBZ(NMESH),BZKP(1,1,NMESH),
     &           VOLCUB(1,NMESH),CLS,NACLS,NACLSMAX,RR,RBASIS,EZOA,ATOM,
     &                 RCLS,ICC,GINP,IDECI,LEFTTINVLL(1,1,1,1,IE), 
     &                 RIGHTTINVLL(1,1,1,1,IE),VACFLAG,
     &                 NLBASIS,NRBASIS,FACTL,NATOMIMP,NSYMAT,DSYMLL,
     &                 RATOM,RROT,NSH1,NSH2,IJTABSYM,IJTABSH,ICHECK,
     &                 INVMOD,REFPOT,TREFLL,TSST,MSST,CFCTOR,  
     &                 CFCTORINV,CREL,RC,RREL,SRREL,IRREL,NRREL,DROTQ,
     &                 SYMUNITARY,KMROT,NATYP,NCPA,ICPA,ITCPAMAX,
     &                 CPATOL,NOQ,IQAT,ITOQ,CONC,IPRINT,ICPAFLAG,
     &                 ISPIN,NSPINDD,
     &                 TQDOS,IQDOSRUN,                         ! qdos
     &    DTREFLL,DTMATLL,DGINP,LLY_GRTR(IE,ISPIN),TRACET(IE,ISPIN),LLY) ! LLY Lloyd
C
C           Skip this part if first part of the qdos is running
            IF ( .NOT.(OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) ) THEN
              IF (NCPA.NE.0) THEN
               IF (ICPAFLAG .NE. 0) THEN
                  NCPAFAIL = NCPAFAIL + 1
                  IECPAFAIL(NCPAFAIL)= IE
                END IF
              END IF  ! (NCPA.NE.0)

              DO I1 = 1,NSHELL(0)
                DO LM1=1,LMMAXD
                  DO LM2=1,LMMAXD
                     GMAT0(LM1,LM2)=GMATLL(LM1,LM2,I1)
                  ENDDO
                ENDDO
C               IREC = IE + IELAST* (ISPIN-1) + IELAST*NSPIN* (I1-1) ! <-- before introducing qdos
                IREC = IQ + NQDOS * (IE-1) + NQDOS * IELAST *        ! qdos ruess: (without qdos, IQ=NQ=1)
     &                   (ISPIN-1) + NQDOS * IELAST * NSPIN * (I1-1) ! qdos ruess
                WRITE (70,REC=IREC) GMAT0
              ENDDO
              IF (TEST('gmatasci')) THEN
                 WRITE(*,*) 'Writing out gmat.ascii'
                 DO I1 = 1,NSHELL(0)
                 DO LM1=1,LMMAXD
                 DO LM2=1,LMMAXD
                    WRITE(298347,FMT='(3I5,2E25.16)') 
     &              I1,LM1,LM2,GMATLL(LM1,LM2,I1)
                 ENDDO
                 ENDDO
                 ENDDO
              ENDIF


            IF ( NATOMIMP==1 .and. OPT('KKRFLEX ') ) THEN
                I1=ATOMIMP(1)
                    irec = ielast*(ispin-1)+ ie+1
                      ILM=0
                      GIMP=(0.e0,0.e0) !complex*8
                      DO LM2=1,LMMAXD
                          DO LM1=1,LMMAXD
                            ILM=ILM+1
                            GIMP(ILM)=GMATLL(LM1,LM2,I1)
                          ENDDO
                      ENDDO
                    WRITE(888,REC=irec) GIMP
                    IF ( ( OPT('GPLAIN  ') ) ) THEN
                      WRITE(8888,'(50000E)') GIMP
                    END IF
!                   END DO
!                 END DO
            END IF

            IF ( LCPAIJ ) THEN
               DO I1 = 1,NATYP
                  DO LM2=1,LMMAXD
                     DO LM1=1,LMMAXD
                        GMAT0(LM1,LM2) = TSST(LM1,LM2,I1)
                        W1(LM1,LM2)    = MSST(LM1,LM2,I1)
                     END DO
                  END DO
                  IREC = IE + IELAST* (ISPIN-1) + IELAST*NSPIN* (I1-1)
                  WRITE (71,REC=IREC) GMAT0,W1
               END DO
              END IF  ! ( LCPAIJ )
C
            ENDIF    ! ( .NOT.(OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) )
 200        CONTINUE ! IQ = 1,NQDOS                                       ! qdos ruess


            IF (LLY.NE.0) THEN                                            ! LLY Lloyd

               IF (LLY.NE.2) THEN                                         ! LLY Lloyd
                  CDOS_LLY(IE,ISPIN) =   TRALPHA(IE,ISPIN)                ! LLY Lloyd
     &                - LLY_GRTR(IE,ISPIN) / VOLBZ(1) + LLY_G0TR(IE)      ! LLY Lloyd
               ELSE                                                       ! LLY Lloyd
                  CDOS_LLY(IE,ISPIN) =    TRACET(IE,ISPIN)                ! LLY Lloyd
     &                 + TRALPHAREF(IE)        ! LLY Lloyd
     &                 - LLY_GRTR(IE,ISPIN) / VOLBZ(1) +  LLY_G0TR(IE)    ! LLY Lloyd
               ENDIF                                                      ! LLY Lloyd

               IF (ISPIN.EQ.1)                                            ! LLY Lloyd
     &              CDOSREF_LLY(IE) = TRALPHAREF(IE) - LLY_G0TR(IE)       ! LLY Lloyd

               IF (TEST('GMAT=0  ')) THEN                                 ! LLY Lloyd
                  CDOS_LLY(IE,ISPIN) = TRALPHA(IE,ISPIN)                  ! LLY Lloyd
                  IF (LLY.EQ.2) CDOS_LLY(IE,ISPIN) =                      ! LLY Lloyd
     &                 TRACET(IE,ISPIN) + TRALPHAREF(IE)                  ! LLY Lloyd
               ENDIF                                                      ! LLY Lloyd

               CDOS_LLY(IE,ISPIN) = CDOS_LLY(IE,ISPIN) / PI               ! LLY Lloyd

            ENDIF                                                         ! LLY Lloyd

c ------------------------------------------------------------------------

CMPI     END IF
 360     CONTINUE               ! IE = 1,IELAST

CMPI     IF( MYRANK.EQ.0 ) THEN
         IF( NCPAFAIL .NE. 0 ) THEN
            WRITE(6,*)
            WRITE(6,'(1X,79(''*''),/)')
            WRITE(6,99019) CPATOL, NCPAFAIL,
     &           (IECPAFAIL(IE),DBLE(EZ(IECPAFAIL(IE))),
     &           IE=1,NCPAFAIL)
            WRITE(6,'(1X,79(''*''),/)')
            WRITE(6,*)
         ELSE
            IF( NCPA .NE. 0 ) THEN
               WRITE(6,*)
               WRITE(6,99020)
               WRITE(6,*)
            END IF
         END IF
CMPI     END IF
C         
 370  CONTINUE                  !  ISPIN = 1,NSPIN
      

      ELSE ! NEW SOC SOLVER 

      DO 460 IE=1,IELAST

       
c read in Green function of reference system
       READ (68,REC=IE) GINP
       IF (LLY.NE.0) READ(681,REC=IE) DGINP  ! LLY
   
       ERYD = EZ(IE)
       NMESH = KMESH(IE)
            WRITE (6,'(A,I3,A,2(1X,F10.6),A,I3)') 
     &             ' ************ IE = ',IE,' ENERGY =',EZ(IE),
     &                ' KMESH = ', NMESH

c construct t matrix for reference system (now is always double matrix)
       DO I1 = 1,NREF
        CALL CALCTREF13(ERYD,VREF(I1),RMTREF(I1),LMAX,LM1,
     +                WN1,WN2,                                ! LLY
     &                ALPHAREF(0,I1),DALPHAREF(0,I1),         ! LLY 
     &                LMAXD+1,LMGF0D)
        DO I=1,LM1
         TREFLL(I,I,I1) = WN1(I,I)
         TREFLL(LM1+I,LM1+I,I1) = WN1(I,I)
         DTREFLL(I,I,I1) = WN2(I,I)                 ! LLY
         DTREFLL(LM1+I,LM1+I,I1) = WN2(I,I)         ! LLY
        ENDDO
       ENDDO ! I1
       TRALPHA(IE,NSPINDD)=CZERO
       TRALPHAREF(IE)=CZERO
c read in t matrix

       LREAD = .FALSE.
       INQUIRE(file='nonco_angle.dat',EXIST=LREAD)
       IF (LREAD) OPEN(UNIT=10,FILE='nonco_angle.dat',FORM='FORMATTED')
       THETA = 0.D0
       PHI = 0.D0


       DO I1 = 1,NATYPD

c read in theta and phi for noncolinear
        IF (LREAD) READ(10,*) THETA,PHI

c read in t-matrix from file
        IREC = IE + IELAST*(I1-1)
        READ (69,REC=IREC) TMAT
         

c rotate t-matrix from local to global frame
        CALL ROTATEMATRIX(TMAT,THETA,PHI,LMGF0D,0)

         DO LM1=1,LMMAXD
          DO LM2=1,LMMAXD
           TSST(LM1,LM2,I1)=TMAT(LM1,LM2)
          ENDDO
         ENDDO

       IF (LLY.NE.0) THEN
        READ(691,REC=IREC) TMAT ! LLY
        CALL ROTATEMATRIX(TMAT,THETA,PHI,LMGF0D,0) ! LLY
         DO LM1=1,LMMAXD
          DO LM2=1,LMMAXD
           DTMATLL(LM1,LM2,I1)=TMAT(LM1,LM2) ! LLY
          ENDDO
         ENDDO
         READ(692,REC=IREC) TRALPHA1
         TRALPHA(IE,NSPINDD)=TRALPHA(IE,NSPINDD)+TRALPHA1 ! LLY
         TRALPHA1=CZERO
          DO L1=0,LMAX
           TRALPHA1=TRALPHA1+(2*L1+1)*
     &              DALPHAREF(L1,REFPOT(I1))/ALPHAREF(L1,REFPOT(I1)) ! LLY
          ENDDO
          TRALPHAREF(IE)=TRALPHAREF(IE)+TRALPHA1 ! LLY
       ENDIF ! LLY
       ENDDO ! I1
       CLOSE(10)
       IF (LLY.NE.0) READ(682,FMT='(2E24.16)') LLY_G0TR(IE) ! LLY

c QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS
            IF (OPT('readcpa ').OR.
     &        (OPT('qdos    ').AND.(IQDOSRUN.EQ.1))) THEN     ! qdos ruess: read in cpa t-matrix
               DO ISITE = 1,NAEZ                              ! qdos ruess
                  TQDOS(:,:,ISITE) = CZERO                    ! qdos ruess
                  READ(37,*) TEXT                             ! qdos ruess
                  READ(37,*) TEXT                             ! qdos ruess
 9920             CONTINUE                                    ! qdos ruess
                  READ(37,99014) LM1,LM2,TREAD                ! qdos ruess
99014             FORMAT (2I5,1P,2D22.14)                     ! qdos ruess
                  IF ( (LM1+LM2).NE.0 ) THEN                  ! qdos ruess
                     TQDOS(LM1,LM2,ISITE) = TREAD / CFCTORINV ! qdos ruess
                     IF ( (LM1+LM2).LT.2*LMMAXD ) GOTO 9920   ! qdos ruess
                  END IF                                      ! qdos ruess
               ENDDO                                          ! qdos ruess
            END IF                                            ! qdos ruess
C  Loop over all QDOS points and change volume for KLOOPZ run accordingly
            DO 220 IQ = 1,NQDOS                               ! qdos ruess
            IF (OPT('qdos    ')) BZKP(:,1,1) = QVEC(:,IQ)     ! qdos ruess: Set q-point x,y,z

            CALL KLOOPZ1_QDOS(ERYD,GMATLL,INS,ALAT,IE,IGF,NSHELL,NAEZ,
     &              NOFKS(NMESH),VOLBZ(NMESH),BZKP(1,1,NMESH),
     &           VOLCUB(1,NMESH),CLS,NACLS,NACLSMAX,RR,RBASIS,EZOA,ATOM,
     &              RCLS,ICC,GINP,IDECI,LEFTTINVLL(1,1,1,1,IE),
     &              RIGHTTINVLL(1,1,1,1,IE),VACFLAG,
     &              NLBASIS,NRBASIS,FACTL,NATOMIMP,NSYMAT,DSYMLL,
     &              RATOM,RROT,NSH1,NSH2,IJTABSYM,IJTABSH,ICHECK,
     &              INVMOD,REFPOT,TREFLL,TSST,MSST,CFCTOR,
     &              CFCTORINV,CREL,RC,RREL,SRREL,IRREL,NRREL,DROTQ,
     &              SYMUNITARY,KMROT,NATYP,NCPA,ICPA,ITCPAMAX,
     &              CPATOL,NOQ,IQAT,ITOQ,CONC,IPRINT,ICPAFLAG,
     &              1,NSPINDD,
     &              TQDOS,IQDOSRUN,     ! qdos
     &    DTREFLL,DTMATLL,DGINP,LLY_GRTR(IE,1),TRACET(IE,1),LLY) ! LLY Lloyd
          
!           Skip this part if first part of the qdos is running
            IF ( .NOT.(OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) ) THEN
               DO I1 = 1,NSHELL(0)
                  GMAT0(1:LMMAXD,1:LMMAXD) =GMATLL(1:LMMAXD,1:LMMAXD,I1)
                  !IREC = IE + IELAST*(I1-1)     ! <-- before introducing qdos
                  IREC = IQ + NQDOS * (IE-1) + NQDOS * IELAST * (I1-1) ! qdos ruess
                  WRITE (70,REC=IREC) GMAT0
               ENDDO
               IF (TEST('gmatasci')) THEN
                 WRITE(*,*) 'Writing out gmat.ascii'
                 DO I1 = 1,NSHELL(0)
                 DO LM1=1,LMMAXD
                 DO LM2=1,LMMAXD
                    WRITE(298347,FMT='(3I5,2E25.16)') 
     &              I1,LM1,LM2,GMATLL(LM1,LM2,I1)
                 ENDDO
                 ENDDO
                 ENDDO
               ENDIF

               IF ( NATOMIMP==1 .and. OPT('KKRFLEX ') ) THEN
                  I1=ATOMIMP(1)
                  IREC = IE+1
                  ILM=0
                  GIMP=(0.e0,0.e0) ! complex*8
                  DO LM2=1,LMMAXD
                     DO LM1=1,LMMAXD
                        ILM=ILM+1
                        GIMP(ILM)=GMATLL(LM1,LM2,I1)
                     ENDDO
                  ENDDO
                  WRITE(888,REC=IREC) GIMP
               ENDIF

c       STOP
            ENDIF               ! ( .NOT.(OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) )
 220     CONTINUE               ! IQ = 1,NQ                                         ! qdos ruess

            IF (LLY.NE.0) THEN                                      ! LLY 

               CDOS_LLY(IE,1) =   TRALPHA(IE,1)                     ! LLY
     &                - LLY_GRTR(IE,1) / VOLBZ(1) + 2d0*LLY_G0TR(IE)    ! LLY 

               CDOS_LLY(IE,1) = CDOS_LLY(IE,1) / PI                 ! LLY 

               CDOSREF_LLY(IE) = TRALPHAREF(IE) - LLY_G0TR(IE)      ! LLY 


            ENDIF                                                   ! LLY 

 460  CONTINUE                  ! IE=1,IELAST
      
      ENDIF


C
C    ----------------------------------------------------------------
C    |           END of do loop over spins and energies             |
C    ----------------------------------------------------------------
C
      CLOSE (68)
      IF ( IGF.NE.0 ) CLOSE (88)

      IF (LLY.NE.0) THEN                                                 ! LLY Lloyd
         OPEN (701,FILE='cdosdiff_lly.dat',FORM='FORMATTED')             ! LLY Lloyd
         OPEN (702,FILE='cdosref_lly.dat',FORM='FORMATTED')              ! LLY Lloyd
         DO IE = 1,IELAST                                                ! LLY
           WRITE(702,FMT='(10E25.16)') DREAL(EZ(IE)),CDOSREF_LLY(IE),    ! LLY
     &                TRALPHAREF(IE),LLY_G0TR(IE)
         ENDDO                                                           ! LLY
         IF (.NOT.OPT('NEWSOSOL')) THEN
          DO ISPIN = 1,NSPIN                               ! LLY Lloyd
           DO IE = 1,IELAST                                              ! LLY
            WRITE(701,FMT='(10E25.16)') DREAL(EZ(IE)),CDOS_LLY(IE,ISPIN)
     &      ,TRALPHA(IE,ISPIN),LLY_GRTR(IE,ISPIN)                        ! LLY
           ENDDO                                                         ! LLY
          ENDDO                         
         ELSE                         
          DO IE = 1,IELAST                                               ! LLY
           WRITE(701,FMT='(10E25.16)') DREAL(EZ(IE)),CDOS_LLY(IE,1)     
     &      ,TRALPHA(IE,1),LLY_GRTR(IE,1)                                ! LLY
          ENDDO                                                          ! LLY
         ENDIF ! .NOT.OPT('NEWSOSOL')                                    ! LLY
         CLOSE(701)                                                      ! LLY
         CLOSE(702)                                                      ! LLY
      ENDIF                                                              ! LLY
C
CMPI  IF( MYRANK.EQ.0 ) THEN
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      IF ( ( OPT('XCPL    ') ).AND.( ICC.LE.0 ) ) THEN
!       open(23,FILE='energies.dat',STATUS='unknown',FORM='formatted')
         DO IE = 1,IELAST
            ERYD = EZ(IE)
            DF = WEZ(IE)/DBLE(NSPIN)
!          write(23,*)ERYD
            CALL TBXCCPLJIJ(69,IE,IELAST,ERYD,DF,NCPA,NAEZ,NATYP,NOQ,
     &                      ITOQ,IQAT,NSHELL,NATOMIMP,ATOMIMP,RATOM,
     &                      FACTL,NOFGIJ,NQCALC,IQCALC,IJTABCALC,
     &                      IJTABSYM,IJTABSH,ISH,JSH,
     &                      DSYMLL,IPRINT,NATYPD,NAEZD,NSHELD,LMMAXD)
         END DO
!         close (23)
      END IF
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CMPI  END IF
C
      CLOSE (69)
      CLOSE (70)
      IF ( LCPAIJ ) CLOSE (71)

      CLOSE(37)                                    ! qdos ruess
C     Finished first qdos run. Now re-run the whole kkr1b program to 
C     calculate the GF for every energy (defined in inputcard) and 
C     kpoint (defined in qvec.dat)
      IQDOSRUN=IQDOSRUN+1                          ! qdos ruess
      CLOSE(681)
      CLOSE(682)
      CLOSE(691)
      CLOSE(692)
      IF (IQDOSRUN.EQ.1) GO TO 210                 ! qdos ruess
C ----------------------------------------------------------------------C
CT3E  CALL SYSTEM_CLOCK(MEND)
CT3E  TIME3 = (MEND-MSTART)
CT3E  TIME4 = 0.0D0
CT3E  IF (MYRANK.EQ.0) WRITE (6,FMT=*) 'elapsed time ',
CT3E +  TIME1,' seconds',
CT3E +  TIME2,' seconds',
CT3E +  TIME3,' seconds',
CT3E +  TIME4,' seconds'
CT3E  IF (MYRANK.EQ.0) WRITE (6,FMT=*) 'elapsed time ',
CT3E +  TIME1+TIME2+TIME3+TIME4,' seconds'
CT3E  TIME1=TIME1*NROFNODES
CT3E  TIME2=TIME2*NROFNODES
CT3E  TIME3=TIME3*NROFNODES
CT3E  TIME4=TIME4*NROFNODES
CT3E  IF (MYRANK.EQ.0) WRITE (6,FMT=*) 'elapsed time ',
CT3E +  TIME1,' seconds',
CT3E +  TIME2,' seconds',
CT3E +  TIME3,' seconds',
CT3E +  TIME4,' seconds'
CT3E  IF (MYRANK.EQ.0) WRITE (6,FMT=*) 'elapsed time ',
CT3E +  TIME1+TIME2+TIME3+TIME4,' seconds'
CMPI  CALL MPI_FINALIZE(IERR)
      STOP
 2100 FORMAT(/,79(1H=),/,5X,' Inversion algorithm used : ',A,/,
     &     79(1H=),/)
99019 FORMAT (/,1X,79('*'),/,
     &             ' tolerance for CPA-cycle:',F15.7,/,
     &             ' CPA not converged for',I3,' energies:',/,
     &             3(' E:',I3,F7.4,:,2X))
99020 FORMAT (/,1X,79('*'),/,25X,'no problems with',
     &         '  CPA-cycle ',/,1X,79('*'),/)


      END
