      MODULE MOD_MAIN1B
      
      implicit none
      
      contains
      
      
      
      subroutine main1b()
      
      use mod_timing
      use mod_types, only: t_tgmat,t_inc,t_lloyd,t_cpa,init_t_cpa,t_imp
#ifdef CPP_MPI
     &                     ,t_mpi_c_grid, save_t_mpi_c_grid,
     &                      get_ntot_pT_ioff_pT_2D, 
     &                      init_params_t_imp, init_t_imp, 
     &                      bcast_t_imp_scalars, bcast_t_imp_arrays
      use mpi
#endif
      use mod_mympi, only: myrank, master
#ifdef CPP_MPI
     &                     ,find_dims_2d,
     &                     distribute_linear_on_tasks, MPIadapt
#endif
      use mod_wunfiles
      use mod_tbxccpljij, only: tbxccpljij
      use mod_tbxccpljijdij, only: tbxccpljijdij
      use mod_version_info

      IMPLICIT NONE
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
      INTEGER IELAST,NPOL,INS,LMAX,NATYP,NREF,NSPIN,NSRA
      INTEGER KMROT,INVMOD,ICC,IGF,IPRINT,IDECI
      INTEGER MAXMESH,NAEZ,I,ISPIN,I1,IE,IREC,L,LM1,LM2,L1,RECLENGTH
      INTEGER NATOMIMP,NSYMAT,NCPA,NCPAFAIL,ICPAFLAG,ITCPAMAX
      INTEGER NMESH,NOFGIJ,NQCALC
      INTEGER ITMPDIR,ILTMP
      INTEGER IQ,NQDOS ! qdos ruess:number of qdos points
      INTEGER IX,ISITE ! qdos ruess
      INTEGER IQDOSRUN ! qdos ruess: counter to organise qdos run
      INTEGER NCLS,NACLSMAX,IC,LRECGRF1
      DOUBLE COMPLEX TREAD ! qdos ruess
C     RFCTOR=A/(2*PI) conversion factor to p.u.
      DOUBLE PRECISION CPATOL,ALAT,RFCTOR,THETA,PHI
      DOUBLE PRECISION PI,THETAS(NATYPD),PHIS(NATYPD)
      DOUBLE PRECISION RCLSIMP(3,NATOMIMPD)
      DOUBLE COMPLEX ERYD,CFCTOR,CFCTORINV
      LOGICAL OPT,TEST,LCPAIJ
      CHARACTER*80 TMPDIR
#ifndef CPP_MPI
      CHARACTER*80 TEXT !qdos ruess
#endif
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
      DOUBLE COMPLEX, allocatable :: DTMTRX(:,:) ! for GREENIMP

C     .. effective (site-dependent) Delta_t^(-1) matrix ..
      DOUBLE COMPLEX WN2(LMGF0D,LMGF0D)                  ! LLY
      DOUBLE COMPLEX MSST(LMMAXD,LMMAXD,NATYPD)
      double complex, allocatable :: LEFTTINVLL(:,:,:,:,:),
     +                               RIGHTTINVLL(:,:,:,:,:)
      DOUBLE COMPLEX W1(LMMAXD,LMMAXD),WN1(LMGF0D,LMGF0D)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD)
      DOUBLE COMPLEX FACTL(LMMAXD,LMMAXD)
      DOUBLE PRECISION RATOM(3,NSHELD),RROT(48,3,NSHELD)
      INTEGER NSH1(NSHELD),NSH2(NSHELD)
      INTEGER IJTABCALC(NOFGIJD),IJTABCALC_I(NOFGIJD)!,
!      &        ISH(NSHELD,NOFGIJD),JSH(NSHELD,NOFGIJD)
      INTEGER, ALLOCATABLE :: ISH(:,:),JSH(:,:)
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
      DOUBLE PRECISION, allocatable :: BZKP(:,:,:)
      DOUBLE PRECISION VOLCUB(KPOIBZ,MAXMSHD),VOLBZ(MAXMSHD)
      INTEGER ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD)
      CHARACTER*35 INVALG(0:2)
#ifdef CPP_MPI
      integer :: ie_start
#endif
      integer :: ie_num, ie_end, ierr

      ! for OPERATOR option
      logical :: lexist, operator_imp
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
C     .. Data statements
      DATA INVALG /'FULL MATRIX                        ',
     &             'BANDED MATRIX (slab)               ',
     &             'BANDED + CORNERS MATRIX (supercell)' /
      DATA IPRINT / 0 /
C     ..
C
      if(t_inc%i_write>0) IPRINT=1

      ! allocatable arrays
      allocate(LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD),
     &         RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD),
     &         BZKP(3,KPOIBZ,MAXMSHD), stat=ierr)
      if(ierr/=0) stop '[main1b] Error allocating leftinvll etc.'
      allocate(ISH(NSHELD,NOFGIJD),JSH(NSHELD,NOFGIJD), stat=ierr)
      if(ierr/=0) stop '[main1b] Error allocating ISH and JSH'

Consistency check
      IF ( (KREL.LT.0) .OR. (KREL.GT.1) )
     &     STOP ' set KREL=0/1 (non/fully) relativistic mode in inc.p'
      IF ( (KREL.EQ.1) .AND. (NSPIND.EQ.2) ) 
     &   STOP ' set NSPIND = 1 for KREL = 1 in inc.p'

      PI = 4.D0*DATAN(1.D0)
C
C ======================================================================
C =             read in variables from unformatted files               =
C ======================================================================
C
C -------------------------------------------------------------- input1b
      call get_params_1b(t_params,NATYPD,NACLSD,IELAST,NPOL,
     &       NCLSD,NREFD,NEMBD,NAEZD,NSRA,INS,NAEZ,NATYP,
     &       NSPIN,LMAX,NCLS,NREF,LLY,KREL,ATOM,CLS,NACLS,REFPOT,EZ,
     &       ITMPDIR,ILTMP,ALAT,RCLS,IEMXD,RMTREF,VREF,TMPDIR,
     &       NSHELD,NPRINCD,KPOIBZ,ATOMIMP,NATOMIMPD,
     &       ICC,IGF,NLBASIS,NRBASIS,NCPA,ICPA,ITCPAMAX,CPATOL,NRD,
     &       IDECI,RBASIS,RR,EZOA,NSHELL,KMROT,KAOEZ,ISH,JSH,NSH1,NSH2,
     &       NOQ,IQAT,NOFGIJD,NOFGIJ,NATOMIMP,CONC,KMESH,MAXMESH,NSYMAT,
     &       NQCALC,RATOM,RROT,DROTQ,IJTABCALC,IJTABCALC_I,IJTABSYM,
     &       IJTABSH,IQCALC,
     &       DSYMLL,INVMOD,ICHECK,SYMUNITARY,RC,CREL,RREL,SRREL,
     &       NRREL,IRREL,LEFTTINVLL,RIGHTTINVLL,VACFLAG,NOFKS,VOLBZ,
     &       BZKP,VOLCUB,WEZ,NEMBD1,LMMAXD,NSYMAXD,NSPINDD,MAXMSHD,
     &       RCLSIMP)
     
     
      IF (TEST('gmatasci'))
     &           OPEN(298347,FILE='gmat.ascii',FORM='formatted')

C ======================================================================
C =                     End read in variables                          =
C ======================================================================
C
C
C ************************************************************************** !fswrt
C     open file to store the output for the (external) Fermi-surface program !fswrt
C      ---> this file is already partly filled with data by main0. More data !fswrt
C           will be stored in                                                !fswrt
      IF (OPT('FERMIOUT') .and. myrank==master) then                         !fswrt
         OPEN(6801,FILE='TBkkr_container.txt', FORM='formatted',             !fswrt
     &             POSITION='append')                                        !fswrt
      ENDIF                                                                  !fswrt

!     open file for WRTGREEN option (writes green_host file for
!     GMATLL_GES creation in zulapi part) file is filled in ROTGLL called in kloopz
      IF (OPT('WRTGREEN') .and. myrank==master) THEN
         OPEN(58,FILE='green_host', FORM='formatted')
      ENDIF
C
C
C
C
C     If qdos option is used set IQDOSRUN so that in a first run the 
C     (t(E)-t_ref(E))^-1 matrix (tmat.qdos) and the gref matrix can be 
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
      
      if ((opt('qdos    ')) .and. (opt('deci-out'))) then
         stop 'ERROR: qdos and deci-out cannot be used simultaniously'
      elseif (opt('qdos    ')) then
#ifdef CPP_MPI
         OPEN (37,ACCESS='direct', RECL=WLENGTH*4,
     &         FILE='tmat.qdos',FORM='unformatted')
#else
         OPEN(37, File='tmat.qdos', form='formatted')      ! qdos ruess
#endif
      elseif (opt('deci-out')) then
         open(37, file='decifile', form='formatted', position='append') ! ruess: needed in case of deci-out option to prepare decifile
      end if
         
         
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
      if(t_inc%i_write>0) WRITE(1337,2100) INVALG(INVMOD)
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


      if (t_tgmat%gref_to_file) then
         CALL OPENDAFILE(68,'gref',4,LRECGRF1,TMPDIR,ITMPDIR,ILTMP)
      end if
      if (t_tgmat%tmat_to_file) then
         CALL OPENDAFILE(69,'tmat',4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)
      end if
      if (t_tgmat%gmat_to_file) then
         CALL OPENDAFILE(70,'gmat',4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)
      end if
      IF (LLY.NE.0) THEN                                               ! LLY Lloyd
         if (t_lloyd%dgref_to_file) then
         CALL OPENDAFILE(681,'dgrefde',7,LRECGRF1,TMPDIR,ITMPDIR,ILTMP) ! LLY Lloyd: derivative of Gref
         end if
         if (t_lloyd%g0tr_to_file) then
         OPEN(682,FILE='lly_g0tr_ie.ascii',FORM='FORMATTED')           ! LLY Lloyd: trace eq.5.27 PhD Thiess
         end if
         if (t_lloyd%dtmat_to_file) then
         CALL OPENDAFILE(691,'dtmatde',7,LRECTMT,TMPDIR,ITMPDIR,ILTMP) ! LLY Lloyd: derivative of t-matrix
         end if
         if (t_lloyd%tralpha_to_file) then
         CALL OPENDAFILE(692,'tralpha',7,LRECTRA,TMPDIR,ITMPDIR,ILTMP) ! LLY Lloyd: Tr[alpha^{-1} dalpha/dE]
         end if
      ENDIF                                                            ! LLY Lloyd

C
      LCPAIJ = .FALSE.
      IF ( ( NCPA.NE.0 ).AND.( NSHELL(0).GT.NATYP ) ) LCPAIJ = .TRUE.
C
      IF ( LCPAIJ ) THEN
       if(t_cpa%dmatproj_to_file) then
           OPEN (71,ACCESS='direct',RECL=2*LRECTMT,
     +           FILE='dmatproj.unformatted',FORM='unformatted')
       else
#ifndef CPP_MPI
           call init_t_cpa(t_inc,t_cpa,IELAST)
#else
           call init_t_cpa(t_inc,t_cpa,t_mpi_c_grid%ntot2)
#endif
       end if!t_cpa%dmatproj_to_file
      END IF!LCPAIJ
C
      IF ( IGF.NE.0 ) THEN

          IF ( ( OPT('GPLAIN  ') ) ) THEN
            OPEN (8888,FILE='kkrflex_green.dat')
          END IF

          IF ( ( OPT('KKRFLEX ') ) ) THEN
              RECLENGTH = WLENGTH*2*NATOMIMP*LMMAXD*NATOMIMP*LMMAXD
              ! sometimes (lmax=2) the record length might be too small to store the parameters, then reclength needs to be bigger
              IF(RECLENGTH<8*IELAST+6) THEN
                 stop '[main1b] record length for kkrflex_green is
     &too small to store parameters, use either more atoms in the 
     &cluster, a higher lmax or less energy points'
              END IF
              OPEN (888,ACCESS='direct',
     &                   RECL=RECLENGTH,
     &                   FILE='kkrflex_green',FORM='unformatted')
          END IF

C==== the following write-out has been disabled, because it was assumed to be    !no-green
C====  obsolete with the implementation of the MPI-communicated arrays. If I am  !no-green
C====  wrong and the write-out is needed in subsequent parts, construct a        !no-green
C====  test-option around it so that it is only written out in this case.        !no-green
C         OPEN (88,ACCESS='direct',RECL=LRECGREEN,                               !no-green
C    &             FILE='green',FORM='unformatted')                              !no-green
          IREC=1

          IF ( ( OPT('KKRFLEX ') ) ) THEN
            if(myrank==master) then
              WRITE(888,REC=IREC) IELAST,NSPIN,NATOMIMP,
     &                NATOMIMP, LMMAXD, KORBIT,
     &             (EZ(IE),IE=1,IELAST),(WEZ(IE),IE=1,IELAST)
              IF ( ( OPT('GPLAIN  ') ) ) THEN
                WRITE(8888,'(5I5,50000F14.7)') IELAST,NSPIN,NATOMIMP,
     &               NATOMIMP, (LMAX+1)**2,
     &               (EZ(IE),IE=1,IELAST),(WEZ(IE),IE=1,IELAST)
              END IF
            end if
#ifdef CPP_MPI
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
          END IF
c         IF ( (.not. OPT('KKRFLEX ') ) ) THEN                !no-green
c           WRITE(88,REC=IREC) IELAST,NSPIN,                  !no-green
c    &         (EZ(IE),IE=1,IELAST),(WEZ(IE),IE=1,IELAST),    !no-green
c    &         NATOMIMPD*LMMAXD                               !no-green
c         END IF                                              !no-green
      END IF

C Value of NQDOS changes to a read-in value if option qdos is applied, otherwise:
      NQDOS = 1                                         ! qdos ruess
      IF (OPT('qdos    ').AND.(IQDOSRUN.EQ.1)) THEN     ! qdos ruess
C        Read BZ path for qdos calculation:
         OPEN(67,FILE='qvec.dat')                       ! qdos ruess
         READ(67,*) NQDOS                               ! qdos ruess
         DEALLOCATE(QVEC, stat=ierr)                    ! qdos ruess: deallocate in first run allocated array to change it
         if(ierr/=0) stop '[main1b] Error deallocating qvec' 
         ALLOCATE(QVEC(3,NQDOS), stat=ierr)             ! qdos ruess
         if(ierr/=0) stop '[main1b] Error reallocating qvec'
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
C Call the k loop just once with one k point to write out the tmat.qdos file
         ALLOCATE(QVEC(3,NQDOS), stat=ierr)             ! qdos ruess
         if(ierr/=0) stop '[main1b] Error allocating qvec'
         QVEC(1:3,1) = 0.D0                             ! qdos ruess
         KMESH(1:IELAST) = 1                            ! qdos ruess
         NOFKS(1) = 1                                   ! qdos ruess
         VOLCUB(1,1) = VOLBZ(1)                         ! qdos ruess
      END IF                                            ! qdos ruess

      NCPAFAIL = 0

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
C
!EEEEEEEEEEEEEEEEEEEEEEEEE  NO-SOC EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#ifdef CPP_MPI         
         ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
         ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
         
         DO 360 ie_num=1,ie_end
            IE = ie_start+ie_num
            
            !start timing measurement for this ie, needed for MPIadapt
            if(MPIadapt.and.t_mpi_c_grid%myrank_ie==0) then
              call timing_start('time_1b_ie')
            end if
#else
         DO 360 IE = 1,IELAST
            ie_num = ie
            ie_end = ielast
#endif
C
            ! write energy into green_host file
            IF (OPT('WRTGREEN') .and. myrank==master) THEN
               WRITE(58,'(2e17.9)') EZ(IE)
            ENDIF


            if (t_tgmat%gref_to_file) then
               READ (68,REC=IE) GINP
            else
               ginp(:,:,:) = t_tgmat%gref(:,:,:,ie_num)
            end if
            if(t_lloyd%dgref_to_file) then
               IF (LLY.NE.0) READ (681,REC=IE) DGINP   ! LLY Lloyd
            else
               IF (LLY.NE.0) DGINP(:,:,:) = t_lloyd%dgref(:,:,:,ie_num)
            endif

            ERYD = EZ(IE)
            NMESH = KMESH(IE)
            if(t_inc%i_write>0) WRITE (1337,'(A,I3,A,2(1X,F10.6),A,I3)') 
     &             ' ************ IE = ',IE,' ENERGY =',EZ(IE),
     &                ' KMESH = ', NMESH
C ********************************************************** I1 = 1,NREF
C -> calculate t(ll') of the reference system (on clusters)
C
#ifdef CPP_TIMING
        call timing_start('main1b - calctref13')
#endif
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
#ifdef CPP_TIMING
        call timing_pause('main1b - calctref13')
#endif    
C
C **********************************************************************
            TRALPHA(IE,ISPIN) = CZERO                                    ! LLY
            TRALPHAREF(IE) = CZERO                                       ! LLY
            DO I1 = 1,NATYP

               if (t_tgmat%tmat_to_file) then
                  IREC = IE + IELAST* (ISPIN-1) + IELAST*NSPIN* (I1-1)
                  READ (69,REC=IREC) TMAT
               else
                  irec = ie_num+ie_end*(ISPIN-1)+ie_end*NSPIN*(I1-1) 
                  tmat(:,:) = t_tgmat%tmat(:,:,irec)
               end if
               TSST(1:LMMAXD,1:LMMAXD,I1)=TMAT(1:LMMAXD,1:LMMAXD)

               IF (LLY.NE.0) THEN                                         ! LLY
                  if(t_lloyd%dtmat_to_file) then
                     IREC =IE+IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1)
                     READ (691,REC=IREC) TMAT                                ! LLY dt/dE
                  else
                     irec = ie_num+ie_end*(ISPIN-1)+ie_end*NSPIN*(I1-1)
                     TMAT(:,:) = t_lloyd%dtmat(:,:,irec)
                  end if
                  DTMATLL(1:LMMAXD,1:LMMAXD,I1) =TMAT(1:LMMAXD,1:LMMAXD)  ! LLY
                  if(t_lloyd%dtmat_to_file) then
                     IREC=IE+IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1)
                     READ (692,REC=IREC) TRALPHA1                            ! LLY
                  else
                     irec = ie_num+ie_end*(ISPIN-1)+ie_end*NSPIN*(I1-1)
                     TRALPHA1 = t_lloyd%tralpha(irec)
                  end if
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

            END DO  !i1 = 1,natyp
            
            
            if(t_lloyd%g0tr_to_file) then
               IF (LLY.NE.0.AND.ISPIN.EQ.1)                         ! LLY
     &                      READ(682,FMT='(2E24.16)') LLY_G0TR(IE)  ! LLY
            else
               IF (LLY.NE.0.AND.ISPIN.EQ.1)                         ! LLY
     &                     LLY_G0TR(IE) = t_lloyd%g0tr(ie_num)
            end if
C ------------------------------------------------------------------------
C
C --> setting up of Delta_t moved to < KLOOPZ1 >
C
            IF (OPT('readcpa ').OR.
     &        (OPT('qdos    ').AND.(IQDOSRUN.EQ.1))) THEN     ! qdos ruess: read in cpa t-matrix
               DO ISITE = 1,NAEZ                              ! qdos ruess
                  TQDOS(:,:,ISITE) = CZERO                    ! qdos ruess
#ifdef CPP_MPI
                  do lm1=1,lmmaxd
                  do lm2=1,lmmaxd
                  irec = LM2+(LM1-1)*LMMAXD+LMMAXD**2*(isite-1)+
     &                   LMMAXD**2*naez*(ie-1)+
     &                   LMMAXD**2*ielast*naez*(ispin-1)
                  read(37,rec=irec) tread
                  IF ( (LM1+LM2).NE.0 ) THEN                  ! qdos ruess
                     TQDOS(LM1,LM2,ISITE) = TREAD / CFCTORINV ! qdos ruess
                  END IF                                      ! qdos ruess
                  end do
                  end do

#else
                  READ(37,*) TEXT                             ! qdos ruess
                  READ(37,*) TEXT                             ! qdos ruess
 9921             CONTINUE                                    ! qdos ruess
                  READ(37,99013) LM1,LM2,TREAD                ! qdos ruess
99013             FORMAT (2I5,1P,2D22.14)                     ! qdos ruess
                  IF ( (LM1+LM2).NE.0 ) THEN                  ! qdos ruess
                     TQDOS(LM1,LM2,ISITE) = TREAD / CFCTORINV ! qdos ruess
                     IF ( (LM1+LM2).LT.2*LMMAXD ) GOTO 9921   ! qdos ruess
                  END IF                                      ! qdos ruess
#endif
               ENDDO                                          ! qdos ruess
            END IF                                            ! qdos ruess
C  Loop over all QDOS points and change volume for KLOOPZ run accordingly
            DO 200 IQ = 1,NQDOS                               ! qdos ruess
            IF (OPT('qdos    ')) BZKP(:,1,1) = QVEC(:,IQ)     ! qdos ruess: Set q-point x,y,z
C
#ifdef CPP_TIMING
        call timing_start('main1b - kloopz')
#endif
            CALL KLOOPZ1_QDOS(ERYD,GMATLL,INS,ALAT,IE,IGF,NSHELL,NAEZ,
     &                 NOFKS(NMESH),VOLBZ(NMESH),BZKP(1,1,NMESH),
     &                 VOLCUB(1,NMESH),CLS,NACLS,NACLSMAX,NCLS,RR,
     &                 RBASIS,EZOA,ATOM,RCLS,ICC,GINP,IDECI,
     &                 LEFTTINVLL(1,1,1,1,IE),RIGHTTINVLL(1,1,1,1,IE),
     &                 VACFLAG,NLBASIS,NRBASIS,FACTL,NATOMIMP,NSYMAT,
     &                 DSYMLL,RATOM,RROT,NSH1,NSH2,IJTABSYM,IJTABSH,
     &                 ICHECK,INVMOD,REFPOT,TREFLL,TSST,MSST,CFCTOR,  
     &                 CFCTORINV,CREL,RC,RREL,SRREL,IRREL,NRREL,DROTQ,
     &                 SYMUNITARY,KMROT,NATYP,NCPA,ICPA,ITCPAMAX,
     &                 CPATOL,NOQ,IQAT,ITOQ,CONC,IPRINT,ICPAFLAG,
     &                 ISPIN,NSPINDD,
     &                 TQDOS,IQDOSRUN,                         ! qdos
     &    DTREFLL,DTMATLL,DGINP,LLY_GRTR(IE,ISPIN),TRACET(IE,ISPIN),LLY) ! LLY Lloyd
     
#ifdef CPP_TIMING
        call timing_pause('main1b - kloopz')
#endif
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
                GMAT0(:,:)=GMATLL(:,:,I1)
                IREC = IQ + NQDOS * (IE-1) + NQDOS * IELAST *        ! qdos ruess: (without qdos, IQ=NQ=1)
     &                   (ISPIN-1) + NQDOS * IELAST * NSPIN * (I1-1) ! qdos ruess
                if (t_tgmat%gmat_to_file) then
                   WRITE (70,REC=IREC) GMAT0
                else
                   IREC = IQ + NQDOS * (ie_num-1) + NQDOS * 
     &                    ie_end * (ISPIN-1) + 
     &                    NQDOS * ie_end * NSPIN * (I1-1)
                   t_tgmat%gmat(:,:,irec) = gmat0
                end if
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


              ! writeout of host green function for impurity code for single-atom cluster (not captured in rotgll)
              IF ( NATOMIMP==1 ) THEN
                I1=ATOMIMP(1)
                IF ( OPT('KKRFLEX ') ) THEN
                  irec = ielast*(ispin-1)+ ie+1
                  ILM=0
                  GIMP=(0.e0,0.e0) !complex*8
                  DO LM2=1,LMMAXD
                      DO LM1=1,LMMAXD
                        ILM=ILM+1
                        GIMP(ILM)=GMATLL(LM1,LM2,I1)
                      ENDDO
                  ENDDO
                  irec = ielast*(ispin-1)+ ie+1
                  WRITE(888,REC=irec) GIMP
                  IF ( OPT('GPLAIN  ') ) THEN
                    WRITE(8888,'(50000E)') GIMP
                  END IF
                ENDIF ! KKRFLEX
                IF (OPT('WRTGREEN') .and. myrank==master) THEN
                  DO LM2=1,LMMAXD
                    DO LM1=1,LMMAXD
                      ILM=ILM+1
                      GIMP(ILM)=GMATLL(LM1,LM2,I1)
                    ENDDO
                ENDDO
                irec = ielast*(ispin-1)+ ie+1
                WRITE(888,REC=irec) GIMP
                IF ( OPT('GPLAIN  ') ) THEN
                  WRITE(8888,'(50000E14.7)') GIMP
                END IF
              ENDIF ! KKRFLEX
              IF (OPT('WRTGREEN') .and. myrank==master) THEN
                DO LM2=1,LMMAXD
                  DO LM1=1,LMMAXD
                    ! writeout of green_host for WRTGREEN option
                    WRITE(58,'((2I5),(2e17.9))') LM2, LM1,
     &                                     GMATLL(LM1,LM2,I1)
                    END DO
                  END DO
                ENDIF ! WRTGREEN
              END IF !( NATOMIMP==1 )
              
              IF ( LCPAIJ ) THEN
               if(t_cpa%dmatproj_to_file) then
                 DO I1 = 1,NATYP
                    GMAT0(:,:) = TSST(:,:,I1)
                    W1(:,:)    = MSST(:,:,I1)
                    IREC = IE + IELAST*(ISPIN-1) + IELAST*NSPIN*(I1-1)
                    WRITE (71,REC=IREC) GMAT0,W1
                 END DO
               else!t_cpa%dmatproj_to_file
                 irec = ie_num + ie_end*(ISPIN-1)
                 t_cpa%dmatts(:,:,:,irec) = TSST(:,:,:)
                 t_cpa%dtilts(:,:,:,irec) = MSST(:,:,:)
               end if!t_cpa%dmatproj_to_file
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

#ifdef CPP_MPI
            !stop timing measurement for this ie, needed for MPIadapt
            if(MPIadapt.and.t_mpi_c_grid%myrank_ie==0) then
              call timing_stop('time_1b_ie', save_out=timings_1b(ie) )
            end if
#endif

 360     CONTINUE               ! IE = 1,IELAST
#ifdef CPP_TIMING
        if(.not.OPT('GREENIMP')) then
           if(t_inc%i_time>0) call timing_stop('main1b - calctref13')
           if(t_inc%i_time>0) call timing_pause('main1b_fourier')
           if(t_inc%i_time>0) call timing_stop('main1b_inversion')
           if(t_inc%i_time>0) call timing_stop('main1b - kloopz')
        endif
#endif
!EEEEEEEEEEEEEEEEEEEEEEEEE  NO-SOC EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

         IF( NCPAFAIL .NE. 0 ) THEN
            if(t_inc%i_write>0) then
            WRITE(1337,*)
            WRITE(1337,'(1X,79(''*''),/)')
            WRITE(1337,99019) CPATOL, NCPAFAIL,
     &           (IECPAFAIL(IE),DBLE(EZ(IECPAFAIL(IE))),
     &           IE=1,NCPAFAIL)
            WRITE(1337,'(1X,79(''*''),/)')
            WRITE(1337,*)
            endif
         ELSE
            IF( NCPA .NE. 0 ) THEN
               if(t_inc%i_write>0) then
               WRITE(1337,*)
               WRITE(1337,99020)
               WRITE(1337,*)
               endif
            END IF
         END IF
C         
 370  CONTINUE                  !  ISPIN = 1,NSPIN
      

      ELSE ! NEW SOC SOLVER 

      ! nonco angles
      call read_angles(t_params,NATYP,THETAS,PHIS)

!EEEEEEEEEEEEEEEEEEEEEEEEEEEE  SOC  EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#ifdef CPP_MPI
      ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
      ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)

      DO 460 ie_num=1,ie_end
         IE = ie_start+ie_num
   
         !start timing measurement for this ie, needed for MPIadapt
         if(MPIadapt.and.t_mpi_c_grid%myrank_ie==0) then
           call timing_start('time_1b_ie')
         end if

#else
      ie_end = IELAST
      DO 460 IE=1,IELAST
         ie_num = ie
#endif

       ! write energy into green_host file
       IF (OPT('WRTGREEN') .and. myrank==master) THEN
          WRITE(58,'(2e17.9)') EZ(IE)
       ENDIF
       
c read in Green function of reference system
       if (t_tgmat%gref_to_file) then
          READ (68,REC=IE) GINP
       else
          ginp(:,:,:) = t_tgmat%gref(:,:,:,ie_num)
       end if
       if(t_lloyd%dgref_to_file) then
          IF (LLY.NE.0) READ (681,REC=IE) DGINP   ! LLY Lloyd
       else
          IF (LLY.NE.0) DGINP(:,:,:) = t_lloyd%dgref(:,:,:,ie_num)
       endif
       ERYD = EZ(IE)
       NMESH = KMESH(IE)
       if(t_inc%i_write>0) WRITE (1337,'(A,I3,A,2(1X,F10.6),A,I3)') 
     &             ' ************ IE = ',IE,' ENERGY =',EZ(IE),
     &                ' KMESH = ', NMESH

c construct t matrix for reference system (now is always double matrix)
#ifdef CPP_TIMING
       call timing_start('main1b - calctref13')
#endif
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
#ifdef CPP_TIMING
       call timing_pause('main1b - calctref13')
#endif
       TRALPHA(IE,NSPINDD)=CZERO
       TRALPHAREF(IE)=CZERO
c read in t matrix

       DO I1 = 1,NATYPD

! c read in theta and phi for noncolinear
        THETA = THETAS(I1)
        PHI   = PHIS(I1)

c read in t-matrix from file
        IREC = IE + IELAST*(I1-1)
        if (t_tgmat%tmat_to_file) then
           READ (69,REC=IREC) TMAT
        else
           irec = ie_num + ie_end * (i1-1)
           tmat(:,:) = t_tgmat%tmat(:,:,irec)
        end if
      

c rotate t-matrix from local to global frame
        CALL ROTATEMATRIX(TMAT,THETA,PHI,LMGF0D,0)

         DO LM1=1,LMMAXD
          DO LM2=1,LMMAXD
           TSST(LM1,LM2,I1)=TMAT(LM1,LM2)
          ENDDO
         ENDDO

       IF (LLY.NE.0) THEN
        IREC = IE + IELAST*(I1-1)
        if(t_lloyd%dtmat_to_file) then
           READ(691,REC=IREC) TMAT ! LLY
        else
           irec = ie_num + ie_end * (i1-1)
           TMAT(:,:) = t_lloyd%dtmat(:,:,irec)
        end if
        CALL ROTATEMATRIX(TMAT,THETA,PHI,LMGF0D,0) ! LLY
         DO LM1=1,LMMAXD
          DO LM2=1,LMMAXD
           DTMATLL(LM1,LM2,I1)=TMAT(LM1,LM2) ! LLY
          ENDDO
         ENDDO
        IREC = IE + IELAST*(I1-1)
        if(t_lloyd%dtmat_to_file) then
           READ(692,REC=IREC) TRALPHA1
        else
           irec = ie_num + ie_end * (i1-1)
           TRALPHA1 = t_lloyd%tralpha(irec)
        end if
         
         TRALPHA(IE,NSPINDD)=TRALPHA(IE,NSPINDD)+TRALPHA1 ! LLY
         TRALPHA1=CZERO
          DO L1=0,LMAX
           TRALPHA1=TRALPHA1+(2*L1+1)*
     &              DALPHAREF(L1,REFPOT(I1))/ALPHAREF(L1,REFPOT(I1)) ! LLY
          ENDDO
          TRALPHAREF(IE)=TRALPHAREF(IE)+TRALPHA1 ! LLY
       ENDIF ! LLY
       ENDDO ! I1
       if(t_lloyd%g0tr_to_file) then
          IF (LLY.NE.0) READ(682,FMT='(2E24.16)') LLY_G0TR(IE)  ! LLY
       else
          IF (LLY.NE.0) LLY_G0TR(IE) = t_lloyd%g0tr(ie_num)
       end if

c QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS
            IF (OPT('readcpa ').OR.
     &        (OPT('qdos    ').AND.(IQDOSRUN.EQ.1))) THEN     ! qdos ruess: read in cpa t-matrix
               DO ISITE = 1,NAEZ                              ! qdos ruess
                  TQDOS(:,:,ISITE) = CZERO                    ! qdos ruess
#ifdef CPP_MPI
                  do lm1=1,lmmaxd
                  do lm2=1,lmmaxd
                  irec = LM2+(LM1-1)*LMMAXD+LMMAXD**2*(isite-1)+
     &                   LMMAXD**2*naez*(ie-1)
                  read(37,rec=irec) tread
                  IF ( (LM1+LM2).NE.0 ) THEN                  ! qdos ruess
                     TQDOS(LM1,LM2,ISITE) = TREAD / CFCTORINV ! qdos ruess
                  END IF                                      ! qdos ruess
                  end do
                  end do

#else
                  READ(37,*) TEXT                             ! qdos ruess
                  READ(37,*) TEXT                             ! qdos ruess
 9920             CONTINUE                                    ! qdos ruess
                  READ(37,99014) LM1,LM2,TREAD                ! qdos ruess
99014             FORMAT (2I5,1P,2D22.14)                     ! qdos ruess
                  IF ( (LM1+LM2).NE.0 ) THEN                  ! qdos ruess
                     TQDOS(LM1,LM2,ISITE) = TREAD / CFCTORINV ! qdos ruess
                     IF ( (LM1+LM2).LT.2*LMMAXD ) GOTO 9920   ! qdos ruess
                  END IF                                      ! qdos ruess
#endif
               ENDDO                                          ! qdos ruess
            END IF                                            ! qdos ruess
C  Loop over all QDOS points and change volume for KLOOPZ run accordingly
            DO 220 IQ = 1,NQDOS                               ! qdos ruess
            IF (OPT('qdos    ')) BZKP(:,1,1) = QVEC(:,IQ)     ! qdos ruess: Set q-point x,y,z

#ifdef CPP_TIMING
            call timing_start('main1b - kloopz')
#endif
            CALL KLOOPZ1_QDOS(ERYD,GMATLL,INS,ALAT,IE,IGF,NSHELL,NAEZ,
     &              NOFKS(NMESH),VOLBZ(NMESH),BZKP(1,1,NMESH),
     &              VOLCUB(1,NMESH),CLS,NACLS,NACLSMAX,NCLS,RR,RBASIS,
     &              EZOA,ATOM,RCLS,ICC,GINP,IDECI,
     &              LEFTTINVLL(1,1,1,1,IE),RIGHTTINVLL(1,1,1,1,IE),
     &              VACFLAG,NLBASIS,NRBASIS,FACTL,NATOMIMP,NSYMAT,
     &              DSYMLL,RATOM,RROT,NSH1,NSH2,IJTABSYM,IJTABSH,
     &              ICHECK,INVMOD,REFPOT,TREFLL,TSST,MSST,CFCTOR,
     &              CFCTORINV,CREL,RC,RREL,SRREL,IRREL,NRREL,DROTQ,
     &              SYMUNITARY,KMROT,NATYP,NCPA,ICPA,ITCPAMAX,
     &              CPATOL,NOQ,IQAT,ITOQ,CONC,IPRINT,ICPAFLAG,
     &              1,NSPINDD,
     &              TQDOS,IQDOSRUN,     ! qdos
     &    DTREFLL,DTMATLL,DGINP,LLY_GRTR(IE,1),TRACET(IE,1),LLY) ! LLY Lloyd     
#ifdef CPP_TIMING
            call timing_pause('main1b - kloopz')
#endif
          
!           Skip this part if first part of the qdos is running
            IF ( .NOT.(OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) ) THEN
               IF (NCPA.NE.0) THEN
                IF (ICPAFLAG .NE. 0) THEN
                   NCPAFAIL = NCPAFAIL + 1
                   IECPAFAIL(NCPAFAIL)= IE
                 END IF
               END IF  ! (NCPA.NE.0)

               DO I1 = 1,NSHELL(0)
                  GMAT0(1:LMMAXD,1:LMMAXD) =GMATLL(1:LMMAXD,1:LMMAXD,I1)
                  IREC = IQ + NQDOS * (IE-1) + NQDOS * IELAST * (I1-1) ! qdos ruess
                  if (t_tgmat%gmat_to_file) then
                     WRITE (70,REC=IREC) GMAT0
                     ! human readable writeout if test option is hit
                     if(test('fileverb')) then
                        write(707070,'(i9,200000F15.7)') irec, gmat0
                     end if
                  else
                     IREC = IQ + NQDOS * (ie_num-1) + NQDOS * 
     &                    ie_end * (I1-1)
                     t_tgmat%gmat(:,:,irec) = gmat0
                  end if
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

               IF ( NATOMIMP==1 ) THEN
                  I1=ATOMIMP(1)
                  IF(OPT('KKRFLEX')) THEN
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
                  IF (OPT('WRTGREEN') .and. myrank==master) THEN
                    DO LM2=1,LMMAXD
                      DO LM1=1,LMMAXD
                        ! writeout of green_host for WRTGREEN option
                        WRITE(58,'((2I5),(2e17.9))') LM2, LM1,
     &                                         GMATLL(LM1,LM2,I1)
                      END DO
                    END DO
                  ENDIF ! WRTGREEN
               ENDIF

               IF ( LCPAIJ ) THEN
                if(t_cpa%dmatproj_to_file) then
                  DO I1 = 1,NATYP
                     DO LM2=1,LMMAXD
                        DO LM1=1,LMMAXD
                           GMAT0(LM1,LM2) = TSST(LM1,LM2,I1)
                           W1(LM1,LM2)    = MSST(LM1,LM2,I1)
                        END DO
                     END DO
                     IREC = IE + IELAST*(I1-1)
                     WRITE (71,REC=IREC) GMAT0,W1
                  END DO!I1
                else!t_cpa%dmatproj_to_file
                  irec = ie_num
                  t_cpa%dmatts(:,:,:,irec) = TSST(:,:,:)
                  t_cpa%dtilts(:,:,:,irec) = MSST(:,:,:)
                end if!t_cpa%dmatproj_to_file
               END IF  ! ( LCPAIJ )

            ENDIF               ! ( .NOT.(OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) )
 220     CONTINUE               ! IQ = 1,NQ                                         ! qdos ruess

         IF (LLY.NE.0) THEN                                      ! LLY 

            CDOS_LLY(IE,1) =   TRALPHA(IE,1)                     ! LLY
     &             - LLY_GRTR(IE,1) / VOLBZ(1) + 2d0*LLY_G0TR(IE)    ! LLY 

            CDOS_LLY(IE,1) = CDOS_LLY(IE,1) / PI                 ! LLY 

            CDOSREF_LLY(IE) = TRALPHAREF(IE) - LLY_G0TR(IE)      ! LLY 


         ENDIF                                                   ! LLY 
         

#ifdef CPP_MPI
         !stop timing measurement for this ie, needed for MPIadapt
         if(MPIadapt.and.t_mpi_c_grid%myrank_ie==0) then
            call timing_stop('time_1b_ie', save_out=timings_1b(ie) )
         end if
#endif

 460  CONTINUE                  ! IE=1,IELAST
#ifdef CPP_TIMING
        if(.not.OPT('GREENIMP')) then
           if(t_inc%i_time>0) call timing_stop('main1b - calctref13')
           if(t_inc%i_time>0) call timing_pause('main1b_fourier')
           if(t_inc%i_time>0) call timing_stop('main1b_inversion')
           if(t_inc%i_time>0) call timing_stop('main1b - kloopz')
        endif
#endif
!EEEEEEEEEEEEEEEEEEEEEEEEEEEE  SOC  EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

         IF( NCPAFAIL .NE. 0 ) THEN
            if(t_inc%i_write>0) then
            WRITE(1337,*)
            WRITE(1337,'(1X,79(''*''),/)')
            WRITE(1337,99019) CPATOL, NCPAFAIL,
     &           (IECPAFAIL(IE),DBLE(EZ(IECPAFAIL(IE))),
     &           IE=1,NCPAFAIL)
            WRITE(1337,'(1X,79(''*''),/)')
            WRITE(1337,*)
            endif
         ELSE
            IF( NCPA .NE. 0 ) THEN
               if(t_inc%i_write>0) then
               WRITE(1337,*)
               WRITE(1337,99020)
               WRITE(1337,*)
               endif
            END IF
         END IF

      ENDIF ! NEWSOSOL


C
C    ----------------------------------------------------------------
C    |           END of do loop over spins and energies             |
C    ----------------------------------------------------------------

      ! close green_host file after write out
      IF (OPT('WRTGREEN') .and. myrank==master) THEN
         CLOSE(58)
      ENDIF
C
      CLOSE (68)
c     IF ( IGF.NE.0 ) CLOSE (88)  !no-green

      IF (LLY.NE.0) THEN                                                 ! LLY Lloyd
         if(t_lloyd%cdos_diff_lly_to_file) then
            OPEN (701,FILE='cdosdiff_lly.dat',FORM='FORMATTED')          ! LLY Lloyd
         end if
         
         IF (.NOT.OPT('NEWSOSOL')) THEN
          DO ISPIN = 1,NSPIN                                             ! LLY
#ifdef CPP_MPI
         ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
         ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
           DO ie_num=1,ie_end
            IE = ie_start+ie_num
#else
           DO IE = 1,IELAST                                              ! LLY
            ie_num = ie
#endif
            if(t_lloyd%cdos_diff_lly_to_file) then
               WRITE(701,FMT='(10E25.16)') EZ(IE),
     &          CDOS_LLY(IE,ISPIN),
     &          TRALPHA(IE,ISPIN),LLY_GRTR(IE,ISPIN)                     ! LLY
            else
               t_lloyd%cdos(ie_num,ISPIN) = CDOS_LLY(IE,ISPIN)
            end if
           ENDDO                                                         ! LLY
          ENDDO                                                          ! LLY
         ELSE                         
#ifdef CPP_MPI
         ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
         ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
           DO ie_num=1,ie_end
            IE = ie_start+ie_num
#else
           DO IE = 1,IELAST                                              ! LLY
            ie_num = ie
#endif
            if(t_lloyd%cdos_diff_lly_to_file) then
              WRITE(701,FMT='(10E25.16)') EZ(IE),CDOS_LLY(IE,1),     
     &                                    TRALPHA(IE,1),LLY_GRTR(IE,1)   ! LLY
            else
               t_lloyd%cdos(ie_num,1) = CDOS_LLY(IE,1)
            end if
           
          ENDDO                                                          ! LLY
         ENDIF ! .NOT.OPT('NEWSOSOL')                                    ! LLY
         if(t_lloyd%cdos_diff_lly_to_file) CLOSE(701)                    ! LLY
      ENDIF                                                              ! LLY
C
      IF ( ( OPT('XCPL    ') ).AND.( ICC.LE.0 ) ) THEN
#ifdef CPP_TIMING
         call timing_start('main1b - tbxccpl')
#endif
         IF(NQDOS/=1) stop 'QDOS option not compatible with XCPL'
         IF (.NOT.OPT('NEWSOSOL')) THEN
            CALL TBXCCPLJIJ(69,IELAST,EZ,WEZ,NSPINDD,NCPA,NAEZ,NATYP,NOQ
     &                     ,ITOQ,IQAT,NSHELL,NATOMIMP,ATOMIMP,RATOM,
     &                      NOFGIJ,NQCALC,IQCALC,IJTABCALC,
     &                      IJTABSYM,IJTABSH,ISH,JSH,
     &                      DSYMLL,IPRINT,NATYPD,NSHELD,LMMAXD,
     &                      NPOL)
         ELSE !.NOT.OPT('NEWSOSOL'))
           CALL TBXCCPLJIJDIJ(naezd,natypd,lmmaxd,lmgf0d,natomimpd,
     &                        nsymaxd,iemxd,thetas,phis,
     &                natomimp,atomimp,nofgijD,iqat,rclsimp,
     &                ijtabcalc,ijtabcalc_I,ijtabsh,
     &                ijtabsym,
     &                ielast,ez,wez,NPOL,dsymll,noq,itoq,ncpa)
         END IF
#ifdef CPP_TIMING
         call timing_stop('main1b - tbxccpl')
#endif
      END IF
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
      CLOSE (69)
      CLOSE (70)
      IF (OPT('FERMIOUT') .and. myrank==master) CLOSE(6801)  !fswrt
      IF ( LCPAIJ .and. t_cpa%dmatproj_to_file ) CLOSE (71)

      CLOSE(37)                                    ! qdos ruess
C     Finished first qdos run. Now re-run the whole kkr1b program to 
C     calculate the GF for every energy (defined in inputcard) and 
C     kpoint (defined in qvec.dat)
      IQDOSRUN=IQDOSRUN+1                          ! qdos ruess
      if(t_lloyd%dgref_to_file)   CLOSE(681)
      if(t_lloyd%g0tr_to_file)    CLOSE(682)
      if(t_lloyd%dtmat_to_file)   CLOSE(691)
      if(t_lloyd%tralpha_to_file) CLOSE(692)
      IF (IQDOSRUN.EQ.1) GO TO 210                 ! qdos ruess



      IF(OPT('OPERATOR') .and. myrank==master) THEN
        ! check if impurity files are present (otherwise no imp.
        ! wavefunctions can be calculated)
        operator_imp = .true.
        inquire(file='potential_imp', exist=lexist)
        if (.not.lexist) operator_imp = .false.
        inquire(file='shapefun_imp', exist=lexist)
        if (.not.lexist) operator_imp = .false.
        inquire(file='scoef', exist=lexist)
        if (.not.lexist) operator_imp = .false.
      ELSE
        operator_imp = .false.
      ENDIF
#ifdef CPP_MPI
      call MPI_Bcast(operator_imp, 1, MPI_LOGICAL, master, 
     &               MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error broadcasting operator_imp'
#endif

      ! Do stuff for GREENIMP option (previously done in zulapi code)
      ! run for OPERATOR option to precalculate impurity wavefunctions
      ! that are then stored in t_imp (used only if potential_imp,
      ! scoef, shapefun_imp)
      IF (OPT('GREENIMP') .or. operator_imp) THEN

        ! consistency checks
        if(.not.(ielast==1.or.ielast==3)) stop 
     + 'Error: GREENIMP option only possible with 1 () or 3 () energy
     +  points in contour'
        if(ielast==1 .and.dabs(dimag(ez(1)))>10**-10) stop 
     + 'Error: T>0 for GREENIMP (DTMTRX writeout, IELAST==3)'
        if(ielast==3 .and.dabs(dimag(ez(1)))<10**-10) stop 
     + 'Error: T==0 for GREENIMP (GMATLL_GES writeout, IELAST==3)'
        ! end consistency checks

#ifdef CPP_MPI
        ! init arrays and communicate parameters of t_imp for all ranks
        ! that are not the master
        call bcast_t_imp_scalars(t_imp)
        if(myrank/=master) call init_t_imp(t_inc,t_imp)
        call bcast_t_imp_arrays(t_imp, t_inc)
#endif
        
        do ie=1,ielast ! big ie loop (use only for GMATLL output)

          IF (OPT('GREENIMP')) THEN
             IF (IELAST.EQ.3 .and. ie==1 .and. myrank==master) then
                OPEN(UNIT=59,FILE='GMATLL_GES',FORM='FORMATTED')
                OPEN(UNIT=60,FILE='green_host',FORM='FORMATTED')
             ENDIF
          ENDIF ! OPT('GREENIMP')
          ALLOCATE(DTMTRX(LMMAXD*t_imp%NATOMIMP,LMMAXD*t_imp%NATOMIMP), 
     &             stat=ierr)
          if(ierr/=0) stop 'Error allocating DTMTRX'
          DTMTRX=CZERO
 
         ! find DTMTRX (written out for IELAST==1), parallelized with
         ! mpi over atoms
          CALL TMATIMP_NEWSOLVER(t_params%NSPIN,
     &       t_params%LMAX,t_params%R,t_params%ZAT,t_params%IELAST,
     +       EZ(IE),NSRA-1,t_params%IRWS,t_params%IPAN,
     &       t_params%IRCUT,t_params%IRMIN,CVLIGHT,
     &       t_params%CLEB(1,1),t_params%ICLEB,t_params%IEND,
     +       t_params%NPAN_LOG,t_params%NPAN_EQ,t_params%NCHEB,
     &       t_params%R_LOG,t_params%VINS,t_params%VISP,t_imp%NATOMIMP,
     +       t_params%RCLSIMP,t_params%ATOMIMP(1:t_imp%NATOMIMP),
     +       t_imp%IHOST,t_imp%HOSTIMP(1:t_imp%NATOMIMP),
     +       t_imp%RIMP(1:t_params%IRMD,1:t_imp%NATOMIMP),
     +       t_imp%ZIMP(1:t_imp%NATOMIMP),
     +       t_imp%IRWSIMP(1:t_imp%NATOMIMP),
     +       t_imp%IPANIMP(1:t_imp%NATOMIMP),
     &       t_imp%IRCUTIMP(0:t_params%IPAND,1:t_imp%NATOMIMP),
     +       t_imp%IRMINIMP(1:t_imp%NATOMIMP),t_imp%VINSIMP,
     +       t_imp%VISPIMP,DTMTRX)

         ! compute GMATLL_GES, on master rank only
          IF (OPT('GREENIMP').and.IELAST.EQ.3.and.myrank==master) THEN
           CALL GREENIMP(t_imp%NATOMIMP,DTMTRX,t_params%EZ(IE))
          ENDIF
          
          DEALLOCATE(DTMTRX, stat=ierr)
          if(ierr/=0) stop 'Error deallocating DTMTRX'

        end do ! ie-loop

        ! done with GREENIMP option, stopping now
        if(.not. OPT('OPERATOR')) then
          if(myrank==master) write(*,*) 'done with GREENIMP, stop here!'
#ifdef CPP_MPI
          call MPI_FINALIZE(ierr)
#endif
          stop
        end if

      ENDIF ! GREENIMP .or. OPERATOR

! ------------------------------------------------------------------------
! determine the spin operator, torque operator and spin flux operator
! used in FScode do compute spin expectation values etc. within Boltzmann
! formalism
! ------------------------------------------------------------------------
      if (OPT('OPERATOR')) then
#ifdef CPP_TIMING
        call timing_start('main1b - operator')
#endif

        call operators_for_FScode(KORBIT, operator_imp)

#ifdef CPP_TIMING
        call timing_stop('main1b - operator')
#endif
      end if ! OPERATOR


C ----------------------------------------------------------------------
C
      if( test('rhoqtest') .and. (myrank==master) ) then
        open(9999, file='params.txt')
        write(9999,*) 2*LMMAXD, t_params%natyp
        write(9999,*) t_params%naez, t_params%nclsd, t_params%nrd, 
     &                t_params%nembd1, t_params%lmax
        write(9999,*) t_params%alat, naclsmax
        close(9999)

        open(9999, file='host.txt')
        write(9999,*) t_params%rbasis(1:3,1:t_params%natyp)
        write(9999,*) t_params%rcls(1:3,1:t_params%nclsd,
     &                              1:t_params%nclsd), 
     &                t_params%rr(1:3,0:t_params%nrd), 
     &                t_params%atom(1:t_params%nclsd,
     &                              1:t_params%naez+t_params%nembd1)
        write(9999,*) t_params%cls(1:t_params%naez+t_params%nembd1), 
     &                t_params%ezoa(1:t_params%nclsd,
     &                              1:t_params%naez+t_params%nembd1), 
     &                t_params%nacls(1:t_params%nclsd)
        close(9999)
      endif


      ! clean up allocations before exiting
      deallocate(LEFTTINVLL, RIGHTTINVLL, BZKP, stat=ierr)
      if(ierr/=0) stop '[main1b] Error deallocating leftinvll etc.'
      if(opt('qdos    ')) then
        deallocate(qvec, stat=ierr)
        if(ierr/=0) stop '[main1b] Error deallocating qvec after main1b'
      end if
      deallocate(ISH, JSH, stat=ierr)
      if(ierr/=0) stop '[main1b] Error deallocating ISH and JSH'



      if(t_inc%i_write>0) 
     &   WRITE (1337,'(79(1H=),/,30X,"< KKR1b finished >",/,79(1H=),/)')


 2100 FORMAT(/,79(1H=),/,5X,' Inversion algorithm used : ',A,/,
     &     79(1H=),/)
99019 FORMAT (/,1X,79('*'),/,
     &             ' tolerance for CPA-cycle:',F15.7,/,
     &             ' CPA not converged for',I3,' energies:',/,
     &             3(' E:',I3,F7.4,:,2X))
99020 FORMAT (/,1X,79('*'),/,25X,'no problems with',
     &         '  CPA-cycle ',/,1X,79('*'),/)



      END SUBROUTINE !main1b
      
      END MODULE
