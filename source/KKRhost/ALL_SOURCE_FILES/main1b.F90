!-------------------------------------------------------------------------------
! MODULE: MOD_MAIN1B
!> @brief Module concerning gmat
!> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
!> and many others ...
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
module MOD_MAIN1B

   use Profiling
   use Constants
   use global_variables

   implicit none

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: main1b
   !> @brief Main subroutine regarding the claculation of the gmat
   !> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
   !> and many others ...
   !----------------------------------------------------------------------------
   subroutine main1b(NR,IRM,LLY,ICC,IGF,INS,NPOL,LMAX,NREF,NSRA,NCLS,NCPA,NEMB,  &
      NAEZ,NATYP,NCLSD,NSPIN,KMROT,LMMAXD,IELAST,INVMOD,NSYMAT,NEMBD1,LMGF0D,    &
      NOFGIJ,NQCALC,NSPINDD,NLBASIS,NRBASIS,MAXMESH,NATOMIMP,ITCPAMAX,ALAT,      &
      CPATOL,NOQ,CLS,IQAT,NSH1,NSH2,ICPA,NACLS,NSHELL,REFPOT,ATOMIMP,EZOA,ATOM,  &
      KAOEZ,ICHECK,CONC,RMTREF,RR,RATOM,RBASIS,RROT,RCLS,SYMUNITARY,KMESH,IQCALC,&
      IJTABSH,IJTABSYM,IJTABCALC,IJTABCALC_I,NRREL,IRREL,VREF,RCLSIMP,RC,RREL,   &
      CREL,SRREL,DROTQ,DSYMLL,VACFLAG)

      use mod_types, only: t_tgmat,t_inc,t_lloyd,t_cpa,init_t_cpa,t_imp
#ifdef CPP_MPI
      use mod_types, only: t_mpi_c_grid, save_t_mpi_c_grid,          &
                           get_ntot_pT_ioff_pT_2D,init_params_t_imp, &
                           init_t_imp,bcast_t_imp_scalars, bcast_t_imp_arrays
      use mpi
#endif
      use mod_mympi, only: myrank, master
#ifdef CPP_MPI
      use mod_mympi, only: find_dims_2d,distribute_linear_on_tasks, MPIadapt
#endif
      use mod_timing
      use mod_wunfiles
      use mod_tbxccpljij, only: tbxccpljij
      use mod_version_info
      use global_variables
      use mod_tbxccpljijdij, only: tbxccpljijdij


      implicit none
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! For KREL = 1 (relativistic mode)
      !
      !  NPOTD = 2 * NATYP
      !  LMMAXD = 2 * (LMAX+1)^2
      !  NSPIND = 1
      !  LMGF0D = (LMAX+1)^2 dimension of the reference system Green
      !          function, set up in the spin-independent non-relativstic
      !          (l,m_l)-representation
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !    ..
      ! .. Input variables
      integer, intent(inout) :: NR           !< Number of real space vectors rr
      integer, intent(inout) :: IRM          !< Maximum number of radial points
      integer, intent(inout) :: LLY          !< LLY <> 0 : apply Lloyd's formula
      integer, intent(inout) :: ICC          !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
      integer, intent(inout) :: IGF          !< Do not print or print (0/1) the KKRFLEX_* files
      integer, intent(inout) :: INS          !< 0 (MT), 1(ASA), 2(Full Potential)
      integer, intent(inout) :: NPOL         !< Number of Matsubara Poles (EMESHT)
      integer, intent(inout) :: LMAX         !< Maximum l component in wave function expansion
      integer, intent(inout) :: NREF         !< Number of diff. ref. potentials
      integer, intent(inout) :: NSRA
      integer, intent(inout) :: NCLS         !< Number of reference clusters
      integer, intent(inout) :: NCPA         !< NCPA = 0/1 CPA flag
      integer, intent(inout) :: NEMB         !< Number of 'embedding' positions
      integer, intent(inout) :: NAEZ         !< Number of atoms in unit cell
      integer, intent(inout) :: NATYP        !< Number of kinds of atoms in unit cell
      integer, intent(inout) :: NCLSD        !< Maximum number of different TB-clusters
      integer, intent(inout) :: NSPIN        !< Counter for spin directions
      integer, intent(inout) :: KMROT        !< 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
      integer, intent(inout) :: LMMAXD       !< (KREL+KORBIT+1)(LMAX+1)^2
      integer, intent(inout) :: IELAST
      integer, intent(inout) :: INVMOD       !< Inversion scheme
      integer, intent(inout) :: NSYMAT
      integer, intent(inout) :: NEMBD1       !< NEMB+1
      integer, intent(inout) :: LMGF0D       !< (LMAX+1)**2
      integer, intent(inout) :: NOFGIJ       !< number of GF pairs IJ to be calculated as determined from IJTABCALC<>0
      integer, intent(inout) :: NQCALC
      integer, intent(inout) :: NSPINDD      !< NSPIND-KORBIT
      integer, intent(inout) :: NLBASIS      !< Number of basis layers of left host (repeated units)
      integer, intent(inout) :: NRBASIS      !< Number of basis layers of right host (repeated units)
      integer, intent(inout) :: MAXMESH
      integer, intent(inout) :: NATOMIMP     !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      integer, intent(inout) :: ITCPAMAX     !< Max. number of CPA iterations
      double precision, intent(inout) :: ALAT      !< Lattice constant in a.u.
      double precision, intent(inout) :: CPATOL    !< Convergency tolerance for CPA-cycle
      integer, dimension(NAEZ), intent(inout)      :: NOQ         !< Number of diff. atom types located
      integer, dimension(NAEZ+NEMB), intent(inout) :: CLS         !< Cluster around atomic sites
      integer, dimension(NATYP), intent(inout)     :: IQAT        !< The site on which an atom is located on a given site
      integer, dimension(NSHELD), intent(inout)    :: NSH1        !< Corresponding index of the sites I/J in (NSH1/2) in the unit cell in a shell
      integer, dimension(NSHELD), intent(inout)    :: NSH2        !< Corresponding index of the sites I/J in (NSH1/2) in the unit cell in a shell
      integer, dimension(NAEZ), intent(inout)      :: ICPA        !< ICPA = 0/1 site-dependent CPA flag
      integer, dimension(NCLSD), intent(inout)     :: NACLS       !< Number of atoms in cluster
      integer, dimension(0:NSHELD), intent(inout)  :: NSHELL      !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
      integer, dimension(NAEZ+NEMB), intent(inout) :: REFPOT      !< Ref. pot. card  at position
      integer, dimension(NATOMIMP), intent(inout)  :: ATOMIMP
      integer, dimension(NACLSD,NAEZ+NEMB), intent(inout)            :: EZOA  !< EZ of atom at site in cluster
      integer, dimension(NACLSD,NAEZ+NEMB), intent(inout)            :: ATOM  !< Atom at site in cluster
      integer, dimension(NATYP,NAEZ+NEMB), intent(inout)             :: KAOEZ !< Kind of atom at site in elem. cell
      integer, dimension(NAEZ/NPRINCD,NAEZ/NPRINCD), intent(inout)   :: ICHECK
      double precision, dimension(NATYP), intent(inout)              :: CONC     !< Concentration of a given atom
      double precision, dimension(NREF), intent(inout)               :: RMTREF   !< Muffin-tin radius of reference system
      double precision, dimension(3,0:NR), intent(inout)             :: RR       !< Set of real space vectors (in a.u.)
      double precision, dimension(3,NSHELD), intent(inout)           :: RATOM
      double precision, dimension(3,NAEZ+NEMB), intent(inout)        :: RBASIS   !< Position of atoms in the unit cell in units of bravais vectors
      double precision, dimension(48,3,NSHELD), intent(inout)        :: RROT
      double precision, dimension(3,NACLSD,NCLSD), intent(inout)     :: RCLS  !< Real space position of atom in cluster
      logical, dimension(2), intent(inout) :: VACFLAG
      logical, dimension(NSYMAXD), intent(inout) :: SYMUNITARY !< unitary/antiunitary symmetry flag
      integer, dimension(IEMXD), intent(inout)  :: KMESH
      integer, dimension(NAEZ), intent(inout)   :: IQCALC
      integer, dimension(NOFGIJ), intent(inout) :: IJTABSH    !< Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
      integer, dimension(NOFGIJ), intent(inout) :: IJTABSYM   !< Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
      integer, dimension(NOFGIJ), intent(inout) :: IJTABCALC  !< Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
      integer, dimension(NOFGIJ), intent(inout) :: IJTABCALC_I
      integer, dimension(2,LMMAXD), intent(inout)        :: NRREL
      integer, dimension(2,2,LMMAXD), intent(inout)      :: IRREL
      double precision, dimension(NREF), intent(inout)         :: VREF
      double precision, dimension(3,NATOMIMPD), intent(inout)  :: RCLSIMP

      double complex, dimension(LMMAXD,LMMAXD), intent(inout) :: RC    !< NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
      double complex, dimension(LMMAXD,LMMAXD), intent(inout) :: RREL  !< Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
      double complex, dimension(LMMAXD,LMMAXD), intent(inout) :: CREL  !< Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
      double complex, dimension(2,2,LMMAXD), intent(inout)              :: SRREL
      double complex, dimension(LMMAXD,LMMAXD,NAEZ), intent(inout)      :: DROTQ !< Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity
      double complex, dimension(LMMAXD,LMMAXD,NSYMAXD), intent(inout)   :: DSYMLL
      ! ..
      ! .. Local variables
      integer :: L
      integer :: I
      integer :: I1
      integer :: IE
      integer :: IQ
      integer :: IX
      integer :: IC
      integer :: L1
      integer :: ILM
      integer :: LM1
      integer :: LM2
      integer :: IREC
      integer :: ILTMP
      integer :: NMESH
      integer :: NQDOS ! qdos ruess:number of qdos points
      integer :: ISITE ! qdos ruess
      integer :: IDECI
      integer :: ISPIN
      integer :: IPRINT
      integer :: ITMPDIR
      integer :: LRECGRF
      integer :: LRECTMT
      integer :: LRECTRA   ! LLY Lloyd
      integer :: IQDOSRUN  ! qdos ruess: counter to organise qdos run
      integer :: NACLSMAX
      integer :: LRECGRF1
      integer :: NCPAFAIL
      integer :: ICPAFLAG
      integer :: RECLENGTH
      integer :: LRECGREEN
      double precision :: PHI
      double precision :: THETA
      double precision :: RFCTOR !< rfctor=a/(2*pi) conversion factor to p.u.
      double complex :: ERYD
      double complex :: TREAD ! qdos ruess
      double complex :: CFCTOR
      double complex :: TRALPHA1 ! LLY Lloyd
      double complex :: CFCTORINV
      logical :: OPT
      logical :: TEST
      logical :: LCPAIJ

      character(len=80) :: TMPDIR
#ifndef CPP_MPI
      character(len=80) :: TEXT !qdos ruess
#endif
      ! .. Local arrays
      integer, dimension(MAXMSHD)   :: NOFKS
      integer, dimension(IEMXD)     :: IECPAFAIL
      integer, dimension(NATYP,NAEZ) :: ITOQ
      double precision, dimension(NATYP)     :: PHIS
      double precision, dimension(MAXMSHD)   :: VOLBZ
      double precision, dimension(NATYP)     :: THETAS
      double precision, dimension(KPOIBZ,MAXMSHD) :: VOLCUB
      double complex, dimension(IEMXD) :: EZ
      double complex, dimension(IEMXD) :: WEZ
      double complex, dimension(IEMXD) :: LLY_G0TR             !< LLY Lloyd  Trace[ X ], Eq.5.27 PhD Thiess
      double complex, dimension(IEMXD) :: TRALPHAREF           ! LLY Lloyd
      double complex, dimension(IEMXD) :: CDOSREF_LLY          ! LLY Lloyd
      double complex, dimension(LMMAXD,LMMAXD)  :: W1
      double complex, dimension(LMGF0D,LMGF0D)  :: WN1
      double complex, dimension(LMGF0D,LMGF0D)  :: WN2         ! LLY
      double complex, dimension(LMMAXD,LMMAXD)  :: TMAT
      double complex, dimension(LMMAXD,LMMAXD)  :: GMAT0
      double complex, dimension(LMMAXD,LMMAXD)  :: FACTL
      double complex, dimension(IEMXD,NSPIND)   :: TRACET      !< Tr[ (t-tref)^-1 d(t-tref)/dE ]  ! LLY Lloyd
      double complex, dimension(IEMXD,NSPIND)   :: TRALPHA
      double complex, dimension(IEMXD,NSPIND)   :: LLY_GRTR    !< LLY Lloyd  Trace[ M^-1 dM/dE ], Eq.5.38 PhD Thiess
      double complex, dimension(IEMXD,NSPIND)   :: CDOS_LLY
      double complex, dimension(0:LMAX,NREF)    :: ALPHAREF
      double complex, dimension(0:LMAX,NREF)    :: DALPHAREF   !< LLY Lloyd Alpha matrix and deriv.
      double complex, dimension(LMMAXD,LMMAXD,NATYP)  :: MSST
      double complex, dimension(LMMAXD,LMMAXD,NATYP)  :: TSST
      double complex, dimension(LMMAXD,LMMAXD,NAEZ)   :: TQDOS  ! qdos ruess
      double complex, dimension(LMMAXD,LMMAXD,NREF)   :: TREFLL
      double complex, dimension(LMMAXD,LMMAXD,NSHELD) :: GMATLL   !< GMATLL = diagonal elements of the G matrix (system)
      double complex, dimension(LMMAXD,LMMAXD,NREF)   :: DTREFLL  !< LLY Lloyd dtref/dE
      double complex, dimension(LMMAXD,LMMAXD,NAEZ)   :: DTMATLL  !< LLY Lloyd  dt/dE
      complex*8, dimension(LMMAXD*LMMAXD) :: GIMP !<  Cluster GF (ref. syst.)
      character(len=35), dimension(0:2) :: INVALG

      ! .. Allocatable local arrays
      integer, dimension(:,:), allocatable :: ISH
      integer, dimension(:,:), allocatable :: JSH
      double precision, dimension(:,:), allocatable   :: QVEC     !< qdos ruess, q-vectors for qdos
      double precision, dimension(:,:,:), allocatable :: BZKP
      double complex, dimension(:,:), allocatable     :: DTMTRX   !< For GREENIMP
      double complex, dimension(:,:,:), allocatable   :: GINP     !< Cluster GF (ref syst.) GINP(NACLSD*LMGF0D,LMGF0D,NCLSD)
      double complex, dimension(:,:,:), allocatable   :: DGINP    !< LLY Lloyd Energy derivative of GINP DGINP(NACLSD*LMGF0D,LMGF0D,NCLSD)
      double complex, dimension(:,:,:,:,:), allocatable :: LEFTTINVLL
      double complex, dimension(:,:,:,:,:), allocatable :: RIGHTTINVLL

#ifdef CPP_MPI
      integer :: ie_start
#endif
      integer :: ie_num, ie_end, ierr, i_stat, i_all
      !-------------------------------------------------------------------------
      !     for conductivity calculation
      !     INTEGER NCPAIRD
      !     PARAMETER(NCPAIRD=10)
      !     INTEGER IATCONDL(NCPAIRD),IATCONDR(NCPAIRD),NCONDPAIR
      !-------------------------------------------------------------------------
      !     .. Intrinsic Functions ..
      intrinsic :: ATAN

      !.. Set the parameters
      LRECTRA=WLENGTH*4 ! LLY Lloyd
      LRECGRF=WLENGTH*4*NACLSD*LMGF0D*LMGF0D*NCLSD ! 4 words = 16 bytes / complex number (in ifort 4; in gfort 16) word/byte distiction moved to subroutine opendafile to be the same for all unformatted files
      LRECTMT=WLENGTH*4*LMMAXD*LMMAXD
      LRECGREEN=WLENGTH*2*NATOMIMP*LMMAXD*NATOMIMP*LMMAXD
      !     ..
      !     .. Data statements
      DATA INVALG /'FULL MATRIX                        ',   &
                  'BANDED MATRIX (slab)               ',    &
                  'BANDED + CORNERS MATRIX (supercell)' /
      DATA IPRINT / 0 /
      !     ..
      if(t_inc%i_write>0) IPRINT=1

      ! allocatable arrays
      allocate(BZKP(3,KPOIBZ,MAXMSHD), stat=i_stat)
      call memocc(i_stat,product(shape(BZKP))*kind(BZKP),'BZKP','main1b')
      allocate(LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), stat=i_stat)
      call memocc(i_stat,product(shape(LEFTTINVLL))*kind(LEFTTINVLL),'LEFTTINVLL','main1b')
      allocate(RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD), stat=i_stat)
      call memocc(i_stat,product(shape(RIGHTTINVLL))*kind(RIGHTTINVLL),'RIGHTTINVLL','main1b')
      allocate(ISH(NSHELD,NOFGIJ), stat=i_stat)
      call memocc(i_stat,product(shape(ISH))*kind(ISH),'ISH','main1b')
      allocate(JSH(NSHELD,NOFGIJ), stat=i_stat)
      call memocc(i_stat,product(shape(JSH))*kind(JSH),'JSH','main1b')
      !Consistency check
      if ( (KREL.LT.0) .OR. (KREL.GT.1) ) stop ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
      if ( (KREL.EQ.1) .AND. (NSPIND.EQ.2) ) stop ' set NSPIND = 1 for KREL = 1 in the inputcard'

      !-------------------------------------------------------------------------
      ! This routine previously used to read from unformatted files created by
      ! the main0 module, now  instead of unformatted files take parameters from
      ! types defined in wunfiles.F90
      !-------------------------------------------------------------------------
      call get_params_1b(t_params,NATYP,NACLSD,IELAST,NPOL,NCLSD,NREF,NEMB,   &
         NAEZ,NSRA,INS,NSPIN,LMAX,NCLS,LLY,KREL,ATOM,CLS,NACLS,REFPOT,EZ,     &
         ITMPDIR,ILTMP,ALAT,RCLS,IEMXD,RMTREF,VREF,TMPDIR,NSHELD,NPRINCD,     &
         KPOIBZ,ATOMIMP,NATOMIMPD,ICC,IGF,NLBASIS,NRBASIS,NCPA,ICPA,          &
         ITCPAMAX,CPATOL,NR,IDECI,RBASIS,RR,EZOA,NSHELL,KMROT,KAOEZ,ISH,      &
         JSH,NSH1,NSH2,NOQ,IQAT,NOFGIJ,NATOMIMP,CONC,KMESH,MAXMESH,NSYMAT,    &
         NQCALC,RATOM,RROT,DROTQ,IJTABCALC,IJTABCALC_I,IJTABSYM,IJTABSH,      &
         IQCALC,DSYMLL,INVMOD,ICHECK,SYMUNITARY,RC,CREL,RREL,SRREL,NRREL,     &
         IRREL,LEFTTINVLL,RIGHTTINVLL,VACFLAG,NOFKS,VOLBZ,BZKP,VOLCUB,WEZ,    &
         NEMBD1,LMMAXD,NSYMAXD,NSPINDD,MAXMSHD,RCLSIMP)


      if (TEST('gmatasci')) open(298347,FILE='gmat.ascii',FORM='formatted')
      !-------------------------------------------------------------------------
      ! End of reading the variables
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------    !fswrt
      ! open file to store the output for the (external) Fermi-surface program      !fswrt
      ! this file is already partly filled with data by main0. More data            !fswrt
      !       will be stored in                                                     !fswrt
      !-------------------------------------------------------------------------    !fswrt
      if (OPT('FERMIOUT') .and. myrank==master) then                                !fswrt
         open(6801,FILE='TBkkr_container.txt', FORM='formatted',POSITION='append')  !fswrt
      endif                                                                         !fswrt

      !-------------------------------------------------------------------------    !fswrt
      ! open file for WRTGREEN option (writes green_host file for
      ! GMATLL_GES creation in zulapi part) file is filled in ROTGLL called in kloopz
      !-------------------------------------------------------------------------    !fswrt
      if (OPT('WRTGREEN') .and. myrank==master) then
         open(58,FILE='green_host', FORM='formatted')
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! If qdos option is used set IQDOSRUN so that in a first run the
      ! (t(E)-t_ref(E))^-1 matrix (tmat.qdos) and the gref matrix can be
      ! written out for one k point, in a second run these matrices are
      ! read in to continue the calculation with the k points specified by
      ! the user in the qvec.dat file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( OPT('qdos    ') ) then                       ! qdos ruess
         IQDOSRUN=0                                     ! qdos ruess
      else                                              ! qdos ruess
         IQDOSRUN=-1                                    ! qdos ruess
      endif                                             ! qdos ruess
      ! Jump back here to continue with second run if qdos option is selected
210   continue                                          ! qdos ruess
      ! Reset GMATLL for calculation in second run
      if ( IQDOSRUN.eq.1 ) then                         ! qdos ruess
         do I1 = 1,NSHELL(0)                            ! qdos ruess
            GMATLL(1:LMMAXD,1:LMMAXD,I1) = CZERO        ! qdos ruess
         enddo                                          ! qdos ruess
      endif                                             ! qdos ruess

      if ((opt('qdos    ')) .and. (opt('deci-out'))) then
         stop 'ERROR: qdos and deci-out cannot be used simultaniously'
      elseif (opt('qdos    ')) then
#ifdef CPP_MPI
         open (37,ACCESS='direct', RECL=WLENGTH*4,FILE='tmat.qdos',FORM='unformatted')
#else
         open(37, File='tmat.qdos', form='formatted')      ! qdos ruess
#endif
      elseif (opt('deci-out')) then
         open(37, file='decifile', form='formatted', position='append') ! ruess: needed in case of deci-out option to prepare decifile
      end if

      do I=1,NAEZ
         do L=1,NOQ(I)
            ITOQ(L,I) = KAOEZ(L,I)
         end do
      end do
      RFCTOR = ALAT/(2*PI)           ! = ALAT/(2*PI)
      CFCTOR = CONE*RFCTOR
      CFCTORINV = CONE/RFCTOR

      call SETFACTL(FACTL,LMAX,KREL,LMMAXD)
      !
      if(t_inc%i_write>0) write(1337,2100) INVALG(INVMOD)

      NACLSMAX = 1
      do IC = 1,NCLS
         if (NACLS(IC).GT.NACLSMAX) NACLSMAX = NACLS(IC)
      enddo
      LRECGRF1 = WLENGTH*4*NACLSMAX*LMGF0D*LMGF0D*NCLS

      if (.not.allocated(GINP)) then
         allocate(GINP(NACLSMAX*LMGF0D,LMGF0D,NCLS),stat=i_stat)
         call memocc(i_stat,product(shape(GINP))*kind(GINP),'GINP','main1b')
      endif
      if (.not.allocated(DGINP)) then
         allocate(DGINP(NACLSMAX*LMGF0D,LMGF0D,NCLS),stat=i_stat)
         call memocc(i_stat,product(shape(DGINP))*kind(DGINP),'DGINP','main1b')
      endif

      if (t_tgmat%gref_to_file) then
         call OPENDAFILE(68,'gref',4,LRECGRF1,TMPDIR,ITMPDIR,ILTMP)
      end if
      if (t_tgmat%tmat_to_file) then
         call OPENDAFILE(69,'tmat',4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)
      end if
      if (t_tgmat%gmat_to_file) then
         call OPENDAFILE(70,'gmat',4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)
      end if
      if (LLY.NE.0) then                                               ! LLY Lloyd
         if (t_lloyd%dgref_to_file) then
            call OPENDAFILE(681,'dgrefde',7,LRECGRF1,TMPDIR,ITMPDIR,ILTMP) ! LLY Lloyd: derivative of Gref
         end if
         if (t_lloyd%g0tr_to_file) then
            open(682,FILE='lly_g0tr_ie.ascii',FORM='FORMATTED')           ! LLY Lloyd: trace eq.5.27 PhD Thiess
         end if
         if (t_lloyd%dtmat_to_file) then
            call OPENDAFILE(691,'dtmatde',7,LRECTMT,TMPDIR,ITMPDIR,ILTMP) ! LLY Lloyd: derivative of t-matrix
         end if
         if (t_lloyd%tralpha_to_file) then
            call OPENDAFILE(692,'tralpha',7,LRECTRA,TMPDIR,ITMPDIR,ILTMP) ! LLY Lloyd: Tr[alpha^{-1} dalpha/dE]
         end if
      endif                                                            ! LLY Lloyd

      LCPAIJ = .FALSE.
      if ( ( NCPA.NE.0 ).AND.( NSHELL(0).GT.NATYP ) ) LCPAIJ = .TRUE.
      !
      if ( LCPAIJ ) then
         if(t_cpa%dmatproj_to_file) then
            open (71,ACCESS='direct',RECL=2*LRECTMT,FILE='dmatproj.unformatted',FORM='unformatted')
         else
#ifndef CPP_MPI
            call init_t_cpa(t_inc,t_cpa,IELAST)
#else
            call init_t_cpa(t_inc,t_cpa,t_mpi_c_grid%ntot2)
#endif
         end if !t_cpa%dmatproj_to_file
      end if !LCPAIJ
      !
      if ( IGF.NE.0 ) then

         if ( ( OPT('GPLAIN  ') ) ) then
            open (8888,FILE='kkrflex_green.dat')
         end if

         if ( ( OPT('KKRFLEX ') ) ) then
            RECLENGTH = WLENGTH*2*NATOMIMP*LMMAXD*NATOMIMP*LMMAXD
            ! sometimes (lmax=2) the record length might be too small to store the parameters, then reclength needs to be bigger
            if(RECLENGTH<8*IELAST+6) then
               stop '[main1b] record length for kkrflex_green is too small to store parameters, use either more atoms in the cluster, a higher lmax or less energy points'
            end if
            open (888,ACCESS='direct',RECL=RECLENGTH,FILE='kkrflex_green',FORM='unformatted')
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! the following write-out has been disabled, because it was assumed to be    !no-green
         !  obsolete with the implementation of the MPI-communicated arrays. If I am  !no-green
         !  wrong and the write-out is needed in subsequent parts, construct a        !no-green
         !  test-option around it so that it is only written out in this case.        !no-green
         !     OPEN (88,ACCESS='direct',RECL=LRECGREEN,                               !no-green
         !&             FILE='green',FORM='unformatted')                              !no-green
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IREC=1

         if ( ( OPT('KKRFLEX ') ) ) then
            if(myrank==master) then
               write(888,REC=IREC) IELAST,NSPIN,NATOMIMP,NATOMIMP,LMMAXD,  &
               KORBIT,(EZ(IE),IE=1,IELAST),(WEZ(IE),IE=1,IELAST)
               if ( ( OPT('GPLAIN  ') ) ) then
                  !WRITE(8888,'(I5,50000F)') IELAST,NSPIN,NATOMIMP,NATOMIMP,&
                  WRITE(8888,*) IELAST,NSPIN,NATOMIMP,NATOMIMP,&
                  (LMAX+1)**2,(EZ(IE),IE=1,IELAST),(WEZ(IE),IE=1,IELAST)
               end if
            end if
#ifdef CPP_MPI
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
         end if
         !         IF ( (.not. OPT('KKRFLEX ') ) ) THEN                !no-green
         !           WRITE(88,REC=IREC) IELAST,NSPIN,                  !no-green
         !    &         (EZ(IE),IE=1,IELAST),(WEZ(IE),IE=1,IELAST),    !no-green
         !    &         NATOMIMPD*LMMAXD                               !no-green
         !         END IF                                              !no-green
      end if

      ! Value of NQDOS changes to a read-in value if option qdos is applied, otherwise:
      NQDOS = 1                                         ! qdos ruess
      if (OPT('qdos    ').AND.(IQDOSRUN.EQ.1)) then     ! qdos ruess
         ! Read BZ path for qdos calculation:
         open(67,FILE='qvec.dat')                       ! qdos ruess
         read(67,*) NQDOS                               ! qdos ruess
         i_all=-product(shape(QVEC))*kind(QVEC)
         deallocate(QVEC, stat=i_stat)                  ! qdos ruess: deallocate in first run allocated array to change it
         call memocc(i_stat,i_all,'QVEC','main1b')

         allocate(QVEC(3,NQDOS), stat=i_stat)           ! qdos ruess
         call memocc(i_stat,product(shape(QVEC))*kind(QVEC),'QVEC','main1b')
         do IQ = 1,NQDOS                                ! qdos ruess
            read(67,*) (QVEC(IX,IQ),IX=1,3)             ! qdos ruess
         enddo                                          ! qdos ruess
         close(67)                                      ! qdos ruess
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Prepare k-mesh information to be appropriate for qdos calculation.
         ! The idea is that subr. KLOOPZ1 is called for only one point at a time,
         ! with weight equal to the full BZ; in this way we avoid changing the
         ! calling list or the contents of kloopz1.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         KMESH(1:IELAST) = 1                             ! qdos ruess
         NOFKS(1) = 1                                    ! qdos ruess
         VOLCUB(1,1) = VOLBZ(1)                          ! qdos ruess
         NSYMAT = 1
      else if (OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) then ! qdos ruess
         ! Call the k loop just once with one k point to write out the tmat.qdos file
         allocate(QVEC(3,NQDOS), stat=i_stat)            ! qdos ruess
         call memocc(i_stat,product(shape(QVEC))*kind(QVEC),'QVEC','main1b')
         if(ierr/=0) stop '[main1b] Error allocating qvec'
         QVEC(1:3,1) = 0.D0                              ! qdos ruess
         KMESH(1:IELAST) = 1                             ! qdos ruess
         NOFKS(1) = 1                                    ! qdos ruess
         VOLCUB(1,1) = VOLBZ(1)                          ! qdos ruess
      end if                                             ! qdos ruess

      NCPAFAIL = 0

      ! Initialize trace for Lloyd formula                    ! LLY Lloyd
      LLY_GRTR(:,:) = CZERO ! 1:IEMXD,1:NSPIND          ! LLY Lloyd

      if (.not.OPT('NEWSOSOL')) then
         !----------------------------------------------------------------------
         ! BEGIN do loop over spins and energies
         !----------------------------------------------------------------------
         do 370 ISPIN = 1,NSPIN
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! NO-SOC
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CPP_MPI
            ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
            ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)

            do 360 ie_num=1,ie_end
               IE = ie_start+ie_num

               !start timing measurement for this ie, needed for MPIadapt
               if(MPIadapt.and.t_mpi_c_grid%myrank_ie==0) then
                  call timing_start('time_1b_ie')
               end if
#else
            do 360 IE = 1,IELAST
               ie_num = ie
               ie_end = ielast
#endif

               ! write energy into green_host file
               if (OPT('WRTGREEN') .and. myrank==master) then
                  write(58,'(2e17.9)') EZ(IE)
               endif

               if (t_tgmat%gref_to_file) then
                  read (68,REC=IE) GINP
               else
                  ginp(:,:,:) = t_tgmat%gref(:,:,:,ie_num)
               end if
               if(t_lloyd%dgref_to_file) then
                  if (LLY.NE.0) read (681,REC=IE) DGINP   ! LLY Lloyd
               else
                  if (LLY.NE.0) DGINP(:,:,:) = t_lloyd%dgref(:,:,:,ie_num)
               endif

               ERYD = EZ(IE)
               NMESH = KMESH(IE)
               if(t_inc%i_write>0) write (1337,'(A,I3,A,2(1X,F10.6),A,I3)')&
                  ' ************ IE = ',IE,' ENERGY =',EZ(IE),' KMESH = ', NMESH
               !----------------------------------------------------------------
               !  I1 = 1,NREF
               ! calculate t(ll') of the reference system (on clusters)
               !----------------------------------------------------------------
#ifdef CPP_TIMING
               call timing_start('main1b - calctref13')
#endif
               if ( KREL.EQ.0 ) then
                  do I1 = 1,NREF
                     call CALCTREF13(ERYD,VREF(I1),RMTREF(I1),LMAX,LM1, &  ! LLY Lloyd
                        TREFLL(1,1,I1),DTREFLL(1,1,I1),ALPHAREF(0,I1),  &  ! LLY Lloyd
                        DALPHAREF(0,I1),LMAX+1,LMMAXD)                    ! LLY Lloyd
                  end do
               else
                  do I1 = 1,NREF
                     call CALCTREF13(ERYD,VREF(I1),RMTREF(I1),LMAX,LM1, &        ! LLY Lloyd
                        WN1,WN2,ALPHAREF(0,I1),DALPHAREF(0,I1),LMAX+1,LMGF0D)   ! LLY Lloyd
                        !-------------------------------------------------------
                        ! add second spin-block for relativistic calculation and transform
                        ! from NREL to REL representation
                        !-------------------------------------------------------
                     call CINIT(LMMAXD*LMMAXD,W1)
                     if (LMMAXD.NE.LM1*2) stop 'LMMAXD <> LM1*2 '
                     do I=1,LM1
                        W1(I,I) = WN1(I,I)
                        W1(LM1+I,LM1+I) = WN1(I,I)
                     end do
                     call CHANGEREP(W1,'RLM>REL',TREFLL(1,1,I1),LMMAXD,&
                        LMMAXD,RC,CREL,RREL,'TREFLL',0)
                  end do
               end if
#ifdef CPP_TIMING
               call timing_pause('main1b - calctref13')
#endif
               !----------------------------------------------------------------
               TRALPHA(IE,ISPIN) = CZERO                                    ! LLY
               TRALPHAREF(IE) = CZERO                                       ! LLY
               do I1 = 1,NATYP

                  if (t_tgmat%tmat_to_file) then
                     IREC = IE + IELAST* (ISPIN-1) + IELAST*NSPIN* (I1-1)
                     read (69,REC=IREC) TMAT
                  else
                     irec = ie_num+ie_end*(ISPIN-1)+ie_end*NSPIN*(I1-1)
                     tmat(:,:) = t_tgmat%tmat(:,:,irec)
                  end if
                  TSST(1:LMMAXD,1:LMMAXD,I1)=TMAT(1:LMMAXD,1:LMMAXD)

                  if (LLY.ne.0) then                                             ! LLY
                     if(t_lloyd%dtmat_to_file) then
                        IREC =IE+IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1)
                        read (691,REC=IREC) TMAT                                 ! LLY dt/dE
                     else
                        irec = ie_num+ie_end*(ISPIN-1)+ie_end*NSPIN*(I1-1)
                        TMAT(:,:) = t_lloyd%dtmat(:,:,irec)
                     end if
                     DTMATLL(1:LMMAXD,1:LMMAXD,I1) =TMAT(1:LMMAXD,1:LMMAXD)      ! LLY
                     if(t_lloyd%dtmat_to_file) then
                        IREC=IE+IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1)
                        read (692,REC=IREC) TRALPHA1                             ! LLY
                     else
                        irec = ie_num+ie_end*(ISPIN-1)+ie_end*NSPIN*(I1-1)
                        TRALPHA1 = t_lloyd%tralpha(irec)
                     end if
                     TRALPHA(IE,ISPIN) = TRALPHA(IE,ISPIN) + TRALPHA1            ! LLY Tr[ alpha^{-1} dalpha/dE]

                     if (ISPIN.eq.1) then  ! Ref. system is spin-independent     ! LLY
                        TRALPHA1 = CZERO                                         ! LLY
                        do L1 = 0,LMAX                                           ! LLY
                           TRALPHA1 = TRALPHA1 + (2*L1+1) * &                    ! LLY
                              DALPHAREF(L1,REFPOT(I1))/ALPHAREF(L1,REFPOT(I1))   ! LLY
                        enddo
                        TRALPHAREF(IE) = TRALPHAREF(IE) + TRALPHA1               ! LLY Tr[ alpharef^{-1} dalpharef/dE
                     endif
                  endif                                                          ! LLY

               end do  !i1 = 1,natyp


               if(t_lloyd%g0tr_to_file) then
                  if (LLY.ne.0.and.ISPIN.eq.1) read(682,FMT='(2E24.16)') LLY_G0TR(IE)  ! LLY
               else
                  if (LLY.ne.0.and.ISPIN.eq.1) LLY_G0TR(IE) = t_lloyd%g0tr(ie_num)        ! LLY
               end if
               ! ------------------------------------------------------------------
               ! Setting up of Delta_t moved to < KLOOPZ1 >
               ! ------------------------------------------------------------------
               if (OPT('readcpa ').or.(OPT('qdos    ').and.(IQDOSRUN.eq.1))) then     ! qdos ruess: read in cpa t-matrix
                  do ISITE = 1,NAEZ                              ! qdos ruess
                     TQDOS(:,:,ISITE) = CZERO                    ! qdos ruess
#ifdef CPP_MPI
                     do lm1=1,lmmaxd
                        do lm2=1,lmmaxd
                           irec = LM2+(LM1-1)*LMMAXD+LMMAXD**2*(isite-1)+&
                              LMMAXD**2*naez*(ie-1)+LMMAXD**2*ielast*naez*(ispin-1)
                           read(37,rec=irec) tread
                           if ( (LM1+LM2).ne.0 ) then                  ! qdos ruess
                              TQDOS(LM1,LM2,ISITE) = TREAD / CFCTORINV ! qdos ruess
                           end if                                      ! qdos ruess
                        end do
                     end do

#else
                     read(37,*) TEXT                             ! qdos ruess
                     read(37,*) TEXT                             ! qdos ruess
 9921                continue                                    ! qdos ruess
                     read(37,99013) LM1,LM2,TREAD                ! qdos ruess
                     99013 format ('(2I5,1P,2D22.14)')                 ! qdos ruess
                     if ( (LM1+LM2).ne.0 ) then                  ! qdos ruess
                        TQDOS(LM1,LM2,ISITE) = TREAD / CFCTORINV ! qdos ruess
                        if ( (LM1+LM2).lt.2*LMMAXD ) goto 9921   ! qdos ruess
                     end if                                      ! qdos ruess
#endif
                  enddo                                          ! qdos ruess
               end if                                            ! qdos ruess
               !  Loop over all QDOS points and change volume for KLOOPZ run accordingly
               do 200 IQ = 1,NQDOS                               ! qdos ruess
                  if (OPT('qdos    ')) BZKP(:,1,1) = QVEC(:,IQ)     ! qdos ruess: Set q-point x,y,z
!
#ifdef CPP_TIMING
                  call timing_start('main1b - kloopz')
#endif
                  call KLOOPZ1_QDOS(NR,NEMBD1,LMMAXD,LMGF0D,LMAX,NREF,ERYD,GMATLL,INS,ALAT,IE,IGF,NSHELL,NAEZ,  &
                     NOFKS(NMESH),VOLBZ(NMESH),BZKP(1,1,NMESH),               &
                     VOLCUB(1,NMESH),CLS,NACLS,NACLSMAX,NCLS,RR,              &
                     RBASIS,EZOA,ATOM,RCLS,ICC,GINP,IDECI,                    &
                     LEFTTINVLL(1,1,1,1,IE),RIGHTTINVLL(1,1,1,1,IE),          &
                     VACFLAG,NLBASIS,NRBASIS,FACTL,NATOMIMP,NSYMAT,           &
                     DSYMLL,RATOM,RROT,NSH1,NSH2,IJTABSYM,IJTABSH,            &
                     ICHECK,INVMOD,REFPOT,TREFLL,TSST,MSST,CFCTOR,            &
                     CFCTORINV,CREL,RC,RREL,SRREL,IRREL,NRREL,DROTQ,          &
                     SYMUNITARY,KMROT,NATYP,NCPA,ICPA,ITCPAMAX,               &
                     CPATOL,NOQ,IQAT,ITOQ,CONC,IPRINT,ICPAFLAG,               &
                     ISPIN,NSPINDD,                                           &
                     TQDOS,IQDOSRUN,                                          &  ! qdos
                     DTREFLL,DTMATLL,DGINP,LLY_GRTR(IE,ISPIN),                &  ! LLY Lloyd
                     TRACET(IE,ISPIN),LLY)                                       ! LLY Lloyd

#ifdef CPP_TIMING
                  call timing_pause('main1b - kloopz')
#endif
                  !-------------------------------------------------------------
                  ! Skip this part if first part of the qdos is running
                  !-------------------------------------------------------------
                  if ( .not.(OPT('qdos    ').and.(IQDOSRUN.eq.0)) ) then
                     if (NCPA.ne.0) then
                        if (ICPAFLAG .ne. 0) then
                           NCPAFAIL = NCPAFAIL + 1
                           IECPAFAIL(NCPAFAIL)= IE
                        end if
                     end if  ! (NCPA.NE.0)

                     do I1 = 1,NSHELL(0)
                        GMAT0(:,:)=GMATLL(:,:,I1)
                        IREC = IQ + NQDOS * (IE-1) + NQDOS * IELAST *        &   ! qdos ruess: (without qdos, IQ=NQ=1)
                           (ISPIN-1) + NQDOS * IELAST * NSPIN * (I1-1)           ! qdos ruess
                        if (t_tgmat%gmat_to_file) then
                           write (70,REC=IREC) GMAT0
                        else
                           IREC = IQ + NQDOS * (ie_num-1) + NQDOS *  &
                              ie_end * (ISPIN-1) + NQDOS * ie_end * NSPIN * (I1-1)
                           t_tgmat%gmat(:,:,irec) = gmat0
                        end if
                     enddo
                     if (TEST('gmatasci')) then
                        write(*,*) 'Writing out gmat.ascii'
                        do I1 = 1,NSHELL(0)
                           do LM1=1,LMMAXD
                              do LM2=1,LMMAXD
                                 write(298347,FMT='(3I5,2E25.16)')&
                                 I1,LM1,LM2,GMATLL(LM1,LM2,I1)
                              enddo
                           enddo
                        enddo
                     endif

                     ! writeout of host green function for impurity code for single-atom cluster (not captured in rotgll)
                     if ( NATOMIMP==1 ) then
                        I1=ATOMIMP(1)
                        if ( OPT('KKRFLEX ') ) then
                           irec = ielast*(ispin-1)+ ie+1
                           ILM=0
                           GIMP=(0.e0,0.e0) !complex*8
                           do LM2=1,LMMAXD
                              do LM1=1,LMMAXD
                                 ILM=ILM+1
                                 GIMP(ILM)=GMATLL(LM1,LM2,I1)
                              enddo
                           enddo
                           irec = ielast*(ispin-1)+ ie+1
                           write(888,REC=irec) GIMP
                           if ( OPT('GPLAIN  ') ) then
                              !write(8888,'(50000E)') GIMP
                              write(8888,*) GIMP
                           end if
                        endif ! KKRFLEX
                        if (OPT('WRTGREEN') .and. myrank==master) then
                           do LM2=1,LMMAXD
                              do LM1=1,LMMAXD
                                 ! writeout of green_host for WRTGREEN option
                                 write(58,'((2I5),(2e17.9))') LM2, LM1,GMATLL(LM1,LM2,I1)
                              end do
                           end do
                        endif ! WRTGREEN
                     end if !( NATOMIMP==1 )

                     if ( LCPAIJ ) then
                        if(t_cpa%dmatproj_to_file) then
                           do I1 = 1,NATYP
                              GMAT0(:,:) = TSST(:,:,I1)
                              W1(:,:)    = MSST(:,:,I1)
                              IREC = IE + IELAST*(ISPIN-1) + IELAST*NSPIN*(I1-1)
                              write (71,REC=IREC) GMAT0,W1
                           end do
                        else !t_cpa%dmatproj_to_file
                           irec = ie_num + ie_end*(ISPIN-1)
                           t_cpa%dmatts(:,:,:,irec) = TSST(:,:,:)
                           t_cpa%dtilts(:,:,:,irec) = MSST(:,:,:)
                        end if!t_cpa%dmatproj_to_file
                     end if  ! ( LCPAIJ )

                  endif    ! ( .NOT.(OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) )
 200           continue ! IQ = 1,NQDOS                                       ! qdos ruess


               if (LLY.ne.0) then                                                      ! LLY Lloyd

                  if (LLY.ne.2) then                                                   ! LLY Lloyd
                     CDOS_LLY(IE,ISPIN) =   TRALPHA(IE,ISPIN)        &                 ! LLY Lloyd
                     - LLY_GRTR(IE,ISPIN) / VOLBZ(1) + LLY_G0TR(IE)                    ! LLY Lloyd
                  else                                                                 ! LLY Lloyd
                     CDOS_LLY(IE,ISPIN) =    TRACET(IE,ISPIN)+ TRALPHAREF(IE) &        ! LLY Lloyd
                     - LLY_GRTR(IE,ISPIN) / VOLBZ(1) +  LLY_G0TR(IE)                   ! LLY Lloyd
                  endif                                                                ! LLY Lloyd

                  if (ISPIN.EQ.1) CDOSREF_LLY(IE) = TRALPHAREF(IE) - LLY_G0TR(IE)      ! LLY Lloyd

                  if (TEST('GMAT=0  ')) then                                           ! LLY Lloyd
                     CDOS_LLY(IE,ISPIN) = TRALPHA(IE,ISPIN)                            ! LLY Lloyd
                     if (LLY.EQ.2) then                                                ! LLY Lloyd
                        CDOS_LLY(IE,ISPIN)=TRACET(IE,ISPIN) + TRALPHAREF(IE)           ! LLY Lloyd
                     endif                                                             ! LLY Lloyd
                  endif                                                                ! LLY Lloyd

                  CDOS_LLY(IE,ISPIN) = CDOS_LLY(IE,ISPIN) / PI                         ! LLY Lloyd

               endif                                                                   ! LLY Lloyd

               ! ---------------------------------------------------------------

#ifdef CPP_MPI
               !stop timing measurement for this ie, needed for MPIadapt
               if(MPIadapt.and.t_mpi_c_grid%myrank_ie==0) then
                  call timing_stop('time_1b_ie', save_out=timings_1b(ie) )
               end if
#endif

 360        continue               ! IE = 1,IELAST
#ifdef CPP_TIMING
            if(.not.OPT('GREENIMP')) then
               if(t_inc%i_time>0) call timing_stop('main1b - calctref13')
               if(t_inc%i_time>0) call timing_pause('main1b_fourier')
               if(t_inc%i_time>0) call timing_stop('main1b_inversion')
               if(t_inc%i_time>0) call timing_stop('main1b - kloopz')
            endif
#endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! NO-SOC
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if( NCPAFAIL .NE. 0 ) then
               if(t_inc%i_write>0) then
                  write(1337,*)
                  write(1337,'(1X,79(''*''),/)')
                  write(1337,99019) CPATOL, NCPAFAIL,    &
                  (IECPAFAIL(IE),DBLE(EZ(IECPAFAIL(IE))),IE=1,NCPAFAIL)
                  write(1337,'(1X,79(''*''),/)')
                  write(1337,*)
               endif
            else
               if( NCPA .ne. 0 ) then
                  if(t_inc%i_write>0) then
                     write(1337,*)
                     write(1337,99020)
                     write(1337,*)
                  endif
               end if
            end if
!
 370     continue                  !  ISPIN = 1,NSPIN


      else ! NEW SOC SOLVER

         ! nonco angles
         call read_angles(t_params,NATYP,THETAS,PHIS)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! SOC
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CPP_MPI
         ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
         ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)

         do 460 ie_num=1,ie_end
            IE = ie_start+ie_num

            !start timing measurement for this ie, needed for MPIadapt
            if(MPIadapt.and.t_mpi_c_grid%myrank_ie==0) then
               call timing_start('time_1b_ie')
            end if

#else
         ie_end = IELAST
         do 460 IE=1,IELAST
            ie_num = ie
#endif

            ! write energy into green_host file
            if (OPT('WRTGREEN') .and. myrank==master) then
               write(58,'(2e17.9)') EZ(IE)
            endif

            ! read in Green function of reference system
            if (t_tgmat%gref_to_file) then
               READ (68,REC=IE) GINP
            else
               ginp(:,:,:) = t_tgmat%gref(:,:,:,ie_num)
            end if
            if(t_lloyd%dgref_to_file) then
               if (LLY.NE.0) read (681,REC=IE) DGINP   ! LLY Lloyd
            else
               if (LLY.NE.0) DGINP(:,:,:) = t_lloyd%dgref(:,:,:,ie_num)
            endif
            ERYD = EZ(IE)
            NMESH = KMESH(IE)
            if(t_inc%i_write>0) then
               write (1337,'(A,I3,A,2(1X,F10.6),A,I3)') &
                  ' ************ IE = ',IE,' ENERGY =',EZ(IE),' KMESH = ', NMESH
            endif

            ! construct t matrix for reference system (now is always double matrix)
#ifdef CPP_TIMING
            call timing_start('main1b - calctref13')
#endif
            do I1 = 1,NREF
               call CALCTREF13(ERYD,VREF(I1),RMTREF(I1),LMAX,LM1,WN1,WN2,  &  ! LLY
                  ALPHAREF(0,I1),DALPHAREF(0,I1),LMAX+1,LMGF0D)              ! LLY
               do I=1,LM1
                  TREFLL(I,I,I1) = WN1(I,I)
                  TREFLL(LM1+I,LM1+I,I1) = WN1(I,I)
                  DTREFLL(I,I,I1) = WN2(I,I)                 ! LLY
                  DTREFLL(LM1+I,LM1+I,I1) = WN2(I,I)         ! LLY
               enddo
            enddo ! I1
#ifdef CPP_TIMING
            call timing_pause('main1b - calctref13')
#endif
            TRALPHA(IE,NSPINDD)=CZERO
            TRALPHAREF(IE)=CZERO
            ! read in t matrix

            do I1 = 1,NATYP

               ! read in theta and phi for noncolinear
               THETA = THETAS(I1)
               PHI   = PHIS(I1)

               ! read in t-matrix from file
               IREC = IE + IELAST*(I1-1)
               if (t_tgmat%tmat_to_file) then
                  read (69,REC=IREC) TMAT
               else
                  irec = ie_num + ie_end * (i1-1)
                  tmat(:,:) = t_tgmat%tmat(:,:,irec)
               end if

               ! rotate t-matrix from local to global frame
               call ROTATEMATRIX(TMAT,THETA,PHI,LMGF0D,0)

               do LM1=1,LMMAXD
                  do LM2=1,LMMAXD
                     TSST(LM1,LM2,I1)=TMAT(LM1,LM2)
                  enddo
               enddo

               if (LLY.NE.0) then
                  IREC = IE + IELAST*(I1-1)
                  if(t_lloyd%dtmat_to_file) then
                     read(691,REC=IREC) TMAT ! LLY
                  else
                     irec = ie_num + ie_end * (i1-1)
                     TMAT(:,:) = t_lloyd%dtmat(:,:,irec)
                  end if
                  call ROTATEMATRIX(TMAT,THETA,PHI,LMGF0D,0) ! LLY
                  do LM1=1,LMMAXD
                     do LM2=1,LMMAXD
                        DTMATLL(LM1,LM2,I1)=TMAT(LM1,LM2) ! LLY
                     enddo
                  enddo
                  IREC = IE + IELAST*(I1-1)
                  if(t_lloyd%dtmat_to_file) then
                     read(692,REC=IREC) TRALPHA1
                  else
                     irec = ie_num + ie_end * (i1-1)
                     TRALPHA1 = t_lloyd%tralpha(irec)
                  end if

                  TRALPHA(IE,NSPINDD)=TRALPHA(IE,NSPINDD)+TRALPHA1 ! LLY
                  TRALPHA1=CZERO
                  do L1=0,LMAX
                     TRALPHA1=TRALPHA1+(2*L1+1)*   &
                     DALPHAREF(L1,REFPOT(I1))/ALPHAREF(L1,REFPOT(I1)) ! LLY
                  enddo
                  TRALPHAREF(IE)=TRALPHAREF(IE)+TRALPHA1 ! LLY
               endif ! LLY
            enddo ! I1
            if(t_lloyd%g0tr_to_file) then
               if (LLY.NE.0) read(682,FMT='(2E24.16)') LLY_G0TR(IE)  ! LLY
            else
               if (LLY.NE.0) LLY_G0TR(IE) = t_lloyd%g0tr(ie_num)
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS QDOS
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (OPT('readcpa ').or.(OPT('qdos    ').and.(IQDOSRUN.eq.1))) then     ! qdos ruess: read in cpa t-matrix
               do ISITE = 1,NAEZ                              ! qdos ruess
                  TQDOS(:,:,ISITE) = CZERO                    ! qdos ruess
#ifdef CPP_MPI
                  do lm1=1,lmmaxd
                     do lm2=1,lmmaxd
                        irec = LM2+(LM1-1)*LMMAXD+LMMAXD**2*(isite-1)+  &
                        LMMAXD**2*naez*(ie-1)
                        read(37,rec=irec) tread
                        if ( (LM1+LM2).NE.0 ) then                  ! qdos ruess
                           TQDOS(LM1,LM2,ISITE) = TREAD / CFCTORINV ! qdos ruess
                        end if                                      ! qdos ruess
                     end do
                  end do

#else
                  reaD(37,*) TEXT                             ! qdos ruess
                  reaD(37,*) TEXT                             ! qdos ruess
 9920             continue                                    ! qdos ruess
                  read(37,99014) LM1,LM2,TREAD                ! qdos ruess
99014             format (2I5,1P,2D22.14)                     ! qdos ruess
                  if ( (LM1+LM2).NE.0 ) then                  ! qdos ruess
                     TQDOS(LM1,LM2,ISITE) = TREAD / CFCTORINV ! qdos ruess
                     if ( (LM1+LM2).LT.2*LMMAXD ) goto 9920   ! qdos ruess
                  end if                                      ! qdos ruess
#endif
               enddo                                          ! qdos ruess
            end if                                            ! qdos ruess
            !-------------------------------------------------------------------
            !  Loop over all QDOS points and change volume for KLOOPZ run accordingly
            !-------------------------------------------------------------------
            do 220 IQ = 1,NQDOS                               ! qdos ruess
               if (OPT('qdos    ')) BZKP(:,1,1) = QVEC(:,IQ)     ! qdos ruess: Set q-point x,y,z

#ifdef CPP_TIMING
               call timing_start('main1b - kloopz')
#endif
               call KLOOPZ1_QDOS(NR,NEMBD1,LMMAXD,LMGF0D,LMAX,NREF,ERYD,GMATLL,INS,ALAT,IE,IGF,NSHELL,NAEZ,  &
                  NOFKS(NMESH),VOLBZ(NMESH),BZKP(1,1,NMESH),               &
                  VOLCUB(1,NMESH),CLS,NACLS,NACLSMAX,NCLS,RR,              &
                  RBASIS,EZOA,ATOM,RCLS,ICC,GINP,IDECI,                    &
                  LEFTTINVLL(1,1,1,1,IE),RIGHTTINVLL(1,1,1,1,IE),          &
                  VACFLAG,NLBASIS,NRBASIS,FACTL,NATOMIMP,NSYMAT,           &
                  DSYMLL,RATOM,RROT,NSH1,NSH2,IJTABSYM,IJTABSH,            &
                  ICHECK,INVMOD,REFPOT,TREFLL,TSST,MSST,CFCTOR,            &
                  CFCTORINV,CREL,RC,RREL,SRREL,IRREL,NRREL,DROTQ,          &
                  SYMUNITARY,KMROT,NATYP,NCPA,ICPA,ITCPAMAX,               &
                  CPATOL,NOQ,IQAT,ITOQ,CONC,IPRINT,ICPAFLAG,               &
                  1,NSPINDD,                                               &
                  TQDOS,IQDOSRUN,                                          &  ! qdos
                  DTREFLL,DTMATLL,DGINP,LLY_GRTR(IE,ISPIN),                &  ! LLY Lloyd
                  TRACET(IE,ISPIN),LLY)                                       ! LLY Lloyd
#ifdef CPP_TIMING
               call timing_pause('main1b - kloopz')
#endif

               ! Skip this part if first part of the qdos is running
               if ( .NOT.(OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) ) then
                  if (NCPA.NE.0) then
                     if (ICPAFLAG .NE. 0) then
                        NCPAFAIL = NCPAFAIL + 1
                        IECPAFAIL(NCPAFAIL)= IE
                     end if
                  end if  ! (NCPA.NE.0)
                  do I1 = 1,NSHELL(0)
                     GMAT0(1:LMMAXD,1:LMMAXD) =GMATLL(1:LMMAXD,1:LMMAXD,I1)
                     IREC = IQ + NQDOS * (IE-1) + NQDOS * IELAST * (I1-1) ! qdos ruess
                     if (t_tgmat%gmat_to_file) then
                        write (70,REC=IREC) GMAT0
                        ! human readable writeout if test option is hit
                        if(test('fileverb')) then
                           write(707070,'(i9,200000F15.7)') irec, gmat0
                        end if
                     else
                        IREC = IQ + NQDOS * (ie_num-1) + NQDOS * &
                        ie_end * (I1-1)
                        t_tgmat%gmat(:,:,irec) = gmat0
                     end if
                  enddo
                  if (TEST('gmatasci')) then
                     write(*,*) 'Writing out gmat.ascii'
                     do I1 = 1,NSHELL(0)
                        do LM1=1,LMMAXD
                           do LM2=1,LMMAXD
                              write(298347,FMT='(3I5,2E25.16)')   &
                              I1,LM1,LM2,GMATLL(LM1,LM2,I1)
                           enddo
                        enddo
                     enddo
                  endif

                  if ( NATOMIMP==1 ) THEN
                     I1=ATOMIMP(1)
                     if(OPT('KKRFLEX')) THEN
                       IREC = IE+1
                       ILM=0
                       GIMP=(0.e0,0.e0) ! complex*8
                       do LM2=1,LMMAXD
                          do LM1=1,LMMAXD
                             ILM=ILM+1
                             GIMP(ILM)=GMATLL(LM1,LM2,I1)
                          enddo
                       enddo
                       write(888,REC=IREC) GIMP
                     endif
                     if (OPT('WRTGREEN') .and. myrank==master) THEN
                       do LM2=1,LMMAXD
                         do LM1=1,LMMAXD
                           ! writeout of green_host for WRTGREEN option
                           write(58,'((2I5),(2e17.9))') LM2, LM1, GMATLL(LM1,LM2,I1)
                         enddo
                       enddo
                     endif ! WRTGREEN
                  endif

                  if ( LCPAIJ ) then
                     if(t_cpa%dmatproj_to_file) then
                        do I1 = 1,NATYP
                           do LM2=1,LMMAXD
                              do LM1=1,LMMAXD
                                 GMAT0(LM1,LM2) = TSST(LM1,LM2,I1)
                                 W1(LM1,LM2)    = MSST(LM1,LM2,I1)
                              end do
                           end do
                           IREC = IE + IELAST*(I1-1)
                           write (71,REC=IREC) GMAT0,W1
                        end do !I1
                     else !t_cpa%dmatproj_to_file
                        irec = ie_num
                        t_cpa%dmatts(:,:,:,irec) = TSST(:,:,:)
                        t_cpa%dtilts(:,:,:,irec) = MSST(:,:,:)
                     end if !t_cpa%dmatproj_to_file
                  end if  ! ( LCPAIJ )

               endif               ! ( .NOT.(OPT('qdos    ').AND.(IQDOSRUN.EQ.0)) )
 220        continue               ! IQ = 1,NQ                                         ! qdos ruess

            if (LLY.NE.0) then                                       ! LLY
               CDOS_LLY(IE,1) =   TRALPHA(IE,1) &                    ! LLY
               - LLY_GRTR(IE,1) / VOLBZ(1) + 2d0*LLY_G0TR(IE)        ! LLY
               CDOS_LLY(IE,1) = CDOS_LLY(IE,1) / PI                  ! LLY

               CDOSREF_LLY(IE) = TRALPHAREF(IE) - LLY_G0TR(IE)       ! LLY

            endif                                                    ! LLY

#ifdef CPP_MPI
            !stop timing measurement for this ie, needed for MPIadapt
            if(MPIadapt.and.t_mpi_c_grid%myrank_ie==0) then
               call timing_stop('time_1b_ie', save_out=timings_1b(ie) )
            end if
#endif

 460     continue                  ! IE=1,IELAST
#ifdef CPP_TIMING
         if(.not.OPT('GREENIMP')) then
            if(t_inc%i_time>0) call timing_stop('main1b - calctref13')
            if(t_inc%i_time>0) call timing_pause('main1b_fourier')
            if(t_inc%i_time>0) call timing_stop('main1b_inversion')
            if(t_inc%i_time>0) call timing_stop('main1b - kloopz')
         endif
#endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! SOC
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if( NCPAFAIL .NE. 0 ) then
            if(t_inc%i_write>0) then
               write(1337,*)
               write(1337,'(1X,79(''*''),/)')
               write(1337,99019) CPATOL, NCPAFAIL,    &
                  (IECPAFAIL(IE),DBLE(EZ(IECPAFAIL(IE))),IE=1,NCPAFAIL)
               write(1337,'(1X,79(''*''),/)')
               write(1337,*)
            endif
         else
            if( NCPA .NE. 0 ) then
               if(t_inc%i_write>0) then
                  write(1337,*)
                  write(1337,99020)
                  write(1337,*)
               endif
            end if
         end if

      endif ! NEWSOSOL

      !-------------------------------------------------------------------------
      ! END of do loop over spins and energies
      !-------------------------------------------------------------------------

      ! close green_host file after write out
      if (OPT('WRTGREEN') .and. myrank==master) then
         close(58)
      endif

      close (68)
!     IF ( IGF.NE.0 ) CLOSE (88)  !no-green

      if (LLY.NE.0) then                                                 ! LLY Lloyd
         if (t_lloyd%cdos_diff_lly_to_file) then
            open (701,FILE='cdosdiff_lly.dat',FORM='FORMATTED')          ! LLY Lloyd
         end if

         if (.not.OPT('NEWSOSOL')) then
            do ISPIN = 1,NSPIN                                             ! LLY
#ifdef CPP_MPI
               ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
               ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
               do ie_num=1,ie_end
                  IE = ie_start+ie_num
#else
               do IE = 1,IELAST                                         ! LLY
                  ie_num = ie
#endif
                  if (t_lloyd%cdos_diff_lly_to_file) then
                     write(701,FMT='(10E25.16)') EZ(IE),CDOS_LLY(IE,ISPIN),   &
                     TRALPHA(IE,ISPIN),LLY_GRTR(IE,ISPIN)                  ! LLY
                  else
                     t_lloyd%cdos(ie_num,ISPIN) = CDOS_LLY(IE,ISPIN)
                  end if
               enddo ! IE                                                  ! LLY
            enddo                                                          ! LLY
         else
#ifdef CPP_MPI
            ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
            ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
            do 879 ie_num=1,ie_end
               IE = ie_start+ie_num
#else
            do 879 IE = 1,IELAST                                            ! LLY
               ie_num = ie
#endif
               if(t_lloyd%cdos_diff_lly_to_file) then
                 write(701,FMT='(10E25.16)') EZ(IE),CDOS_LLY(IE,1),  &
                 TRALPHA(IE,1),LLY_GRTR(IE,1)   ! LLY
              else
                 t_lloyd%cdos(ie_num,1) = CDOS_LLY(IE,1)
              end if

            879 continue                                                 ! LLY
         endif ! .NOT.OPT('NEWSOSOL')                                    ! LLY
         if(t_lloyd%cdos_diff_lly_to_file) close(701)                    ! LLY
      endif                                                              ! LLY

      if ( ( OPT('XCPL    ') ).AND.( ICC.LE.0 ) ) then
#ifdef CPP_TIMING
         call timing_start('main1b - tbxccpl')
#endif
         if(NQDOS/=1) stop 'QDOS option not compatible with XCPL'
         if (.NOT.OPT('NEWSOSOL')) then
            call TBXCCPLJIJ(69,IELAST,EZ,WEZ,NSPINDD,NCPA,NAEZ,NATYP,NOQ,  &
               ITOQ,IQAT,NSHELL,NATOMIMP,ATOMIMP,RATOM,                    &
               NOFGIJ,NQCALC,IQCALC,IJTABCALC,                             &
               IJTABSYM,IJTABSH,ISH,JSH,                                   &
               DSYMLL,IPRINT,NATYP,NSHELD,LMMAXD,                         &
               NPOL)
         else !.NOT.OPT('NEWSOSOL'))
            call TBXCCPLJIJDIJ(naez,natyp,lmmaxd,lmgf0d,natomimpd, &
               iemxd,thetas,phis,                            &
               natomimp,atomimp,nofgij,iqat,rclsimp,                &
               ijtabcalc,ijtabcalc_I,ijtabsh,                        &
               ijtabsym,                                             &
               ielast,ez,wez,NPOL,dsymll,noq,itoq,ncpa)
         end if
#ifdef CPP_TIMING
         call timing_stop('main1b - tbxccpl')
#endif
      end if
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
      close (69)
      close (70)
      if (OPT('FERMIOUT') .and. myrank==master) close(6801)  !fswrt
      if ( LCPAIJ .and. t_cpa%dmatproj_to_file ) close (71)

      close(37)                                    ! qdos ruess
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finished first qdos run. Now re-run the whole kkr1b program to
      ! calculate the GF for every energy (defined in inputcard) and
      ! kpoint (defined in qvec.dat)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IQDOSRUN=IQDOSRUN+1                          ! qdos ruess
      if(t_lloyd%dgref_to_file)   close(681)
      if(t_lloyd%g0tr_to_file)    close(682)
      if(t_lloyd%dtmat_to_file)   close(691)
      if(t_lloyd%tralpha_to_file) close(692)
      if (IQDOSRUN.EQ.1) goto 210                 ! qdos ruess

      !Do stuff for GREENIMP option (previously done in zulapi code)
      if (OPT('GREENIMP')) then

        ! consistency checks
        if(.not.(ielast==1.or.ielast==3)) stop 'Error: GREENIMP option only possible with 1 () or 3 () energy points in contour'
        if(ielast==1 .and.dabs(dimag(ez(1)))>10**(-10)) stop 'Error: T>0 for GREENIMP (DTMTRX writeout, IELAST==3)'
        if(ielast==3 .and.dabs(dimag(ez(1)))<10**(-10)) stop 'Error: T==0 for GREENIMP (GMATLL_GES writeout, IELAST==3)'
        ! end consistency checks

#ifdef CPP_MPI
        ! init arrays and communicate parameters of t_imp for all ranks
        ! that are not the master
        call bcast_t_imp_scalars(t_imp)
        if(myrank/=master) call init_t_imp(t_inc,t_imp)
        call bcast_t_imp_arrays(t_imp, t_inc)
#endif

        do ie=1,ielast ! big ie loop (use only for GMATLL output)
            if (IELAST.EQ.3 .and. ie==1 .and. myrank==master) then
               open(UNIT=59,FILE='GMATLL_GES',FORM='FORMATTED')
               open(UNIT=60,FILE='green_host',FORM='FORMATTED')
            endif
            allocate(DTMTRX(LMMAXD*t_imp%NATOMIMP,LMMAXD*t_imp%NATOMIMP),stat=i_stat)
            call memocc(i_stat,product(shape(DTMTRX))*kind(DTMTRX),'DTMTRX','main1b')
            DTMTRX=CZERO

            ! find DTMTRX (written out for IELAST==1), parallelized with
            ! mpi over atoms
            call TMATIMP_NEWSOLVER(IRM,NSRA-1,t_params%LMAX,t_params%IEND, &
               t_params%IRID,t_params%LPOT,t_params%NATYP,t_params%NCLEB,  &
               t_params%IPAND,t_params%IRNSD,t_params%NFUND,t_imp%IHOST,   &
               t_params%NTOTD,t_params%NSPIN,t_params%LMPOT,t_params%NCHEB,&
               t_params%LMMAXD,t_params%KORBIT,t_params%NSPOTD,            &
               t_params%IELAST,t_params%IRMIND,t_params%NPAN_EQ,           &
               t_params%NPAN_LOG,t_imp%NATOMIMP,CVLIGHT,t_params%R_LOG,    &
               t_params%IPAN,t_params%IRWS,t_params%IRMIN,                 &
               t_imp%HOSTIMP(1:t_imp%NATOMIMP),                            &
               t_imp%IPANIMP(1:t_imp%NATOMIMP),                            &
               t_imp%IRWSIMP(1:t_imp%NATOMIMP),                            &
               t_params%ATOMIMP(1:t_imp%NATOMIMP),                         &
               t_imp%IRMINIMP(1:t_imp%NATOMIMP),t_params%ICLEB,            &
               t_params%IRCUT,                                             &
               t_imp%IRCUTIMP(0:t_params%IPAND,1:t_imp%NATOMIMP),          &
               t_params%ZAT,t_imp%ZIMP(1:t_imp%NATOMIMP),       &
               t_params%R,t_params%CLEB(1,1),                              &
               t_imp%RIMP(1:t_params%IRM,1:t_imp%NATOMIMP),               &
               t_params%RCLSIMP,EZ(IE),t_params%VISP,t_imp%VISPIMP,        &
               t_params%VINS,t_imp%VINSIMP,DTMTRX)

            ! compute GMATLL_GES, on master rank only
            if (IELAST.EQ.3 .and. myrank==master) then
               CALL GREENIMP(t_imp%NATOMIMP,DTMTRX,t_params%EZ(IE))
            endif

            i_all=-product(shape(DTMTRX))*kind(DTMTRX)
            deallocate(DTMTRX, stat=i_stat)
            call memocc(i_stat,i_all,'DTMTRX','main1b')

         end do ! ie-loop

         ! done with GREENIMP option, stopping now
#ifdef CPP_MPI
         call MPI_FINALIZE(ierr)
#endif
         if(myrank==master) write(*,*) 'done with GREENIMP, stop here!'
         stop

      endif ! GREENIMP

! ----------------------------------------------------------------------
!
      if( test('rhoqtest') .and. (myrank==master) ) then
         open(9999, file='params.txt')
         write(9999,*) 2*LMMAXD, t_params%natyp
         write(9999,*) t_params%naez, t_params%nclsd, t_params%nr,  &
                       t_params%nembd1, t_params%lmax
         write(9999,*) t_params%alat, naclsmax
         close(9999)

         open(9999, file='host.txt')
         write(9999,*) t_params%rbasis(1:3,1:t_params%natyp)
         write(9999,*) t_params%rcls(1:3,1:t_params%nclsd,1:t_params%nclsd),  &
                       t_params%rr(1:3,0:t_params%nr),                       &
                       t_params%atom(1:t_params%nclsd,1:t_params%naez+t_params%nembd1)
         write(9999,*) t_params%cls(1:t_params%naez+t_params%nembd1),          &
                       t_params%ezoa(1:t_params%nclsd,1:t_params%naez+t_params%nembd1),&
                       t_params%nacls(1:t_params%nclsd)
         close(9999)
      endif

      ! clean up allocations
      i_all=-product(shape(LEFTTINVLL))*kind(LEFTTINVLL)
      deallocate(LEFTTINVLL,stat=i_stat)
      call memocc(i_stat,i_all,'LEFTTINVLL','main1b')
      i_all=-product(shape(RIGHTTINVLL))*kind(RIGHTTINVLL)
      deallocate(RIGHTTINVLL,stat=i_stat)
      call memocc(i_stat,i_all,'RIGHTTINVLL','main1b')
      i_all=-product(shape(BZKP))*kind(BZKP)
      deallocate(BZKP,stat=i_stat)
      call memocc(i_stat,i_all,'BZKP','main1b')
      if(opt('qdos    ')) then
         i_all=-product(shape(qvec))*kind(qvec)
         deallocate(qvec, stat=i_stat)
         call memocc(i_stat,i_all,'qvec','main1b')
      end if
      i_all=-product(shape(ISH))*kind(ISH)
      deallocate(ISH,stat=i_stat)
      call memocc(i_stat,i_all,'ISH','main1b')
      i_all=-product(shape(JSH))*kind(JSH)
      deallocate(JSH,stat=i_stat)
      call memocc(i_stat,i_all,'JSH','main1b')

      if(t_inc%i_write>0) write (1337,'(79(1H=),/,30X,"< KKR1b finished >",/,79(1H=),/)')

2100 FORMAT('(/,79(1H=),/,5X," Inversion algorithm used : ",A,/,79(1H=),/)')
99019 FORMAT('(/,1X,79(*),/," tolerance for CPA-cycle:",F15.7,/," CPA not converged for",I3," energies:",/,3(" E:",I3,F7.4,:,2X))')
99020 FORMAT('(/,1X,79(*),/,25X,"no problems with","  CPA-cycle ",/,1X,79(*),/)')

   end subroutine main1b

end module MOD_MAIN1B
