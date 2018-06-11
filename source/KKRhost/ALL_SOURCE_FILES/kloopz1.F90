!-------------------------------------------------------------------------------
! SUBROUTINE: KLOOPZ1_QDOS
!> @note
!> - Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine KLOOPZ1_QDOS(NR,NEMBD1,LMMAXD,LMGF0D,LMAX,NREF,ERYD,GMATLL,INS,ALAT,IE,IGF,  &
   NSHELL,NAEZ,NOFKS,VOLBZ,BZKP,VOLCUB,CLS,NACLS,     &
   NACLSMAX,NCLS,RR,RBASIS,EZOA,ATOM,RCLS,ICC,        &
   GINP,IDECI,LEFTTINVLL,RIGHTTINVLL,VACFLAG,         &
   NLBASIS,NRBASIS,FACTL,NATOMIMP,NSYMAT,DSYMLL,      &
   RATOM,RROT,NSH1,NSH2,IJTABSYM,IJTABSH,ICHECK,      &
   INVMOD,REFPOT,TREFLL,TSST,MSST,CFCTOR,             &
   CFCTORINV,CREL,RC,RREL,SRREL,IRREL,NRREL,DROTQ,    &
   SYMUNITARY,KMROT,NATYP,NCPA,ICPA,ITCPAMAX,         &
   CPATOL,NOQ,IQAT,ITOQ,CONC,IPRINT,ICPAFLAG,         &
   ISPIN,NSPIN,                                       &
   TQDOS,IQDOSRUN,                                    &  !qdos ruess
   DTREFLL,DTMATLL,DGINP,LLY_GRTR,TRACET,LLY)            ! LLY Lloyd

   use mod_types, only: t_inc
   use mod_mympi, only: myrank, master
   use global_variables
   use Constants
   use Profiling
      Use mod_datatypes, Only: dp

   implicit none
   !
   ! .. Parameters
   integer :: LINMAX
   parameter (LINMAX=1)
   real (kind=dp) :: TOLMSSQ
   parameter ( TOLMSSQ=1.0D-6 )
   ! .. Input variables
   integer, intent(in) :: NR        !< Number of real space vectors rr
   integer, intent(in) :: IE
   integer, intent(in) :: LLY       !< LLY <> 0 => use Lloyd formula
   integer, intent(in) :: IGF       !< Do not print or print (0/1) the KKRFLEX_* files
   integer, intent(in) :: INS       !< 0 (MT), 1(ASA), 2(Full Potential)
   integer, intent(in) :: ICC       !< Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
   integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
   integer, intent(in) :: NCPA      !< NCPA = 0/1 CPA flag
   integer, intent(in) :: NCLS      !< Number of reference clusters
   integer, intent(in) :: NREF      !< Number of diff. ref. potentials
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: KMROT     !< 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
   integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NOFKS
   integer, intent(in) :: IDECI
   integer, intent(in) :: ISPIN
   integer, intent(in) :: NSPIN     !< Counter for spin directions
   integer, intent(in) :: NEMBD1    !< NEMB+1
   integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
   integer, intent(in) :: LMGF0D    !< (LMAX+1)**2
   integer, intent(in) :: NSYMAT
   integer, intent(in) :: INVMOD    !< Inversion scheme
   integer, intent(in) :: IPRINT
   integer, intent(in) :: NLBASIS   !< Number of basis layers of left host (repeated units)
   integer, intent(in) :: NRBASIS   !< Number of basis layers of right host (repeated units)
   integer, intent(in) :: NATOMIMP  !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
   integer, intent(in) :: NACLSMAX
   integer, intent(in) :: IQDOSRUN  !< qdos ruess: counts qdos run
   real (kind=dp), intent(in) :: ALAT   !< Lattice constant in a.u.
   real (kind=dp), intent(in) :: VOLBZ
   real (kind=dp), intent(in) :: CPATOL !< Convergency tolerance for CPA-cycle
   complex (kind=dp), intent(in) :: ERYD
   complex (kind=dp), intent(in) :: CFCTOR
   complex (kind=dp), intent(in) :: CFCTORINV
   ! .. Input arrays
   integer, dimension(*), intent(in)         :: CLS      !< Cluster around atomic sites
   integer, dimension(NAEZ), intent(in)      :: NOQ      !< Number of diff. atom types located
   integer, dimension(*), intent(in)         :: NSH1     !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
   integer, dimension(*), intent(in)         :: NSH2     !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
   integer, dimension(NAEZ), intent(in)      :: ICPA     !< ICPA = 0/1 site-dependent CPA flag
   integer, dimension(NATYP), intent(in)     :: IQAT     !< The site on which an atom is located on a given site
   integer, dimension(*), intent(in)         :: NACLS    !< Number of atoms in cluster
   integer, dimension(0:NSHELD), intent(in)  :: NSHELL   !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
   integer, dimension(*), intent(in)         :: REFPOT   !< Ref. pot. card  at position
   integer, dimension(*), intent(in)         :: IJTABSH  !< Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
   integer, dimension(*), intent(in)         :: IJTABSYM !< Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
   integer, dimension(NACLSD,*), intent(in)                    :: ATOM  !< Atom at site in cluster
   integer, dimension(NACLSD,*), intent(in)                    :: EZOA  !< EZ of atom at site in cluster
   integer, dimension(NATYP,NAEZ), intent(in)                  :: ITOQ
   integer, dimension(2,LMMAXD), intent(in)                    :: NRREL
   integer, dimension(NAEZ/NPRINCD,NAEZ/NPRINCD), intent(in)   :: ICHECK
   integer, dimension(2,2,LMMAXD), intent(in) :: IRREL
   real (kind=dp), dimension(NATYP), intent(in)  :: CONC        !< Concentration of a given atom
   real (kind=dp), dimension(KPOIBZ), intent(in) :: VOLCUB
   real (kind=dp), dimension(3,0:NR), intent(in)    :: RR       !< Set of real space vectors (in a.u.)
   real (kind=dp), dimension(3,KPOIBZ), intent(in)  :: BZKP
   real (kind=dp), dimension(3,*), intent(in)       :: RATOM
   real (kind=dp), dimension(3,*), intent(in)       :: RBASIS   !< Position of atoms in the unit cell in units of bravais vectors
   real (kind=dp), dimension(3,NACLSD,*), intent(in)   :: RCLS  !< Real space position of atom in cluster
   real (kind=dp), dimension(48,3,*), intent(in)       :: RROT
   complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in) :: RC     !< NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
   complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in) :: CREL   !< Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
   complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in) :: RREL   !< Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
   complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in) :: FACTL
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NATYP), intent(in)  :: TSST
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NATYP), intent(in)  :: MSST
   complex (kind=dp), dimension(2,2,LMMAXD), intent(in)           :: SRREL
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NAEZ), intent(in)   :: TQDOS  ! qdos : Read-in inverse t-matrix
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NAEZ), intent(in)   :: DROTQ   !< Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NREF), intent(in)   :: TREFLL
   complex (kind=dp), dimension(LMMAXD,LMMAXD,*), intent(in)      :: DSYMLL
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NREF), intent(in)   :: DTREFLL !< LLY Lloyd dtref/dE
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NAEZ), intent(in)   :: DTMATLL  ! LLY  dt/dE (should be av.-tmatrix in CPA)
   complex (kind=dp), dimension(LMGF0D*NACLSMAX,LMGF0D,NCLS), intent(in) :: GINP !< Cluster GF (ref syst.)
   complex (kind=dp), dimension(LMGF0D*NACLSMAX,LMGF0D,NCLS), intent(in) :: DGINP !< LLY Lloyd Energy derivative of GINP
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NEMBD1,NSPIN), intent(in) :: LEFTTINVLL
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NEMBD1,NSPIN), intent(in) :: RIGHTTINVLL
   logical, dimension(2), intent(in) :: VACFLAG
   logical, dimension(*), intent(in) :: SYMUNITARY !< unitary/antiunitary symmetry flag
   ! .. Output variables
   integer, intent(out) :: ICPAFLAG
   ! .. In/Out variables
   integer, intent(inout) :: ITCPAMAX  !< Max. number of CPA iterations
   complex (kind=dp), intent(inout) :: TRACET   !< \f$Tr\left[ (t-tref)^{-1} \frac{d(t-tref)}{dE} \right]\f$
   complex (kind=dp), intent(inout) :: LLY_GRTR !< Trace Eq.5.38 PhD Thiess (k-integrated)! LLY Lloyd
   complex (kind=dp), dimension(LMMAXD,LMMAXD,*), intent(inout) :: GMATLL  !< GMATLL = diagonal elements of the G matrix (system)
   ! .. Local Scalars
   integer :: i_stat,i_all
   integer :: IH,LM1,LM2,NS,NSDIA,ICALL,IREC
   integer :: IQ,JQ,IT,I,J,IQTAU,ICPASTART,ITCPA,IU,NSMAX
   real (kind=dp) :: CPAERRL,CPAERR,CPACORR,CPACHNG
   complex (kind=dp) :: EZ,CNSYMAT,TAUVBZ
   logical :: LDIA
   character(len=4) :: STR4
   character(len=10) :: STR10
   ! .. Local Arrays ..
   integer, dimension(NAEZ)      :: NKMQ
   integer, dimension(NAEZ)      :: NLINQ
   integer, dimension(NSYMAXD)   :: ISUMG
   integer, dimension(LINMAX)    :: IKM1LIN
   integer, dimension(LINMAX)    :: IKM2LIN
   integer, dimension(NSYMAXD,NAEZ)    :: ISUMQ
   integer, dimension(NSYMAXD,NATYP)   :: ISUMT
   complex (kind=dp), dimension(LMMAXD,LMMAXD)  :: GLL
   complex (kind=dp), dimension(LINMAX,NATYP)   :: TAUTLIN
   ! .. Effective (site-dependent) Delta_t^(-1) matrix ..
   complex (kind=dp), dimension(LMMAXD,LMMAXD) :: XC
   complex (kind=dp), dimension(LMMAXD,LMMAXD) :: W1
   complex (kind=dp), dimension(LMMAXD,LMMAXD) :: W2
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NAEZ) :: MSSQ
   ! .. Local allocatable arrays
   complex (kind=dp), dimension(:,:,:), allocatable :: DMSSQ
   complex (kind=dp), dimension(:,:,:), allocatable :: TAUDELQ
   complex (kind=dp), dimension(:,:,:), allocatable :: TAUDELT
   complex (kind=dp), dimension(:,:,:,:), allocatable :: GS

#ifdef CPP_MPI
   integer :: irec0
#endif
   ! .. External Functions/Subroutines
   logical :: TEST,OPT
   external :: TEST,OPT,KKRMAT01,ROTGLL,ROTATE,SYMETRMAT,ZCOPY,ZGEMM,ZSCAL
   ! .. Intrinsic Functions
   intrinsic :: DBLE
   ! .. Data statement
   data ICALL / 0 /
   ! .. Save statement
   save :: ICALL,ISUMG,CNSYMAT

   !----------------------------------------------------------------------------
   if ( TEST('flow    ') ) then
      WRITE (1337,*)'>>> KLOOPZ1: invert delta_t and do Fourier transformation'
   endif
   ICALL = ICALL + 1
   !
   ! Reinitialise the ICALL counter for second run of kloopz1                    ! qdos ruess
   if ( (OPT('qdos    ').AND.(IQDOSRUN.EQ.1)) ) ICALL = 1                        ! qdos ruess
   !
   if (OPT('FERMIOUT') .and. myrank==master) then                                ! fswrt
      write(6801,'(A)') 'energy(ie):'                                            ! fswrt
      write(6801,'(2ES25.16)') ERYD                                              ! fswrt
   end if                                                                        ! fswrt
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! The arrays ISUM are used in the symmetrisation routine SYMETRMAT
   ! Symmetrising single-site : same matrix for each symmetry
   ! Symmetrising G matrix    : pick G(ISYM) for symmetry ISYM
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if ( ICALL.EQ.1 ) then
      CNSYMAT = CONE/DBLE(NSYMAT)
      do IU = 1,NSYMAXD
         do IT = 1,NATYP
            ISUMT(IU,IT) = IT
         end do
         do IQ = 1,NAEZ
            ISUMQ(IU,IQ) = IQ
         end do
         ISUMG(IU) = IU
      end do
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! TSST in the LOCAL frame is used to set up
   ! MSST = (TSST-TREF)^(-1) in the LOCAL frame to be used in <CPAMILLSX> and < PROJTAU > below
   !
   ! MSSQ = the inverse of the effective (on-site) Delta_t matrix in the GLOBAL frame;
   !        the Average T-matrix Approximation (ATA) (ICPASTART=1) is used
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ICPASTART = 1
   call MSSINIT(NCPA,ICPASTART,TSST,MSST,MSSQ,TREFLL,DROTQ, &
      REFPOT,IQAT,ITOQ,NOQ,CONC,                            &
      KMROT,NATYP,NAEZ,LMMAXD)  ! nref was taken out of calling list 1.2.2012
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !    VIRTUAL ATOMS:
   !    Be careful! in case of  OPT('VIRATOMS')==1 MSSQ is the Tmatrix
   !    not the inverse T-matrix!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !----------------------------------------------------------------------------
   ! Output now:
   ! MSSQ is the (Delta_t)^(-1) matrix in the GLOBAL frame
   !      and refers to a site IQ occupied by NOQ(IQ) atoms having
   !      the occupancies CONC(ITOQ(1..NOQ(IQ))
   !
   ! TSST is the t-matrix in the LOCAL frame
   !      for each of the NATYP atoms
   ! MSST is the Delta_t^(-1) matrix in the LOCAL frame
   !      for each of the NATYP atoms
   !----------------------------------------------------------------------------
   !
   !     Write gref to the TBkkr-container-file
   if (OPT('FERMIOUT') .and. myrank==master) then                                ! fswrt
      write(6801,'(A)') 'GINP(ie):'                                              ! fswrt
      do I=1,NCLS                                                                ! fswrt
         do LM2=1,LMGF0D                                                         ! fswrt
            do LM1=1,NACLSMAX*LMGF0D                                             ! fswrt
               write(6801,'(2ES25.16)')  GINP(LM1,LM2,I)                         ! fswrt
            end do                                                               ! fswrt
         end do                                                                  ! fswrt
      end do                                                                     ! fswrt
   end if                                                                        ! fswrt
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! BEGIN CPA - LOOP  (if required)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !----------------------------------------------------------------------------
   !  ikm1lin,ikm2lin,nlinq --> dummy settings for <PROJTAU>
   !----------------------------------------------------------------------------
   IKM1LIN(1) = 1
   IKM2LIN(1) = 1
   do IQ = 1,NAEZ
      NKMQ(IQ) = LMMAXD
      NLINQ(IQ) = 1
   end do
   !
   CPAERRL = 1.0D+6
   ITCPA = 0
   EZ = ERYD
   !
   if ( NCPA.NE.0 ) then
      allocate(DMSSQ(LMMAXD,LMMAXD,NAEZ),stat=i_stat)
      call memocc(i_stat,product(shape(DMSSQ))*kind(DMSSQ),'DMSSQ','kloopz1')
      DMSSQ(:,:,:) = CZERO
      call CINIT(LMMAXD*LMMAXD*NAEZ,DMSSQ)
   end if
   !----------------------------------------------------------------------------
   ! < GLL2K >  incorporated now
   !----------------------------------------------------------------------------
   TAUVBZ = 1.D0/VOLBZ
   !----------------------------------------------------------------------------
   ! Convert inverted delta_t-matrices to p.u.
   !----------------------------------------------------------------------------
   do I = 1,NAEZ
      !
      if ( .not. OPT('VIRATOMS') ) then
         call ZSCAL(LMMAXD*LMMAXD,CFCTOR,MSSQ(1,1,I),1)
      else
         call ZSCAL(LMMAXD*LMMAXD,CFCTORINV,MSSQ(1,1,I),1)
      end if   !( .not. OPT('VIRATOMS') )
      !
      if (.NOT.OPT('NEWSOSOL') ) then
         if (KMROT.EQ.0) then
            do LM2 = 1,LMMAXD
               do LM1 = 1,LM2
                  MSSQ(LM1,LM2,I) = 0.5D0 *( MSSQ(LM1,LM2,I) + MSSQ(LM2,LM1,I) )
                  MSSQ(LM2,LM1,I) = MSSQ(LM1,LM2,I)
               end do
            end do
         end if
      endif
   end do
   !
   do I = 1,NATYP
      call ZSCAL(LMMAXD*LMMAXD,CFCTOR,MSST(1,1,I),1)
   end do
   !----------------------------------------------------------------------------
   ! Symmetrise the delta_t^(-1) matrices
   !----------------------------------------------------------------------------
   if (.NOT.OPT('NEWSOSOL')) then
      if ( ( KREL.EQ.1 ).OR.( INS.NE.0 ) ) then
         do IQ = 1,NAEZ
            call SYMETRMAT(NSYMAT,CNSYMAT,DSYMLL,SYMUNITARY,MSSQ, &
               ISUMQ(1,IQ),MSSQ(1,1,IQ),LMMAXD)
            if ( KMROT.EQ.0 ) then
               do I = 1,NOQ(IQ)
                  IT = ITOQ(I,IQ)
                  call SYMETRMAT(NSYMAT,CNSYMAT,DSYMLL,SYMUNITARY,MSST, &
                     ISUMT(1,IT),MSST(1,1,IT),LMMAXD)
               end do
            end if
         end do
      end if
   endif
   !
   EZ = EZ*CFCTOR*CFCTOR
   NSDIA = MAX(NAEZ,NATYP)
   !
   allocate(TAUDELQ(LMMAXD,LMMAXD,NSHELL(0)),stat=i_stat)
   call memocc(i_stat,product(shape(TAUDELQ))*kind(TAUDELQ),'TAUDELQ','kloopz1')
   TAUDELQ(:,:,:) = CZERO
   !----------------------------------------------------------------------------
   ! Jonathan Chico: This possibly can be removed by a do while
   100  continue
   !----------------------------------------------------------------------------
   ITCPA = ITCPA + 1
   !
   allocate(GS(LMMAXD,LMMAXD,NSYMAXD,NSHELL(0)),stat=i_stat)
   call memocc(i_stat,product(shape(GS))*kind(GS),'GS','kloopz1')
   GS(:,:,:,:) = CZERO
   !
   ! copy read-in cpa t-matrix but only after fort.37 was created in first run   !qdos ruess
   if (OPT('readcpa ').OR.(OPT('qdos    ').AND.(IQDOSRUN.EQ.1))) then            !qdos ruess
      MSSQ(:,:,:) = TQDOS(:,:,:) ! lmmaxd,lmmaxd,naez                            !qdos ruess
   endif
   !
   call KKRMAT01(NR,NREF,LMAX,LMGF0D,LMMAXD,BZKP,NOFKS,GS,VOLCUB,     &
      MSSQ,RROT,NSHELL(0),NSDIA,ALAT,NSYMAT,NAEZ,CLS,NACLS,NACLSMAX,    &
      RR,EZOA,ATOM,NSH1,NSH2,GINP,RBASIS,RCLS,LEFTTINVLL(1,1,1,ISPIN),  &
      RIGHTTINVLL(1,1,1,ISPIN),VACFLAG,NLBASIS,NRBASIS,FACTL,ICHECK,    &
      INVMOD,IDECI,SRREL,IRREL,NRREL,DTREFLL,DTMATLL,DGINP,REFPOT,      &
      LLY_GRTR,TRACET,CFCTOR,LLY) ! LLY
   !
   NSMAX = NSHELL(0)
   !----------------------------------------------------------------------------
   ! NS=1,NSMAX
   !----------------------------------------------------------------------------
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! For qdos calculation, do not symmetrize the GF matrix, but call routine
   ! symetrmat anyway because of factor tauvbz:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do NS = 1,NSMAX
      !-------------------------------------------------------------------------
      ! symmetrise GS, get GLL as sum over all symmetry wedges of GS
      ! (isumg(i) = i, see above)
      !-------------------------------------------------------------------------
      call SYMETRMAT(NSYMAT,TAUVBZ,DSYMLL,SYMUNITARY,GS(1,1,1,NS),ISUMG,GLL,LMMAXD)
      !
      if (NS.LE.NATYP) then
         IQTAU = IQAT(NS)
      else
         IQTAU = NS
      end if
      !
      TAUDELQ(1:LMMAXD,1:LMMAXD,IQTAU) = -GLL(1:LMMAXD,1:LMMAXD)
   end do
   !----------------------------------------------------------------------------
   ! NS=1,NSMAX
   !----------------------------------------------------------------------------
   i_all=-product(shape(GS))*kind(GS)
   deallocate(GS, stat=i_stat)
   call memocc(i_stat,i_all,'GS','kloopz1')
   !----------------------------------------------------------------------------
   if ( NCPA.GT.0 ) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do one CPA iteration, the output is a new MSSQ = (Delta_t)^(-1)
      ! in the GLOBAL frame
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ICPAFLAG = 0
      if (OPT('readcpa ').OR.(OPT('qdos    ').AND.IQDOSRUN.GT.0)) then
         ! copy read-in cpa t-matrix in second run
         call CPAMILLSX(ITCPA,CPAERR,CPACORR,CPACHNG,IPRINT,ICPA, &
            NAEZ,NKMQ,NOQ,ITOQ,CONC,MSSQ,MSST,TAUDELQ,DMSSQ,      &
            KMROT,DROTQ,NATYP,NAEZ,LMMAXD)
         MSSQ(:,:,:) = TQDOS(:,:,:) ! lmmaxd,lmmaxd,naez   ! qdos
         ITCPAMAX = 0
         CPAERR = 0.D0
      else
         call CPAMILLSX(ITCPA,CPAERR,CPACORR,CPACHNG,IPRINT,ICPA, &
            NAEZ,NKMQ,NOQ,ITOQ,CONC,MSSQ,MSST,TAUDELQ,DMSSQ,      &
            KMROT,DROTQ,NATYP,NAEZ,LMMAXD)
         !----------------------------------------------------------------------
         ! Symmetrise m-CPA
         !----------------------------------------------------------------------
         if (.NOT.OPT('NEWSOSOL')) then
            do IQ = 1,NAEZ
               if ( ICPA(IQ).NE.0 ) then
                  call SYMETRMAT(NSYMAT,CNSYMAT,DSYMLL,  &
                     SYMUNITARY,MSSQ,ISUMQ(1,IQ),MSSQ(1,1,IQ),LMMAXD)
               endif
            end do
         endif
         !
         if ( IPRINT.GE.1 .and. (t_inc%i_write>0))  then
            write (1337,99004) CPAERR,CPACORR,CPACHNG
         endif
         !
         if ( CPAERR.LE.CPATOL ) then
            if ( IPRINT.GT.0 .and. (t_inc%i_write>0)) then
               write (1337,99001) ITCPA,CPAERR,CPACORR,CPACHNG
            endif
         else if ( ITCPA.GT.ITCPAMAX) then
            if(t_inc%i_write>0) then
               write (1337,99002) ITCPA,CPAERR,CPACORR,CPACHNG
            endif
            ICPAFLAG = 1
         else if ( CPAERR.GT.20*CPAERRL ) then
            if(t_inc%i_write>0) write (1337,99003) ITCPA
            if(t_inc%i_write>0) write (1337,99004) CPAERR,CPACORR,CPACHNG
            ICPAFLAG = 2
         else
            !-------------------------------------------------------------------
            ! Go to the next CPA iteration if not converged
            !-------------------------------------------------------------------
            CPAERRL = CPAERR
            goto 100
         end if              !   ( CPAERR.LE.CPATOL )
      end if                 !  (OPT('readcpa ').OR.OPT('qdos    '))
      i_all=-product(shape(DMSSQ))*kind(DMSSQ)
      deallocate(DMSSQ, stat=i_stat)
      call memocc(i_stat,i_all,'DMSSQ','kloopz1')
      !
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! END CPA - LOOP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !----------------------------------------------------------------------------
   !  The inverse of the Delta_t(CPA)-matrix is calculated, now write
   !  this on a file in case of decimation output.
   !  attention: new format Dec. 2004
   !----------------------------------------------------------------------------
   !
   !  In first qdos run write out t matrix which is read in for the calculation for every k point
   if ( OPT('deci-out').OR.(IQDOSRUN.EQ.0) ) then  ! qdos ruess
      do IH = 1,NAEZ
#ifdef CPP_MPI
         irec0 = LMMAXD**2*(ih-1)+LMMAXD**2*naez*(ie-1)+ &
         LMMAXD**2*t_inc%IELAST*naez*(ispin-1)
         if (OPT('deci-out'))  write (37,99005) IE,ERYD,IH
#else
         write (37,99005) IE,ERYD,IH
#endif
         do LM1 = 1,LMMAXD
            do LM2 = 1,LMMAXD
               if ( LM1.EQ.LM2 ) then
#ifdef CPP_MPI
                  if (.not.OPT('deci-out')) then
                     irec = irec0 + LM2 + LMMAXD*(LM1-1)
                     write(37,rec=irec) MSSQ(LM1,LM2,IH)*CFCTORINV
                  else
                     write(37,99006) LM1,LM2,MSSQ(LM1,LM2,IH)*CFCTORINV
                  end if
#else
                  write(37,99006) LM1,LM2,MSSQ(LM1,LM2,IH)*CFCTORINV
#endif
               else
#ifdef CPP_MPI
                  irec = irec0 + LM2 + LMMAXD*(LM1-1)
                  if ( ABS(MSSQ(LM1,LM2,IH)/MSSQ(LM1,LM1,IH)).GT.TOLMSSQ ) then
                     if (.not.OPT('deci-out')) then
                        write(37,rec=irec) MSSQ(LM1,LM2,IH)*CFCTORINV
                     else
                        write(37,99006)LM1,LM2,MSSQ(LM1,LM2,IH)*CFCTORINV
                     end if
                  else
                     if (.not.OPT('deci-out')) then
                        write(37,rec=irec) (0.d0, 0.d0)
                     end if
                  end if
#else
                  if ( ABS(MSSQ(LM1,LM2,IH)/MSSQ(LM1,LM1,IH)).GT.TOLMSSQ )  then
                     write(37,99006)LM1,LM2,MSSQ(LM1,LM2,IH)*CFCTORINV
                  endif
#endif
               end if
            end do
         end do
      end do
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! calculate the component-projected site-diagonal
   !  TAU-matrices TAUDELT(IT).
   !     - there are NSMAX = NAEZ/NATYP (NCPA=0/1) site-diagonal
   !       elements TAUDELQ(1..NSMAX) in the GLOBAL (crystal) frame of
   !       reference
   !     - they are NSMAX site-diagonal elements projected on atomic types
   !       in the array TAUDELT(1..NSMAX) in the LOCAL frame of reference
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   allocate(TAUDELT(LMMAXD,LMMAXD,NSHELL(0)),stat=i_stat)
   call memocc(i_stat,product(shape(TAUDELT))*kind(TAUDELT),'TAUDELT','kloopz1')
   TAUDELT(:,:,:) = CZERO
   if ( NCPA.GT.0) then
      NSMAX = NATYP
      call PROJTAU(ICPAFLAG,CPACHNG,KMROT,.FALSE.,.FALSE.,9,CZERO,   &
         NATYP,NAEZ,NKMQ,MSST,MSSQ,NLINQ,IQAT,CONC,TAUDELQ,TAUDELT,  &
         TAUTLIN,IKM1LIN,IKM2LIN,DROTQ,NSHELD,NAEZ,LMMAXD,LINMAX)
   else
      NSMAX = NAEZ
      do NS=1,NSMAX
         call ZCOPY(LMMAXD*LMMAXD,TAUDELQ(1,1,NS),1,TAUDELT(1,1,NS),1)
         if ( KMROT.NE.0 ) then
            do J = 1,LMMAXD
               call ZCOPY(LMMAXD,TAUDELT(1,J,NS),1,W1(1,J),1)
            end do
            call ROTATE(W1,'G->L',TAUDELT(1,1,NS),LMMAXD,DROTQ(1,1,NS),LMMAXD)
         end if
      end do
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! off-diagonal elements (if present) are in the same array
   ! TAUDELQ in the range (NSMAX+1,..,NSHELL(0)).
   ! The G_ij's are left UNPROJECTED and in the GLOBAL frame of
   ! reference (remember this for further use!), i.e. in case of CPA,
   ! G_ij(CPA) is stored. The component-projected G-elements can
   ! be obtained through the projection matrices (see below)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do NS = NSMAX+1, NSHELL(0)
      call ZCOPY(LMMAXD*LMMAXD,TAUDELQ(1,1,NS),1,TAUDELT(1,1,NS),1)
   end do
   !----------------------------------------------------------------------------
   ! Switch from TAU to G multiplying with MSST on the site-diagonal,
   ! system-matrix (1..NSMAX) and with MSSQ on the rest of the elements
   !----------------------------------------------------------------------------
   !
   if ( TEST('Gmat    ') .and. (t_inc%i_write>0)) then
      write (1337,'(/,4X,70("-"),/,4X,A,I4)')'system G_ii matrix for i = 1,',NSMAX
   endif
   do NS = 1,NSHELL(0)
      !
      LDIA = (DABS(RATOM(1,NS)**2+RATOM(2,NS)**2+RATOM(3,NS)**2).LT.1D-6)
      !-------------------------------------------------------------------------
      ! GLL = -TAU
      !-------------------------------------------------------------------------
      call ZCOPY(LMMAXD*LMMAXD,TAUDELT(1,1,NS),1,GLL,1)
      call ZSCAL(LMMAXD*LMMAXD,-CONE,GLL,1)
      !-------------------------------------------------------------------------
      ! XC = (Delta_t(I))^(-1) * GLL
      !-------------------------------------------------------------------------
      if ( NS.LE.NSMAX ) then
         !----------------------------------------------------------------------
         ! deal with the site-diagonal GFUN of the system - use MSST(I)
         ! on both sides of TAU multiplication
         !----------------------------------------------------------------------
         do J = 1,LMMAXD
            call ZCOPY(LMMAXD,MSST(1,J,NS),1,W1(1,J),1)
         end do
         call ZCOPY(LMMAXD*LMMAXD,W1,1,W2,1)
      else
         !----------------------------------------------------------------------
         ! deal with the Gij of the cluster - use MSSQ(IQ) and MSSQ(JQ)
         !----------------------------------------------------------------------
         IQ = NSH1(NS)
         JQ = NSH2(NS)
         do J = 1,LMMAXD
            call ZCOPY(LMMAXD,MSSQ(1,J,IQ),1,W1(1,J),1)
         end do
         do J = 1,LMMAXD
            call ZCOPY(LMMAXD,MSSQ(1,J,JQ),1,W2(1,J),1)
         end do
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! NS.LT.NSMAX
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( .not. OPT('VIRATOMS') ) then
         if ( .not. TEST('testgmat') ) then
            CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,W1,   &
               LMMAXD,GLL,LMMAXD,CZERO,XC,LMMAXD)
            !
            if ( LDIA ) then
               !----------------------------------------------------------------
               ! GLL = - (Delta_t(I))^(-1) - (Delta_t(I))^(-1) * GLL * (Delta_t(I))^(-1)
               !----------------------------------------------------------------
               call ZCOPY(LMMAXD*LMMAXD,W2,1,GLL,1)
               call ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,-CONE,XC,LMMAXD,  &
                  W2,LMMAXD,-CONE,GLL,LMMAXD)
            else
               !----------------------------------------------------------------
               ! GLL =  - (Delta_t(I))^(-1)  * GLL * (Delta_t(J))^(-1)
               !----------------------------------------------------------------
               call ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,-CONE,XC,LMMAXD,  &
                  W2,LMMAXD,CZERO,GLL,LMMAXD)
            end if
         end if !( .not. TEST('testgmat') ) THEN
      end if !( .not. OPT('VIRATOMS') ) THEN
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDIA
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !-------------------------------------------------------------------------
      ! GMATLL = GLL/RFCTOR
      !-------------------------------------------------------------------------
      GMATLL(1:LMMAXD,1:LMMAXD,NS) = GLL(1:LMMAXD,1:LMMAXD)*CFCTORINV
      !
      if ( ( NS.LE.NSMAX ).AND.( TEST('Gmat    ') ) ) then
         write(STR4,'(I4)') NS
         STR10 = '   i ='//STR4(1:4)
         call CMATSTR(STR10,10,GMATLL(1,1,NS),LMMAXD,LMMAXD, &
            2*KREL+1,2*KREL+1,0,1d-8,6)
      end if
   end do
   !----------------------------------------------------------------------------
   if ( TEST('Gmat    ') .and. (t_inc%i_write>0)) write (1337,'(/,4X,70("-"))')
   !----------------------------------------------------------------------------
   ! it calculates the rest of the G n n' matrix from the
   ! knowledge of the representative pairs (shells) using the
   ! real space symmetries (added 23.2.2000)
   !----------------------------------------------------------------------------
   if ( ICC.GT.0 ) then
      IREC = 1+IE+t_inc%IELAST*(ISPIN-1) ! added for mpi run
      call ROTGLL(GMATLL,NATOMIMP,IJTABSYM,IJTABSH,  &
         DSYMLL,SYMUNITARY,IGF,RC,CREL,RREL,KREL,LMMAXD,IREC)
   end if
   !       IF ( OPT('VIRATOMS') ) THEN
   !         write(*,*) 'VIRTUAL ATOM OPTION : stop calculation '
   !         stop
   !       END IF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! In the case of NCPA.NE.0 and NSHELL(0).GT.NATYP the projection
   ! matrices DMAT and DTIL which are used to get
   !                     ij            ij    _
   !                    G     =  D  * G    * D
   !                     ab       a    CPA    b
   ! - with a/b the atom of type a/b sitting on site i/j - are calculated
   !   and stored for later use.  the allocated work space for
   !   TSST (DMAT) and MSST (DTIL) is used.
   !   for an atom having occupancy 1, DMAT/DTIL = unit matrix
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (( NCPA.NE.0 ).AND. ( NSHELL(0).GT.NSMAX )) then
      do IT = 1,NATYP
         IQ = IQAT(IT)
         IH = REFPOT(IQ)
         !
         if ( KMROT.NE.0 ) then
            call ROTATE(TSST(1,1,IT),'L->G',W1,LMMAXD,DROTQ(1,1,IQ),LMMAXD)
         else
            call ZCOPY(LMMAXD*LMMAXD,TSST(1,1,IT),1,W1,1)
         end if
         !
         do J = 1,LMMAXD
            call ZAXPY(LMMAXD,-CONE,TREFLL(1,J,IH),1,W1(1,J),1)
         end do
         !
         call GIJDMAT(TAUDELQ(1,1,IQ),W1,MSSQ(1,1,IQ),   &
            TSST(1,1,IT),MSST(1,1,IT),CFCTORINV,IPRINT,  &
            IE,IT,KREL,LMMAXD)
      end do
   end if
   !----------------------------------------------------------------------------
   if ( TEST('flow    ') .and. (t_inc%i_write>0)) write (1337,*) '<<< KLOOPZ1'
   !
   i_all=-product(shape(TAUDELQ))*kind(TAUDELQ)
   deallocate(TAUDELQ, stat=i_stat)
   call memocc(i_stat,i_all,'TAUDELQ','kloopz1')
   i_all=-product(shape(TAUDELT))*kind(TAUDELT)
   deallocate(TAUDELT, stat=i_stat)
   call memocc(i_stat,i_all,'TAUDELT','kloopz1')
   !
   return
   !
   99001 format (' CPA converged after',I3,' iterations   ','ERR:',F9.6,   &
      ' CORR:',F9.6,' CHNG:',F9.6)
   99002 format (10X,10('!'),' CPA-cycle  NOT  converged after ',I3,       &
      ' iterations ',10('!'),/,14X,'ERROR ',F12.8,' CORRECTION ',F15.8,' CHANGE ',F15.8)
   99003 format (' CPA: ERROR increased by more than ','20*TOL for ITCPA=',&
      i4,' >>>> iteration stopped ',5X,10('!'))
   99004 format (' CPA: ERROR ',F12.8,'    CORRECTION ',F15.8,'    CHANGE ',F15.8)
   !
   99005 format (/,80('*')/,'ENERGY ',I5,2D16.8,' SITE ',I3)
   99006 format (2I5,1P,2D22.14)

end subroutine KLOOPZ1_QDOS
