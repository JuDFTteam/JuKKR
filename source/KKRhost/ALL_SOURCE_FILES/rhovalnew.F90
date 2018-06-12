!-------------------------------------------------------------------------------
! SUBROUTINE: RHOVALNEW
!> @note Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine RHOVALNEW(LMPOT,   &
   LDORHOEF,IELAST,NSRA,NSPIN,LMAX,EZ,WEZ,ZAT,SOCSCALE,CLEB,ICLEB,IEND,IFUNM,    &
   LMSP,NCHEB,NPAN_TOT,NPAN_LOG,NPAN_EQ,RMESH,IRWS,RPAN_INTERVALL,IPAN_INTERVALL,&
   RNEW,VINSNEW,THETASNEW,THETA,PHI,I1,IPOT,DEN_out,ESPV,RHO2NS,R2NEF,MUORB,     &
   angles_new,IDOLDAU,LOPT,PHILDAU,WLDAU,DENMATN,NATYP)

#ifdef CPP_OMP
   use omp_lib
#endif
#ifdef CPP_MPI
   use mpi
#endif
#ifdef CPP_TIMING
   use mod_timing
#endif

   use mod_types, only: t_tgmat, t_inc
#ifdef CPP_MPI
   use mod_types, only: gather_tmat,t_mpi_c_grid, save_t_mpi_c_grid, &
   get_ntot_pT_ioff_pT_2D
#endif
   use mod_mympi, only: myrank, master
#ifdef CPP_MPI
   use mod_mympi, only: find_dims_2d,distribute_linear_on_tasks,mympi_main1c_comm_newsosol
#endif
   use mod_save_wavefun, only: t_wavefunctions, read_wavefunc
   use mod_version_info
   use global_variables
   use Constants
   use Profiling
   Use mod_datatypes, Only: dp

   IMPLICIT NONE

   integer, intent(in) :: I1
   integer, intent(in) :: NSRA
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: IEND      !< Number of nonzero gaunt coefficients
   integer, intent(in) :: IPOT
   integer, intent(in) :: IRWS      !< R point at WS radius for a given atom
   integer, intent(in) :: LOPT      !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
   integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NSPIN     !< Counter for spin directions
   integer, intent(in) :: NCHEB     !< Number of Chebychev pannels for the new solver
   integer, intent(in) :: LMPOT     !< (LPOT+1)**2
   integer, intent(in) :: IELAST
   integer, intent(in) :: IDOLDAU   !< flag to perform LDA+U
   integer, intent(in) :: NPAN_EQ   !< Number of intervals from [R_LOG] to muffin-tin radius Used in conjunction with runopt NEWSOSOL
   integer, intent(in) :: NPAN_TOT
   integer, intent(in) :: NPAN_LOG  !< Number of intervals from nucleus to [R_LOG] Used in conjunction with runopt NEWSOSOL
   real (kind=dp), intent(in) :: ZAT       !< Nuclear charge for a given atom
   real (kind=dp), intent(in) :: SOCSCALE  !< Spin-orbit scaling for a given atom
   logical, intent(in) :: LDORHOEF
   integer, dimension(LMXSPD), intent(in)    :: LMSP !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
   integer, dimension(LMXSPD), intent(in)    :: IFUNM
   integer, dimension(0:NTOTD), intent(in)   :: IPAN_INTERVALL
   integer, dimension(NCLEB,4), intent(in)   :: ICLEB
   real (kind=dp), dimension(*), intent(in)   :: CLEB !< GAUNT coefficients (GAUNT)
   real (kind=dp), dimension(IRMD), intent(in) :: RMESH
   real (kind=dp), dimension(MMAXD,MMAXD,NSPIND), intent(in) :: WLDAU !< potential matrix
   complex (kind=dp), dimension(IRMD), intent(in) :: PHILDAU

   ! .. In/Out variables
   real (kind=dp), intent(inout) :: PHI
   real (kind=dp), intent(inout) :: THETA
   real (kind=dp), dimension(NRMAXD), intent(inout)       :: RNEW
   real (kind=dp), dimension(0:NTOTD), intent(inout)      :: RPAN_INTERVALL
   real (kind=dp), dimension(0:LMAX+1,3), intent(inout)   :: MUORB
   real (kind=dp), dimension(NRMAXD,NFUND), intent(inout) :: THETASNEW
   real (kind=dp), dimension(NRMAXD,LMPOT,NSPOTD), intent(inout) :: VINSNEW  !< Non-spherical part of the potential
   complex (kind=dp), dimension(IEMXD), intent(inout) :: EZ
   complex (kind=dp), dimension(IEMXD), intent(inout) :: WEZ
   ! .. Output variables
   real (kind=dp), dimension(2), intent(out)              :: angles_new
   real (kind=dp), dimension(0:LMAX+1,2), intent(out)     :: ESPV
   real (kind=dp), dimension(IRMD,LMPOT,4), intent(out)    :: R2NEF
   real (kind=dp), dimension(IRMD,LMPOT,4), intent(out)    :: RHO2NS
   complex (kind=dp), dimension(0:LMAX+1,IEMXD,2), intent(out) :: DEN_out

   ! .. Local variables
   integer :: LMAXD1
   integer :: IR,IREC,USE_SRATRICK,NVEC,LM1,LM2,IE,IRMDNEW,IMT1,JSPIN,IDIM,IORB, L1
   integer :: i_stat, i_all
   integer :: IQ, NQDOS             ! qdos ruess: number of qdos points
   integer :: LRECGFLLE,IERR        ! lmlm-dos
   integer :: LMLO,LMHI,IS,JS,MMAX  ! LDAU
   integer :: IX,M1       ! qdos ruess

   real (kind=dp) :: THETANEW,PHINEW
   real (kind=dp) :: TOTMOMENT
   real (kind=dp) :: TOTXYMOMENT
   complex (kind=dp) :: EK
   complex (kind=dp) :: DF
   complex (kind=dp) :: ERYD
   complex (kind=dp) :: TEMP1
   complex (kind=dp) :: DENTEMP
   complex (kind=dp) :: GMATPREFACTOR
   integer, dimension(4)         :: LMSHIFT1
   integer, dimension(4)         :: LMSHIFT2
   integer, dimension(2*LMMAXSO) :: JLK_INDEX
   real (kind=dp), dimension(3) :: MOMENT
   real (kind=dp), dimension(3) :: DENORBMOM
   real (kind=dp), dimension(3) :: DENORBMOMNS
   real (kind=dp), dimension(2,4) :: DENORBMOMSP
   real (kind=dp), dimension(0:LMAX,3) :: DENORBMOMLM
   complex (kind=dp), dimension(4)                       :: RHO2
   complex (kind=dp), dimension(4)                       :: RHO2INT
   complex (kind=dp), dimension(2*(LMAX+1))              :: ALPHASPH
   complex (kind=dp), dimension(LMMAXSO,LMMAXSO)         :: GMAT0
   complex (kind=dp), dimension(LMMAXSO,LMMAXSO)         :: GLDAU   ! LDAU
   complex (kind=dp), dimension(LMMAXSO,LMMAXSO)         :: TMATLL
   complex (kind=dp), dimension(LMMAXSO,LMMAXSO)         :: ALPHALL ! LLY
   complex (kind=dp), dimension(LMMAXSO,LMMAXSO)         :: TMATTEMP
   complex (kind=dp), dimension(2,2)                     :: RHO2NS_TEMP
   complex (kind=dp), dimension(LMMAXSO,LMMAXSO,IEMXD)   :: GMATLL
   complex (kind=dp), dimension(MMAXD,MMAXD,2,2)         :: DENMATN ! LDAU

   ! .. Local allocatable arrays
   real (kind=dp), dimension(:,:), allocatable :: QVEC ! qdos ruess: q-vectors for qdos
   real (kind=dp), dimension(:,:,:), allocatable :: VINS
   complex (kind=dp), dimension(:,:), allocatable :: TMATSPH
   complex (kind=dp), dimension(:,:), allocatable :: CDENTEMP
   complex (kind=dp), dimension(:,:), allocatable :: RHOTEMP
   complex (kind=dp), dimension(:,:), allocatable :: RHONEWTEMP
   complex (kind=dp), dimension(:,:,:), allocatable :: HLK
   complex (kind=dp), dimension(:,:,:), allocatable :: JLK
   complex (kind=dp), dimension(:,:,:), allocatable :: HLK2
   complex (kind=dp), dimension(:,:,:), allocatable :: JLK2
   complex (kind=dp), dimension(:,:,:), allocatable :: CDENNS
   complex (kind=dp), dimension(:,:,:), allocatable :: R2NEFC
   complex (kind=dp), dimension(:,:,:), allocatable :: RHO2NSC
   complex (kind=dp), dimension(:,:,:), allocatable :: VNSPLL0
   complex (kind=dp), dimension(:,:,:), allocatable :: R2NEFNEW
   complex (kind=dp), dimension(:,:,:), allocatable :: RHO2NSNEW
   complex (kind=dp), dimension(:,:,:), allocatable :: GFLLE_PART
   complex (kind=dp), dimension(:,:,:,:), allocatable :: RLL
   complex (kind=dp), dimension(:,:,:,:), allocatable :: SLL
   complex (kind=dp), dimension(:,:,:,:), allocatable :: DEN
   complex (kind=dp), dimension(:,:,:,:), allocatable :: CDEN
   complex (kind=dp), dimension(:,:,:,:), allocatable :: GFLLE
   complex (kind=dp), dimension(:,:,:,:), allocatable :: DENLM
   complex (kind=dp), dimension(:,:,:,:), allocatable :: CDENLM
   complex (kind=dp), dimension(:,:,:,:), allocatable :: VNSPLL
   complex (kind=dp), dimension(:,:,:,:), allocatable :: R2ORBC
   complex (kind=dp), dimension(:,:,:,:), allocatable :: VNSPLL1
   complex (kind=dp), dimension(:,:,:,:), allocatable :: RLLLEFT
   complex (kind=dp), dimension(:,:,:,:), allocatable :: SLLLEFT
   complex (kind=dp), dimension(:,:,:,:), allocatable :: R2NEFC_loop
   complex (kind=dp), dimension(:,:,:,:), allocatable :: RHO2NSC_loop

#ifdef CPP_MPI
   complex (kind=dp), dimension(2) :: DENTOT         ! qdos ruess
   ! communication
   complex (kind=dp), dimension(:,:,:,:), allocatable :: workc
#endif
   ! OMP - number of threads, thread id
   integer nth,ith
   integer :: ie_start,ie_end, ie_num
   integer :: i1_myrank ! lmlm-dos, needed for MPI with more than one rank per energy point (nranks_ie>1)
   ! read in wavefunctions
   logical :: rll_was_read_in, sll_was_read_in,rllleft_was_read_in, sllleft_was_read_in
   !      ..
   logical :: TEST,OPT
   external :: TEST,OPT

   ! determine if omp is used
   ith = 0
   nth = 1
#ifdef CPP_OMP
   !$omp parallel shared(nth,ith)
   !$omp single
   nth = omp_get_num_threads()
   if(t_inc%i_write>0) write(1337,*) 'nth =',nth
   !$omp end single
   !$omp end parallel
#endif

   ! .. Parameters
   LMAXD1= LMAX+1

   IRMDNEW= NPAN_TOT*(NCHEB+1)
   IMT1=IPAN_INTERVALL(NPAN_LOG+NPAN_EQ)+1
   allocate(VINS(IRMDNEW,LMPOT,NSPIN),stat=i_stat)
   call memocc(i_stat,product(shape(VINS))*kind(VINS),'VINS','RHOVALNEW')
   VINS=0d0
   do LM1=1,LMPOT
      do IR=1,IRMDNEW
         VINS(IR,LM1,1)=VINSNEW(IR,LM1,IPOT)
         VINS(IR,LM1,NSPIN)=VINSNEW(IR,LM1,IPOT+NSPIN-1)
      enddo
   enddo

   !! set up the non-spherical ll' matrix for potential VLL'
   if (NSRA.EQ.2) then
      USE_SRATRICK=1
   else
      USE_SRATRICK=0
   endif
   allocate(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW),stat=i_stat)
   call memocc(i_stat,product(shape(VNSPLL0))*kind(VNSPLL0),'VNSPLL0','RHOVALNEW')
   VNSPLL0=CZERO
   allocate(VNSPLL1(LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(VNSPLL1))*kind(VNSPLL1),'VNSPLL1','RHOVALNEW')
   VNSPLL0=CZERO
   !
   call VLLMAT(1,NRMAXD,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,VINS,LMPOT,CLEB,ICLEB,IEND,&
      NSPIN,ZAT,RNEW,USE_SRATRICK,NCLEB)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LDAU
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (IDOLDAU.EQ.1) then
      LMLO=LOPT**2+1
      LMHI=(LOPT+1)**2
      do IR=1,IRMDNEW
         VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)=VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)+WLDAU(1:MMAXD,1:MMAXD,1)
      enddo
      LMLO=LMLO+LMMAXD
      LMHI=LMHI+LMMAXD
      do IR=1,IRMDNEW
         VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)=VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)+WLDAU(1:MMAXD,1:MMAXD,2)
      enddo
   endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! LDAU
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! initial allocate
   if (NSRA.EQ.2) then
      allocate(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
      call memocc(i_stat,product(shape(VNSPLL))*kind(VNSPLL),'VNSPLL','RHOVALNEW')
   else
      allocate(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
      call memocc(i_stat,product(shape(VNSPLL))*kind(VNSPLL),'VNSPLL','RHOVALNEW')
   endif

   allocate(HLK(4*(LMAX+1),IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(HLK))*kind(HLK),'HLK','RHOVALNEW')
   allocate(JLK(4*(LMAX+1),IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(JLK))*kind(JLK),'JLK','RHOVALNEW')
   allocate(HLK2(4*(LMAX+1),IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(HLK2))*kind(HLK2),'HLK2','RHOVALNEW')
   allocate(JLK2(4*(LMAX+1),IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(JLK2))*kind(JLK2),'JLK2','RHOVALNEW')
   allocate(TMATSPH(2*(LMAX+1),0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(TMATSPH))*kind(TMATSPH),'TMATSPH','RHOVALNEW')
   allocate(RLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(RLL))*kind(RLL),'RLL','RHOVALNEW')
   allocate(SLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(SLL))*kind(SLL),'SLL','RHOVALNEW')
   allocate(RLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(RLLLEFT))*kind(RLLLEFT),'RLLLEFT','RHOVALNEW')
   allocate(SLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(SLLLEFT))*kind(SLLLEFT),'SLLLEFT','RHOVALNEW')
   allocate(CDEN(IRMDNEW,0:LMAX,4,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(CDEN))*kind(CDEN),'CDEN','RHOVALNEW')
   allocate(CDENLM(IRMDNEW,LMMAXD,4,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(CDENLM))*kind(CDENLM),'CDENLM','RHOVALNEW')
   allocate(CDENNS(IRMDNEW,4,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(CDENNS))*kind(CDENNS),'CDENNS','RHOVALNEW')
   allocate(RHO2NSC(IRMDNEW,LMPOT,4),stat=i_stat)
   call memocc(i_stat,product(shape(RHO2NSC))*kind(RHO2NSC),'RHO2NSC','RHOVALNEW')
   allocate(RHO2NSC_loop(IRMDNEW,LMPOT,4,ielast),stat=i_stat)
   call memocc(i_stat,product(shape(RHO2NSC_loop))*kind(RHO2NSC_loop),'RHO2NSC_loop','RHOVALNEW')
   allocate(RHO2NSNEW(IRMD,LMPOT,4),stat=i_stat)
   call memocc(i_stat,product(shape(RHO2NSNEW))*kind(RHO2NSNEW),'RHO2NSNEW','RHOVALNEW')
   allocate(R2NEFC(IRMDNEW,LMPOT,4),stat=i_stat)
   call memocc(i_stat,product(shape(R2NEFC))*kind(R2NEFC),'R2NEFC','RHOVALNEW')
   allocate(R2NEFC_loop(IRMDNEW,LMPOT,4,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(R2NEFC_loop))*kind(R2NEFC_loop),'R2NEFC_loop','RHOVALNEW')
   allocate(R2NEFNEW(IRMD,LMPOT,4),stat=i_stat)
   call memocc(i_stat,product(shape(R2NEFNEW))*kind(R2NEFNEW),'R2NEFNEW','RHOVALNEW')
   allocate(R2ORBC(IRMDNEW,LMPOT,4,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(R2ORBC))*kind(R2ORBC),'R2ORBC','RHOVALNEW')
   allocate(CDENTEMP(IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(CDENTEMP))*kind(CDENTEMP),'CDENTEMP','RHOVALNEW')
   allocate(GFLLE_PART(LMMAXSO,LMMAXSO,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(GFLLE_PART))*kind(GFLLE_PART),'GFLLE_PART','RHOVALNEW')
   allocate(GFLLE(LMMAXSO,LMMAXSO,IELAST,1),stat=i_stat)
   call memocc(i_stat,product(shape(GFLLE))*kind(GFLLE),'GFLLE','RHOVALNEW')
   allocate(DEN(0:LMAXD1,IEMXD,1,2),DENLM(LMMAXD,IEMXD,1,2),stat=i_stat)
   call memocc(i_stat,product(shape(DEN))*kind(DEN),'DEN','RHOVALNEW')
   RHO2NSC=CZERO
   RHO2NSC_loop=CZERO
   R2NEFC=CZERO
   R2NEFC_loop=CZERO
   R2ORBC=CZERO
   RHO2NS=0.D0  ! fivos 19.7.2014, this was CZERO
   R2NEF=0.D0   ! fivos 19.7.2014, this was CZERO
   RHO2NSNEW=CZERO
   R2NEFNEW=CZERO
   DEN=CZERO
   DENLM=CZERO
   ESPV=0d0
   RHO2INT=CZERO
   DENORBMOM=0d0
   DENORBMOMSP=0d0
   DENORBMOMLM=0d0
   DENORBMOMNS=0d0
   THETANEW=0d0
   PHINEW=0d0
   GFLLE_PART=CZERO
   GFLLE=CZERO
   GLDAU=CZERO
   ! LM shifts for correct density summation
   LMSHIFT1(1)=0                                                   ! qdos ruess
   LMSHIFT1(2)=LMMAXD                                              ! qdos ruess
   LMSHIFT1(3)=0                                                   ! qdos ruess
   LMSHIFT1(4)=LMMAXD                                              ! qdos ruess
   LMSHIFT2(1)=0                                                   ! qdos ruess
   LMSHIFT2(2)=LMMAXD                                              ! qdos ruess
   LMSHIFT2(3)=LMMAXD                                              ! qdos ruess
   LMSHIFT2(4)=0                                                   ! qdos ruess

   !      DO IR=1,3
   !       DO LM1=0,LMAXD1+1
   !        MUORB(LM1,IR)=0d0  !zimmer: initialization shifted to main1c
   !       ENDDO
   !      ENDDO

   NQDOS = 1                                                                     ! qdos ruess
   if (OPT('qdos    ')) then                                                     ! qdos ruess
      !        Read BZ path for qdos calculation:                                ! qdos ruess
      open(67,FILE='qvec.dat',STATUS='old',IOSTAT=IERR,ERR=3000)                 ! qdos ruess
      read(67,*) NQDOS                                                           ! qdos ruess
      allocate(QVEC(3,NQDOS),stat=i_stat)                                        ! qdos ruess
      call memocc(i_stat,product(shape(QVEC))*kind(QVEC),'QVEC','RHOVALNEW')     ! qdos ruess
      do IQ = 1,NQDOS                                                            ! qdos ruess
         read(67,*) (QVEC(IX,IQ),IX=1,3)                                         ! qdos ruess
      enddo                                                                      ! qdos ruess
      close(67)                                                                  ! qdos ruess
      !        Change allocation for GFLLE to be suitabel for qdos run           ! qdos ruess
      i_all=-product(shape(GFLLE))*kind(GFLLE)                                   ! qdos ruess
      deallocate(GFLLE, stat=i_stat)                                             ! qdos ruess
      call memocc(i_stat,i_all,'GFLLE','RHOVALNEW')                              ! qdos ruess
      i_all=-product(shape(DEN))*kind(DEN)                                       ! qdos ruess
      deallocate(DEN, stat=i_stat)                                               ! qdos ruess
      call memocc(i_stat,i_all,'DEN','RHOVALNEW')                                ! qdos ruess
      i_all=-product(shape(DENLM))*kind(DENLM)                                   ! qdos ruess
      deallocate(DENLM, stat=i_stat)                                             ! qdos ruess
      call memocc(i_stat,i_all,'DENLM','RHOVALNEW')                              ! qdos ruess
      !                                                                          ! qdos ruess
      allocate(GFLLE(LMMAXSO,LMMAXSO,IELAST,NQDOS),stat=i_stat)                  ! qdos ruess
      call memocc(i_stat,product(shape(GFLLE))*kind(GFLLE),'GFLLE','RHOVALNEW')  ! qdos ruess
      allocate(DEN(0:LMAXD1,IEMXD,NQDOS,2),stat=i_stat)                          ! qdos ruess
      call memocc(i_stat,product(shape(DEN))*kind(DEN),'DEN','RHOVALNEW')        ! qdos ruess
      allocate(DENLM(LMMAXD,IEMXD,NQDOS,2),stat=i_stat)                          ! qdos ruess
      call memocc(i_stat,product(shape(DENLM))*kind(QVEC),'DENLM','RHOVALNEW')   ! qdos ruess
      3000  if (IERR.NE.0) stop 'ERROR READING ''qvec.dat'''                     ! qdos ruess
   end if  ! OPT('qdos    ')                                                     ! qdos ruess

#ifdef CPP_MPI
   i1_myrank = i1 - t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie)    ! lmlm-dos ruess
#else
   i1_myrank = i1                                                    ! lmlm-dos ruess
#endif
   if ((OPT('lmlm-dos')).AND.(I1_myrank.EQ.1)) then                  ! lmlm-dos ruess
      LRECGFLLE = 4*LMMAXSO*LMMAXSO*IELAST*NQDOS                     ! lmlm-dos ruess
      open(91,ACCESS='direct',RECL=LRECGFLLE,FILE='gflle',&          ! lmlm-dos ruess
         FORM='unformatted',STATUS='replace',ERR=3001,IOSTAT=IERR)   ! lmlm-dos ruess
      3001 if (IERR.NE.0) stop 'ERROR CREATING ''gflle'''            ! lmlm-dos ruess
   endif                                                             ! lmlm-dos ruess

   ! initialize to zero
   DEN=CZERO
   DENLM=CZERO
   ! energy loop
   if(myrank==master.and.t_inc%i_write>0) write(1337,*) 'atom: ',I1
#ifdef CPP_MPI
   ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
   ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
#else
   ie_start = 0 ! offset
   ie_end   = IELAST
#endif

#ifdef CPP_OMP
   ! omp: start parallel region here
   !$omp parallel do default(none)                                         &
   !$omp private(eryd,ie,ir,irec,lm1,lm2,gmatprefactor,nvec)               &
   !$omp private(jlk_index,tmatll,ith)                                     &
   !$omp private(iq,df,ek,tmattemp,gmatll,gmat0,iorb,dentemp)              &
   !$omp private(rho2ns_temp,rho2,temp1,jspin)                             &
   !$omp private(alphasph,alphall,ie_num)                                  &
   !$omp private(rll_was_read_in, sll_was_read_in)                         &
   !$omp private(rllleft_was_read_in, sllleft_was_read_in)                 &
   !!$omp firstprivate(t_inc)                                              &
   !$omp shared(t_inc)                                                     &
   !$omp shared(ldorhoef,nqdos,lmshift1,lmshift2,wez,lmsp,imt1,ifunm)      &
   !$omp shared(r2orbc,r2nefc,cden,cdenlm,cdenns,rho2nsc_loop)             &
   !$omp shared(nspin,nsra,iend,ipot,ielast,npan_tot,ncheb,lmax)           &
   !$omp shared(zat,socscale,ez,rmesh,cleb,rnew,nth,icleb,thetasnew,i1)    &
   !$omp shared(rpan_intervall,vinsnew,ipan_intervall,r2nefc_loop)         &
   !$omp shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0)               &
   !$omp shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,cdentemp)         &
   !$omp shared(tmatsph,den,denlm,gflle,gflle_part,rllleft,sllleft)        &
   !$omp shared(t_tgmat,ie_end, ie_start, t_wavefunctions)                 &
   !$omp shared(LMMAXSO,LMMAXD,LMPOT,NRMAXD,NTOTD,LMAXD1)                  &
   !$omp reduction(+:rho2int,espv) reduction(-:muorb)                      &
   !$omp reduction(-:denorbmom,denorbmomsp,denorbmomlm,denorbmomns)
#endif
   do ie_num=1,ie_end
      IE = ie_start+ie_num

#ifdef CPP_OMP
      ith = omp_get_thread_num()
#else
      ith = 0
#endif

      ERYD=EZ(IE)
      EK=SQRT(ERYD)
      DF=WEZ(IE)/DBLE(NSPIN)
      if (NSRA.EQ.2) then
         EK = SQRT(ERYD+ERYD*ERYD/(CVLIGHT*CVLIGHT))*(1d0+ERYD/(CVLIGHT*CVLIGHT))
      endif
#ifdef CPP_OMP
      !$omp critical
#endif
      if(t_inc%i_write>0) write(1337,*) 'energy:',IE,'',ERYD
#ifdef CPP_OMP
      !$omp end critical
#endif

      if(t_wavefunctions%Nwfsavemax>0) then ! read wavefunctions?
         ! read in wavefunction from memory
         call read_wavefunc(t_wavefunctions,rll, rllleft, sll, sllleft, &
            i1, ie, NSRA, LMMAXSO, IRMDNEW, ith, nth,                   &
            rll_was_read_in, sll_was_read_in,                           &
            rllleft_was_read_in, sllleft_was_read_in)
      end if

      ! recalculate wavefuntions, also include left solution
      ! contruct the spin-orbit coupling hamiltonian and add to potential
      call SPINORBIT_HAM(LMAX,LMMAXD,VINS,RNEW,    &
         ERYD,ZAT,CVLIGHT,SOCSCALE,NSPIN,LMPOT,    &
         THETA,PHI,IPAN_INTERVALL,RPAN_INTERVALL,  &
         NPAN_TOT,NCHEB,IRMDNEW,NRMAXD,            &
         VNSPLL0,VNSPLL1(:,:,:,ith),'1')

      ! extend matrix for the SRA treatment
      VNSPLL(:,:,:,ith)=CZERO
      if (NSRA.EQ.2) then
         if (USE_SRATRICK.EQ.0) then
            call VLLMATSRA(VNSPLL1(:,:,:,ith),VNSPLL(:,:,:,ith),RNEW,LMMAXSO,&
               IRMDNEW,NRMAXD,ERYD,LMAX,0,'Ref=0')
         elseif (USE_SRATRICK.EQ.1) then
            call VLLMATSRA(VNSPLL1(:,:,:,ith),VNSPLL(:,:,:,ith),RNEW,LMMAXSO,&
               IRMDNEW,NRMAXD,ERYD,LMAX,0,'Ref=Vsph')
         endif
      else
         VNSPLL(:,:,:,ith)=VNSPLL1(:,:,:,ith)
      endif

      if( (t_wavefunctions%Nwfsavemax>0 .and..not. rll_was_read_in)  &
         .or. (t_wavefunctions%Nwfsavemax==0)) then ! read/recalc wavefunctions

         ! calculate the source terms in the Lippmann-Schwinger equation
         ! these are spherical hankel and bessel functions
         HLK(:,:,ith)=CZERO
         JLK(:,:,ith)=CZERO
         HLK2(:,:,ith)=CZERO
         JLK2(:,:,ith)=CZERO
         GMATPREFACTOR=CZERO
         JLK_INDEX=0
         call RLLSLLSOURCETERMS(NSRA,NVEC,ERYD,RNEW,IRMDNEW,NRMAXD,LMAX,&
            LMMAXSO,1,JLK_INDEX,HLK(:,:,ith),                           &
            JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),                   &
            GMATPREFACTOR)

         ! using spherical potential as reference
         if (USE_SRATRICK.EQ.1) then
            call CALCSPH(NSRA,IRMDNEW,NRMAXD,LMAX,NSPIN,ZAT,CVLIGHT,ERYD,  &
               LMPOT,LMMAXSO,RNEW,VINS,NCHEB,NPAN_TOT,RPAN_INTERVALL,      &
               JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),HLK2(:,:,ith),          &
               JLK2(:,:,ith),GMATPREFACTOR,TMATSPH(:,ith),                 &
               ALPHASPH,USE_SRATRICK)
         endif

         ! calculate the tmat and wavefunctions
         RLL(:,:,:,ith)=CZERO
         SLL(:,:,:,ith)=CZERO

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Right solutions
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         TMATLL=CZERO
         ! faster calculation of RLL.
         ! no irregular solutions SLL are needed in self-consistent iterations
         ! because the density depends only on RLL, RLLLEFT and SLLLEFT
         if(OPT('RLL-SLL ') .and.  .not. ( OPT('XCPL    ').or.OPT('OPERATOR') ) ) then
            call rll_global_solutions(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),  &
               RLL(:,:,:,ith),TMATLL,                                         &
               NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),                &
               IRMDNEW,NRMAXD,NSRA,JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),       &
               HLK2(:,:,ith),JLK2(:,:,ith),                                   &
               GMATPREFACTOR,'1',USE_SRATRICK,ALPHALL)
         else
            call RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),          &
               RLL(:,:,:,ith),SLL(:,:,:,ith),TMATLL,                    &
               NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),          &
               IRMDNEW,NRMAXD,NSRA,JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith), &
               HLK2(:,:,ith),JLK2(:,:,ith),                             &
               GMATPREFACTOR,'1','1','0',USE_SRATRICK,ALPHALL)
         endif
         if (NSRA.EQ.2) then
            RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
            SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
         endif

      end if ! read/recalc wavefunctions

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Left solutions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if( (t_wavefunctions%Nwfsavemax>0 .and.                  &
         .not. (rllleft_was_read_in.and.sllleft_was_read_in) ) &
         .or. (t_wavefunctions%Nwfsavemax==0)) then
         ! read/recalc wavefunctions left contruct the TRANSPOSE spin-orbit coupling hamiltonian and add to potential
         call SPINORBIT_HAM(LMAX,LMMAXD,VINS,RNEW,ERYD,ZAT, &
            CVLIGHT,SOCSCALE,NSPIN,LMPOT,THETA,PHI,         &
            IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,NCHEB,   &
            IRMDNEW,NRMAXD,VNSPLL0,VNSPLL1(:,:,:,ith),      &
            'transpose')
         ! extend matrix for the SRA treatment
         VNSPLL(:,:,:,ith)=CZERO
         if (NSRA.EQ.2) then
            if (USE_SRATRICK.EQ.0) then
               call VLLMATSRA(VNSPLL1(:,:,:,ith),VNSPLL(:,:,:,ith),RNEW,LMMAXSO,&
                  IRMDNEW,NRMAXD,ERYD,LMAX,0,'Ref=0')
            elseif (USE_SRATRICK.EQ.1) then
               call VLLMATSRA(VNSPLL1(:,:,:,ith),VNSPLL(:,:,:,ith),RNEW,LMMAXSO,&
                  IRMDNEW,NRMAXD,ERYD,LMAX,0,'Ref=Vsph')
            endif
         else
            VNSPLL(:,:,:,ith)=VNSPLL1(:,:,:,ith)
         endif

         ! calculate the source terms in the Lippmann-Schwinger equation
         ! these are spherical hankel and bessel functions
         HLK(:,:,ith)=CZERO
         JLK(:,:,ith)=CZERO
         HLK2(:,:,ith)=CZERO
         JLK2(:,:,ith)=CZERO
         GMATPREFACTOR=CZERO
         JLK_INDEX=0
         call RLLSLLSOURCETERMS(NSRA,NVEC,ERYD,RNEW,IRMDNEW,NRMAXD,LMAX,&
            LMMAXSO,1,JLK_INDEX,HLK(:,:,ith),                           &
            JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),                   &
            GMATPREFACTOR)

         ! using spherical potential as reference
         ! notice that exchange the order of left and right hankel/bessel functions
         if (USE_SRATRICK.EQ.1) then
            call CALCSPH(NSRA,IRMDNEW,NRMAXD,LMAX,NSPIN,ZAT,CVLIGHT,ERYD,  &
               LMPOT,LMMAXSO,RNEW,VINS,NCHEB,NPAN_TOT,RPAN_INTERVALL,      &
               JLK_INDEX,HLK2(:,:,ith),JLK2(:,:,ith),                      &
               HLK(:,:,ith),JLK(:,:,ith),GMATPREFACTOR,                    &
               ALPHASPH,TMATSPH(:,ith),USE_SRATRICK)
         endif

         ! calculate the tmat and wavefunctions
         RLLLEFT(:,:,:,ith)=CZERO
         SLLLEFT(:,:,:,ith)=CZERO

         ! left solutions
         ! notice that exchange the order of left and right hankel/bessel functions
         TMATTEMP=CZERO
         ! faster calculation of RLLLEFT and SLLLEFT.
         if(OPT('RLL-SLL ') .and.  .not. ( OPT('XCPL    ').or.OPT('OPERATOR') ) ) then
            call rll_global_solutions(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),  &
               RLLLEFT(:,:,:,ith),TMATTEMP,                                   &
               NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),                &
               IRMDNEW,NRMAXD,NSRA,JLK_INDEX,HLK2(:,:,ith),JLK2(:,:,ith),     &
               HLK(:,:,ith),JLK(:,:,ith),                                     &
               GMATPREFACTOR,'1',USE_SRATRICK,ALPHALL)
            call sll_global_solutions(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),  &
               SLLLEFT(:,:,:,ith),                                            &
               NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),                &
               IRMDNEW,NRMAXD,NSRA,JLK_INDEX,HLK2(:,:,ith),JLK2(:,:,ith),     &
               HLK(:,:,ith),JLK(:,:,ith),                                     &
               GMATPREFACTOR,'1',USE_SRATRICK,ALPHALL)
         else
            call RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),             &
               RLLLEFT(:,:,:,ith),SLLLEFT(:,:,:,ith),TMATTEMP,             &
               NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),             &
               IRMDNEW,NRMAXD,NSRA,JLK_INDEX,HLK2(:,:,ith),JLK2(:,:,ith),  &
               HLK(:,:,ith),JLK(:,:,ith),                                  &
               GMATPREFACTOR,'1','1','0',USE_SRATRICK,ALPHALL)
         end if
         if (NSRA.EQ.2) then
            RLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=RLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
            SLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=SLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
         endif
      end if ! read/recalc wavefunctions left

      do IQ = 1,NQDOS                                                ! qdos
         ! read in GF
         IREC = IQ + NQDOS * (IE-1) +  NQDOS * IELAST * (I1-1)       ! qdos
#ifdef CPP_OMP
         !$omp critical
#endif
         if (t_tgmat%gmat_to_file) then
            read(69,REC=IREC) GMAT0
         else
            IREC = IQ + NQDOS * (ie_num-1) + NQDOS *ie_end * (I1-1)
            GMAT0(:,:) = t_tgmat%gmat(:,:,irec)
         end if
#ifdef CPP_OMP
         !$omp end critical
#endif

         ! rotate gmat from global frame to local frame
         call ROTATEMATRIX(GMAT0,THETA,PHI,LMMAXD,1)

         do LM1=1,LMMAXSO
            do LM2=1,LMMAXSO
               GMATLL(LM1,LM2,IE)=GMAT0(LM1,LM2)
            enddo
         enddo
         ! calculate density
         call RHOOUTNEW(NSRA,LMAX,GMATLL(1,1,IE),EK,  &
            LMPOT,DF,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,                 &
            IRMDNEW,THETASNEW,IFUNM,IMT1,LMSP,                       &
            RLL(:,:,:,ith),                                          & !SLL(:,:,:,ith), commented out since sll is not used in rhooutnew
            RLLLEFT(:,:,:,ith),SLLLEFT(:,:,:,ith),                   &
            CDEN(:,:,:,ith),CDENLM(:,:,:,ith),                       &
            CDENNS(:,:,ith),RHO2NSC_loop(:,:,:,ie),0,                &
            GFLLE(:,:,IE,IQ),RPAN_INTERVALL,IPAN_INTERVALL)

         do JSPIN=1,4
            do LM1 = 0,LMAX
               CDENTEMP(:,ith)=CZERO
               DENTEMP=CZERO
               do IR=1,IRMDNEW
                  CDENTEMP(IR,ith)=CDEN(IR,LM1,JSPIN,ith)
               enddo
               call INTCHEB_CELL(CDENTEMP(:,ith),DENTEMP,RPAN_INTERVALL,&
                  IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
               RHO2(JSPIN)=DENTEMP
               RHO2INT(JSPIN)=RHO2INT(JSPIN)+RHO2(JSPIN)*DF
               if (JSPIN.LE.2) then
                  DEN(LM1,IE,IQ,JSPIN)=RHO2(JSPIN)
               endif
            enddo

            if (JSPIN.LE.2) then
               do LM1 = 1,LMMAXD
                  CDENTEMP(:,ith)=CZERO
                  DENTEMP=CZERO
                  do IR=1,IRMDNEW
                     CDENTEMP(IR,ith)=CDENLM(IR,LM1,JSPIN,ith)
                  enddo
                  call INTCHEB_CELL(CDENTEMP(:,ith),DENTEMP,RPAN_INTERVALL,&
                     IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
                  DENLM(LM1,IE,IQ,JSPIN)=DENTEMP
               enddo
               CDENTEMP(:,ith)=CZERO
               DENTEMP=CZERO
               do IR=1,IRMDNEW
                  CDENTEMP(IR,ith)=CDENNS(IR,JSPIN,ith)
               enddo
               call INTCHEB_CELL(CDENTEMP(:,ith),DENTEMP,RPAN_INTERVALL,&
                  IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
               DEN(LMAXD1,IE,IQ,JSPIN)=DENTEMP
               RHO2INT(JSPIN)=RHO2INT(JSPIN)+DEN(LMAXD1,IE,IQ,JSPIN)*DF
            endif
         enddo ! JSPIN

         do JSPIN=1,4
            if (JSPIN.LE.2) then
               do LM1=0,LMAXD1
                  ESPV(LM1,JSPIN)=ESPV(LM1,JSPIN)+aimag( ERYD * DEN(LM1,IE,IQ,JSPIN) * DF )
               enddo
            endif
         enddo
      end do   ! IQ = 1,NQDOS

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Get charge at the Fermi energy (IELAST)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (IE.EQ.IELAST.AND.LDORHOEF) then
         call RHOOUTNEW(NSRA,LMAX,GMATLL(1,1,IE),EK,  &
            LMPOT,CONE,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,               &
            IRMDNEW,THETASNEW,IFUNM,IMT1,LMSP,                       &
            RLL(:,:,:,ith),                                          & !SLL(:,:,:,ith), ! commented out since sll is not used in rhooutnew
            RLLLEFT(:,:,:,ith),SLLLEFT(:,:,:,ith),                   &
            CDEN(:,:,:,ith),CDENLM(:,:,:,ith),                       &
            CDENNS(:,:,ith),R2NEFC_loop(:,:,:,ith),0,                &
            GFLLE_PART(:,:,ith),RPAN_INTERVALL,IPAN_INTERVALL)
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Get orbital moment
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do IORB=1,3
         call RHOOUTNEW(NSRA,LMAX,GMATLL(1,1,IE),EK,  &
            LMPOT,CONE,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,               &
            IRMDNEW,THETASNEW,IFUNM,IMT1,LMSP,                       &
            RLL(:,:,:,ith),                                          & !SLL(:,:,:,ith), ! commented out since sll is not used in rhooutnew
            RLLLEFT(:,:,:,ith),SLLLEFT(:,:,:,ith),                   &
            CDEN(:,:,:,ith),CDENLM(:,:,:,ith),                       &
            CDENNS(:,:,ith),R2ORBC(:,:,:,ith),IORB,                  &
            GFLLE_PART(:,:,ith),RPAN_INTERVALL,IPAN_INTERVALL)
         do JSPIN=1,4
            if (JSPIN.LE.2) then
               do LM1=0,LMAX
                  CDENTEMP(:,ith)=CZERO
                  DENTEMP=CZERO
                  do IR=1,IRMDNEW
                     CDENTEMP(IR,ith)=CDEN(IR,LM1,JSPIN,ith)
                  enddo
                  call INTCHEB_CELL(CDENTEMP(:,ith),DENTEMP,RPAN_INTERVALL,&
                     IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
                  !
                  RHO2(JSPIN)=DENTEMP
                  MUORB(LM1,JSPIN)=MUORB(LM1,JSPIN)-aimag(RHO2(JSPIN)*DF)
                  DENORBMOM(IORB)=DENORBMOM(IORB)-aimag(RHO2(JSPIN)*DF)
                  DENORBMOMSP(JSPIN,IORB)=DENORBMOMSP(JSPIN,IORB)-aimag(RHO2(JSPIN)*DF)
                  DENORBMOMLM(LM1,IORB)=DENORBMOMLM(LM1,IORB)-aimag(RHO2(JSPIN)*DF)
                  CDENTEMP(:,ith)=CZERO
                  !
                  do IR=1,IRMDNEW
                     CDENTEMP(IR,ith)=CDENNS(IR,JSPIN,ith)
                  enddo
                  call INTCHEB_CELL(CDENTEMP(:,ith),TEMP1,RPAN_INTERVALL,&
                     IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
                  DENORBMOMNS(IORB)=DENORBMOMNS(IORB)-aimag(TEMP1*DF)
               enddo
            endif
         enddo
      enddo ! IORB
   enddo ! IE loop
#ifdef CPP_OMP
   !$omp end parallel do
#endif

   ! omp: move sum from rhooutnew here after parallel calculation
   do IR=1,IRMDNEW
      do LM1=1,LMPOT
         do JSPIN=1,4
            do IE=1,IELAST
               RHO2NSC(IR,LM1,JSPIN) = RHO2NSC(IR,LM1,JSPIN) +RHO2NSC_loop(IR,LM1,JSPIN,IE)
            enddo
         enddo
      enddo
   enddo
   ! omp: don't forget to do the same with density at fermi energy:
   do ith=0,nth-1
      R2NEFC(:,:,:) = R2NEFC(:,:,:) + R2NEFC_loop(:,:,:,ith)
   enddo

#ifdef CPP_MPI
   if (OPT('qdos    ')) then                                                     ! qdos
      ! first communicate den array to write out qdos files                      ! qdos
      IDIM = (LMAXD1+1)*IEMXD*2*NQDOS                                            ! qdos
      allocate(workc(0:LMAXD1,IEMXD,2,NQDOS),stat=i_stat)                        ! qdos
      call memocc(i_stat,product(shape(workc))*kind(workc),'workc','RHOVALNEW')  ! qdos
      workc = CZERO                                                              ! qdos
      call MPI_REDUCE(DEN, workc, IDIM, MPI_DOUBLE_COMPLEX, MPI_SUM,master, &    ! qdos
         t_mpi_c_grid%myMPI_comm_at, IERR)                                       ! qdos
      call ZCOPY(IDIM,WORKC,1,DEN,1)                                             ! qdos
      i_all=-product(shape(workc))*kind(workc)                                   ! qdos
      deallocate(workc, stat=i_stat)                                             ! qdos
      call memocc(i_stat,i_all,'workc','RHOVALNEW')                              ! qdos
      !                                                                          ! qdos
      if(t_mpi_c_grid%myrank_at==master) then                                    ! qdos
         ie_start = 0                                                            ! qdos
         ie_end   = IELAST                                                       ! qdos
         do ie_num=1,ie_end                                                      ! qdos
            IE = ie_start+ie_num                                                 ! qdos
            do IQ=1,NQDOS                                                        ! qdos
               if ((IQ.EQ.1).AND.(IE_num.EQ.1)) then                             ! qdos
                  if(t_inc%NATYP.ge.100) then                                    ! qdos
                     open(31,                                                  & ! qdos
                        FILE="qdos."//char(48+I1/100)//char(48+mod(I1/10,10))//& ! qdos
                        char(48+mod(I1,10))//"."//char(48+1)//".dat")            ! qdos
                     open(32,                                                  & ! qdos
                        FILE="qdos."//char(48+I1/100)//char(48+mod(I1/10,10))//& ! qdos
                        char(48+mod(I1,10))//"."//char(48+2)//".dat")            ! qdos
                  else                                                           ! qdos
                     open(31,                                                  & ! qdos
                        FILE="qdos."//char(48+mod(I1/10,10))//                 & ! qdos
                        char(48+mod(I1,10))//"."//char(48+1)//".dat")            ! qdos
                     open(32,                                                  & ! qdos
                        FILE="qdos."//char(48+mod(I1/10,10))//                 & ! qdos
                        char(48+mod(I1,10))//"."//char(48+2)//".dat")            ! qdos
                  end if                                                         ! qdos
                  call version_print_header(31)                                  ! qdos
                  write (31,*) ' '                                               ! qdos
                  write (31,8600) '# ISPIN=',1,' I1=',I1                         ! qdos
                  write (31,'(7(A,3X))') '#   Re(E)','Im(E)','k_x','k_y','k_z'&  ! qdos
                     ,'DEN_tot','DEN_s,p,...'                                    ! qdos
                  if(NSPIN.gt.1) then                                            ! qdos
                     call version_print_header(32)                               ! qdos
                     write (32,*) ' '                                            ! qdos
                     write (32,8600) '# ISPIN=',2,' I1=',I1                      ! qdos
                     write (32,'(7(A,3X))') '#   Re(E)','Im(E)',           &     ! qdos
                     'k_x','k_y','k_z','DEN_tot','DEN_s,p,...'                   ! qdos
                  end if                                                         ! qdos
               endif   ! IQ.EQ.1                                                 ! qdos
               do JSPIN =1,2                                                     ! qdos
                  DENTOT(JSPIN) = DCMPLX(0.D0,0.D0, kind=dp)                              ! qdos
                  do L1 = 0,LMAXD1                                               ! qdos
                     DENTOT(JSPIN) = DENTOT(JSPIN) + DEN(L1,IE,IQ,JSPIN)         ! qdos
                  enddo                                                          ! qdos
               enddo                                                             ! qdos
               !    write qdos.nn.s.dat                                          ! qdos
               write(31,9000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),  &        ! qdos
                  -aimag(DENTOT(1))/PI,(-aimag(DEN(L1,IE,IQ,1))/PI,L1=0,LMAXD1)  ! qdos
               write(32,9000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),  &        ! qdos
                  -aimag(DENTOT(2))/PI,(-aimag(DEN(L1,IE,IQ,2))/PI,L1=0,LMAXD1)  ! qdos
               9000        format(5F10.6,40E16.8)                                ! qdos
               !                                                                 ! qdos
               if(test('compqdos')) then                                         ! complex qdos
                  if ((IQ.EQ.1).AND.(IE_num.EQ.1)) then                          ! complex qdos
                     if(t_inc%NATYP.ge.100) then                                 ! complex qdos
                        open(31,                                              &  ! complex qdos
                           FILE="cqdos."//char(48+I1/100)//                   &  ! complex qdos
                           char(48+mod(I1/10,10))//char(48+mod(I1,10))//"."   &  ! complex qdos
                           //char(48+1)//".dat")                                 ! complex qdos
                        open(32,                                              &  ! complex qdos
                           FILE="cqdos."//char(48+I1/100)//                   &  ! complex qdos
                           char(48+mod(I1/10,10))//char(48+mod(I1,10))//"."   &  ! complex qdos
                           //char(48+2)//".dat")                                 ! complex qdos
                     else                                                        ! complex qdos
                        open(31,FILE="cqdos."//char(48+mod(I1/10,10))//    &     ! complex qdos
                           char(48+mod(I1,10))//"."//char(48+1)//".dat")         ! complex qdos
                        open(32, FILE="cqdos."//char(48+mod(I1/10,10))//   &     ! complex qdos
                           char(48+mod(I1,10))//"."//char(48+2)//".dat")         ! complex qdos
                     end if                                                      ! complex qdos
                     call version_print_header(31)                               ! complex qdos
                     write(31,*) ' '                                             ! complex qdos
                     write(31,'(A)') '#   lmax, natyp, nspin, nqdos, ielast:'    ! complex qdos
                     write(31,'(5I9)') lmax, natyp, nspin, nqdos, ielast         ! complex qdos
                     write(31,'(7(A,3X))') '#   Re(E)','Im(E)',   &              ! complex qdos
                        'k_x','k_y','k_z','DEN_tot','DEN_s,p,...'                ! complex qdos
                     if(NSPIN.gt.1) then                                         ! complex qdos
                        call version_print_header(32)                            ! complex qdos
                        write(32,*) ' '                                          ! complex qdos
                        write(32,'(A)') '# lmax, natyp, nspin, nqdos, ielast:'   ! complex qdos
                        write(32,'(5I9)') lmax, natyp, nspin, nqdos, ielast      ! complex qdos
                        write(32,'(7(A,3X))') '#   Re(E)','Im(E)',   &           ! complex qdos
                           'k_x','k_y','k_z','DEN_tot','DEN_s,p,...'             ! complex qdos
                     end if                                                      ! complex qdos
                  endif   ! IQ.EQ.1                                              ! complex qdos
                  do JSPIN =1,2                                                  ! complex qdos
                     DENTOT(JSPIN) = DCMPLX(0.D0,0.D0, kind=dp)                           ! complex qdos
                     do L1 = 0,LMAXD1                                            ! complex qdos
                        DENTOT(JSPIN) = DENTOT(JSPIN) + DEN(L1,IE,IQ,JSPIN)      ! complex qdos
                     enddo                                                       ! complex qdos
                  enddo                                                          ! complex qdos
                  !    write qdos.nn.s.dat                                       ! complex qdos
                  write(31,9002) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),  &     ! complex qdos
                     DENTOT(1), (DEN(L1,IE,IQ,1),L1=0,LMAXD1)                    ! complex qdos
                  write(32,9002) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),  &     ! complex qdos
                     DENTOT(2), (DEN(L1,IE,IQ,2),L1=0,LMAXD1)                    ! complex qdos
                  9002          format(6F10.6,80E16.8)                           ! complex qdos
               end if                                                            ! complex qdos
               ! qdos
            enddo !IQ                                                            ! qdos
         enddo !IE                                                               ! qdos
      end if !myrank_at==master                                                  ! qdos
   endif      ! OPT('qdos    ')                                                  ! qdos
#endif

#ifdef CPP_MPI
   ! do communication only when compiled with MPI
#ifdef CPP_TIMING
   call timing_start('main1c - communication')
#endif
   ! reset NQDOS to avoid endless communication
   if( .not.OPT('lmdos    ')) then
      NQDOS=1
   else
      if(myrank==master) write(*,*) 'lmlm-dos option, communcation might take a while!', IEMXD,NQDOS
   endif
   ! set these arrays to zero to avoid double counting in cases where extra ranks are used
   if ( t_mpi_c_grid%myrank_ie>(t_mpi_c_grid%dims(1)-1) ) then
      den         = CZERO
      denlm       = CZERO
      gflle       = CZERO
      r2nefc      = CZERO
      rho2nsc     = CZERO
      rho2int     = CZERO
      muorb       = 0.0d0
      espv        = 0.0d0
      denorbmom   = 0.0d0
      denorbmomsp = 0.0d0
      denorbmomlm = 0.0d0
      denorbmomns = 0.0d0
   endif
   call mympi_main1c_comm_newsosol(IRMDNEW,LMPOT,LMAX,LMAXD1,  &
      LMMAXD,LMMAXSO,IELAST,NQDOS,                       &
      den,denlm,gflle,rho2nsc,r2nefc,                          &
      rho2int,espv,muorb,denorbmom,                            &
      denorbmomsp,denorbmomlm,denorbmomns,                     &
      t_mpi_c_grid%myMPI_comm_at)
#ifdef CPP_TIMING
   call timing_pause('main1c - communication')
#endif

   !MPI: do these writeout/data collection steps only on master and broadcast important results afterwards
   if(t_mpi_c_grid%myrank_at==master) then
#endif
! CPP_MPI
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDAU
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (IDOLDAU.EQ.1) then
         ! calculate WLDAU
         do IE=1,IELAST
            do LM1=1,LMMAXSO
               do LM2=1,LMMAXSO
                  GLDAU(LM1,LM2)=GLDAU(LM1,LM2)+GFLLE(LM1,LM2,IE,1)*WEZ(IE)/DBLE(NSPIN)
               enddo
            enddo
         enddo
         ! calculate occupation matrix
         MMAX=2*LOPT+1
         do IS=1,2
            do JS=1,2
               LMLO=LOPT**2+1+(IS-1)*LMMAXD
               LMHI=(LOPT+1)**2+(JS-1)*LMMAXD
               LM2=LOPT**2+1+(JS-1)*LMMAXD
               do M1=1,MMAX
                  LM1=LMLO-1+M1
                  DENMATN(1:MMAX,M1,JS,IS)=(1.0/(2.0*CI))*&
                     (GLDAU(LM2:LMHI,LM1)-CONJG(GLDAU(LM1,LM2:LMHI)))
               enddo
            enddo
         enddo
      endif ! LDAU
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDAU
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(.not.OPT('qdos    ')) then
         ! omp: moved write-out of dos files out of parallel energy loop
         ! Write out lm-dos:                                                     ! lm-dos
         if (OPT('lmdos    ')) then                                              ! qdos ruess
            do IE=1,IELAST                                                       ! lm-dos
               IQ=1                                                              ! lm-dos
               if (IE.EQ.1) then                                                 ! lm-dos
                  if(t_inc%NATYP.ge.100) then                                    ! lm-dos
                     open(29,                               &                    ! lm-dos
                        FILE="lmdos."//char(48+I1/100)//    &                    ! lm-dos
                        char(48+mod(I1/10,10))//            &                    ! lm-dos
                        char(48+mod(I1,10))//"."//char(48+1)//".dat")            ! lm-dos
                     open(30,                               &                    ! lm-dos
                        FILE="lmdos."//char(48+I1/100)//    &                    ! lm-dos
                        char(48+mod(I1/10,10))//            &                    ! lm-dos
                        char(48+mod(I1,10))//"."//char(48+2)//".dat")            ! lm-dos
                  else                                                           ! lm-dos
                     open(29,                                  &                 ! lm-dos
                        FILE="lmdos."//char(48+I1/10)//        &                 ! lm-dos
                        char(48+mod(I1,10))//"."//char(48+1)//".dat")            ! lm-dos
                     open(30,FILE="lmdos."//char(48+I1/10)//   &                 ! lm-dos
                        char(48+mod(I1,10))//"."//char(48+2)//".dat")            ! lm-dos
                  end if                                                         ! lm-dos
                  call version_print_header(29)                                  ! lm-dos
                  write (29,*) ' '                                               ! lm-dos
                  write (29,8600) '# ISPIN=',1,' I1=',I1                         ! lm-dos
                  call version_print_header(30)                                  ! lm-dos
                  write (30,*) ' '                                               ! lm-dos
                  write (30,8600) '# ISPIN=',2,' I1=',I1                         ! lm-dos
               endif !IE==1                                                      ! lm-dos
               write(29,9001) EZ(IE),(-aimag(DENLM(L1,IE,IQ,1))/PI,L1=1,LMMAXD)  ! lm-dos
               write(30,9001) EZ(IE),(-aimag(DENLM(L1,IE,IQ,2))/PI,L1=1,LMMAXD)  ! lm-dos
               9001    format(30E12.4)                                           ! lm-dos
               8600    format (a8,I3,a4,I5)                                      ! lm-dos/qdos ruess
            enddo !IE
         endif
      end if ! .not. OPT('qdos    ')

      !      write gflle to file                                                 ! lmlm-dos
      if (OPT('lmlm-dos')) then                                                  ! lmlm-dos
         if(t_inc%i_write>0) then                                                ! lmlm-dos
            write(1337,*) 'gflle:',shape(gflle),shape(gflle_part),LRECGFLLE      ! lmlm-dos
         endif                                                                   ! lmlm-dos
         write(91,REC=I1) GFLLE                                                  ! lmlm-dos
      endif                                                                      ! lmlm-dos
      !
      allocate(RHOTEMP(IRMDNEW,LMPOT),stat=i_stat)
      call memocc(i_stat,product(shape(RHOTEMP))*kind(RHOTEMP),'RHOTEMP','RHOVALNEW')
      allocate(RHONEWTEMP(IRWS,LMPOT),stat=i_stat)
      call memocc(i_stat,product(shape(RHONEWTEMP))*kind(RHONEWTEMP),'RHONEWTEMP','RHOVALNEW')
      !
      do JSPIN=1,4
         RHOTEMP=CZERO
         RHONEWTEMP=CZERO
         do LM1=1,LMPOT
            do IR=1,IRMDNEW
               RHOTEMP(IR,LM1)=RHO2NSC(IR,LM1,JSPIN)
            enddo
         enddo
         call CHEB2OLDGRID(IRWS,IRMDNEW,LMPOT,RMESH,NCHEB,NPAN_TOT,&
            RPAN_INTERVALL,IPAN_INTERVALL,RHOTEMP,RHONEWTEMP,IRMD)
         do LM1=1,LMPOT
            do IR=1,IRWS
               RHO2NSNEW(IR,LM1,JSPIN)=RHONEWTEMP(IR,LM1)
            enddo
         enddo

         RHOTEMP=CZERO
         RHONEWTEMP=CZERO
         do LM1=1,LMPOT
            do IR=1,IRMDNEW
               RHOTEMP(IR,LM1)=R2NEFC(IR,LM1,JSPIN)
            enddo
         enddo
         call CHEB2OLDGRID(IRWS,IRMDNEW,LMPOT,RMESH,NCHEB,NPAN_TOT,&
            RPAN_INTERVALL,IPAN_INTERVALL,RHOTEMP,RHONEWTEMP,IRMD)
         do LM1=1,LMPOT
            do IR=1,IRWS
               R2NEFNEW(IR,LM1,JSPIN)=RHONEWTEMP(IR,LM1)
            enddo
         enddo
      enddo
      i_all=-product(shape(RHOTEMP))*kind(RHOTEMP)
      deallocate(RHOTEMP, stat=i_stat)
      call memocc(i_stat,i_all,'RHOTEMP','RHOVALNEW')
      i_all=-product(shape(RHONEWTEMP))*kind(RHONEWTEMP)
      deallocate(RHONEWTEMP, stat=i_stat)
      call memocc(i_stat,i_all,'RHONEWTEMP','RHOVALNEW')
      ! calculate new THETA and PHI for non-colinear
      if (.NOT.TEST('FIXMOM  ')) then
         RHO2NS_TEMP(1,1)=RHO2INT(1)
         RHO2NS_TEMP(2,2)=RHO2INT(2)
         RHO2NS_TEMP(1,2)=RHO2INT(3)
         RHO2NS_TEMP(2,1)=RHO2INT(4)
         !
         call ROTATEMATRIX(RHO2NS_TEMP,THETA,PHI,1,0)
         !
         RHO2INT(1)=RHO2NS_TEMP(1,1)
         RHO2INT(2)=RHO2NS_TEMP(2,2)
         RHO2INT(3)=RHO2NS_TEMP(1,2)
         RHO2INT(4)=RHO2NS_TEMP(2,1)
         !
         MOMENT(1)=aimag(RHO2INT(3)+RHO2INT(4))
         MOMENT(2)=-real(RHO2INT(3)-RHO2INT(4))
         MOMENT(3)=aimag(-RHO2INT(1)+RHO2INT(2))
         !
         TOTMOMENT=SQRT(MOMENT(1)**2+MOMENT(2)**2+MOMENT(3)**2)
         TOTXYMOMENT=SQRT(MOMENT(1)**2+MOMENT(2)**2)
         !
         if (ABS(TOTXYMOMENT).GT.1d-05) then
            if (ABS(MOMENT(3)).LT.1d-05) then
               THETANEW=PI/2d0
            else
               THETANEW=ACOS(MOMENT(3)/TOTMOMENT)
            endif
            if (TOTXYMOMENT.LT.1d-05) then
               PHINEW=0d0
            else
               PHINEW=DATAN2(MOMENT(2),MOMENT(1))
            endif
         endif
         !
         if(t_inc%i_write>0) then
            write(1337,*) 'moment',myrank,MOMENT(1),MOMENT(2),MOMENT(3)
            write(1337,*) THETANEW/(2.0D0*PI)*360.0D0,PHINEW/(2.0D0*PI)*360.0D0
         endif
         !only on master different from zero:
         angles_new(1) = THETANEW
         angles_new(2) = PHINEW
         call ROTATEVECTOR(RHO2NSNEW,RHO2NS,IRWS,LMPOT,THETANEW,PHINEW,&
            THETA,PHI,IRMD)
         call ROTATEVECTOR(R2NEFNEW,R2NEF,IRWS,LMPOT,THETANEW,PHINEW,  &
            THETA,PHI,IRMD)
      else
         RHO2NS(:,:,:)=aimag(RHO2NSNEW(:,:,:))
         R2NEF(:,:,:)=aimag(R2NEFNEW(:,:,:))
      endif
      !
      IDIM = IRMD*LMPOT
      call DSCAL(IDIM,2.D0,RHO2NS(1,1,1),1)
      call DAXPY(IDIM,-0.5D0,RHO2NS(1,1,1),1,RHO2NS(1,1,2),1)
      call DAXPY(IDIM,1.0D0,RHO2NS(1,1,2),1,RHO2NS(1,1,1),1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Do the same at the Fermi energy
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call DSCAL(IDIM,2.D0,R2NEF(1,1,1),1)
      call DAXPY(IDIM,-0.5D0,R2NEF(1,1,1),1,R2NEF(1,1,2),1)
      call DAXPY(IDIM,1.0D0,R2NEF(1,1,2),1,R2NEF(1,1,1),1)
      !
      do LM1=0,LMAXD1
         do IE=1,IEMXD
            do JSPIN=1,NSPIN
               DEN_out(LM1,IE,JSPIN) =  DEN(LM1,IE,1,JSPIN)
            enddo
         enddo
      enddo
      !
#ifdef CPP_MPI
   endif !(myrank==master)

   ! communicate den_out to all processors with the same atom number
   IDIM = (LMAX+2)*IEMXD*2
   call MPI_Bcast(den_out, idim, MPI_DOUBLE_COMPLEX, master,&
      t_mpi_c_grid%myMPI_comm_at, ierr)
   if(ierr/=MPI_SUCCESS) stop 'error bcast den_out in rhovalnew'
   IDIM = 2
   call MPI_Bcast(angles_new, idim, MPI_DOUBLE_PRECISION, master,&
      t_mpi_c_grid%myMPI_comm_at, ierr)
   if(ierr/=MPI_SUCCESS) stop 'error bcast angles_new in rhovalnew'
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Deallocate arrays
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   i_all=-product(shape(VINS))*kind(VINS)
   deallocate(VINS, stat=i_stat)
   call memocc(i_stat,i_all,'VINS','RHOVALNEW')
   i_all=-product(shape(VNSPLL0))*kind(VNSPLL0)
   deallocate(VNSPLL0, stat=i_stat)
   call memocc(i_stat,i_all,'VNSPLL0','RHOVALNEW')
   i_all=-product(shape(VNSPLL1))*kind(VNSPLL1)
   deallocate(VNSPLL1, stat=i_stat)
   call memocc(i_stat,i_all,'VNSPLL1','RHOVALNEW')
   i_all=-product(shape(VNSPLL))*kind(VNSPLL)
   deallocate(VNSPLL, stat=i_stat)
   call memocc(i_stat,i_all,'VNSPLL','RHOVALNEW')
   i_all=-product(shape(HLK))*kind(HLK)
   deallocate(HLK, stat=i_stat)
   call memocc(i_stat,i_all,'HLK','RHOVALNEW')
   i_all=-product(shape(JLK))*kind(JLK)
   deallocate(JLK, stat=i_stat)
   call memocc(i_stat,i_all,'JLK','RHOVALNEW')
   i_all=-product(shape(HLK2))*kind(HLK2)
   deallocate(HLK2, stat=i_stat)
   call memocc(i_stat,i_all,'HLK2','RHOVALNEW')
   i_all=-product(shape(JLK2))*kind(JLK2)
   deallocate(JLK2, stat=i_stat)
   call memocc(i_stat,i_all,'JLK2','RHOVALNEW')
   i_all=-product(shape(TMATSPH))*kind(TMATSPH)
   deallocate(TMATSPH, stat=i_stat)
   call memocc(i_stat,i_all,'TMATSPH','RHOVALNEW')
   i_all=-product(shape(RLL))*kind(RLL)
   deallocate(RLL, stat=i_stat)
   call memocc(i_stat,i_all,'RLL','RHOVALNEW')
   i_all=-product(shape(SLL))*kind(SLL)
   deallocate(SLL, stat=i_stat)
   call memocc(i_stat,i_all,'SLL','RHOVALNEW')
   i_all=-product(shape(RLLLEFT))*kind(RLLLEFT)
   deallocate(RLLLEFT, stat=i_stat)
   call memocc(i_stat,i_all,'RLLLEFT','RHOVALNEW')
   i_all=-product(shape(SLLLEFT))*kind(SLLLEFT)
   deallocate(SLLLEFT, stat=i_stat)
   call memocc(i_stat,i_all,'SLLLEFT','RHOVALNEW')
   i_all=-product(shape(CDEN))*kind(CDEN)
   deallocate(CDEN, stat=i_stat)
   call memocc(i_stat,i_all,'CDEN','RHOVALNEW')
   i_all=-product(shape(CDENLM))*kind(CDENLM)
   deallocate(CDENLM, stat=i_stat)
   call memocc(i_stat,i_all,'CDENLM','RHOVALNEW')
   i_all=-product(shape(CDENNS))*kind(CDENNS)
   deallocate(CDENNS, stat=i_stat)
   call memocc(i_stat,i_all,'CDENNS','RHOVALNEW')
   i_all=-product(shape(RHO2NSC))*kind(RHO2NSC)
   deallocate(RHO2NSC, stat=i_stat)
   call memocc(i_stat,i_all,'RHO2NSC','RHOVALNEW')
   i_all=-product(shape(RHO2NSC_loop))*kind(RHO2NSC_loop)
   deallocate(RHO2NSC_loop, stat=i_stat)
   call memocc(i_stat,i_all,'RHO2NSC_loop','RHOVALNEW')
   i_all=-product(shape(RHO2NSNEW))*kind(RHO2NSNEW)
   deallocate(RHO2NSNEW, stat=i_stat)
   call memocc(i_stat,i_all,'RHO2NSNEW','RHOVALNEW')
   i_all=-product(shape(R2NEFC))*kind(R2NEFC)
   deallocate(R2NEFC, stat=i_stat)
   call memocc(i_stat,i_all,'R2NEFC','RHOVALNEW')
   i_all=-product(shape(R2NEFC_loop))*kind(R2NEFC_loop)
   deallocate(R2NEFC_loop, stat=i_stat)
   call memocc(i_stat,i_all,'R2NEFC_loop','RHOVALNEW')
   i_all=-product(shape(R2NEFNEW))*kind(R2NEFNEW)
   deallocate(R2NEFNEW, stat=i_stat)
   call memocc(i_stat,i_all,'R2NEFNEW','RHOVALNEW')
   i_all=-product(shape(R2ORBC))*kind(R2ORBC)
   deallocate(R2ORBC, stat=i_stat)
   call memocc(i_stat,i_all,'R2ORBC','RHOVALNEW')
   i_all=-product(shape(CDENTEMP))*kind(CDENTEMP)
   deallocate(CDENTEMP, stat=i_stat)
   call memocc(i_stat,i_all,'CDENTEMP','RHOVALNEW')
   i_all=-product(shape(GFLLE_PART))*kind(GFLLE_PART)
   deallocate(GFLLE_PART, stat=i_stat)
   call memocc(i_stat,i_all,'GFLLE_PART','RHOVALNEW')
   i_all=-product(shape(GFLLE))*kind(GFLLE)
   deallocate(GFLLE, stat=i_stat)
   call memocc(i_stat,i_all,'GFLLE','RHOVALNEW')
   i_all=-product(shape(DEN))*kind(DEN)
   deallocate(DEN, stat=i_stat)
   call memocc(i_stat,i_all,'DEN','RHOVALNEW')
   i_all=-product(shape(DENLM))*kind(DENLM)
   deallocate(DENLM, stat=i_stat)
   call memocc(i_stat,i_all,'DENLM','RHOVALNEW')

end subroutine RHOVALNEW
