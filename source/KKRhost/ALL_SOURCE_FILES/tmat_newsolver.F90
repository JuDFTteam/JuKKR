!-------------------------------------------------------------------------------
! SUBROUTINE: TMAT_NEWSOLVER
!> @brief Calculation of the T-Matrix
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,ZAT,SOCSCALE,EZ,NSRA,CLEB,ICLEB,  &
   IEND,NCHEB,NPAN_TOT,RPAN_INTERVALL,IPAN_INTERVALL,RNEW,VINSNEW,THETA,PHI,  &
   I1,IPOT,LMPOT,LLY,DELTAE,IDOLDAU,LOPT,WLDAU,t_dtmatJij_at)

#ifdef CPP_OMP
   use omp_lib        ! necessary for omp functions
#endif
#ifdef CPP_MPI
   use mpi
#endif
   use mod_mympi, only: myrank, nranks, master
#ifdef CPP_MPI
   use mod_mympi, only: distribute_linear_on_tasks, MPIadapt
   use mod_timing
#endif
   use mod_types, only: t_tgmat,t_inc,t_mpi_c_grid,init_tgmat, &
                        t_lloyd,init_tlloyd, type_dtmatJijDij, &
                        init_t_dtmatJij_at
   use Constants
   use Profiling
   use mod_wunfiles, only: t_params
   use mod_jijhelp, only: calc_dtmatJij
   use mod_save_wavefun, only: t_wavefunctions, find_isave_wavefun,save_wavefunc
   use global_variables
   use mod_DataTypes

   implicit none

   integer, intent(in) :: I1
   integer, intent(in) :: LLY       !< LLY <> 0: apply Lloyds formula
   integer, intent(in) :: LOPT      !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: NSRA
   integer, intent(in) :: IEND      !< Number of nonzero gaunt coefficients
   integer, intent(in) :: IPOT
   integer, intent(in) :: NCHEB     !< Number of Chebychev pannels for the new solver
   integer, intent(in) :: NSPIN     !< Counter for spin directions
   integer, intent(in) :: LMPOT     !< (LPOT+1)**2
   integer, intent(in) :: IELAST
   integer, intent(in) :: NPAN_TOT
   integer, intent(in) :: IDOLDAU   !< flag to perform LDA+U
   real (kind=dp), intent(in) :: ZAT       !< Nuclear charge for a given atom
   real (kind=dp), intent(in) :: PHI
   real (kind=dp), intent(in) :: THETA
   real (kind=dp), intent(in) :: SOCSCALE  !< Spin-orbit scaling for a given atom
   complex (kind=dp), intent(in) :: DELTAE      !< Energy difference for numerical derivative
   integer, dimension(0:NTOTD), intent(in) :: IPAN_INTERVALL
   integer, dimension(NCLEB,4), intent(in) :: ICLEB
   real (kind=dp), dimension(NCLEB), intent(in) :: CLEB !< GAUNT coefficients (GAUNT)
   real (kind=dp), dimension(NRMAXD), intent(in) :: RNEW
   real (kind=dp), dimension(0:NTOTD), intent(in) :: RPAN_INTERVALL
   real (kind=dp), dimension(MMAXD,MMAXD,NSPIND), intent(in) :: WLDAU !< potential matrix
   real (kind=dp), dimension(NRMAXD,LMPOT,NSPOTD), intent(in) :: VINSNEW
   complex (kind=dp), dimension(IEMXD), intent(in) :: EZ
   ! .. In/Out variables
   type(type_dtmatJijDij), intent(inout) :: t_dtmatJij_at

   ! .. Local variables
   integer :: IR,IREC,USE_SRATRICK,NVEC,LM1,LM2,IE,IRMDNEW
   integer :: i_stat, i_all, lmsize
   complex (kind=dp) :: ERYD
   complex (kind=dp), dimension(2*(LMAX+1)) :: ALPHASPH
   ! .. Local allocatable arrays
   integer, dimension(:), allocatable :: JLK_INDEX
   real (kind=dp), dimension(:,:,:), allocatable :: VINS  !< Non-spherical part of the potential
   complex (kind=dp), dimension(:,:), allocatable     :: AUX      ! LLY
   complex (kind=dp), dimension(:,:), allocatable     :: TMAT0
   complex (kind=dp), dimension(:,:), allocatable     :: ALPHA0   ! LLY
   complex (kind=dp), dimension(:,:), allocatable     :: TMATLL
   complex (kind=dp), dimension(:,:), allocatable     :: DTMATLL
   complex (kind=dp), dimension(:,:), allocatable     :: TMATSPH
   complex (kind=dp), dimension(:,:), allocatable     :: ALPHALL  ! LLY
   complex (kind=dp), dimension(:,:), allocatable     :: DALPHALL ! LLY
   complex (kind=dp), dimension(:,:,:), allocatable   :: HLK
   complex (kind=dp), dimension(:,:,:), allocatable   :: JLK
   complex (kind=dp), dimension(:,:,:), allocatable   :: HLK2
   complex (kind=dp), dimension(:,:,:), allocatable   :: JLK2
   complex (kind=dp), dimension(:,:,:), allocatable   :: VNSPLL0
   complex (kind=dp), dimension(:,:,:,:), allocatable :: RLL
   complex (kind=dp), dimension(:,:,:,:), allocatable :: SLL
   complex (kind=dp), dimension(:,:,:,:), allocatable :: VNSPLL
   complex (kind=dp), dimension(:,:,:,:), allocatable :: VNSPLL1
   complex (kind=dp), dimension(:,:,:,:), allocatable :: RLLLEFT
   complex (kind=dp), dimension(:,:,:,:), allocatable :: SLLLEFT

   ! .. LDAU local variables
   integer :: LMLO,LMHI
   ! .. LLoyd local variables
   integer :: IDERIV,SIGNDE
   complex (kind=dp) :: TRALPHA
   complex (kind=dp) :: GMATPREFACTOR
   integer, dimension(:), allocatable :: IPIV    ! LLY
   ! .. OMP local variables
   integer :: nth,ith ! total number of threads and thread id

#ifdef CPP_MPI
   integer, dimension(0:nranks-1) :: ntot_pT, ioff_pT
#endif
   integer :: ie_end, ie_num, ie_start, ierr
      
   !rhoqtest
   logical, external :: test, opt
   integer :: mu0, nscoef

   ! BdG
   character(len=100) :: filename

   lmsize = lmmaxd/2

   ! .. Allocation of local arrays
   allocate(AUX(LMMAXSO,LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(AUX))*kind(AUX),'AUX','tmat_newsolver')
   AUX=CZERO
   allocate(IPIV(LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(IPIV))*kind(IPIV),'IPIV','tmat_newsolver')
   IPIV=0
   allocate(TMAT0(LMMAXSO,LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(TMAT0))*kind(TMAT0),'TMAT0','tmat_newsolver')
   TMAT0=CZERO
   allocate(TMATLL(LMMAXSO,LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(TMATLL))*kind(TMATLL),'TMATLL','tmat_newsolver')
   TMATLL=CZERO
   allocate(ALPHA0(LMMAXSO,LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(ALPHA0))*kind(ALPHA0),'ALPHA0','tmat_newsolver')
   ALPHA0=CZERO
   allocate(DTMATLL(LMMAXSO,LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(DTMATLL))*kind(DTMATLL),'DTMATLL','tmat_newsolver')
   DTMATLL=CZERO
   allocate(ALPHALL(LMMAXSO,LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(ALPHALL))*kind(ALPHALL),'ALPHALL','tmat_newsolver')
   ALPHALL=CZERO
   allocate(DALPHALL(LMMAXSO,LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(DALPHALL))*kind(DALPHALL),'DALPHALL','tmat_newsolver')
   DALPHALL=CZERO
   allocate(JLK_INDEX(2*LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(JLK_INDEX))*kind(JLK_INDEX),'JLK_INDEX','tmat_newsolver')
   JLK_INDEX=0

#ifdef CPP_OMP
   ! determine if omp parallelisation is used (compiled with -openmp flag and OMP_NUM_THREADS>1)
   !$omp parallel shared(nth,ith)
   !$omp single
   nth = omp_get_num_threads()
   if(t_inc%i_write>0) write(1337,*) 'nth =',nth
   !$omp end single
   !$omp end parallel
#else
   nth = 1
   ith = 0
#endif

   IRMDNEW= NPAN_TOT*(NCHEB+1)
   allocate(VINS(IRMDNEW,LMPOT,NSPIN),stat=i_stat)
   call memocc(i_stat,product(shape(VINS))*kind(VINS),'VINS','tmat_newsolver')
   VINS=0.0d0

   do LM1=1,LMPOT
      do IR=1,IRMDNEW
         VINS(IR,LM1,1)=VINSNEW(IR,LM1,IPOT)
         VINS(IR,LM1,NSPIN)=VINSNEW(IR,LM1,IPOT+NSPIN-1)
      enddo
   enddo
   ! set up the non-spherical ll' matrix for potential VLL'
   if (NSRA.EQ.2) then
      USE_SRATRICK=1
   elseif (NSRA.EQ.1) then
      USE_SRATRICK=0
   endif

   allocate(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW),stat=i_stat)
   call memocc(i_stat,product(shape(VNSPLL0))*kind(VNSPLL0),'VNSPLL0','tmat_newsolver')
   VNSPLL0=CZERO
   allocate(VNSPLL1(LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(VNSPLL1))*kind(VNSPLL1),'VNSPLL1','tmat_newsolver')
   VNSPLL1=CZERO

   VNSPLL0=CZERO
   call VLLMAT(1,NRMAXD,IRMDNEW,lmsize,LMMAXSO,VNSPLL0,VINS,LMPOT,CLEB,ICLEB,&
      IEND,NSPIN,ZAT,RNEW,USE_SRATRICK,NCLEB)
   ! LDAU
   if (IDOLDAU.EQ.1) then
      LMLO=LOPT**2+1
      LMHI=(LOPT+1)**2
      do IR=1,IRMDNEW
         VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)=VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)+  &
            WLDAU(1:MMAXD,1:MMAXD,1)
      enddo
      LMLO=LMLO+lmsize
      LMHI=LMHI+lmsize
      do IR=1,IRMDNEW
         VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)=VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)+  &
            WLDAU(1:MMAXD,1:MMAXD,2)
      enddo
   endif
   ! LDAU

   ! initial allocate
   ! potential array
   if (NSRA.EQ.2) THEN
      allocate(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
      call memocc(i_stat,product(shape(VNSPLL))*kind(VNSPLL),'VNSPLL','tmat_newsolver')
      VNSPLL=CZERO
   else
      allocate(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
      call memocc(i_stat,product(shape(VNSPLL))*kind(VNSPLL),'VNSPLL','tmat_newsolver')
      VNSPLL=CZERO
   endif

   ! source terms (bessel and hankel functions)
   allocate(HLK(1:4*(LMAX+1),IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(HLK))*kind(HLK),'HLK','tmat_newsolver')
   HLK=CZERO
   allocate(JLK(1:4*(LMAX+1),IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(JLK))*kind(JLK),'JLK','tmat_newsolver')
   JLK=CZERO
   allocate(HLK2(1:4*(LMAX+1),IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(HLK2))*kind(HLK2),'HLK2','tmat_newsolver')
   HLK2=CZERO
   allocate(JLK2(1:4*(LMAX+1),IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(JLK2))*kind(JLK2),'JLK2','tmat_newsolver')
   JLK2=CZERO

   ! Spherical part of tmatrix (used with SRATRICK)
   allocate(TMATSPH(2*(LMAX+1),0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(TMATSPH))*kind(TMATSPH),'TMATSPH','tmat_newsolver')
   TMATSPH=CZERO

   ! Regular and irregular wavefunctions
   allocate(RLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(RLL))*kind(RLL),'RLL','tmat_newsolver')
   RLL=CZERO
   allocate(SLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
   call memocc(i_stat,product(shape(SLL))*kind(SLL),'SLL','tmat_newsolver')
   SLL=CZERO

   ! Left regular and irregular wavefunctions (used here only in case of XCPL or saving of left wavefunctions)
   if( opt('XCPL    ') .or. (t_wavefunctions%save_rllleft .or.t_wavefunctions%save_sllleft .or. test('rhoqtest')) ) then
      allocate(RLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
      call memocc(i_stat,product(shape(RLLLEFT))*kind(RLLLEFT),'RLLLEFT','tmat_newsolver')
      RLLLEFT=CZERO
      allocate(SLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),stat=i_stat)
      call memocc(i_stat,product(shape(SLLLEFT))*kind(SLLLEFT),'SLLLEFT','tmat_newsolver')
      SLLLEFT=CZERO

   else
      allocate(RLLLEFT(1,1,1,0:nth-1),stat=i_stat)
      call memocc(i_stat,product(shape(RLLLEFT))*kind(RLLLEFT),'RLLLEFT','tmat_newsolver')
      RLLLEFT=CZERO
      allocate(SLLLEFT(1,1,1,0:nth-1),stat=i_stat)
      call memocc(i_stat,product(shape(SLLLEFT))*kind(SLLLEFT),'SLLLEFT','tmat_newsolver')
      SLLLEFT=CZERO
   end if ! ( opt('XCPL    ') .or. ... )

   ! Energy loop
   if(myrank==master.and.(t_inc%i_write>0)) WRITE(1337,*) 'atom: ',I1,' NSRA:',NSRA

#ifdef CPP_MPI
   call distribute_linear_on_tasks(t_mpi_c_grid%nranks_at,  &
      t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at+(i1-1), & ! print this info only for first atom at master
      master,IELAST,ntot_pT,ioff_pT,.true.,.true.)

   ie_start = ioff_pT(t_mpi_c_grid%myrank_at)
   ie_end   = ntot_pT(t_mpi_c_grid%myrank_at)

   t_mpi_c_grid%ntot2=ie_end   !t_mpi_c_grid%dims(1)
   if (.not.(allocated(t_mpi_c_grid%ntot_pT2).or.allocated(t_mpi_c_grid%ioff_pT2))) then
      allocate(t_mpi_c_grid%ntot_pT2(0:t_mpi_c_grid%nranks_at-1),stat=i_stat)
      call memocc(i_stat,product(shape(t_mpi_c_grid%ntot_pT2))*kind(t_mpi_c_grid%ntot_pT2),'t_mpi_c_grid%ntot_pT2','tmat_newsolver')
      t_mpi_c_grid%ntot_pT2=0
      allocate(t_mpi_c_grid%ioff_pT2(0:t_mpi_c_grid%nranks_at-1),stat=i_stat)
      call memocc(i_stat,product(shape(t_mpi_c_grid%ioff_pT2))*kind(t_mpi_c_grid%ioff_pT2),'t_mpi_c_grid%ioff_pT2','tmat_newsolver')
      t_mpi_c_grid%ioff_pT2=0
   endif
   t_mpi_c_grid%ntot_pT2 = ntot_pT
   t_mpi_c_grid%ioff_pT2 = ioff_pT
   ! Now initialize arrays for tmat, gmat, and gref
   call init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
   if(lly.ne.0) call init_tlloyd(t_inc,t_lloyd,t_mpi_c_grid)
       
   if(test('rhoqtest')) then
      if(ielast/=3) stop 'Error: wrong energy contour for rhoqtest'
      ie_start=1
      ie_end=1
   end if
#else
   if(.not.(allocated(t_mpi_c_grid%ntot_pT2).or.allocated(t_mpi_c_grid%ioff_pT2))) then
      allocate(t_mpi_c_grid%ntot_pT2(1),stat=i_stat)
      call memocc(i_stat,product(shape(t_mpi_c_grid%ntot_pT2))*kind(t_mpi_c_grid%ntot_pT2),'t_mpi_c_grid%ntot_pT2','tmat_newsolver')
      t_mpi_c_grid%ntot_pT2=0
      allocate(t_mpi_c_grid%ioff_pT2(1),stat=i_stat)
      call memocc(i_stat,product(shape(t_mpi_c_grid%ioff_pT2))*kind(t_mpi_c_grid%ioff_pT2),'t_mpi_c_grid%ioff_pT2','tmat_newsolver')
      t_mpi_c_grid%ioff_pT2=0
   endif
   t_mpi_c_grid%ntot2      =IELAST
   t_mpi_c_grid%ntot_pT2   = IELAST
   t_mpi_c_grid%ioff_pT2   = 0
   ! now initialize arrays for tmat, gmat, and gref
   call init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
   if(lly.ne.0) call init_tlloyd(t_inc,t_lloyd,t_mpi_c_grid)

   ie_start = 0
   ie_end   = IELAST
#endif

   ! For Jij-tensor calculation: allocate array to hold additional t-matrices
   call init_t_dtmatJij_at(t_inc, t_mpi_c_grid, t_dtmatJij_at)

   ! Initialize wfsave
   if(t_inc%i_iteration==0) then
      call find_isave_wavefun(t_wavefunctions)
      ! Reset Nwfsavemax to 0 if test option 'STOP1B  ' is found
      ! to prevent unnessesary storing of wavefunctions
      if(test('STOP1B  ') .and. .not. opt('OPERATOR')) then
         t_wavefunctions%Nwfsavemax = 0
      endif
   endif

#ifdef CPP_OMP
   !$omp parallel do default(none)                                            &
   !$omp private(eryd,ie,ir,nvec,lm1,lm2,gmatprefactor)                       &
   !$omp private(jlk_index,tmatll,ith,irec, ie_num)                           &
   !$omp private(tralpha, aux, ideriv, ipiv)                                  &
   !$omp private(alpha0)                                                      &
   !$omp private(alphall)                                                     &
   !$omp private(tmat0)                                                       &
   !$omp private(alphasph)                                                    &
   !$omp private(dtmatll)                                                     &
   !$omp private(dalphall)                                                    &
   !$omp shared(t_inc)                                                        &
   !!$omp firstprivate(t_inc)                                                 &
   !$omp shared(nspin,nsra,lmax,lmsize,iend,ipot,ielast,npan_tot,ncheb)       &
   !$omp shared(zat,socscale,ez,cleb,rnew,nth,LMPOT,NRMAXD,LMMAXSO,NTOTD)     &
   !$omp shared(rpan_intervall,vinsnew,ipan_intervall,NCLEB)                  &
   !$omp shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0)                  &
   !$omp shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,rllleft,sllleft)     &
   !$omp shared(tmatsph, ie_end,t_tgmat,t_lloyd, ie_start, t_dtmatjij_at)     &
   !$omp shared(lly,deltae,i1,t_mpi_c_grid, t_wavefunctions, icleb)           &
   !$omp shared(mu0, nscoef)
#endif

   do ie_num=1,ie_end

      IE = ie_start+ie_num

#ifdef CPP_MPI
      !start timing measurement for this pair of ie and i1, needed for MPIadapt
      call timing_start('time_1a_ieiatom')
#endif

      ! get current thread
      if (nth>=1) then
#ifdef CPP_OMP
         ith = omp_get_thread_num()
#endif
      else
         ith = 0
      endif

      ! In case of Lloyds formula the derivative of t is needed.
      ! Then calculate t at E+dE, E-dE and average for t, subtract for dt/dE
      TMATLL=CZERO
      ALPHALL=CZERO ! LLY
      DTMATLL=CZERO ! LLY
      DALPHALL=CZERO ! LLY
      IDERIV=0
      if (LLY.NE.0) IDERIV=1
      do SIGNDE=-IDERIV,IDERIV,2
         ERYD = EZ(IE)+real(SIGNDE, kind=dp)*DELTAE/2.D0 ! LLY
#ifdef CPP_OMP
         !$omp critical
#endif
         if(t_inc%i_write>0) WRITE(1337,*) 'energy:',IE,'',ERYD
#ifdef CPP_OMP
         if(ie==1.and.(t_inc%i_write>0)) write(1337,*) 'nested omp?',omp_get_nested()
         !$omp end critical
#endif

         ! Contruct the spin-orbit coupling hamiltonian and add to potential
         call SPINORBIT_HAM(LMAX,lmsize,VINS,RNEW,ERYD,ZAT,CVLIGHT,SOCSCALE,  &
            NSPIN,LMPOT,THETA,PHI,IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,    &
            NCHEB,IRMDNEW,NRMAXD,VNSPLL0(:,:,:),VNSPLL1(:,:,:,ith),'1')
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

         ! Calculate the source terms in the Lippmann-Schwinger equation
         ! these are spherical hankel and bessel functions
         HLK(:,:,ith)=CZERO
         JLK(:,:,ith)=CZERO
         HLK2(:,:,ith)=CZERO
         JLK2(:,:,ith)=CZERO
         GMATPREFACTOR=CZERO
         call RLLSLLSOURCETERMS(NSRA,NVEC,ERYD,RNEW,IRMDNEW,NRMAXD,LMAX,LMMAXSO, &
            1,JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),   &
            GMATPREFACTOR)

         if (test('BdG_dev ')) then
           write(filename, '(A,I0.3,A,I0.3,A)') 'rll_source_jlk_atom_',i1,'_energ_',ie,'.dat'
           open(888888, file=trim(filename), form='formatted')
           write(888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=',4*(LMAX+1),' IRMDNEW=', IRMDNEW, ' ; ERYD=', ERYD
           write(888888, '(2ES21.9)') jlk(:,:,ith)
           close(888888)
           write(filename, '(A,I0.3,A,I0.3,A)') 'rll_source_hlk_atom_',i1,'_energ_',ie,'.dat'
           open(888888, file=trim(filename), form='formatted')
           write(888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=',4*(LMAX+1),' IRMDNEW=', IRMDNEW, ' ; ERYD=', ERYD
           write(888888, '(2ES21.9)') hlk(:,:,ith)
           close(888888)
           write(filename, '(A,I0.3,A,I0.3,A)') 'rll_source_jlk2_atom_',i1,'_energ_',ie,'.dat'
           open(888888, file=trim(filename), form='formatted')
           write(888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=',4*(LMAX+1),' IRMDNEW=', IRMDNEW, ' ; ERYD=', ERYD
           write(888888, '(2ES21.9)') jlk2(:,:,ith)
           close(888888)
           write(filename, '(A,I0.3,A,I0.3,A)') 'rll_source_hlk2_atom_',i1,'_energ_',ie,'.dat'
           open(888888, file=trim(filename), form='formatted')
           write(888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=',4*(LMAX+1),' IRMDNEW=', IRMDNEW, ' ; ERYD=', ERYD
           write(888888, '(2ES21.9)') hlk2(:,:,ith)
           close(888888)
         end if

         ! Using spherical potential as reference
         if (USE_SRATRICK.EQ.1) then
            TMATSPH(:,ith)=CZERO
            call CALCSPH(NSRA,IRMDNEW,NRMAXD,LMAX,NSPIN,ZAT,ERYD,LMPOT, &
               LMMAXSO,RNEW,VINS,NCHEB,NPAN_TOT,RPAN_INTERVALL,JLK_INDEX,        &
               HLK(:,:,ith),JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),            &
               GMATPREFACTOR,TMATSPH(:,ith),ALPHASPH,USE_SRATRICK)
         endif
         ! Calculate the tmat and wavefunctions
         RLL(:,:,:,ith)=CZERO
         SLL(:,:,:,ith)=CZERO

         ! Right solutions
         TMAT0=CZERO
         ALPHA0=CZERO ! LLY
         ! faster calculation of RLL.
         ! no irregular solutions are needed in self-consistent iterations
         ! because the t-matrix depends only on RLL
         if( OPT('RLL-SLL ') .and. .not.(OPT('XCPL    ').or.OPT('OPERATOR')) ) then
            call rll_global_solutions(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),  &
               RLL(:,:,:,ith),TMAT0(:,:),NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO, &
               4*(LMAX+1),IRMDNEW,NSRA,JLK_INDEX,HLK(:,:,ith),         &
               JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),GMATPREFACTOR,'1',    &
               USE_SRATRICK,ALPHA0(:,:))
         else
            call RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),RLL(:,:,:,ith), &
               SLL(:,:,:,ith),TMAT0(:,:),NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO, &
               4*(LMAX+1),IRMDNEW,NSRA,JLK_INDEX,HLK(:,:,ith),         &
               JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),GMATPREFACTOR,'1','1',&
               '0',USE_SRATRICK,ALPHA0(:,:))
         end if

         if (NSRA.EQ.2) then
            RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
            SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
         endif
         if (test('BdG_dev ')) then
           write(filename, '(A,I0.3,A,I0.3,A)') 'rll_atom_',i1,'_energ_',ie,'.dat'
           open(888888, file=trim(filename), form='formatted')
           write(888888, '(A,I9,A,I9,A,I9)') '# dimension: lmmaxso*nvec=',nvec*lmmaxso,' lmmaxso=',lmmaxso,' irmdnew=', irmdnew 
           write(888888, '(2ES21.9)') rll(:,:,:,ith)
           close(888888)
           write(filename, '(A,I0.3,A,I0.3,A)') 'sll_atom_',i1,'_energ_',ie,'.dat'
           open(888888, file=trim(filename), form='formatted')
           write(888888, '(A,I9,A,I9,A,I9)') '# dimension: lmmaxso*nvec=',nvec*lmmaxso,' lmmaxso=',lmmaxso,' irmdnew=', irmdnew 
           write(888888, '(2ES21.9)') sll(:,:,:,ith)
           close(888888)
         end if

         ! add spherical contribution of tmatrix
         if (USE_SRATRICK.EQ.1) then
            do LM1=1,LMMAXSO
               TMAT0(LM1,LM1)=TMAT0(LM1,LM1)+TMATSPH(JLK_INDEX(LM1),ith)
            enddo
            do LM2=1,LMMAXSO
               do LM1=1,LMMAXSO
                  ALPHA0(LM1,LM2)=ALPHASPH(JLK_INDEX(LM1))*ALPHA0(LM1,LM2)        ! LLY
               enddo
            enddo
         endif
         do LM1=1,LMMAXSO
            do LM2=1,LMMAXSO
               TMATLL(LM1,LM2)=TMATLL(LM1,LM2)+TMAT0(LM1,LM2)
            enddo
         enddo
         do LM1=1,LMMAXSO
            do LM2=1,LMMAXSO
               ALPHALL(LM1,LM2)=ALPHALL(LM1,LM2)+ALPHA0(LM1,LM2)
            enddo
         enddo
         if (LLY.NE.0) then
            do LM1=1,LMMAXSO
               do LM2=1,LMMAXSO
                  DTMATLL(LM1,LM2)=DTMATLL(LM1,LM2)+real(SIGNDE, kind=dp)*TMAT0(LM1,LM2) ! LLY
                  DALPHALL(LM1,LM2)=DALPHALL(LM1,LM2)+real(SIGNDE, kind=dp)*ALPHA0(LM1,LM2) ! LLY
               enddo
            enddo
         endif
      enddo ! signde=-ideriv,ideriv,2 ! lly

      ! Average values of t-matrix and alpha at e+de and e-de
      do LM1=1,LMMAXSO
         do LM2=1,LMMAXSO
            TMATLL(LM1,LM2)=TMATLL(LM1,LM2)/real(1+IDERIV, kind=dp) ! LLY
            ALPHALL(LM1,LM2)=ALPHALL(LM1,LM2)/real(1+IDERIV, kind=dp) ! LLY
            if (LLY.NE.0) then
               ! Contruct derivative of t-matrix and alpha
               DTMATLL(LM1,LM2)=DTMATLL(LM1,LM2)/DELTAE ! LLY
               DALPHALL(LM1,LM2)=DALPHALL(LM1,LM2)/DELTAE ! LLY
            endif
         enddo
      enddo
      if (LLY.NE.0) then
         ! calculate Tr[alpha^-1*dalpha/de] for LLoyd's formula
         ALPHA0=CZERO ! LLY
         AUX=CZERO ! LLY
         call ZGEINV1(ALPHALL,ALPHA0,AUX,IPIV,LMMAXSO)
         call ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,ALPHA0,LMMAXSO,DALPHALL,&
            LMMAXSO,CZERO,AUX,LMMAXSO) ! LLY
         ! Trace of AUX
         TRALPHA=CZERO ! LLY
         do LM1=1,LMMAXSO
            TRALPHA=TRALPHA+AUX(LM1,LM1) ! LLY
         enddo
      endif ! LLY

      if(test('rhoqtest') .and. ie==2) then
         ! read in mu0 atom index
         open(9999,file='mu0')
         read(9999,*) mu0, nscoef
         close(9999)
      end if

      ! Calculate additional t-matrices for Jij-tensor calculation
      if (t_dtmatJij_at%calculate .or.( t_wavefunctions%isave_wavefun(i1, ie)>0 .and.   &
          (t_wavefunctions%save_rllleft .or.t_wavefunctions%save_sllleft) )             &
         .or. ((test('rhoqtest') .and. ie==2).and.(i1==mu0)) ) then                     !rhoqtest
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculate the left-hand side solution this needs to be done for the
         ! calculation of t-matrices for Jij tensor or if wavefunctions should be saved
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Contruct the spin-orbit coupling hamiltonian and add to potential
         call SPINORBIT_HAM(LMAX,lmsize,VINS,RNEW,ERYD,ZAT,CVLIGHT,SOCSCALE,NSPIN,  &
            LMPOT,THETA,PHI,IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW,  &
            NRMAXD,VNSPLL0(:,:,:),VNSPLL1(:,:,:,ith),'transpose')

         ! Extend matrix for the SRA treatment
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

         ! Calculate the source terms in the Lippmann-Schwinger equation
         ! these are spherical hankel and bessel functions
         HLK(:,:,ith)=CZERO
         JLK(:,:,ith)=CZERO
         HLK2(:,:,ith)=CZERO
         JLK2(:,:,ith)=CZERO
         GMATPREFACTOR=CZERO
         jlk_index = 0
         call RLLSLLSOURCETERMS(NSRA,NVEC,ERYD,RNEW,IRMDNEW,NRMAXD,LMAX,LMMAXSO,1,  &
            JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),GMATPREFACTOR)

         ! Using spherical potential as reference
         ! notice that exchange the order of left and right hankel/bessel functions
         if (USE_SRATRICK.EQ.1) then
            TMATSPH(:,ith)=CZERO
            call CALCSPH(NSRA,IRMDNEW,NRMAXD,LMAX,NSPIN,ZAT,ERYD,LMPOT, &
               LMMAXSO,RNEW,VINS,NCHEB,NPAN_TOT,RPAN_INTERVALL,JLK_INDEX,        &
               HLK2(:,:,ith),JLK2(:,:,ith),HLK(:,:,ith),JLK(:,:,ith),            &
               GMATPREFACTOR,ALPHASPH,TMATSPH(:,ith),USE_SRATRICK)
         endif

         ! Calculate the tmat and wavefunctions
         RLLLEFT(:,:,:,ith)=CZERO
         SLLLEFT(:,:,:,ith)=CZERO

         ! Left solutions
         ! notice that exchange the order of left and right hankel/bessel functions
         TMAT0=CZERO
         ALPHA0=CZERO ! LLY
         ! faster calculation of RLL.
         ! no left solutions are needed in self-consistent iterations
         ! because the t-matrix depends only on RLL
         if( OPT('RLL-SLL ') .and. .not.(OPT('XCPL    ').or.OPT('OPERATOR')) ) then
            ! do nothing
         else
            call RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),RLLLEFT(:,:,:,ith),      &
               SLLLEFT(:,:,:,ith),TMAT0,NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),&
               IRMDNEW,NSRA,JLK_INDEX,HLK2(:,:,ith),JLK2(:,:,ith),HLK(:,:,ith), &
               JLK(:,:,ith),GMATPREFACTOR,'1','1','0',USE_SRATRICK,ALPHA0)
         end if
         if (NSRA.EQ.2) then
            RLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=RLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
            SLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=SLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
         endif

         if(test('rhoqtest')) then
#ifdef CPP_OMP
          write(*,*) 'rhoqtest does not work in OMP version!!'
          write(*,*) 'please use hybrid compilation mode'
          stop
#else
          open(9999, file='params.txt')
          write(9999,*) lmmaxso, t_params%natyp
          write(9999,*) t_params%naez, t_params%nclsd, t_params%nr, t_params%nembd1-1, t_params%lmax
          write(9999,*) t_params%alat
          close(9999)

          open(9999, file='host.txt')
          write(9999,*) t_params%rbasis(1:3,1:t_params%natyp)
          write(9999,*) t_params%rcls(1:3,1:t_params%nclsd,1:t_params%nclsd), t_params%rr(1:3,0:t_params%nr), t_params%atom(1:t_params%nclsd,1:t_params%naez+t_params%nembd1-1)
          write(9999,*) t_params%cls(1:t_params%naez+t_params%nembd1-1),t_params%ezoa(1:t_params%nclsd,1:t_params%naez+t_params%nembd1-1), t_params%nacls(1:t_params%nclsd)
          close(9999)

          open(9999, file='wavefunctions.txt')
          write(9999,'(100I9)') ntotd, npan_tot, ncheb, nsra, irmdnew
          write(9999,'(1000E26.17)') rnew(1:irmdnew)
          do ir=1,irmdnew
            do lm1=1,nsra*lmmaxso
              do lm2=1,lmmaxso
                 write(9999,'(20000E16.7)') Rll(lm1, lm2, ir, ith),Rllleft(lm1, lm2, ir, ith)
              end do
            end do
          enddo
          do lm1=0,npan_tot
            write(9999,'(E16.7,I9)') rpan_intervall(lm1), ipan_intervall(lm1)
          enddo
          close(9999)
#endif
         end if ! test('rhoqtest')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculate the left-hand side solution
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end if !t_dtmatJij_at%calculate .or. t_wavefunctions%Nwfsavemax>0

      ! save_wavefuncions
      if(t_wavefunctions%Nwfsavemax>0) then
         ! here all four (left, right, regular and irregular) are stored, the memory demand cound be reduced by a factor 2 if only the right solution would be computed here and saved and the left solution would be calculated later in main1c
         call save_wavefunc(t_wavefunctions, rll, rllleft,sll, sllleft, i1, ie, &
            NSRA, LMMAXSO, IRMDNEW, ith)
      end if

      if(t_dtmatJij_at%calculate) then
         call calc_dtmatJij(lmsize,LMMAXSO,LMPOT,NTOTD,NRMAXD,NSRA,IRMDNEW,NSPIN,  &
            VINS,RLLLEFT(:,:,:,ith),RLL(:,:,:,ith),RPAN_INTERVALL,IPAN_INTERVALL,   &
            NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,NCLEB,RNEW,t_dtmatJij_at%dtmat_xyz(:,:,:,ie_num))

      end if!t_dtmatJij_at%calculate

      ! writeout
#ifdef CPP_OMP
      !$omp critical
#endif
      if (t_tgmat%tmat_to_file) then
         IREC = IE + IELAST*(I1-1)
#ifndef CPP_OMP        
         if(test('rhoqtest')) then
           
           if(ie_num==1.and.i1==1) then
             write(*,*)                      ! status bar
             write(*,*) 'rhoq: write-out t-mat', ie_end, t_params%natyp
             write(*, '("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
             write(*,FMT=190, advance='no') !beginning of statusbar
           endif
           
           if(myrank==master) then
             if(t_params%NATYP*ie_end>=50) then
               if(mod( I1+t_params%natyp*(ie_num-1), (t_params%NATYP*ie_end/50))==0 ) write(6,FMT=200, advance='no')
             else
               write(6,FMT=200, advance='no')
             end if
           end if

!          write(*,*) 'rotating with', theta, phi
!          lmGF0D= (LMAXD+1)**2
!          caLL ROTATEMATRIX(TMATLL,THETA,PHI,LMGF0D,0)

         end if    ! test('rhoqtest') 
#endif
         write(69,REC=IREC) TMATLL(:,:)
         ! human readable writeout if test option is hit
         if(test('fileverb')) then
            write(696969,'(i9,20000F15.7)') irec, TMATLL(:,:)
         end if
      else
#ifdef CPP_MPI
         irec = ie_num +ie_end*(I1-t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie)-1)
#else
         irec = ie_num + ie_end * (i1-1)
#endif
         t_tgmat%tmat(:,:,irec) = TMATLL(:,:)
      end if
      if (LLY.NE.0) then
         if(t_lloyd%dtmat_to_file) then
            IREC = IE + IELAST*(I1-1)
            write(691,REC=IREC) DTMATLL(:,:)    ! LLY
            if(test('fileverb')) then
               write(691691691,'(i9,20000F15.7)') irec, DTMATLL(:,:)
            end if
         else
            irec = ie_num + ie_end*(i1-1)
            t_lloyd%dtmat(:,:,irec) = DTMATLL(:,:)
         end if
         if(t_lloyd%tralpha_to_file) then
            IREC = IE + IELAST*(I1-1)
            write(692,REC=IREC) TRALPHA                              ! LLY
            if(test('fileverb')) then
               write(692692692,'(i9,20000F15.7)') irec, TRALPHA
            end if
         else
            irec = ie_num + ie_end*(i1-1)
            t_lloyd%tralpha(irec) = TRALPHA
         end if
      endif
#ifdef CPP_OMP
   !$omp end critical
#endif

#ifdef CPP_MPI
      !stop timing measurement for this pair of ie and i1, needed for MPIadapt
      if(MPIadapt) call timing_stop('time_1a_ieiatom',save_out=timings_1a(ie, i1))
#endif

   enddo ! IE loop
#ifdef CPP_OMP
      !$omp end parallel do
#endif

190     FORMAT('                 |')   ! status bar
200     FORMAT('|')                    ! status bar
        if(test('rhoqtest').and.i1==t_params%natyp.and.myrank==master) write(6,*) ! status bar
        ! finished kpts status bar

   i_all=-product(shape(AUX))*kind(AUX)
   deallocate(AUX,stat=i_stat)
   call memocc(i_stat,i_all,'AUX','tmat_newsolver')
   i_all=-product(shape(ALPHA0))*kind(ALPHA0)
   deallocate(ALPHA0,stat=i_stat)
   call memocc(i_stat,i_all,'ALPHA0','tmat_newsolver')
   i_all=-product(shape(ALPHALL))*kind(ALPHALL)
   deallocate(ALPHALL,stat=i_stat)
   call memocc(i_stat,i_all,'ALPHALL','tmat_newsolver')
   i_all=-product(shape(DALPHALL))*kind(DALPHALL)
   deallocate(DALPHALL,stat=i_stat)
   call memocc(i_stat,i_all,'DALPHALL','tmat_newsolver')
   i_all=-product(shape(TMAT0))*kind(TMAT0)
   deallocate(TMAT0,stat=i_stat)
   call memocc(i_stat,i_all,'TMAT0','tmat_newsolver')
   i_all=-product(shape(DTMATLL))*kind(DTMATLL)
   deallocate(DTMATLL,stat=i_stat)
   call memocc(i_stat,i_all,'DTMATLL','tmat_newsolver')
   i_all=-product(shape(TMATLL))*kind(TMATLL)
   deallocate(TMATLL,stat=i_stat)
   call memocc(i_stat,i_all,'TMATLL','tmat_newsolver')
   i_all=-product(shape(JLK_INDEX))*kind(JLK_INDEX)
   deallocate(JLK_INDEX,stat=i_stat)
   call memocc(i_stat,i_all,'JLK_INDEX','tmat_newsolver')
   i_all=-product(shape(VINS))*kind(VINS)
   deallocate(VINS,stat=i_stat)
   call memocc(i_stat,i_all,'VINS','tmat_newsolver')
   i_all=-product(shape(VNSPLL0))*kind(VNSPLL0)
   deallocate(VNSPLL0,stat=i_stat)
   call memocc(i_stat,i_all,'VNSPLL0','tmat_newsolver')
   i_all=-product(shape(VNSPLL1))*kind(VNSPLL1)
   deallocate(VNSPLL1,stat=i_stat)
   call memocc(i_stat,i_all,'VNSPLL1','tmat_newsolver')
   i_all=-product(shape(VNSPLL))*kind(VNSPLL)
   deallocate(VNSPLL,stat=i_stat)
   call memocc(i_stat,i_all,'VNSPLL','tmat_newsolver')
   i_all=-product(shape(HLK))*kind(HLK)
   deallocate(HLK,stat=i_stat)
   call memocc(i_stat,i_all,'HLK','tmat_newsolver')
   i_all=-product(shape(JLK))*kind(JLK)
   deallocate(JLK,stat=i_stat)
   call memocc(i_stat,i_all,'JLK','tmat_newsolver')
   i_all=-product(shape(HLK2))*kind(HLK2)
   deallocate(HLK2,stat=i_stat)
   call memocc(i_stat,i_all,'HLK2','tmat_newsolver')
   i_all=-product(shape(JLK2))*kind(JLK2)
   deallocate(JLK2,stat=i_stat)
   call memocc(i_stat,i_all,'JLK2','tmat_newsolver')
   i_all=-product(shape(TMATSPH))*kind(TMATSPH)
   deallocate(TMATSPH,stat=i_stat)
   call memocc(i_stat,i_all,'TMATSPH','tmat_newsolver')
   i_all=-product(shape(RLL))*kind(RLL)
   deallocate(RLL,stat=i_stat)
   call memocc(i_stat,i_all,'RLL','tmat_newsolver')
   i_all=-product(shape(SLL))*kind(SLL)
   deallocate(SLL,stat=i_stat)
   call memocc(i_stat,i_all,'SLL','tmat_newsolver')
   i_all=-product(shape(RLLLEFT))*kind(RLLLEFT)
   deallocate(RLLLEFT,stat=i_stat)
   call memocc(i_stat,i_all,'RLLLEFT','tmat_newsolver')
   i_all=-product(shape(SLLLEFT))*kind(SLLLEFT)
   deallocate(SLLLEFT,stat=i_stat)
   call memocc(i_stat,i_all,'SLLLEFT','tmat_newsolver')

end subroutine TMAT_NEWSOLVER
