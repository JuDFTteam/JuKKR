!-------------------------------------------------------------------------------
!> @brief Performs k-space integration, determines scattering path operator
!> \f$\tau = \left(g\left(\mathbf{k},e\right)-t^{-1}\right)^{-1}\f$
!>  and Greens function of the real system -> \f$GS\f$
!
!> @details Modifications according to H. Hoehler ( July 2002)
!> define Fourier transformation as
!>
!> \f$ G\left(\mu,\mu'\right)_{L,L'}=\frac{1}{2}\left[\sum_n G^{n,0}\left(\mu,\mu'\right)_{L,L'}\exp\left(-iKR^n\right)+ \sum_n G^{n,0}\left(\mu,\mu'\right)_{L,L'}\exp\left(-iK\left(-R\right)^n\right)\right]  \f$
!>
!> this operation has to be done to satisfy the point symmetry;
!> the application of the Fourier transformation is just an
!> approximation for the tb system, since the translational invariance
!> is not satisfied --> force it by R, -R
!
!> @note
!> - New version 10.99: up -> left , down -> right, for decimation
!> - Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine KKRMAT01(NR,LMAX,NREF,LMGF0D,LMMAXD,BZKP,NOFKS,GS,VOLCUB,TINVLL,RROT, &
   NSHELL,NSDIA,ALAT,NSYMAT,                          &
   NAEZ,CLS,NACLS,NACLSMAX,RR,EZOA,ATOM,              &
   NSH1,NSH2,GINP,RBASIS,RCLS,                        &
   TINVBUP,TINVBDOWN,VACFLAG,NLBASIS,NRBASIS,         &
   FACTL,ICHECK,INVMOD,IDECI,SRREL,IRREL,NRREL,       &
   DTREFLL,DTMATLL,DGINP,REFPOT,LLY_GRTR,TRACET,CFCTOR,LLY)   ! LLY
#ifdef CPP_MPI
   use mpi

   use mod_types, only: t_mpi_c_grid
   use mod_mympi, only: myrank, nranks, master,distribute_linear_on_tasks
#else
   use mod_mympi, only: myrank, nranks, master
#endif
   use mod_types, only: t_inc
#ifdef CPP_TIMING
   use mod_timing
#endif
#ifdef CPP_HYBRID
      use omp_lib
#endif
      use mod_rhoqtools, only: rhoq_find_kmask, rhoq_saveG, rhoq_write_tau0, rhoq_read_mu0_scoef

   use global_variables
   use Constants
   use Profiling
      Use mod_datatypes, Only: dp

   implicit none
   ! .. Input variables
   integer, intent(in) :: NR        !< Number of real space vectors rr
   integer, intent(in) :: LLY       !< LLY <> 0 --> use Lloyds formula
   integer, intent(in) :: NREF      !< Number of diff. ref. potentials
   integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: NOFKS
   integer, intent(in) :: NSDIA
   integer, intent(in) :: IDECI
   integer, intent(in) :: LMGF0D    !< (LMAX+1)**2
   integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
   integer, intent(in) :: NSHELL    !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
   integer, intent(in) :: NSYMAT
   integer, intent(in) :: INVMOD    !< Inversion scheme
   integer, intent(in) :: NLBASIS   !< Number of basis layers of left host (repeated units)
   integer, intent(in) :: NRBASIS   !< Number of basis layers of right host (repeated units)
   integer, intent(in) :: NACLSMAX
   real (kind=dp), intent(in) :: ALAT         !< Lattice constant in a.u.
   integer, dimension(*), intent(in)                           :: CLS     !< Cluster around atomic sites
   integer, dimension(*), intent(in)                           :: NSH1    !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
   integer, dimension(*), intent(in)                           :: NSH2    !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
   integer, dimension(*), intent(in)                           :: NACLS   !< Number of atoms in cluster
   integer, dimension(*), intent(in)                           :: REFPOT  !< Ref. pot. card  at position ! REFPOT(NAEZD+NEMBD)
   integer, dimension(NACLSD,*), intent(in)                    :: ATOM   !< Atom at site in cluster
   integer, dimension(NACLSD,*), intent(in)                    :: EZOA   !< EZ of atom at site in cluster
   integer, dimension(2,LMMAXD), intent(in)                    :: NRREL
   integer, dimension(NAEZ/NPRINCD,NAEZ/NPRINCD), intent(in)   :: ICHECK
   integer, dimension(2,2,LMMAXD), intent(in)                  :: IRREL
   real (kind=dp), dimension(*), intent(in)            :: VOLCUB
   real (kind=dp), dimension(3,0:NR), intent(in)       :: RR       !< Set of real space vectors (in a.u.)
   real (kind=dp), dimension(3,*), intent(in)          :: BZKP
   real (kind=dp), dimension(3,*), intent(in)          :: RBASIS   !< Position of atoms in the unit cell in units of bravais vectors
   real (kind=dp), dimension(48,3,*), intent(in)       :: RROT
   real (kind=dp), dimension(3,NACLSD,*), intent(in)   :: RCLS  !< Real space position of atom in cluster
   complex (kind=dp), dimension(LMMAXD,LMMAXD), intent(in)              :: FACTL
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NAEZ), intent(in)         :: TINVLL
   complex (kind=dp), dimension(LMMAXD,LMMAXD,*), intent(in)            :: TINVBUP
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NREF), intent(in)         :: DTREFLL ! LLY dtref/dE
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NAEZ), intent(in)         :: DTMATLL ! LLY  dt/dE (should be av.-tmatrix in CPA)
   complex (kind=dp), dimension(LMMAXD,LMMAXD,*), intent(in)            :: TINVBDOWN
   complex (kind=dp), dimension(LMGF0D*NACLSMAX,LMGF0D,*), intent(in)   :: GINP  ! Gref
   complex (kind=dp), dimension(LMGF0D*NACLSMAX,LMGF0D,*), intent(in)   :: DGINP ! LLY dGref/dE
   complex (kind=dp), dimension(2,2,LMMAXD), intent(in)                 :: SRREL
   logical, dimension(2), intent(in) :: VACFLAG
   ! .. In/Out variables
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NSYMAXD,*), intent(inout) :: GS
   ! .. Output variables
   complex (kind=dp), intent(out) :: LLY_GRTR ! Trace Eq.5.38 PhD Thiess  (integrated) ! LLY Lloyd
   ! .. Local variables
   integer :: ALM
   integer :: NDIM
   integer :: ALMGF0
   integer :: i_stat, i_all
   integer :: IKM1,IKM2,IS,N1,N2,J1,J2,I2
   integer :: IQ1,IQ2,IOFF1,IOFF2,JOFF1,JOFF2
   integer :: I,I1,ILM,ISYM,IU,J,JLM,IL1,KPT,LM,LM1,LM2,NS,IL2,JL1,JL2
   real (kind=dp) :: ZKTR
   complex (kind=dp) :: LLY_GRTR_K ! Trace Eq.5.38 PhD Thiess  (k-dependent) ! LLY Lloyd
   complex (kind=dp) :: CSUM1,CSUM2,TRACE,TRACET  ! LLY Lloyd
   complex (kind=dp) :: CARG,CITPI,CFCTOR
   real (kind=dp), dimension(3) :: KP
   real (kind=dp), dimension(6) :: BZKPK
   real (kind=dp), dimension(3,0:NR) :: RRM
   complex (kind=dp), dimension(LMMAXD,LMMAXD) :: G
   complex (kind=dp), dimension(LMMAXD,LMMAXD) :: GAUX1 ! LLY
   complex (kind=dp), dimension(LMMAXD,LMMAXD) :: GAUX2 ! LLY
   complex (kind=dp), dimension(LMMAXD,LMMAXD) :: GAUX3 ! LLY
   complex (kind=dp), dimension(NSYMAXD,NSHELD) :: ETAIKR
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NAEZ) :: T_AUX   ! LLY auxiliary array for t-matrix manipulation
   ! .. Local allocatable arrays
   complex (kind=dp), dimension(:,:), allocatable :: GLLKE,GLLKEM,GLLKEN
   complex (kind=dp), dimension(:,:), allocatable :: DGLLKE,DGLLKEM,DGLLKEN,GREFLLKE ! LLY
   complex (kind=dp), dimension(:,:), allocatable :: GLLKE0V,GLLKE0V2,GLLKETV ! for VIRTUAL ATOMS
   complex (kind=dp), dimension(:,:), allocatable :: GLLKETV_new ! for VIRTUAL ATOMS
   complex (kind=dp), dimension(:,:), allocatable :: GLLKE0,GLLKE0M
   ! .. Parameters
   complex (kind=dp) :: CMI
   parameter (CMI=(0D0,-1D0))

#ifdef CPP_MPI
   integer :: ntot1
   integer, dimension(0:nranks-1) :: ntot_pT, ioff_pT
#endif
   integer :: k_start, k_end
   ! ..
   logical :: TEST,OPT
   ! .. External subroutines ..
   external :: CINIT,DLKE0,OPT,TEST,GTDYSON
   ! .. Intrinsic functions ..
   intrinsic :: ATAN,EXP
   !     ..
!#ifdef CPP_MPI
   complex (kind=dp), dimension(LMMAXD,LMMAXD,NSYMAXD) :: WORK
   integer :: IERR,IWORK
!#endif
   integer :: mu, nscoef, imin, ie
   integer, allocatable :: iatomimp(:)
   
   real (kind=dp), allocatable :: rhoq_kmask(:,:) ! only in reduced number of kpts
   integer, allocatable :: kmask(:) ! logical array over all kpts (determine if kpt=1,nofks is in reduced set)
   integer :: mythread

   !      NDIM=LMGF0D*NAEZ
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Array sizes definitions
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   NDIM   = LMMAXD*NAEZ
   ALM    = NAEZ*LMMAXD
   ALMGF0 = NAEZ*LMGF0D

   if ( TEST('flow     ') .and. (t_inc%i_write>0)) write(1337,*) '>>> kkrmat1: loop over k-points'
   !
   CITPI = CMI*8.D0*ATAN(1.D0)    ! = -i*2*PI
   !
   do NS = 1,NSHELL
      do IU = 1,NSYMAXD
         call CINIT(LMMAXD*LMMAXD,GS(1,1,IU,NS))
      end do
   end do

   LLY_GRTR = CZERO ! LLY Lloyd
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Array allocations
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   allocate(GLLKE(ALM,ALM),stat=i_stat)
   call memocc(i_stat,product(shape(GLLKE))*kind(GLLKE),'GLLKE','kkrmat01')
   ! LLY Lloyd
   if (LLY.NE.0) then
      allocate(DGLLKE(ALM,ALM),stat=i_stat)
      call memocc(i_stat,product(shape(DGLLKE))*kind(DGLLKE),'DGLLKE','kkrmat01')
      allocate(GREFLLKE(ALM,ALM),stat=i_stat)
      call memocc(i_stat,product(shape(GREFLLKE))*kind(GREFLLKE),'GREFLLKE','kkrmat01')
   endif

   if ( OPT('VIRATOMS') ) then
      allocate(GLLKE0V(ALM,ALM),stat=i_stat)
      call memocc(i_stat,product(shape(GLLKE0V))*kind(GLLKE0V),'GLLKE0V','kkrmat01')
      allocate(GLLKE0V2(ALM,ALM),stat=i_stat)
      call memocc(i_stat,product(shape(GLLKE0V2))*kind(GLLKE0V2),'GLLKE0V2','kkrmat01')
      allocate(GLLKETV(ALM,LMMAXD),stat=i_stat)
      call memocc(i_stat,product(shape(GLLKETV))*kind(GLLKETV),'GLLKETV','kkrmat01')
      allocate(GLLKETV_new(LMMAXD,ALM),stat=i_stat)
      call memocc(i_stat,product(shape(GLLKETV_new))*kind(GLLKETV_new),'GLLKETV_new','kkrmat01')
   end if !( OPT('VIRATOMS') ) THEN
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! End of array allocations
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  K-points loop
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(test('rhoqtest')) then
         call rhoq_read_mu0_scoef(iatomimp, mu, nscoef, imin)
         call rhoq_find_kmask(nofks, k_end, bzkp(1:3,1:nofks), kmask, rhoq_kmask)
   end if !test('rhoqtest')

#ifdef CPP_MPI
   ! MPI:
   if(.not.opt('qdos    ')) then
   
     if(test('rhoqtest')) then
       ntot1 = k_end
     else
       ntot1 = NOFKS
     endif

     if(myrank==master) write(1337,*) 'kkrmat k loop:', NOFKS, t_mpi_c_grid%nranks_ie
     call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,master,ntot1,ntot_pT,ioff_pT,.true.)

     k_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
     k_end   = ioff_pT(t_mpi_c_grid%myrank_ie)+ntot_pT(t_mpi_c_grid%myrank_ie)
     t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie)
     t_mpi_c_grid%ntot_pT1 = ntot_pT
     t_mpi_c_grid%ioff_pT1 = ioff_pT
      
   else !.not.opt('qdos    ')
      
     k_start = 1
     k_end = NOFKS
      
   end if !.not.opt('qdos    ')
#else
   k_start = 1
   if(.not.test('rhoqtest')) k_end = NOFKS 
#endif

   ! k-loop not needed for GREENIMP-case
   if(opt('GREENIMP')) then
      if(myrank==master) write(*,*) 'Skipping kloop in kkrmat'
      k_start = 1
      k_end = 0
   end if

   ! Print header of statusbar for k-loop
   if(t_inc%i_write>0) then
      write(1337,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') &
      0, 20, 40, 60, 80, 100
      write(1337,FMT=190, advance='no') ! Beginning of statusbar
   endif


#ifdef CPP_HYBRID
   !$omp parallel default(shared)
   !$omp& private(kpt, ns, i, j, isym, carg, i1, zktr, i2, iq1, iq2, ioff1)
   !$omp& private(ioff2, joff1, joff2, ikm1, ikm2, csum1, is, n1, n2)
   !$omp& private(j1, csum2, j2, il1, il2, lm1, lm2, gaux1, gaux2)
   !$omp& private(jl1, jl2, gaux3, ilm, jlm, mythread )
   !$omp& reduction(+:trace)
   mythread = omp_get_thread_num()
#else
   mythread = 0
#endif

   ! kpts loop
   do KPT = k_start,k_end
      GLLKE(:,:) = CZERO
      if (LLY.NE.0) DGLLKE(:,:) = CZERO
      KP(1:3) = BZKP(1:3,KPT)

      ! overwrite kpt in case of rhoqtest (take only reduced set of kpts)
      if(test('rhoqtest')) kp(1:3) = rhoq_kmask(1:3,kpt)

      ETAIKR(1:NSYMAT,1:NSHELL) = VOLCUB(KPT)
      !-------------------------------------------------------------------------
      ! First NAEZ/NATYP elements of GS() are site-diagonal
      !-------------------------------------------------------------------------
#ifdef CPP_HYBRID
      !$omp do
#endif
      do NS = NSDIA+1,NSHELL
         I = NSH1(NS)
         J = NSH2(NS)
         do ISYM  = 1,NSYMAT
            CARG = CZERO
            do I1 = 1,3
               ZKTR = RROT(ISYM,I1,NS) - RBASIS(I1,J) + RBASIS(I1,I)
               ZKTR = KP(I1)*ZKTR
               CARG =  CARG + ZKTR
            end do
            ETAIKR(ISYM,NS) = ETAIKR(ISYM,NS) * EXP(CARG*CITPI)
         end do
      end do
#ifdef CPP_HYBRID
      !$omp end do
#endif

      BZKPK(1:3) = KP(1:3)
      BZKPK(4:6) = 0.D0

      !-------------------------------------------------------------------------
      ! Fourier transformation
      !-------------------------------------------------------------------------
#ifdef CPP_TIMING
      if(mythread==0 .and. t_inc%i_time>0) call timing_start('main1b - fourier')
#endif

      RRM(1:3,1:NR) = -RR(1:3,1:NR)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  KREL .EQ. 0/1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (KREL.EQ.0) then

         !----------------------------------------------------------------------
         if(mythread==0) then
           allocate(GLLKEN(ALMGF0,ALMGF0),stat=i_stat)
           call memocc(i_stat,product(shape(GLLKEN))*kind(GLLKEN),'GLLKEN','kkrmat01')
           GLLKEN(:,:) = CZERO
         end if
         !----------------------------------------------------------------------
         call DLKE0(GLLKEN,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RR,EZOA,ATOM,BZKPK,RCLS,GINP)
         !----------------------------------------------------------------------
         if(mythread==0) then
           allocate(GLLKEM(ALMGF0,ALMGF0),stat=i_stat)
           call memocc(i_stat,product(shape(GLLKEM))*kind(GLLKEM),'GLLKEM','kkrmat01')
           GLLKEM(:,:) = CZERO
         end if
         !----------------------------------------------------------------------
         call DLKE0(GLLKEM,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RRM,EZOA,ATOM,BZKPK,RCLS,GINP)
         !----------------------------------------------------------------------
         ! LLY Lloyd
         ! Fourier for dGref/dE for Lloyds formula, repeat the above allocation
         ! and Fourier transform for the derivatives.
         !----------------------------------------------------------------------
         IF (LLY.NE.0) THEN
            !-------------------------------------------------------------------
            if(mythread==0) then
              allocate(DGLLKEN(ALMGF0,ALMGF0),stat=i_stat)
              call memocc(i_stat,product(shape(DGLLKEN))*kind(DGLLKEN),'DGLLKEN','kkrmat01')
              DGLLKEN(:,:) = CZERO
            end if
            !-------------------------------------------------------------------
            call DLKE0(DGLLKEN,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RR,EZOA,ATOM,BZKPK,RCLS,DGINP)
            !-------------------------------------------------------------------
            if(mythread==0) then
              allocate(DGLLKEM(ALMGF0,ALMGF0),stat=i_stat)
              call memocc(i_stat,product(shape(DGLLKEM))*kind(DGLLKEM),'DGLLKEM','kkrmat01')
              DGLLKEM(:,:) = CZERO
            end if
            !-------------------------------------------------------------------
            call DLKE0(DGLLKEM,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RRM,EZOA,ATOM,BZKPK,RCLS,DGINP)
         endif
         !----------------------------------------------------------------------
         ! LLY Lloyd
         !--------------------------------------------------------------------
         if (.NOT.OPT('NEWSOSOL')) then
            do I2=1,ALM
               do I1=1,ALM
                  GLLKE(I1,I2)= (GLLKEN(I1,I2) + GLLKEM(I2,I1))*0.5D0
                  if (LLY.NE.0)  then
                     DGLLKE(I1,I2) =( DGLLKEN(I1,I2) + DGLLKEM(I2,I1) )*0.5D0 ! LLY Lloyd
                  endif
               enddo
            enddo
         else                ! (.NOT.OPT('NEWSOSOL'))
            do I2=1,ALMGF0
               do I1=1,ALMGF0
                  GLLKEN(I1,I2)=(GLLKEN(I1,I2) + GLLKEM(I2,I1))*0.5D0
                  if (LLY.NE.0)  then
                     DGLLKEN(I1,I2) =( DGLLKEN(I1,I2) + DGLLKEM(I2,I1) )*0.5D0      ! LLY Lloyd
                  endif
               enddo
            enddo
            ! bigger GLLKE matrix and rearrange with atom block
#ifdef CPP_HYBRID
            !$omp do
#endif
            do IQ1=1,NAEZ
               do IQ2=1,NAEZ
                  IOFF1 = LMMAXD*(IQ1-1)
                  JOFF1 = LMGF0D*(IQ1-1)
                  IOFF2 = LMMAXD*(IQ2-1)
                  JOFF2 = LMGF0D*(IQ2-1)
                  do LM1=1,LMGF0D
                     do LM2=1,LMGF0D
                        GLLKE(IOFF1+LM1,IOFF2+LM2) =                    &
                           GLLKEN(JOFF1+LM1,JOFF2+LM2)
                        GLLKE(IOFF1+LM1+LMGF0D,IOFF2+LM2+LMGF0D) =      &
                           GLLKEN(JOFF1+LM1,JOFF2+LM2)
                        if (LLY.NE.0) then                                       ! LLY Lloyd
                           DGLLKE(IOFF1+LM1,IOFF2+LM2) =                &        ! LLY Lloyd
                              DGLLKEN(JOFF1+LM1,JOFF2+LM2)                       ! LLY Lloyd
                           DGLLKE(IOFF1+LM1+LMGF0D,IOFF2+LM2+LMGF0D)=   &        ! LLY Lloyd
                              DGLLKEN(JOFF1+LM1,JOFF2+LM2)                       ! LLY Lloyd
                        endif ! (LLY.NE.0)                                       ! LLY Lloyd
                     enddo
                  enddo
               enddo
            enddo
#ifdef CPP_HYBRID
            !$omp end do
#endif
         endif               ! (.NOT.OPT('NEWSOSOL'))
         !----------------------------------------------------------------------
         if(mythread==0) then
           i_all=-product(shape(GLLKEM))*kind(GLLKEM)
           deallocate(GLLKEM, stat=i_stat)
           call memocc(i_stat,i_all,'GLLKEM','kkrmat01')
           i_all=-product(shape(GLLKEM))*kind(GLLKEM)
           deallocate(GLLKEN, stat=i_stat)
           call memocc(i_stat,i_all,'GLLKEN','kkrmat01')
           if (LLY.NE.0) then
              i_all=-product(shape(DGLLKEM))*kind(DGLLKEM)
              deallocate(DGLLKEM, stat=i_stat)
              call memocc(i_stat,i_all,'DGLLKEM','kkrmat01')
              i_all=-product(shape(DGLLKEN))*kind(DGLLKEN)
              deallocate(DGLLKEN, stat=i_stat)
              call memocc(i_stat,i_all,'DGLLKEN','kkrmat01')
           endif
         end if
         !----------------------------------------------------------------------
         ! LLY Lloyd At this point DGLLKE contains the Fourier transform of the dGref/dE
         !----------------------------------------------------------------------
      else                   !  (KREL.EQ.0)
         !----------------------------------------------------------------------
         ! LLY Lloyd Not implementing Lloyds formula for KREL=1 (Dirac ASA)
         !----------------------------------------------------------------------
         if(mythread==0) then
           allocate(GLLKE0(ALMGF0,ALMGF0),stat=i_stat)
           call memocc(i_stat,product(shape(GLLKE0))*kind(GLLKE0),'GLLKE0','kkrmat01')
           allocate(GLLKE0M(ALMGF0,ALMGF0),stat=i_stat)
           call memocc(i_stat,product(shape(GLLKE0M))*kind(GLLKE0M),'GLLKE0M','kkrmat01')
         end if
         !----------------------------------------------------------------------
         call DLKE0(GLLKE0,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RR,EZOA,ATOM,BZKPK,RCLS,GINP)
         call DLKE0(GLLKE0M,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RRM,EZOA,ATOM,BZKPK,RCLS,GINP)
         !
         do I2=1,ALMGF0
            do I1=1,ALMGF0
               GLLKE0(I1,I2)=(GLLKE0(I1,I2) + GLLKE0M(I2,I1))*0.5D0
            enddo
         enddo
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Double the GLLKE0 matrix and transform to the REL representation
         !    ==> GLLKE
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CPP_HYBRID
         !$omp do
#endif
         do IQ1=1,NAEZ
            do IQ2=1,NAEZ
               IOFF1 = LMMAXD*(IQ1-1)
               JOFF1 = LMGF0D*(IQ1-1)
               IOFF2 = LMMAXD*(IQ2-1)
               JOFF2 = LMGF0D*(IQ2-1)
               ! ---------------------------------------------------------------
               do IKM2 = 1,LMMAXD
                  do IKM1 = 1,LMMAXD
                     CSUM1 = CZERO
                     do IS = 1,2
                        N1 = NRREL(IS,IKM1)
                        N2 = NRREL(IS,IKM2)
                        do I1 = 1,N1
                           J1 = IRREL(I1,IS,IKM1) + JOFF1
                           CSUM2 = CZERO
                           do I2 = 1,N2
                              J2 = IRREL(I2,IS,IKM2) + JOFF2
                              CSUM2 = CSUM2 +GLLKE0(J1,J2)*SRREL(I2,IS,IKM2)
                           end do
                           CSUM1 = CSUM1 + CONJG(SRREL(I1,IS,IKM1))*CSUM2
                        end do
                     end do
                     GLLKE(IOFF1+IKM1,IOFF2+IKM2) = CSUM1
                  end do
               end do
               !----------------------------------------------------------------
            end do
         end do
#ifdef CPP_HYBRID
         !$omp end do
#endif
         !-----------------------------------------------------------------------
         if(mythread==0) then
           i_all=-product(shape(GLLKE0))*kind(GLLKE0)
           deallocate(GLLKE0, stat=i_stat)
           call memocc(i_stat,i_all,'GLLKE0','kkrmat01')
           i_all=-product(shape(GLLKE0))*kind(GLLKE0)
           deallocate(GLLKE0M, stat=i_stat)
           call memocc(i_stat,i_all,'GLLKE0M','kkrmat01')
         end if
         !----------------------------------------------------------------------
      end if !  (KREL.EQ.0)
#ifdef CPP_TIMING
      if(mythread==0 .and. t_inc%i_time>0) call timing_pause('main1b - fourier')
#endif
      !
      if (LLY.NE.0) GREFLLKE(1:ALM,1:ALM) = GLLKE(1:ALM,1:ALM) ! LLY Save k-dependent Gref
      !
      if ( IDECI.EQ.1 ) then
         call DECIMATE(GLLKE,NAEZ,TINVBUP,TINVBDOWN,VACFLAG,FACTL,NLBASIS,NRBASIS,&
         ALM,NDIM,LMMAXD)
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Construct the matrix M=[-(t)^-1 + G^r] and store it
      ! in the same matrix GLLKE where G^r was stored.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( .not. OPT('VIRATOMS') ) then
         do I1=1,NAEZ
            do LM1 = 1,LMMAXD
               do LM2 = 1,LMMAXD
                  IL1 = LMMAXD*(I1-1)+LM1
                  IL2 = LMMAXD*(I1-1)+LM2
                  GLLKE(IL1,IL2)= GLLKE(IL1,IL2) - TINVLL(LM1,LM2,I1)
               enddo
            enddo
         enddo
         !
#ifdef CPP_TIMING
         if(mythread==0 .and. t_inc%i_time>0) call timing_start('main1b - inversion')
#endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Perform the inversion of matrix M
         ! the output is the scattering path operator TAU stored in GLLKE
         ! Actually -TAU, because TAU = (Deltat^-1 - Gref)^-1
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (LLY.NE.0) then ! If LLY, full inversion is needed
            call INVERSION(GLLKE,0,ICHECK) ! LLY
         else
            call INVERSION(GLLKE,INVMOD,ICHECK)
         endif
#ifdef CPP_TIMING
         if(mythread==0 .and. t_inc%i_time>0) call timing_pause('main1b - inversion')
#endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! LLY Lloyd
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (LLY.NE.0) then
            !-------------------------------------------------------------------
            ! LLY  Prepare quantities for Lloyds formula.
            ! LLY  Needed is Trace[ (1-Gref * Deltat)^-1 * d(1-Gref * Deltat)/dE ] (PhD Thiess Eq.5.38)
            ! LLY  where Deltat = t-tref. This is re-written as:
            ! LLY  -Trace[ Tau * ( dGref/dE + Gref * (dt/dE - dtref/dE) Deltat^-1 ) ]
            ! LLY  where Tau is the scattering path operator Tau = (Deltat^-1 - Gref)^-1
            ! LLY  (negative of array GLLKE) and (t-tref)^-1 is in array TINVLL.
            ! LLY  The quantities Gref, dGref/dE, dt/dE have been prepared by main1a.
            ! LLY  Quantity dtref/dE is in array DTREFLL
            !-------------------------------------------------------------------
            !
            ! First set up (dt/dE - dtref/dE) Deltat^-1, store in array t_aux
            do I1 = 1,NAEZ
               ! GAUX1 = dt/dE-dtref/dE
               GAUX1(1:LMMAXD,1:LMMAXD) = (1.D0/CFCTOR) *   &
                  (  DTMATLL(1:LMMAXD,1:LMMAXD,I1) -        &
                     DTREFLL(1:LMMAXD,1:LMMAXD,REFPOT(I1)) )
               GAUX2(1:LMMAXD,1:LMMAXD) =TINVLL(1:LMMAXD,1:LMMAXD,I1)
               ! T_AUX = (dt/dE-dtref/dE)* Deltat^-1
               call ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,   &
                  GAUX1,LMMAXD,GAUX2,LMMAXD,CZERO,GAUX3,LMMAXD)
               T_AUX(1:LMMAXD,1:LMMAXD,I1) = GAUX3(1:LMMAXD,1:LMMAXD)
            enddo
            !-------------------------------------------------------------------
            ! Now perform dGref/dE + Gref * t_aux
            ! (Gref is ALM*ALM ; t_aux site-diagonal LMMAXD*LMMAXD)
            !-------------------------------------------------------------------
            do J1 = 1,NAEZ               ! Loop over columns of Gref
               JL1 = LMMAXD*(J1-1) + 1
               JL2 = LMMAXD*(J1-1) + LMMAXD
               GAUX3(1:LMMAXD,1:LMMAXD) = T_AUX(1:LMMAXD,1:LMMAXD,J1)
               do I1 = 1,NAEZ            ! Loop over rows of Gref
                  IL1 = LMMAXD*(I1-1) + 1
                  IL2 = LMMAXD*(I1-1) + LMMAXD
                  ! Copy to small matrices
                  GAUX1(1:LMMAXD,1:LMMAXD) =GREFLLKE(IL1:IL2,JL1:JL2)
                  GAUX2(1:LMMAXD,1:LMMAXD) = DGLLKE(IL1:IL2,JL1:JL2)
                  ! GAUX2 = GAUX2 + GAUX1 * T_AUX
                  call ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,GAUX1,   &
                     LMMAXD,GAUX3,LMMAXD,CONE,GAUX2,LMMAXD)
                  ! Copy back to large matrix, use again array DGLLKE
                  ! (the I1-J1 block is not needed any more)
                  DGLLKE(IL1:IL2,JL1:JL2) = GAUX2(1:LMMAXD,1:LMMAXD)
               enddo
            enddo

            ! full matrix multiple
            !            ALLOCATE(GLLKE0(ALM,ALM))
            !            GLLKE0=CZERO
            !            DO I1=1,NAEZ
            !               DO LM1 = 1,LMMAXD
            !                  DO LM2 = 1,LMMAXD
            !                     IL1 = LMMAXD*(I1-1)+LM1
            !                     IL2 = LMMAXD*(I1-1)+LM2
            !                     GLLKE0(IL1,IL2)= T_AUX(LM1,LM2,I1)
            !                  ENDDO
            !               ENDDO
            !            ENDDO
            !            CALL ZGEMM('N','N',ALM,ALM,ALM,CONE,GREFLLKE,
            !     &                 ALM,GLLKE0,ALM,CONE,DGLLKE,ALM)
            !            DEALLOCATE(GLLKE0)

            ! Now array DGLLKE contains
            ! ( dGref/dE + Gref * (dt/dE - dtref/dE) Deltat^-1 )
            ! Build trace of tau * DGLLKE, -tau is conained in GLLKE.
            TRACE = CZERO
            do IL1 = 1,ALM
               do IL2 = 1,ALM
                  TRACE = TRACE + GLLKE(IL1,IL2) * DGLLKE(IL2,IL1)
               enddo
            enddo
            LLY_GRTR_K = TRACE
         endif ! (LLY.NE.0)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! LLY Lloyd
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else                   !  .not. OPT('VIRATOMS')
         ! LLY Lloyd formula not built in yet for viratoms
         GLLKE0V(1:ALM,1:ALM) = GLLKE(1:ALM,1:ALM)
         !
         do I1 = 1,NAEZ
            IL1 = (I1-1)*LMMAXD + 1
            ! GLLKETV = -GLLKE0V * TINVLL,
            ! where TINVLL contains (t-tref) and not 1/(t-tref) in case of opt VIRATOMS
            ! tref=0 for each vir. atom.
            call ZGEMM('N','N',NDIM,LMMAXD,LMMAXD,-CONE, &
               GLLKE0V(1,IL1),ALM,TINVLL(1,1,I1),LMGF0D, &
               CZERO,GLLKETV(1,1),ALM)
            call ZCOPY(ALM*LMMAXD,GLLKETV(1,1),1,GLLKE0V2(1,IL1),1)
         end do
         !----------------------------------------------------------------------
         ! Solve (1-gt)G=g instead of [Gref - t^-1]^-1 for viratoms
         ! because for a virtual atom t=0, t^-1 undefined.
         !----------------------------------------------------------------------
         call GTDYSON(GLLKE0V2,GLLKE,NDIM,ALM,ALM)
      end if                 !  .not. OPT('VIRATOMS')
      !-------------------------------------------------------------------------
      ! Global sum on array gs
      !-------------------------------------------------------------------------
      ! no omp at this loop because of rhoq output
      do NS = 1,NSHELL
         I = NSH1(NS)
         J = NSH2(NS)
         ILM = LMMAXD*(I-1) + 1
         JLM = LMMAXD*(J-1)

         do LM = 1,LMMAXD
            call ZCOPY(LMMAXD,GLLKE(ILM,JLM+LM),1,G(1,LM),1)
         end do
         !----------------------------------------------------------------------
         do ISYM = 1,NSYMAT
            do LM2=1,LMMAXD
               do LM1=1,LMMAXD
                  GS(LM1,LM2,ISYM,NS) = GS(LM1,LM2,ISYM,NS)+ ETAIKR(ISYM,NS) * G(LM1,LM2)
               end do
            end do
            if(test('rhoqtest')) then
               call rhoq_saveG(nscoef,rhoq_kmask,kpt,nofks,k_end,kp,i,j,mu,imin,iatomimp,lmmaxd,G)
            end if
         end do ! isym
         ! ----------------------------------------------------------------------
      end do ! ns
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LLY Lloyd Integration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (LLY.NE.0) then
         LLY_GRTR = LLY_GRTR + LLY_GRTR_K * VOLCUB(KPT) * NSYMAT
      endif
      !
      ! Update statusbar
      if(mythread==0) then
#ifdef CPP_MPI
        if( (((k_end-k_start)/50)==0 .or. &
        mod(KPT-k_start,(k_end-k_start)/50)==0) .and. &
        t_inc%i_write>0 ) write(1337,FMT=200, advance='no')
#else
        if( ((NOFKS/50)==0 .or.mod(KPT,NOFKS/50)==0) .and.t_inc%i_write>0 ) write(1337,FMT=200, advance='no')
#endif
      end if ! mythread==0
   !
   end do ! KPT = 1,NOFKS   end K-points loop
   190 format('                 |')      ! status bar
   200 format('|')                       ! status bar
   if(t_inc%i_write>0) write(1337,*)      ! finalize status bar
   !
#ifdef CPP_HYBRID
   !$omp end parallel
#endif
         
   if(test('rhoqtest')) then
#ifdef CPP_TIMING
     call timing_start('main1b - kkrmat01 - writeout_rhoq')
#endif
      call rhoq_write_tau0(nofks,nshell,nsh1,nsh2,nsymat,nscoef,mu,iatomimp,kmask,lmmaxd,bzkp,imin)
#ifdef CPP_TIMING
     call timing_stop('main1b - kkrmat01 - writeout_rhoq')
#endif
   end if !test('rhoqtest')
   !
   TRACET = CZERO
   if (LLY.EQ.2) then
      ! Add trace of (t-tref)^-1 * d(t-tref)/dE. Remember that in this case
      ! Tr(alpha^-1 d alpha/dE) should be subtracted and
      ! Tr(alpha_ref^-1 d alpha_ref/dE) should be added.
      do I1 = 1,NAEZ
         GAUX1(1:LMMAXD,1:LMMAXD) = CFCTOR*     & !(1.D0/CFCTOR) *
            (  DTMATLL(1:LMMAXD,1:LMMAXD,I1) -  &
               DTREFLL(1:LMMAXD,1:LMMAXD,REFPOT(I1)) )
         do LM1 = 1,LMMAXD
            do LM2 = 1,LMMAXD
               TRACET = TRACET + GAUX1(LM1,LM2) * TINVLL(LM2,LM1,I1)
            enddo
         enddo
      enddo
   endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! deallocate arrays
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   i_all=-product(shape(GLLKE))*kind(GLLKE)
   deallocate(GLLKE, stat=i_stat)
   call memocc(i_stat,i_all,'GLLKE','kkrmat01')
   ! LLY Lloyd
   if (LLY.NE.0) then
      i_all=-product(shape(DGLLKE))*kind(DGLLKE)
      deallocate(DGLLKE, stat=i_stat)
      call memocc(i_stat,i_all,'DGLLKE','kkrmat01')
      i_all=-product(shape(GREFLLKE))*kind(GREFLLKE)
      deallocate(GREFLLKE, stat=i_stat)
      call memocc(i_stat,i_all,'GREFLLKE','kkrmat01')
   endif
   !----------------------------------------------------------------------------
   !
#ifdef CPP_MPI
   IF(.not.opt('qdos    ')) then
     do NS = 1,NSHELL
        IWORK = LMMAXD*LMMAXD*NSYMAXD
        WORK = CZERO
        call MPI_ALLREDUCE(GS(1,1,1,NS),WORK,IWORK,  &
           MPI_DOUBLE_COMPLEX,MPI_SUM,               &
           t_mpi_c_grid%mympi_comm_ie,IERR)
        call ZCOPY(IWORK,WORK,1,GS(1,1,1,NS),1)
     end do
     !
     if (LLY.NE.0) then
        IWORK = 1
        WORK = CZERO
        call MPI_ALLREDUCE(LLY_GRTR,WORK,IWORK,  &
           MPI_DOUBLE_COMPLEX,MPI_SUM,                  &
           t_mpi_c_grid%mympi_comm_ie,IERR)
        call ZCOPY(IWORK,WORK(1,1,1),1,LLY_GRTR,1)
     endif
     !
     if(lly.eq.2) then
        IWORK = 1
        WORK = CZERO
        call MPI_ALLREDUCE(TRACET,WORK,IWORK, &
           MPI_DOUBLE_COMPLEX,MPI_SUM,               &
           t_mpi_c_grid%mympi_comm_ie,IERR)
        call ZCOPY(IWORK,WORK,1,TRACET,1)
     endif
   end if! .not.opt('qdos    ')
#endif

   if ( TEST('flow    ') .and. (t_inc%i_write>0)) write(1337,*) '<<< KKRMAT1'

end subroutine KKRMAT01

!-------------------------------------------------------------------------------
! SUBROUTINE: GTDYSON
!> @brief Solve the Dyson equation \f$(1-g t)  G = g\f$
!> @note
!> - Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine GTDYSON(GTMAT,GMAT,NDIM,LMGF0D,NGD)

   use Constants
      use mod_dataTypes, only: dp

   implicit none
   ! .. Input variables
   integer, intent(in) :: NGD
   integer, intent(in) :: NDIM
   integer, intent(in) :: LMGF0D !< (LMAX+1)**2
   ! .. In/Out variables
   complex (kind=dp), dimension(NGD,LMGF0D), intent(inout)  :: GMAT
   complex (kind=dp), dimension(NGD,NGD), intent(inout)     :: GTMAT
   ! .. Local variables
   integer :: I,INFO
   integer, dimension(NGD) :: IPVT
   ! .. External subroutines
   external :: ZGETRF,ZGETRS
   !
   do I = 1,NDIM
      GTMAT(I,I) = CONE + GTMAT(I,I) ! GTMAT= 1 - G * T
   enddo
   !----------------------------------------------------------------------------
   ! SOLVE THE SYSTEM OF LINEAR EQUATIONS
   !----------------------------------------------------------------------------
   call ZGETRF(NDIM,NDIM,GTMAT,NGD,IPVT,INFO)
   call ZGETRS('N',NDIM,LMGF0D,GTMAT,NGD,IPVT,GMAT,NGD,INFO)
   return

end subroutine GTDYSON
