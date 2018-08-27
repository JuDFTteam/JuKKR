module mod_tmatimp_newsolver

contains

!----------------------------------------------------------------------------
! SUBROUTINE: TMATIMP_NEWSOLVER
!> @brief Calculate and write down impurity tmatrix and delta matrix
!> first calculate t-matrix for the host corresponding to imp. cluster
!> @author N. H. Long
!> @date 05.2013, Juelich
!> @note
!> - Adapted to new routines (mainly changed interfaces) to work in KKRcode
!> also added MPI parallelization. Philipp Rüssmann, Juelich, 09.2017
!> - Jonathan Chico Feb. 2018: Removed inc.p dependencies and rewrote to Fortran90
!----------------------------------------------------------------------------
subroutine TMATIMP_NEWSOLVER(IRM,KSRA,LMAX,IEND,IRID,LPOT,NATYP,NCLEB,IPAND,IRNSD,NFUND,  &
   IHOST,NTOTD,NSPIN,LMPOT,NCHEB,LMMAXD,KORBIT,NSPOTD,IELAST,IRMIND,NPAN_EQ,NPAN_LOG,      &
   NATOMIMP,C,R_LOG,IPAN,IRMIN,HOSTIMP,IPANIMP,IRWSIMP,ATOMIMP,IRMINIMP,       &
   ICLEB,IRCUT,IRCUTIMP,ZAT,ZIMP,R,CLEB,RIMP,RCLSIMP,E,VM2ZIMP,VINSIMP,   &
   DTMTRX, LMMAXSO)

#ifdef CPP_MPI
   use mpi
   use mod_mympi, only: myrank, master, nranks,distribute_linear_on_tasks
#else
   use mod_mympi, only: myrank, master, nranks
#endif
   use mod_types, only: t_inc, t_imp
   use mod_create_newmesh
   use mod_version_info
   use mod_wunfiles, only: t_params
   use Constants
   use mod_Profiling
   Use mod_datatypes, Only: dp

   use mod_calcsph
   use mod_interpolate_poten
   use mod_intcheb_cell
   use mod_rllsll
   use mod_rllsllsourceterms
   use mod_spinorbit_ham
   use mod_rotatespinframe, only: rotatematrix
   use mod_vllmat
   use mod_vllmatsra

   implicit none

   ! .. Input variables
   integer, intent(in) :: IRM       !< Maximum number of radial points
   integer, intent(in) :: KSRA
   integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
   integer, intent(in) :: IEND
   integer, intent(in) :: IRID
   integer, intent(in) :: LPOT      !< Maximum l component in potential expansion
   integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
   integer, intent(in) :: NCLEB     !< Number of Clebsch-Gordon coefficients
   integer, intent(in) :: IPAND     !< Number of panels in non-spherical part
   integer, intent(in) :: IRNSD
   integer, intent(in) :: NFUND     !< Shape functions parameters in non-spherical part
   integer, intent(in) :: NTOTD
   integer, intent(in) :: IHOST
   integer, intent(in) :: NSPIN     !< Counter for spin directions
   integer, intent(in) :: LMPOT     !< (LPOT+1)**2
   integer, intent(in) :: NCHEB     !< Number of Chebychev pannels for the new solver
   integer, intent(in) :: KORBIT    !< Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
   integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
   integer, intent(in) :: LMMAXSO
   integer, intent(in) :: NSPOTD    !< Number of potentials for storing non-sph. potentials
   integer, intent(in) :: IELAST
   integer, intent(in) :: IRMIND    !< IRM-IRNSD
   integer, intent(in) :: NPAN_EQ   !< Number of intervals from [R_LOG] to muffin-tin radius Used in conjunction with runopt NEWSOSOL
   integer, intent(in) :: NPAN_LOG  !< Number of intervals from nucleus to [R_LOG] Used in conjunction with runopt NEWSOSOL
   integer, intent(in) :: NATOMIMP  !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
   real (kind=dp), intent(in) :: C
   real (kind=dp), intent(in) :: R_LOG !< Radius up to which log-rule is used for interval width. Used in conjunction with runopt NEWSOSOL
   integer, dimension(NATYP), intent(in)    :: IPAN  !< Number of panels in non-MT-region
   integer, dimension(NATYP), intent(in)    :: IRMIN !< Max R for spherical treatment
   integer, dimension(NATYP), intent(in)    :: HOSTIMP
   integer, dimension(NATOMIMP), intent(in) :: IPANIMP
   integer, dimension(NATOMIMP), intent(in) :: IRWSIMP
   integer, dimension(NATOMIMP), intent(in) :: ATOMIMP
   integer, dimension(NATOMIMP), intent(in) :: IRMINIMP
   integer, dimension(NCLEB,4), intent(in)            :: ICLEB    !< Pointer array
   integer, dimension(0:IPAND,NATYP), intent(in)      :: IRCUT    !< R points of panel borders
   integer, dimension(0:IPAND,NATOMIMP), intent(in)   :: IRCUTIMP
   real (kind=dp), dimension(NATYP), intent(in)     :: ZAT      !< Nuclear charge
   real (kind=dp), dimension(NATOMIMP), intent(in)  :: ZIMP
   real (kind=dp), dimension(IRM,NATYP), intent(in)    :: R        !< Radial mesh ( in units a Bohr)
   real (kind=dp), dimension(NCLEB,2), intent(in)      :: CLEB     !< GAUNT coefficients (GAUNT)
   real (kind=dp), dimension(IRM,NATOMIMP), intent(in) :: RIMP
   real (kind=dp), dimension(3,NATOMIMP), intent(in)   :: RCLSIMP
   ! .. In/Out variables
   complex (kind=dp), intent(inout) :: E
   real (kind=dp), dimension(IRM,NSPIN*NATOMIMP), intent(inout) :: VM2ZIMP
   real (kind=dp), dimension(IRMIND:IRM,LMPOT,NSPIN*NATOMIMP), intent(inout) :: VINSIMP
   complex (kind=dp), dimension((KORBIT+1)*LMMAXD*NATOMIMP,(KORBIT+1)*LMMAXD*NATOMIMP), intent(inout) :: DTMTRX
   ! .. Local variables
   integer :: ipot
   integer :: I1,IR,NSRA,USE_SRATRICK,NVEC,LM1,LM2,ISPIN,I2,IL1,IL2,IRMDNEWD
   integer :: i_stat, i_all
#ifdef CPP_MPI
   integer :: ierr
#endif
   real (kind=dp) :: THETA,PHI
   complex (kind=dp) :: GMATPREFACTOR
   integer, dimension(NATYP) :: NPAN_TOT
   integer, dimension(NATYP) :: NPAN_INST
   integer, dimension(NATYP) :: NPAN_EQ_AT
   integer, dimension(NATYP) :: NPAN_LOG_AT
   real (kind=dp), dimension(NATOMIMP)  :: PHIimp
   real (kind=dp), dimension(NATYP)     :: PHIhost
   real (kind=dp), dimension(NATOMIMP)  :: THETAimp
   real (kind=dp), dimension(NATYP)     :: THETAhost
   complex (kind=dp), dimension(2*(LMAX+1)) :: TMATSPH
   complex (kind=dp), dimension(2*(LMAX+1)) :: dummy_alpha
   complex (kind=dp), dimension((KORBIT+1)*LMMAXD,(KORBIT+1)*LMMAXD,IHOST) :: TMATLL
   complex (kind=dp), dimension((KORBIT+1)*LMMAXD,(KORBIT+1)*LMMAXD)       :: dummy_alphaget
   ! .. Allocatable variables
   integer, dimension(:), allocatable :: irmdnew
   integer, dimension(:), allocatable :: JLK_INDEX
   integer, dimension(:,:), allocatable :: IPAN_INTERVALL
   real (kind=dp), dimension(:,:), allocatable :: RNEW,RPAN_INTERVALL
   real (kind=dp), dimension(:,:,:), allocatable :: VINSNEW
   complex (kind=dp), dimension(:), allocatable :: DELTATMP
   complex (kind=dp), dimension(:,:), allocatable :: HLK,JLK,HLK2,JLK2
   complex (kind=dp), dimension(:,:), allocatable :: RADIALHOST,RADIALIMP
   complex (kind=dp), dimension(:,:), allocatable :: VLLIMP,DELTAV,DELTAIMP
   complex (kind=dp), dimension(:,:,:), allocatable :: VNSIMP
   complex (kind=dp), dimension(:,:,:), allocatable :: RLL,SLL
   complex (kind=dp), dimension(:,:,:), allocatable :: DELTABG,DELTASM
   complex (kind=dp), dimension(:,:,:), allocatable :: TMATLLIMP,DELTAMTR
   complex (kind=dp), dimension(:,:,:), allocatable :: VNSPLL0,VNSPLL,VNSPLL1
   complex (kind=dp), dimension(:,:,:,:), allocatable :: RLLHOST
   complex (kind=dp), dimension(:,:,:,:), allocatable :: VNSHOST
   ! .. Parallelization variables
   integer :: i1_start, i1_end, i1_start_imp, i1_end_imp
#ifdef CPP_MPI
   integer, dimension(0:nranks-1) :: ntot_pT, ioff_pT
   complex (kind=dp), dimension(:,:,:),allocatable :: temp
   complex (kind=dp), dimension(:,:,:,:),allocatable :: temp2 ! needed for MPI communication
#endif
   logical, external :: OPT
   logical, external :: TEST

   if(myrank==master) write(6,*) 'in tmatimp'
   if (KSRA.GE.1) then
      NSRA=2
   else
      NSRA=1
   endif

   allocate(JLK_INDEX(2*LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(JLK_INDEX))*kind(JLK_INDEX),'JLK_INDEX','tmatimp_newsolver')
   allocate(DELTAV(LMMAXSO,LMMAXSO),stat=i_stat)
   call memocc(i_stat,product(shape(DELTAV))*kind(DELTAV),'DELTAV','tmatimp_newsolver')
   allocate(DELTAMTR(LMMAXSO,LMMAXSO,NATOMIMP),stat=i_stat)
   call memocc(i_stat,product(shape(DELTAMTR))*kind(DELTAMTR),'DELTAMTR','tmatimp_newsolver')
   allocate(TMATLLIMP(LMMAXSO,LMMAXSO,NATOMIMP),stat=i_stat)
   call memocc(i_stat,product(shape(TMATLLIMP))*kind(TMATLLIMP),'TMATLLIMP','tmatimp_newsolver')
   allocate(RLLHOST(NSRA*LMMAXSO,LMMAXSO,IHOST,NTOTD*(NCHEB+1)),stat=i_stat)
   call memocc(i_stat,product(shape(RLLHOST))*kind(RLLHOST),'RLLHOST','tmatimp_newsolver')
   allocate(VNSHOST(NSRA*LMMAXSO,NSRA*LMMAXSO,IHOST,NTOTD*(NCHEB+1)),stat=i_stat)
   call memocc(i_stat,product(shape(VNSHOST))*kind(VNSHOST),'VNSHOST','tmatimp_newsolver')
   TMATLL    = CZERO
   RLLHOST   = CZERO
   VNSHOST   = CZERO
   TMATLLIMP = CZERO
   DELTAMTR  = CZERO
   DELTAV    = CZERO
   TMATSPH   = CZERO

   if(myrank==master) then
      ! read angles from nonco_ange files
      open(UNIT=12,FILE='nonco_angle.dat',FORM='FORMATTED')
      do I1=1,NATYP
         read(12,*) THETAhost(i1),PHIhost(i1)
         THETAhost(i1)=THETAhost(i1)/360.0D0*2.0D0*PI
         PHIhost(i1)=PHIhost(i1)/360.0D0*2.0D0*PI
      enddo
      close(12)
      open(UNIT=13,FILE='nonco_angle_imp.dat',FORM='FORMATTED')
      do I1=1,NATOMIMP
         read(13,*) THETAimp(i1),PHIimp(i1)
         THETAimp(i1)=THETAimp(i1)/360.0D0*2.0D0*PI
         PHIimp(i1)=PHIimp(i1)/360.0D0*2.0D0*PI
      enddo
      close(13)
   endif

#ifdef CPP_MPI
   ! broadcast read-in values to all ranks (MPI_COMM_WORLD since
   ! atom dimension is solely used without energy parallelization)
   call MPI_Bcast(THETAhost, NATYP, MPI_DOUBLE_PRECISION, master,MPI_COMM_WORLD,ierr)
   if(ierr/=0) stop 'Error MPI_Bcast THETAhost in tmatimp'
   call MPI_Bcast(PHIhost, NATYP, MPI_DOUBLE_PRECISION, master,MPI_COMM_WORLD,ierr)
   if(ierr/=0) stop 'Error MPI_Bcast PHIhost in tmatimp'
   call MPI_Bcast(THETAimp, NATOMIMP, MPI_DOUBLE_PRECISION, master,MPI_COMM_WORLD,ierr)
   if(ierr/=0) stop 'Error MPI_Bcast THETAimp in tmatimp'
   call MPI_Bcast(PHIimp, NATOMIMP, MPI_DOUBLE_PRECISION, master,MPI_COMM_WORLD,ierr)
   if(ierr/=0) stop 'Error MPI_Bcast PHIimp in tmatimp'

   ! set start/end of parallel atom loops
   if(t_inc%i_write>0) write(1337, *) 'Parallelization host atoms:'
   call distribute_linear_on_tasks(nranks, myrank, master, IHOST,ntot_pT,ioff_pT,.true.)
   i1_start = ioff_pT(myrank) + 1
   i1_end   = ioff_pT(myrank) + ntot_pT(myrank)
   if(t_inc%i_write>0) write(1337, *) 'Parallelization imp. atoms:'
   call distribute_linear_on_tasks(nranks, myrank, master, NATOMIMP,ntot_pT,ioff_pT,.true.)
   i1_start_imp = ioff_pT(myrank) + 1
   i1_end_imp   = ioff_pT(myrank) + ntot_pT(myrank)
#else
   i1_start = 1
   i1_end = IHOST
   i1_start_imp = 1
   i1_end_imp = NATOMIMP
#endif

   if(OPT('GREENIMP')) then

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! START calculate tmat and radial wavefunctions of host atoms
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      ! create new mesh before loop starts
      ! data for the new mesh
      allocate(irmdnew(natyp), stat=i_stat)
      call memocc(i_stat,product(shape(irmdnew))*kind(irmdnew),'irmdnew','tmatimp_newsolver')
      IRMDNEWD = 0
      do I1=1,NATYP
         NPAN_INST(I1)= IPAN(I1)-1
         NPAN_TOT(I1)= NPAN_LOG+NPAN_EQ+NPAN_INST(I1)
         if(NPAN_TOT(I1)*(NCHEB+1)>IRMDNEWD) then
            IRMDNEWD = NPAN_TOT(I1)*(NCHEB+1)
         endif
         IRMDNEW(I1) = NPAN_TOT(I1)*(NCHEB+1)
      end do
      ! new mesh
      allocate(RNEW(IRMDNEWD,NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(RNEW))*kind(RNEW),'RNEW','tmatimp_newsolver')
      allocate(RPAN_INTERVALL(0:NTOTD,NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(RPAN_INTERVALL))*kind(RPAN_INTERVALL),'RPAN_INTERVALL','tmatimp_newsolver')
      allocate(IPAN_INTERVALL(0:NTOTD,NATYP),stat=i_stat)
      call memocc(i_stat,product(shape(IPAN_INTERVALL))*kind(IPAN_INTERVALL),'IPAN_INTERVALL','tmatimp_newsolver')
      allocate(VINSNEW(IRMDNEWD,LMPOT,NSPOTD),stat=i_stat) !NSPIND*max(NATYP,NATOMIMP)))
      call memocc(i_stat,product(shape(VINSNEW))*kind(VINSNEW),'VINSNEW','tmatimp_newsolver')
   
      ! In second step interpolate potential (gain atom by atom with NATYP==1)
      call CREATE_NEWMESH(NATYP,IRM,IPAND,IRID,NTOTD,NFUND,&
         NCHEB,IRMDNEWD,NSPIN,R(:,:),IRMIN(:),IPAN(:),IRCUT(0:IPAND,:),    &
         R_LOG,NPAN_LOG,NPAN_EQ,NPAN_LOG_AT(:),NPAN_EQ_AT(:),NPAN_TOT(:),  &
         RNEW(:,:),RPAN_INTERVALL(0:NTOTD,:),IPAN_INTERVALL(0:NTOTD,:),1)
   
      ! calculate tmat and radial wavefunctions of host atoms
      ! parallelized with MPI over atoms
      do I2=i1_start, i1_end
         I1=HOSTIMP(I2)
   
         THETA = THETAhost(i1)
         PHI = PHIhost(i1)
         ISPIN=1
         IPOT=NSPIN*(I1-1)+1
         write(6,*) 'HOST',I2,I1
   
         ! set up the non-spherical ll' matrix for potential VLL'
         if (NSRA.EQ.2) then
            USE_SRATRICK=1
         elseif (NSRA.EQ.1) then
            USE_SRATRICK=0
         endif
   
         allocate(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(VNSPLL0))*kind(VNSPLL0),'VNSPLL0','tmatimp_newsolver')
         VNSPLL0=CZERO
         ! output potential onto which SOC is added
         allocate(VNSPLL1(LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(VNSPLL1))*kind(VNSPLL1),'VNSPLL1','tmatimp_newsolver')
         VNSPLL1=CZERO
   
         call VLLMAT(1,IRMDNEW(I1),IRMDNEW(I1),LMMAXD,LMMAXSO,VNSPLL0,     &
            VINSNEW(1:IRMDNEW(I1),1:LMPOT,IPOT:IPOT+NSPIN-1),LMPOT,        &
            CLEB,ICLEB,IEND,NSPIN,                                         &
            ZAT(I1),RNEW(1:IRMDNEW(I1),I1),USE_SRATRICK,NCLEB)
         ! contruct the spin-orbit coupling hamiltonian and add to potential
         call SPINORBIT_HAM(LMAX,LMMAXD,                                &
            VINSNEW(1:IRMDNEW(I1),1:LMPOT,IPOT:IPOT+NSPIN-1),          &
            RNEW(1:IRMDNEW(I1),I1),E,ZAT(I1),C,t_params%SOCSCALE(I1),   &
            NSPIN,LMPOT,THETA,PHI,IPAN_INTERVALL(0:NTOTD,I1),          &
            RPAN_INTERVALL(0:NTOTD,I1),NPAN_TOT(I1),NCHEB,              &
            IRMDNEW(I1),IRMDNEW(I1),VNSPLL0,                            &
            VNSPLL1,'1')
         ! extend matrix for the SRA treatment
         if (NSRA.EQ.2) then
            allocate(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW(I1)),stat=i_stat)
            call memocc(i_stat,product(shape(VNSPLL))*kind(VNSPLL),'VNSPLL','tmatimp_newsolver')
            if (USE_SRATRICK.EQ.0) then
               call VLLMATSRA(VNSPLL1,VNSPLL,RNEW(1:IRMDNEW(I1),I1),LMMAXSO,  &
                  IRMDNEW(I1),IRMDNEW(I1),E,LMAX,0,'Ref=0')
            elseif (USE_SRATRICK.EQ.1) then
               call VLLMATSRA(VNSPLL1,VNSPLL,RNEW(1:IRMDNEW(I1),I1),LMMAXSO,  &
                  IRMDNEW(I1),IRMDNEW(I1),E,LMAX,0,'Ref=Vsph')
            endif
         else
            allocate(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
            call memocc(i_stat,product(shape(VNSPLL))*kind(VNSPLL),'VNSPLL','tmatimp_newsolver')
            VNSPLL(:,:,:)=VNSPLL1(:,:,:)
         endif
   
         ! calculate the source terms in the Lippmann-Schwinger equation
         ! these are spherical hankel and bessel functions
         allocate(HLK(1:4*(LMAX+1),IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(HLK))*kind(HLK),'HLK','tmatimp_newsolver')
         allocate(JLK(1:4*(LMAX+1),IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(JLK))*kind(JLK),'JLK','tmatimp_newsolver')
         allocate(HLK2(1:4*(LMAX+1),IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(HLK2))*kind(HLK2),'HLK2','tmatimp_newsolver')
         allocate(JLK2(1:4*(LMAX+1),IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(JLK2))*kind(JLK2),'JLK2','tmatimp_newsolver')
         HLK=CZERO
         JLK=CZERO
         HLK2=CZERO
         JLK2=CZERO
         GMATPREFACTOR=CZERO
         call RLLSLLSOURCETERMS(NSRA,NVEC,E,RNEW(1:IRMDNEW(I1),I1),  &
            IRMDNEW(I1),IRMDNEW(I1),                                 &
            LMAX,LMMAXSO,1,JLK_INDEX,HLK,JLK,HLK2,                   &
            JLK2,GMATPREFACTOR)
         ! using spherical potential as reference
         if (USE_SRATRICK.EQ.1) then
            call CALCSPH(NSRA,IRMDNEW(I1),IRMDNEW(I1),LMAX,NSPIN,ZAT(I1),E, &
               LMPOT,LMMAXSO,RNEW(1:IRMDNEW(I1),I1),                         &
               VINSNEW(1:IRMDNEW(I1),1:LMPOT,IPOT:IPOT+NSPIN-1),             &
               NCHEB,NPAN_TOT(I1),                                            &
               RPAN_INTERVALL(0:NTOTD,I1),JLK_INDEX,HLK,JLK,HLK2,             &
               JLK2,GMATPREFACTOR,TMATSPH,dummy_alpha,USE_SRATRICK)
         endif
   
         ! calculate the tmat and wavefunctions
         allocate(RLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(RLL))*kind(RLL),'RLL','tmatimp_newsolver')
         allocate(SLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(SLL))*kind(SLL),'SLL','tmatimp_newsolver')
         RLL=CZERO
         SLL=CZERO
   
         ! right solutions
         call RLLSLL(RPAN_INTERVALL(0:NTOTD,I1),RNEW(1:IRMDNEW(I1),I1), &
            VNSPLL,RLL,SLL,TMATLL(:,:,I2),NCHEB,                        &
            NPAN_TOT(I1),LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW(I1),   &
            NSRA,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR, &
            '1','1','0',USE_SRATRICK,dummy_alphaget)
         if (NSRA.EQ.2) then
            do IR=1,IRMDNEW(I1)
               do LM1=1,LMMAXSO
                  do LM2=1,LMMAXSO
                     RLL(LM1+LMMAXSO,LM2,IR)=RLL(LM1+LMMAXSO,LM2,IR)/C
                     SLL(LM1+LMMAXSO,LM2,IR)=SLL(LM1+LMMAXSO,LM2,IR)/C
                  enddo
               enddo
            enddo
         endif
         ! save radial wavefunction for a host
         do IR=1,IRMDNEW(I1)
            do LM1=1,NVEC*LMMAXSO
               do LM2=1,LMMAXSO
                  RLLHOST(LM1,LM2,I2,IR)=RLL(LM1,LM2,IR)
               enddo
            enddo
         enddo
         ! add spherical contribution of tmatrix
   
         if (USE_SRATRICK.EQ.1) then
            do LM1=1,(KORBIT+1)*LMMAXD
               TMATLL(LM1,LM1,I2)=TMATLL(LM1,LM1,I2)+TMATSPH(JLK_INDEX(LM1))
            enddo
         endif
   
         ! rotate tmatrix and radial wavefunction to global frame
         call ROTATEMATRIX(TMATLL(1,1,I2),THETA,PHI,LMMAXD,0)
   
         ! create SRA potential for host
         ! set up the non-spherical ll' matrix for potential VLL'
         VNSPLL0=CZERO
         VNSPLL1=CZERO
         call VLLMAT(1,IRMDNEW(I1),IRMDNEW(I1),LMMAXD,LMMAXSO,VNSPLL0,  &
            VINSNEW(1:IRMDNEW(I1),1:LMPOT,IPOT:IPOT+NSPIN-1), LMPOT,         &
            CLEB,ICLEB,IEND,NSPIN,ZAT(I1),RNEW(1:IRMDNEW(I1),I1),0,NCLEB)
   
         ! contruct the spin-orbit coupling hamiltonian and add to potential
         call SPINORBIT_HAM(LMAX,LMMAXD,                                &
            VINSNEW(1:IRMDNEW(I1),1:LMPOT,IPOT:IPOT+NSPIN-1),          &
            RNEW(1:IRMDNEW(I1),I1),E,ZAT(I1),C,t_params%SOCSCALE(I1),   &
            NSPIN,LMPOT,THETA,PHI,IPAN_INTERVALL(0:NTOTD,I1),          &
            RPAN_INTERVALL(0:NTOTD,I1),NPAN_TOT(I1),NCHEB,              &
            IRMDNEW(I1),IRMDNEW(I1),VNSPLL0,                            &
            VNSPLL1,'1')
   
         ! save potential for a host
         do IR=1,IRMDNEW(I1)
            do LM1=1,LMMAXSO
               do LM2=1,LMMAXSO
                  VNSHOST(LM1,LM2,I2,IR)=VNSPLL1(LM1,LM2,IR)
                  if (NSRA.EQ.2) then
                     VNSHOST(LM1+LMMAXSO,LM2+LMMAXSO,I2,IR)=VNSPLL1(LM1,LM2,IR)
                  endif
               enddo
            enddo
         enddo
   
         i_all=-product(shape(VNSPLL0))*kind(VNSPLL0)
         deallocate(VNSPLL0,stat=i_stat)
         call memocc(i_stat,i_all,'VNSPLL0','tmatimp_newsolver')
         i_all=-product(shape(VNSPLL1))*kind(VNSPLL1)
         deallocate(VNSPLL1,stat=i_stat)
         call memocc(i_stat,i_all,'VNSPLL1','tmatimp_newsolver')
         i_all=-product(shape(VNSPLL))*kind(VNSPLL)
         deallocate(VNSPLL,stat=i_stat)
         call memocc(i_stat,i_all,'VNSPLL','tmatimp_newsolver')
         i_all=-product(shape(HLK))*kind(HLK)
         deallocate(HLK,stat=i_stat)
         call memocc(i_stat,i_all,'HLK','tmatimp_newsolver')
         i_all=-product(shape(JLK))*kind(JLK)
         deallocate(JLK,stat=i_stat)
         call memocc(i_stat,i_all,'JLK','tmatimp_newsolver')
         i_all=-product(shape(HLK2))*kind(HLK2)
         deallocate(HLK2,stat=i_stat)
         call memocc(i_stat,i_all,'HLK2','tmatimp_newsolver')
         i_all=-product(shape(JLK2))*kind(JLK2)
         deallocate(JLK2,stat=i_stat)
         call memocc(i_stat,i_all,'JLK2','tmatimp_newsolver')
         i_all=-product(shape(SLL))*kind(SLL)
         deallocate(SLL,stat=i_stat)
         call memocc(i_stat,i_all,'SLL','tmatimp_newsolver')
         i_all=-product(shape(RLL))*kind(RLL)
         deallocate(RLL,stat=i_stat)
         call memocc(i_stat,i_all,'RLL','tmatimp_newsolver')
      enddo ! I2
   
      i_all=-product(shape(RNEW))*kind(RNEW)
      deallocate(RNEW,stat=i_stat)
      call memocc(i_stat,i_all,'RNEW','tmatimp_newsolver')
      i_all=-product(shape(VINSNEW))*kind(VINSNEW)
      deallocate(VINSNEW,stat=i_stat)
      call memocc(i_stat,i_all,'VINSNEW','tmatimp_newsolver')
      i_all=-product(shape(RPAN_INTERVALL))*kind(RPAN_INTERVALL)
      deallocate(RPAN_INTERVALL,stat=i_stat)
      call memocc(i_stat,i_all,'RPAN_INTERVALL','tmatimp_newsolver')
      i_all=-product(shape(IPAN_INTERVALL))*kind(IPAN_INTERVALL)
      deallocate(IPAN_INTERVALL,stat=i_stat)
      call memocc(i_stat,i_all,'IPAN_INTERVALL','tmatimp_newsolver')
   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
#ifdef CPP_MPI
      ! collect results and write out only on master
      ! communicate VNSHOST, RLLHOST, TMATLL
      i1 = NSRA*LMMAXSO*LMMAXSO*IHOST*NTOTD*(NCHEB+1)
      ! Allocation of temp2 for RLLHOST
      allocate(temp2(NSRA*LMMAXSO,LMMAXSO,IHOST,NTOTD*(NCHEB+1)),stat=i_stat)
      call memocc(i_stat,product(shape(temp2))*kind(temp2),'temp2','tmatimp_newsolver')
      temp2 = CZERO
   
      call MPI_ALLREDUCE(RLLHOST,temp2,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(ierr/=0) stop 'Error in MPI_Allreduce for RLLHOST in tmatimp'
      RLLHOST = temp2
      ! Deallocation of temp2 for RLLHOST
      i_all=-product(shape(temp2))*kind(temp2)
      deallocate(temp2,stat=i_stat)
      call memocc(i_stat,i_all,'temp2','tmatimp_newsolver')
   
      i1 = NSRA*LMMAXSO*NSRA*LMMAXSO*IHOST*NTOTD*(NCHEB+1)
      ! Allocation of temp2 for VNSHOST
      allocate(temp2(NSRA*LMMAXSO,NSRA*LMMAXSO,IHOST,NTOTD*(NCHEB+1)),stat=i_stat)
      call memocc(i_stat,product(shape(temp2))*kind(temp2),'temp2','tmatimp_newsolver')
      temp2 = CZERO
   
      call MPI_ALLREDUCE(VNSHOST,temp2,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(ierr/=0) stop 'Error in MPI_Allreduce for VNSHOST in tmatimp'
      VNSHOST = temp2
      ! Deallocation of temp2 for RLLHOST
      i_all=-product(shape(temp2))*kind(temp2)
      deallocate(temp2,stat=i_stat)
      call memocc(i_stat,i_all,'temp2','tmatimp_newsolver')
   
      i1 = (KORBIT+1)*LMMAXD*(KORBIT+1)*LMMAXD*IHOST
      ! Allocation of temp for TMATLL
      allocate(temp(LMMAXSO,LMMAXSO,NATOMIMP),stat=i_stat)
      call memocc(i_stat,product(shape(temp))*kind(temp),'temp','tmatimp_newsolver')
      temp = CZERO
      call MPI_ALLREDUCE(TMATLL,temp,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(ierr/=0) stop 'Error in MPI_Allreduce for TMATLL in tmatimp'
      TMATLL = temp
      ! Deallocation of temp for TMATLL
      i_all=-product(shape(temp))*kind(temp)
      deallocate(temp,stat=i_stat)
      call memocc(i_stat,i_all,'temp','tmatimp_newsolver')
#endif
   
      ! write out DTMTRX file containgin Delta_T and Delta-matrices
      if(myrank==master) then
         if (IELAST.EQ.1) then
            open(UNIT=20,FILE='DTMTRX',FORM='FORMATTED')
            write(20,'(I5)') NATOMIMP
            do I1=1,NATOMIMP
               write(20,'(3e17.9)') (RCLSIMP(I2,I1),I2=1,3)
            enddo
         endif
      end if
   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END  calculate tmat and radial wavefunctions of host atoms
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   elseif(myrank==master) then

     write(*,*) 'skipping host atom loop in tmatimp_newsolver'

   end if ! (OPT('GREENIMP'))

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! START calculate tmat and radial wavefunctions of impurity atoms
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! create new mesh before loop starts
   ! data for the new mesh
   i_all=-product(shape(irmdnew))*kind(irmdnew)
   deallocate(irmdnew,stat=i_stat)
   call memocc(i_stat,i_all,'irmdnew','tmatimp_newsolver')

   allocate(irmdnew(natomimp), stat=i_stat)
   call memocc(i_stat,product(shape(irmdnew))*kind(irmdnew),'irmdnew','tmatimp_newsolver')
   IRMDNEWD = 0
   do I1=1,NATOMIMP
      NPAN_INST(I1)= IPANIMP(I1)-1
      NPAN_TOT(I1)= NPAN_LOG+NPAN_EQ+NPAN_INST(I1)
      if(NPAN_TOT(I1)*(NCHEB+1)>IRMDNEWD) then
         IRMDNEWD = NPAN_TOT(I1)*(NCHEB+1)
      endif
      IRMDNEW(I1) = NPAN_TOT(I1)*(NCHEB+1)
   end do
   ! new mesh
   allocate(RNEW(IRMDNEWD,NATOMIMP),stat=i_stat)
   call memocc(i_stat,product(shape(RNEW))*kind(RNEW),'RNEW','tmatimp_newsolver')
   allocate(RPAN_INTERVALL(0:NTOTD,NATOMIMP),stat=i_stat)
   call memocc(i_stat,product(shape(RPAN_INTERVALL))*kind(RPAN_INTERVALL),'RPAN_INTERVALL','tmatimp_newsolver')
   allocate(IPAN_INTERVALL(0:NTOTD,NATOMIMP),stat=i_stat)
   call memocc(i_stat,product(shape(IPAN_INTERVALL))*kind(IPAN_INTERVALL),'IPAN_INTERVALL','tmatimp_newsolver')
   allocate(VINSNEW(IRMDNEWD,LMPOT,NSPOTD),stat=i_stat)
   call memocc(i_stat,product(shape(VINSNEW))*kind(VINSNEW),'VINSNEW','tmatimp_newsolver')

   ! initialize with zeros
   TMATLLIMP = CZERO
   TMATSPH = CZERO

   call CREATE_NEWMESH(NATOMIMP,IRM,IPAND,IRID,NTOTD,NFUND,&
      NCHEB,IRMDNEWD,NSPIN,RIMP(:,1:NATOMIMP),IRMINIMP(1:NATOMIMP),        &
      IPANIMP(1:NATOMIMP),IRCUTIMP(0:IPAND,1:NATOMIMP),R_LOG,NPAN_LOG,     &
      NPAN_EQ,NPAN_LOG_AT(1:NATOMIMP),NPAN_EQ_AT(1:NATOMIMP),              &
      NPAN_TOT(1:NATOMIMP),RNEW(1:IRMDNEWD,1:NATOMIMP),                    &
      RPAN_INTERVALL(0:NTOTD,1:NATOMIMP),IPAN_INTERVALL(0:NTOTD,1:NATOMIMP),1)

   ! In second step interpolate potential
   call INTERPOLATE_POTEN(LPOT,IRM,IRNSD,NATOMIMP,IPAND,LMPOT,NSPOTD, &
      NTOTD,IRMDNEWD,                                      &
      NSPIN,RIMP(:,1:NATOMIMP),IRMINIMP(1:NATOMIMP),              &
      IRWSIMP(1:NATOMIMP),                                        &
      IRCUTIMP(0:IPAND,1:NATOMIMP),                               &
      VINSIMP(IRMIND:IRM,1:LMPOT,1:NATOMIMP),                   &
      VM2ZIMP(1:IRM,1:NATOMIMP),NPAN_LOG_AT(1:NATOMIMP),         &
      NPAN_EQ_AT(1:NATOMIMP),NPAN_TOT(1:NATOMIMP),                &
      RNEW(1:IRMDNEWD,1:NATOMIMP),                                &
      IPAN_INTERVALL(0:NTOTD,1:NATOMIMP),VINSNEW)

   ! now start loop over atoms
   do I1=i1_start_imp,i1_end_imp
      THETA = THETAimp(i1)
      PHI = PHIimp(i1)
      ISPIN=1
      IPOT=NSPIN*(I1-1)+1
      write(6,*) 'IMP',I1

      allocate(VNSIMP(NSRA*LMMAXSO,NSRA*LMMAXSO,IRMDNEW(I1)),stat=i_stat)
      call memocc(i_stat,product(shape(VNSIMP))*kind(VNSIMP),'VNSIMP','tmatimp_newsolver')
      VNSIMP=CZERO
      ! set up the non-spherical ll' matrix for potential VLL'
      if (NSRA.EQ.2) then
         USE_SRATRICK=1
      elseif (NSRA.EQ.1) then
         USE_SRATRICK=0
      endif
      allocate(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
      call memocc(i_stat,product(shape(VNSPLL0))*kind(VNSPLL0),'VNSPLL0','tmatimp_newsolver')
      VNSPLL0=CZERO
      ! output potential onto which SOC is added
      allocate(VNSPLL1(LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
      call memocc(i_stat,product(shape(VNSPLL1))*kind(VNSPLL1),'VNSPLL1','tmatimp_newsolver')
      VNSPLL1=CZERO

      call VLLMAT(1,IRMDNEW(I1),IRMDNEW(I1),LMMAXD,LMMAXSO,VNSPLL0,  &
         VINSNEW(1:IRMDNEW(I1),1:LMPOT,IPOT:IPOT+NSPIN-1),LMPOT,          &
         CLEB,ICLEB,IEND,NSPIN,                                      &
         ZIMP(I1),RNEW(1:IRMDNEW(I1),I1),USE_SRATRICK,NCLEB)

      ! Contruct the spin-orbit coupling hamiltonian and add to potential
      call SPINORBIT_HAM(LMAX,LMMAXD,                                &
         VINSNEW(1:IRMDNEW(I1),1:LMPOT,IPOT:IPOT+NSPIN-1),          &
         RNEW(1:IRMDNEW(I1),I1),E,ZIMP(I1),C,t_params%SOCSCALE(I1),  &
         NSPIN,LMPOT,THETA,PHI,IPAN_INTERVALL(0:NTOTD,I1),          &
         RPAN_INTERVALL(0:NTOTD,I1),NPAN_TOT(I1),NCHEB,              &
         IRMDNEW(I1),IRMDNEW(I1),VNSPLL0,                            &
         VNSPLL1,'1')

      ! extend matrix for the SRA treatment
      if (NSRA.EQ.2) then
         allocate(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(VNSPLL))*kind(VNSPLL),'VNSPLL','tmatimp_newsolver')
         if (USE_SRATRICK.EQ.0) then
            call VLLMATSRA(VNSPLL1,VNSPLL,                  &
               RNEW(1:IRMDNEW(I1),I1),LMMAXSO,IRMDNEW(I1),  &
               IRMDNEW(I1),E,LMAX,0,'Ref=0')
         elseif (USE_SRATRICK.EQ.1) then
            call VLLMATSRA(VNSPLL1,VNSPLL,                  &
               RNEW(1:IRMDNEW(I1),I1),LMMAXSO,IRMDNEW(I1),  &
               IRMDNEW(I1),E,LMAX,0,'Ref=Vsph')
         endif
      else
         allocate(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(VNSPLL))*kind(VNSPLL),'VNSPLL','tmatimp_newsolver')
         VNSPLL(:,:,:)=VNSPLL1(:,:,:)
      endif

      ! calculate the source terms in the Lippmann-Schwinger equation
      ! these are spherical hankel and bessel functions
      allocate(HLK(1:4*(LMAX+1),IRMDNEW(I1)),stat=i_stat)
      call memocc(i_stat,product(shape(HLK))*kind(HLK),'HLK','tmatimp_newsolver')
      allocate(JLK(1:4*(LMAX+1),IRMDNEW(I1)),stat=i_stat)
      call memocc(i_stat,product(shape(JLK))*kind(JLK),'JLK','tmatimp_newsolver')
      allocate(HLK2(1:4*(LMAX+1),IRMDNEW(I1)),stat=i_stat)
      call memocc(i_stat,product(shape(HLK2))*kind(HLK2),'HLK2','tmatimp_newsolver')
      allocate(JLK2(1:4*(LMAX+1),IRMDNEW(I1)),stat=i_stat)
      call memocc(i_stat,product(shape(JLK2))*kind(JLK2),'JLK2','tmatimp_newsolver')
      HLK=CZERO
      JLK=CZERO
      HLK2=CZERO
      JLK2=CZERO
      GMATPREFACTOR=CZERO
      call RLLSLLSOURCETERMS(NSRA,NVEC,E,RNEW(1:IRMDNEW(I1),I1),  &
         IRMDNEW(I1),IRMDNEW(I1),                                 &
         LMAX,LMMAXSO,1,JLK_INDEX,                                &
         HLK,JLK,HLK2,JLK2,GMATPREFACTOR)
      ! using spherical potential as reference
      if (USE_SRATRICK.EQ.1) then
         call CALCSPH(NSRA,IRMDNEW(I1),IRMDNEW(I1),LMAX,NSPIN,ZIMP(I1),  &
            E,LMPOT,LMMAXSO,RNEW(1:IRMDNEW(I1),I1),                       &
            VINSNEW(1:IRMDNEW(I1),1:LMPOT,IPOT:IPOT+NSPIN-1),             &
            NCHEB,NPAN_TOT(I1),RPAN_INTERVALL(0:NTOTD,I1),                 &
            JLK_INDEX,HLK,JLK,HLK2,JLK2,                                   &
            GMATPREFACTOR,TMATSPH,dummy_alpha,USE_SRATRICK)
      endif

      ! calculate the tmat and wavefunctions
      allocate(RLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
      call memocc(i_stat,product(shape(RLL))*kind(RLL),'RLL','tmatimp_newsolver')
      allocate(SLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
      call memocc(i_stat,product(shape(SLL))*kind(SLL),'SLL','tmatimp_newsolver')
      RLL=CZERO
      SLL=CZERO

      ! Right solutions
      call RLLSLL(RPAN_INTERVALL(0:NTOTD,I1),RNEW(1:IRMDNEW(I1),I1), &
         VNSPLL,RLL,SLL,TMATLLIMP(:,:,I1),NCHEB,                     &
         NPAN_TOT(I1),LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW(I1),   &
         NSRA,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR, &
         '1','1','0',USE_SRATRICK,dummy_alphaget)
      if (NSRA.EQ.2) then
         do IR=1,IRMDNEW(I1)
            do LM1=1,LMMAXSO
               do LM2=1,LMMAXSO
                  RLL(LM1+LMMAXSO,LM2,IR)=RLL(LM1+LMMAXSO,LM2,IR)/C
                  SLL(LM1+LMMAXSO,LM2,IR)=SLL(LM1+LMMAXSO,LM2,IR)/C
               enddo
            enddo
         enddo
      endif

      ! for OPERATOR option save impurity wavefuncitons
      if (OPT('OPERATOR')) then
        t_imp%RLLIMP(:,:,:,i1) = RLL(:,:,:)
      end if

      ! add spherical contribution of tmatrix
      if (USE_SRATRICK.EQ.1) then
         do LM1=1,(KORBIT+1)*LMMAXD
            TMATLLIMP(LM1,LM1,I1)=TMATLLIMP(LM1,LM1,I1)+TMATSPH(JLK_INDEX(LM1))
         enddo
      endif

      ! rotate tmatrix and radial wavefunction to global frame
      call ROTATEMATRIX(TMATLLIMP(:,:,I1),THETA,PHI,LMMAXD,0)

      ! create SRA potential for impurity
      ! set up the non-spherical ll' matrix for potential VLL'
      VNSPLL0=CZERO
      call VLLMAT(1,IRMDNEW(I1),IRMDNEW(I1),LMMAXD,LMMAXSO,VNSPLL0,  &
         VINSNEW(1:IRMDNEW(I1),1:LMPOT,IPOT:IPOT+NSPIN-1),LMPOT,          &
         CLEB,ICLEB,IEND,NSPIN,                                      &
         ZIMP(I1),RNEW(1:IRMDNEW(I1),I1),0,NCLEB)
      !     +             ZIMP(I1),RNEW(:,I1),USE_SRATRICK)

      ! contruct the spin-orbit coupling hamiltonian and add to potential
      call SPINORBIT_HAM(LMAX,LMMAXD,                                &
         VINSNEW(1:IRMDNEW(I1),1:LMPOT,IPOT:IPOT+NSPIN-1),          &
         RNEW(1:IRMDNEW(I1),I1),E,ZIMP(I1),C,t_params%SOCSCALE(I1),  &
         NSPIN,LMPOT,THETA,PHI,IPAN_INTERVALL(0:NTOTD,I1),          &
         RPAN_INTERVALL(0:NTOTD,I1),NPAN_TOT(I1),NCHEB,              &
         IRMDNEW(I1),IRMDNEW(I1),VNSPLL0,                            &
         VNSPLL1,'1')
      do IR=1,IRMDNEW(I1)
         do LM1=1,LMMAXSO
            do LM2=1,LMMAXSO
               VNSIMP(LM1,LM2,IR)=VNSPLL1(LM1,LM2,IR)
               if (NSRA.EQ.2) then
                  VNSIMP(LM1+LMMAXSO,LM2+LMMAXSO,IR)=VNSPLL1(LM1,LM2,IR)
               endif
            enddo
         enddo
      enddo

      ! calculate delta_t_imp matrix written in TMATLLIMP
      do I2=1,IHOST
         if (ATOMIMP(I1).EQ.HOSTIMP(I2)) then
            do LM1=1,LMMAXSO
               do LM2=1,LMMAXSO
                  TMATLLIMP(LM1,LM2,I1)=TMATLLIMP(LM1,LM2,I1)-TMATLL(LM1,LM2,I2)
               enddo
            enddo
            do LM1=1,NSRA*LMMAXSO
               do LM2=1,NSRA*LMMAXSO
                  do IR=1,IRMDNEW(I1)
                     VNSIMP(LM1,LM2,IR)=VNSIMP(LM1,LM2,IR)-VNSHOST(LM1,LM2,I2,IR)
                  enddo
               enddo
            enddo
         endif
      enddo

      ! calculate delta matrix \delta=int{R_imp*\deltaV*R_host}
      if (IELAST.EQ.1) then
         ALLOCATE(DELTABG(LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(DELTABG))*kind(DELTABG),'DELTABG','tmatimp_newsolver')
         ALLOCATE(DELTASM(LMMAXSO,LMMAXSO,IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(DELTASM))*kind(DELTASM),'DELTASM','tmatimp_newsolver')
         DELTABG=CZERO
         DELTASM=CZERO
         allocate(DELTATMP(IRMDNEW(I1)),stat=i_stat)
         call memocc(i_stat,product(shape(DELTATMP))*kind(DELTATMP),'DELTATMP','tmatimp_newsolver')
         allocate(RADIALHOST(LMMAXSO,LMMAXSO),stat=i_stat)
         call memocc(i_stat,product(shape(RADIALHOST))*kind(RADIALHOST),'RADIALHOST','tmatimp_newsolver')
         allocate(RADIALIMP(LMMAXSO,LMMAXSO),stat=i_stat)
         call memocc(i_stat,product(shape(RADIALIMP))*kind(RADIALIMP),'RADIALIMP','tmatimp_newsolver')
         allocate(VLLIMP(LMMAXSO,LMMAXSO),stat=i_stat)
         call memocc(i_stat,product(shape(VLLIMP))*kind(VLLIMP),'VLLIMP','tmatimp_newsolver')
         DELTATMP=CZERO

         ! big component for SRA stored in DELTABG
         do IR=1,IRMDNEW(I1)
            RADIALHOST=CZERO
            RADIALIMP=CZERO
            VLLIMP=CZERO
            DELTAV=CZERO
            do LM1=1,LMMAXSO
               do LM2=1,LMMAXSO
                  do I2=1,IHOST
                     if (ATOMIMP(I1).EQ.HOSTIMP(I2)) then
                        RADIALHOST(LM1,LM2)=RLLHOST(LM1,LM2,I2,IR)
                     endif
                  enddo !I2
                  RADIALIMP(LM1,LM2)=RLL(LM1,LM2,IR)
                  VLLIMP(LM1,LM2)=VNSIMP(LM1,LM2,IR)
               enddo !LM2
            enddo !LM1

            call ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,VLLIMP,     &
               LMMAXSO,RADIALIMP,LMMAXSO,CZERO,DELTAV,LMMAXSO)
            call ZGEMM('C','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST, &
               LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTABG(:,:,IR),LMMAXSO)

            ! small component for SRA stored in DELTASM
            if (NSRA.EQ.2) then
               RADIALHOST=CZERO
               RADIALIMP=CZERO
               VLLIMP=CZERO
               DELTAV=CZERO
               do LM1=1,LMMAXSO
                  do LM2=1,LMMAXSO
                     do I2=1,IHOST
                        if (ATOMIMP(I1).EQ.HOSTIMP(I2)) then
                           RADIALHOST(LM1,LM2)=RLLHOST(LM1+LMMAXSO,LM2,I2,IR)
                        endif
                     enddo
                     RADIALIMP(LM1,LM2)=RLL(LM1+LMMAXSO,LM2,IR)
                     VLLIMP(LM1,LM2)=VNSIMP(LM1+LMMAXSO,LM2+LMMAXSO,IR)
                  enddo
               enddo
               call ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,VLLIMP,     &
                  LMMAXSO,RADIALIMP,LMMAXSO,CZERO,DELTAV,LMMAXSO)
               call ZGEMM('C','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST, &
                  LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTASM(:,:,IR),LMMAXSO)

               ! sum up big and small component stored in DELTABG
               do LM1=1,LMMAXSO
                  do LM2=1,LMMAXSO
                     DELTABG(LM1,LM2,IR)=DELTABG(LM1,LM2,IR)+DELTASM(LM1,LM2,IR)
                  enddo
               enddo

            endif ! NSRA
         enddo ! IR

         ! integrate
         do LM1=1,LMMAXSO
            do LM2=1,LMMAXSO
               do IR=1,IRMDNEW(I1)
                  DELTATMP(IR)=DELTABG(LM1,LM2,IR)
               enddo
               call INTCHEB_CELL(DELTATMP,DELTAMTR(LM1,LM2,I1),&
                  RPAN_INTERVALL(0:NTOTD,I1),                  &
                  IPAN_INTERVALL(0:NTOTD,I1),                  &
                  NPAN_TOT(I1),NCHEB,IRMDNEW(I1))
            enddo
         enddo

         i_all=-product(shape(DELTATMP))*kind(DELTATMP)
         deallocate(DELTATMP,stat=i_stat)
         call memocc(i_stat,i_all,'DELTATMP','tmatimp_newsolver')
         i_all=-product(shape(RADIALHOST))*kind(RADIALHOST)
         deallocate(RADIALHOST,stat=i_stat)
         call memocc(i_stat,i_all,'RADIALHOST','tmatimp_newsolver')
         i_all=-product(shape(RADIALIMP))*kind(RADIALIMP)
         deallocate(RADIALIMP,stat=i_stat)
         call memocc(i_stat,i_all,'RADIALIMP','tmatimp_newsolver')
         i_all=-product(shape(VLLIMP))*kind(VLLIMP)
         deallocate(VLLIMP,stat=i_stat)
         call memocc(i_stat,i_all,'VLLIMP','tmatimp_newsolver')
         i_all=-product(shape(DELTABG))*kind(DELTABG)
         deallocate(DELTABG,stat=i_stat)
         call memocc(i_stat,i_all,'DELTABG','tmatimp_newsolver')
         i_all=-product(shape(DELTASM))*kind(DELTASM)
         deallocate(DELTASM,stat=i_stat)
         call memocc(i_stat,i_all,'DELTASM','tmatimp_newsolver')

      endif ! IELAST.EQ.1

      i_all=-product(shape(VNSPLL0))*kind(VNSPLL0)
      deallocate(VNSPLL0,stat=i_stat)
      call memocc(i_stat,i_all,'VNSPLL0','tmatimp_newsolver')
      i_all=-product(shape(VNSPLL1))*kind(VNSPLL1)
      deallocate(VNSPLL1,stat=i_stat)
      call memocc(i_stat,i_all,'VNSPLL1','tmatimp_newsolver')
      i_all=-product(shape(VNSPLL))*kind(VNSPLL)
      deallocate(VNSPLL,stat=i_stat)
      call memocc(i_stat,i_all,'VNSPLL','tmatimp_newsolver')
      i_all=-product(shape(HLK))*kind(HLK)
      deallocate(HLK,stat=i_stat)
      call memocc(i_stat,i_all,'HLK','tmatimp_newsolver')
      i_all=-product(shape(JLK))*kind(JLK)
      deallocate(JLK,stat=i_stat)
      call memocc(i_stat,i_all,'JLK','tmatimp_newsolver')
      i_all=-product(shape(HLK2))*kind(HLK2)
      deallocate(HLK2,stat=i_stat)
      call memocc(i_stat,i_all,'HLK2','tmatimp_newsolver')
      i_all=-product(shape(JLK2))*kind(JLK2)
      deallocate(JLK2,stat=i_stat)
      call memocc(i_stat,i_all,'JLK2','tmatimp_newsolver')
      i_all=-product(shape(RLL))*kind(RLL)
      deallocate(RLL,stat=i_stat)
      call memocc(i_stat,i_all,'RLL','tmatimp_newsolver')
      i_all=-product(shape(SLL))*kind(SLL)
      deallocate(SLL,stat=i_stat)
      call memocc(i_stat,i_all,'SLL','tmatimp_newsolver')
      i_all=-product(shape(VNSIMP))*kind(VNSIMP)
      deallocate(VNSIMP,stat=i_stat)
      call memocc(i_stat,i_all,'VNSIMP','tmatimp_newsolver')

   enddo ! I1 impurity

   i_all=-product(shape(RNEW))*kind(RNEW)
   deallocate(RNEW,stat=i_stat)
   call memocc(i_stat,i_all,'RNEW','tmatimp_newsolver')
   i_all=-product(shape(VINSNEW))*kind(VINSNEW)
   deallocate(VINSNEW,stat=i_stat)
   call memocc(i_stat,i_all,'VINSNEW','tmatimp_newsolver')
   i_all=-product(shape(RPAN_INTERVALL))*kind(RPAN_INTERVALL)
   deallocate(RPAN_INTERVALL,stat=i_stat)
   call memocc(i_stat,i_all,'RPAN_INTERVALL','tmatimp_newsolver')
   i_all=-product(shape(IPAN_INTERVALL))*kind(IPAN_INTERVALL)
   deallocate(IPAN_INTERVALL,stat=i_stat)
   call memocc(i_stat,i_all,'IPAN_INTERVALL','tmatimp_newsolver')

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! END calculate tmat and radial wavefunctions of impurity atoms
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! final writeout only on master
#ifdef CPP_MPI
   !collect results and write out only on master
   !collect TMATLLIMP, DELTAMTR
   !communicate TMATLLIMP, DELTAMTR
   i1 = LMMAXSO*LMMAXSO*NATOMIMP
   allocate(temp(LMMAXSO,LMMAXSO,NATOMIMP),stat=i_stat)
   call memocc(i_stat,product(shape(temp))*kind(temp),'temp','tmatimp_newsolver')
   temp = CZERO
   call MPI_ALLREDUCE(TMATLLIMP,temp,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
   if(ierr/=0) stop 'Error in MPI_Allreduce for TMATLLIMP in tmatimp'
   TMATLLIMP = temp
   i_all=-product(shape(temp))*kind(temp)
   deallocate(temp,stat=i_stat)
   call memocc(i_stat,i_all,'temp','tmatimp_newsolver')

   i1 = LMMAXSO*LMMAXSO*NATOMIMP
   allocate(temp(LMMAXSO,LMMAXSO,NATOMIMP),stat=i_stat)
   call memocc(i_stat,product(shape(temp))*kind(temp),'temp','tmatimp_newsolver')
   temp = CZERO
   call MPI_ALLREDUCE(DELTAMTR,temp,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
   if(ierr/=0) stop 'Error in MPI_Allreduce for DELTAMTR in tmatimp'
   DELTAMTR = temp
   i_all=-product(shape(temp))*kind(temp)
   deallocate(temp,stat=i_stat)
   call memocc(i_stat,i_all,'temp','tmatimp_newsolver')
#endif

   ! collect results and writeout only for GREENIMP option
   if(OPT('GREENIMP') .and. myrank==master) then

      do I1=1,NATOMIMP
         do LM1=1,LMMAXSO
            do LM2=1,LMMAXSO
               IL1=LMMAXSO*(I1-1)+LM1
               IL2=LMMAXSO*(I1-1)+LM2
               DTMTRX(IL1,IL2)=TMATLLIMP(LM1,LM2,I1)
            enddo
         enddo
      enddo

      ! write down to the file DTMTRX
      if (IELAST.EQ.1) then
         allocate(DELTAIMP((KORBIT+1)*LMMAXD*NATOMIMP,(KORBIT+1)*LMMAXD*NATOMIMP),stat=i_stat)
         call memocc(i_stat,product(shape(DELTAIMP))*kind(DELTAIMP),'DELTAIMP','tmatimp_newsolver')
         DELTAIMP=CZERO
         do I1=1,NATOMIMP
            do LM1=1,LMMAXSO
               do LM2=1,LMMAXSO
                  IL1=LMMAXSO*(I1-1)+LM1
                  IL2=LMMAXSO*(I1-1)+LM2
                  DELTAIMP(IL1,IL2)=DELTAMTR(LM1,LM2,I1)
               enddo
            enddo
         enddo
         do LM1=1,LMMAXSO*NATOMIMP
            do LM2=1,LMMAXSO*NATOMIMP
               write(20,'((2I5),(4e17.9))') LM2,LM1,DTMTRX(LM2,LM1),DELTAIMP(LM2,LM1)
            enddo
         enddo
         i_all=-product(shape(DELTAIMP))*kind(DELTAIMP)
         deallocate(DELTAIMP,stat=i_stat)
         call memocc(i_stat,i_all,'DELTAIMP','tmatimp_newsolver')
         if(myrank==master) write(6,*) 'end of delta t'
      endif ! IELAST.EQ.1

      close(20) ! output file DTMTRX

   end if ! myrank==master

   i_all=-product(shape(VNSHOST))*kind(VNSHOST)
   deallocate(VNSHOST,stat=i_stat)
   call memocc(i_stat,i_all,'VNSHOST','tmatimp_newsolver')
   i_all=-product(shape(RLLHOST))*kind(RLLHOST)
   deallocate(RLLHOST,stat=i_stat)
   call memocc(i_stat,i_all,'RLLHOST','tmatimp_newsolver')
   i_all=-product(shape(TMATLLIMP))*kind(TMATLLIMP)
   deallocate(TMATLLIMP,stat=i_stat)
   call memocc(i_stat,i_all,'TMATLLIMP','tmatimp_newsolver')
   i_all=-product(shape(DELTAMTR))*kind(DELTAMTR)
   deallocate(DELTAMTR,stat=i_stat)
   call memocc(i_stat,i_all,'DELTAMTR','tmatimp_newsolver')
   i_all=-product(shape(DELTAV))*kind(DELTAV)
   deallocate(DELTAV,stat=i_stat)
   call memocc(i_stat,i_all,'DELTAV','tmatimp_newsolver')
   i_all=-product(shape(irmdnew))*kind(irmdnew)
   deallocate(irmdnew,stat=i_stat)
   call memocc(i_stat,i_all,'irmdnew','tmatimp_newsolver')
end subroutine TMATIMP_NEWSOLVER

end module mod_tmatimp_newsolver
