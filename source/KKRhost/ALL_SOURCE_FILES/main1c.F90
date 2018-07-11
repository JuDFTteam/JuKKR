!-------------------------------------------------------------------------------
! MODULE: MOD_MAIN1C
!> @brief Wrapper module for the calculation of the density for the JM-KKR package
!> @details The code uses the information obtained in the main0 module, this is
!> mostly done via the get_params_1c() call, that obtains parameters of the type
!> t_params and passes them to local variables
!> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
!> and many others ...
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
module MOD_MAIN1C

   use Profiling
   use Constants
   use global_variables
   use mod_DataTypes

   implicit none

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: main1c
   !> @brief Main subroutine regarding the calculation of the electronic density
   !> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
   !> and many others ...
   !----------------------------------------------------------------------------
   subroutine main1c()

#ifdef CPP_MPI
      use mpi
#endif
#ifdef CPP_TIMING
      use mod_timing
#endif
      use mod_types, only: t_tgmat, t_inc, t_lloyd
#ifdef CPP_MPI
      use mod_types, only: gather_tmat, t_mpi_c_grid, save_t_mpi_c_grid,   &
                           get_ntot_pT_ioff_pT_2D
#endif
      use mod_mympi, only: myrank, nranks, master
#ifdef CPP_MPI
      use mod_mympi, only: find_dims_2d,distribute_linear_on_tasks,  &
                           mympi_main1c_comm,mympi_main1c_comm_newsosol2
#endif
      use mod_wunfiles
      use mod_version_info

      use mod_main0

      ! .. Parameters
      integer :: NMVECMAX
      parameter (NMVECMAX = 4)
      !     ..
      ! .. Local variables
      integer :: NQDOS
      integer :: LMAXP1
      integer :: NACLS1
      integer :: LMAXD1
      integer :: LRECTMT
      integer :: LMMAXD1
      integer :: ITMPDIR
      integer :: NSPINPOT

      integer :: i1_start, i1_end,i_stat,i_all
      integer :: ie_start, ie_end, ie_num
      integer :: L,I1,IE,IR,IS,LM,IQ,ierr,IPOT
      integer :: ICELL,IPOT1,ISPIN,IHOST,ILTMP
      real (kind=dp) :: DENEF
      real (kind=dp) :: CHRGSEMICORE
      character(len=5)  :: TEXTNS
      character(len=80) :: TMPDIR
      character(len=8)  :: QDOSOPT
      logical :: LMOMVEC
      logical :: LDORHOEF
      logical :: ITERMVDIR
      complex (kind=dp) :: CSUM                    ! LLY Lloyd
      complex (kind=dp) :: EREAD                   ! LLY
      integer, dimension(20,NATYP) :: NKCORE    !< Number of KAPPA values for a given (n,l) core state
      integer, dimension(20,NPOTD) :: KAPCORE   !< The (maximum 2) values of KAPPA
      real (kind=dp), dimension(NATYPD)                    :: EU
      real (kind=dp), dimension(NATYPD)                    :: EDC
      real (kind=dp), dimension(NATYPD)                    :: PHI
      real (kind=dp), dimension(NATYPD)                    :: THETA
      real (kind=dp), dimension(NATYPD)                    :: DENEFAT
      real (kind=dp), dimension(NSPIND)                   :: CHARGE_LLY  ! LLY
      real (kind=dp), dimension(0:LMAXD+1,NPOTD)           :: ESPV
      real (kind=dp), dimension(0:LMAXD+1,2)               :: ESPV1
      real (kind=dp), dimension(0:LMAXD+1,2)               :: DOSTOT
      real (kind=dp), dimension(KREL*20+(1-KREL),NPOTD)   :: ECOREREL !< for a given (n,l) state the core energies corresponding first/second KAPPA value, AVERAGED over \mu's  These values are written out to the  potential file (routine <RITES>), but the read in (routine <STARTB1>) updates the ECORE array
      real (kind=dp), dimension(2,NATYPD)                  :: angles_new
      real (kind=dp), dimension(0:LMAXD+1,NATYPD,2)         :: CHARGE
      real (kind=dp), dimension(MMAXD,MMAXD,NSPIND,NATYPD) :: WLDAUOLD
      complex (kind=dp), dimension(IEMXD)                   :: DF
      complex (kind=dp), dimension(0:LMAXD+1,IEMXD,2)        :: DEN1
      complex (kind=dp), dimension(NATYPD,3,NMVECMAX)        :: MVEVI          ! OUTPUT
      complex (kind=dp), dimension(NATYPD,3,NMVECMAX)        :: MVEVIEF        ! OUTPUT
      complex (kind=dp), dimension(MMAXD,MMAXD,NPOTD)       :: DENMATC
      complex (kind=dp), dimension(MMAXD,MMAXD,2,2,NATYPD)   :: DENMATN
      complex (kind=dp), dimension(0:LMAXD,3,NMVECMAX)       :: MVEVIL1
      complex (kind=dp), dimension(0:LMAXD,3,NMVECMAX)       :: MVEVIL2        ! WORK ARRAYS
      complex (kind=dp), dimension(0:LMAXD,NATYPD,3,NMVECMAX) :: MVEVIL
      character(len=7), dimension(3) :: TEXTS
      character(len=4), dimension(0:6) :: TEXTL
      !-------------------------------------------------------------------------
      !> @note attention: muorb second index means both spins and total
      !-------------------------------------------------------------------------
      real (kind=dp), dimension(IRMD*KREL + (1-KREL),NATYPD) :: RHOORB   !< orbital density
      real (kind=dp), dimension(0:LMAXD+1+1,3,NATYPD) :: MUORB           !< orbital magnetic moment
      ! ----------------------------------------------------------------------
      !  R2NEF (IRMD,LMPOTD,NATYP,2)  ! rho at FERMI energy
      !  RHO2NS(IRMD,LMPOTD,NATYP,2)  ! radial density
      !   nspin=1            : (*,*,*,1) radial charge density
      !   nspin=2 or krel=1  : (*,*,*,1) rho(2) + rho(1) -> charge
      !                               (*,*,*,2) rho(2) - rho(1) -> mag. moment
      !  RHOC(IRMD,NPOTD)              ! core charge density
      ! ----------------------------------------------------------------------
      real (kind=dp), dimension(IRMD,NPOTD)    :: RHOC   !< core charge density
      real (kind=dp), dimension(IRMD,LMPOTD,4) :: RHO2M1
      real (kind=dp), dimension(IRMD,LMPOTD,4) :: RHO2M2
      real (kind=dp), dimension(NRMAXD,LMPOTD,NSPOTD) :: VINSNEW
      !---------------------------------------------------------------

      ! .. Alocatable arrays
      real (kind=dp), dimension(:,:), allocatable      :: QVEC
      real (kind=dp), dimension(:,:,:), allocatable    :: RHO2N1
      real (kind=dp), dimension(:,:,:), allocatable    :: RHO2N2
      real (kind=dp), dimension(:,:,:,:), allocatable  :: R2NEF     !< rho at FERMI energy
      real (kind=dp), dimension(:,:,:,:), allocatable  :: RHO2NS    !< radial density
      complex (kind=dp), dimension(:,:,:,:), allocatable    :: DEN   ! DEN(0:LMAXD1,IEMXD,NPOTD,NQDOS)
      complex (kind=dp), dimension(:,:,:,:), allocatable    :: DENLM ! DENLM(LMMAXD1,IEMXD,NPOTD,NQDOS)
      complex (kind=dp), dimension(:), allocatable ::  CDOS2          ! LLY Lloyd
      complex (kind=dp), dimension(:), allocatable :: CDOS0
      complex (kind=dp), dimension(:), allocatable :: CDOS1
      complex (kind=dp), dimension(:), allocatable :: CDOSAT0
      complex (kind=dp), dimension(:), allocatable :: CDOSAT1
      complex (kind=dp), dimension(:,:), allocatable :: CDOS_LLY

      !-------------------------------------------------------------------------
      ! MPI parameters
      !-------------------------------------------------------------------------
#ifdef CPP_MPI
      integer :: ntot1
      integer :: ihelp
      integer :: idim, nranks_local, myrank_ie_tmp
      complex (kind=dp) :: DENTOT ! qdos
      integer, dimension(0:nranks-1) :: ntot_pT
      integer, dimension(0:nranks-1) :: ioff_pT
      complex (kind=dp), dimension(:,:), allocatable :: WORK
      complex (kind=dp), dimension(:,:,:,:), allocatable :: workc
#endif

      ! .. Intrinsic Functions ..
      intrinsic ATAN,DBLE,aimag,REAL
      ! .. External Functions ..
      logical OPT,TEST
      external OPT,TEST
      !     ..
      !     .. Data statements ..
      DATA TEXTL/' s =',' p =',' d =',' f =',' g =',' h =',' i ='/
      DATA TEXTS/'spin dn','spin up','       '/
      DATA TEXTNS/' ns ='/
      DATA LDORHOEF/.TRUE./
      DATA IHOST / 1 /          ! this is the host program

      ! .. Calculate parameters
      LMAXD1=LMAX+1
      LMMAXD1=LMMAXD+1
      LRECTMT=WLENGTH*4*LMMAXD*LMMAXD


      ! Allocate arrays
      allocate(RHO2N1(IRMD,LMPOTD,NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(RHO2N1))*kind(RHO2N1),'RHO2N1','main1c')
      allocate(RHO2NS(IRMD,LMPOTD,NATYPD,2),stat=i_stat)
      call memocc(i_stat,product(shape(RHO2NS))*kind(RHO2NS),'RHO2NS','main1c')
      allocate(RHO2N2(IRMD,LMPOTD,NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(RHO2N2))*kind(RHO2N2),'RHO2N2','main1c')
      allocate(R2NEF(IRMD,LMPOTD,NATYPD,2),stat=i_stat)
      call memocc(i_stat,product(shape(R2NEF))*kind(R2NEF),'R2NEF','main1c')

      ! Initialze to zero
      VINS        = 0.0d0
      R2NEF       = 0.0d0
      MUORB       = 0.0d0
      RHO2NS      = 0.0d0
      THETAS      = 0.0d0
      VINSNEW     = 0.0d0
      THETASNEW   = 0.0d0
      angles_new  = 0.0d0

      ! Consistency check
      if ( (KREL.lt.0) .or. (KREL.gt.1) ) stop ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
      if ( (KREL.eq.1) .and. (NSPIND.eq.2) ) stop ' set NSPIND = 1 for KREL = 1 in the inputcard'
      !-------------------------------------------------------------------------
      ! This routine previously used to read from unformatted files created by
      ! the main0 module, now  instead of unformatted files take parameters from
      ! types defined in wunfiles.F90
      !-------------------------------------------------------------------------
      call get_params_1c(t_params,KREL,NAEZD,NATYPD,NCLEB,LM2D,NCHEB,IPAND,LMPOTD,  &
         LMAXD,LMXSPD,NFUND,NPOTD,NTOTD,MMAXD,IEMXD,IRMD,NSRA,INS,NSPIN,NACLS1,    &
         ICST,KMROT,IQAT,IDOLDAU,IRWS,IPAN,IRCUT,IEND,ICLEB, LOFLM,JEND,IFUNM1,  &
         LMSP1,NFU,LLMSP,LCORE,NCORE,NTCELL,IRMIN,ITITLE,INTERVX,INTERVY,INTERVZ,&
         LLY,ITMPDIR,ILTMP,NPAN_EQ_AT,IPAN_INTERVALL,NPAN_LOG_AT,NPAN_TOT,NTLDAU,LOPT, &
         ITLDAU,IELAST,IESEMICORE,NPOL,IRSHIFT,JWSREL,ZREL,ITRUNLDAU,QMTET,QMPHI,&
         CONC,ALAT,ZAT,DRDI,RMESH,A,B,CLEB,THETAS,SOCSCALE,RPAN_INTERVALL,CSCL,  &
         RNEW,SOCSCL,THETASNEW,EFERMI,EREFLDAU,UEFF,JEFF,EMIN,EMAX,TK,VINS,VISP, &
         ECORE,DRDIREL,R2DRDIREL,RMREL,VTREL,BTREL,WLDAU,ULDAU,EZ,WEZ,PHILDAU,   &
         TMPDIR,SOLVER,NSPIND,NSPOTD,IRMIND,LMAXD1,NCELLD,IRID,R_LOG, NAEZ, NATYP, LMAX)

      ! Initialization needed due to merging to one executable
      ESPV(:,:)      = 0.d0
      RHO2N1(:,:,:)  = 0.d0
      RHO2N2(:,:,:)  = 0.
      !-------------------------------------------------------------------------
      !                       End read in variables
      !-------------------------------------------------------------------------

      if ( IDOLDAU.eq.1 ) then
         !
         open (67,FILE='ldau.unformatted',FORM='unformatted')
         read (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         close(67)
      endif

      NQDOS=1
      if (OPT('qdos    ')) then                                                  ! qdos
         !        Read BZ path for qdos calculation:                             ! qdos
         open(67,FILE='qvec.dat')                                                ! qdos
         read(67,*) NQDOS                                                        ! qdos
         allocate(QVEC(3,NQDOS),stat=i_stat)                                     ! qdos
         call memocc(i_stat,product(shape(QVEC))*kind(QVEC),'QVEC','main1c')     ! qdos
         do IQ = 1,NQDOS                                                         ! qdos
            read(67,*) (QVEC(I1,IQ),I1=1,3)                                      ! qdos
         enddo                                                                   ! qdos
         close(67)                                                               ! qdos
      end if

      allocate(DEN(0:LMAXD1,IELAST,NQDOS,NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(DEN))*kind(DEN),'DEN','main1c')
      allocate(DENLM(LMMAXD,IELAST,NQDOS,NPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(DENLM))*kind(DENLM),'DENLM','main1c')

      call CINIT(IELAST*(LMAX+2)*NPOTD*NQDOS,DEN)
      call CINIT(IELAST*(LMMAXD)*NPOTD*NQDOS,DENLM)
      DENEF = 0.0D0
      call RINIT(NATYP,DENEFAT)
      !
      ITERMVDIR = OPT('ITERMDIR')
      LMOMVEC = ( ITERMVDIR .or. ( KMROT.ne.0 ) )
      NSPINPOT = KREL*2 + (1-KREL)*NSPIN
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! no need to calculate charge correction if no host program, if decimation
      ! or if no energy contour
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      LDORHOEF = (IHOST.EQ.1) .AND. (.NOT.OPT('DECIMATE')).AND. (NPOL.NE.0)
      !-------------------------------------------------------------------------
      ! LDA+U
      !-------------------------------------------------------------------------
      if ( IDOLDAU.EQ.1 ) CALL CINIT(MMAXD*MMAXD*NPOTD,DENMATC(1,1,1))
      !-------------------------------------------------------------------------
      ! LDA+U
      !-------------------------------------------------------------------------
      if (t_tgmat%gmat_to_file) then
         call OPENDAFILE(69,'gmat',4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)
      endif

      !    write parameters file that contains passed parameters for further treatment of gflle
      if (OPT('lmlm-dos')) then                                                  ! lmlm-dos
         QDOSOPT = 'n'                                                           ! lmlm-dos
         if (OPT('qdos    ')) then                                               ! lmlm-dos qdos
            QDOSOPT = 'y'                                                        ! lmlm-dos qdos
         endif                                                                   ! lmlm-dos qdos
         open(67,FORM='formatted',FILE='parameters.gflle')                       ! lmlm-dos
         DF(:)=WEZ(:)/DBLE(NSPIN)                                                ! lmlm-dos
         write(67,*) IELAST,IEMXD,NATYP,NSPIN,LMAX,QDOSOPT,DF(1:IELAST), &       ! lmlm-dos
            EZ(1:IELAST),KORBIT                                                  ! lmlm-dos
         close(67)                                                               ! lmlm-dos
      endif  ! OPT('lmlm-dos')                                                   ! lmlm-dos

      ! ------------------------------------------------------------------------
      ! LLY Lloyd
      ! ------------------------------------------------------------------------
      if (LLY.ne.0) then                                                         ! LLY Lloyd

         allocate(CDOS0(IELAST), stat=i_stat)
         call memocc(i_stat,product(shape(CDOS0))*kind(CDOS0),'CDOS0','main1c')
         allocate(CDOS1(IELAST), stat=i_stat)
         call memocc(i_stat,product(shape(CDOS1))*kind(CDOS1),'CDOS1','main1c')
         allocate(CDOS2(NATYPD), stat=i_stat)
         call memocc(i_stat,product(shape(CDOS2))*kind(CDOS2),'CDOS2','main1c')
         allocate(CDOSAT0(IELAST), stat=i_stat)
         call memocc(i_stat,product(shape(CDOSAT0))*kind(CDOSAT0),'CDOSAT0','main1c')
         allocate(CDOSAT1(IELAST), stat=i_stat)            
         call memocc(i_stat,product(shape(CDOSAT1))*kind(CDOSAT1),'CDOSAT1','main1c')
         allocate(CDOS_LLY(IELAST,NSPIND), stat=i_stat)            
         call memocc(i_stat,product(shape(CDOS_LLY))*kind(CDOS_LLY),'CDOS_LLY','main1c')

         ! LLY Lloyd
         ! Calculate free-space contribution to dos                              ! LLY Lloyd
         CDOS0(1:IELAST) = CZERO                                                  ! LLY Lloyd
         CDOS1(1:IELAST) = CZERO                                                  ! LLY Lloyd
         CDOS2(1:NAEZ) = CZERO                                                   ! LLY Lloyd
         do I1 = 1,NAEZ                                                          ! LLY Lloyd
            CDOSAT0(1:IELAST) = CZERO                                             ! LLY Lloyd
            CDOSAT1(1:IELAST) = CZERO                                             ! LLY Lloyd
            ICELL = NTCELL(I1)                                                   ! LLY Lloyd
            do IE = 1,IELAST                                                     ! LLY Lloyd
               call RHOVAL0(EZ(IE),DRDI(1,I1),RMESH(1,I1),IPAN(I1),  &           ! LLY Lloyd
                  IRCUT(0,I1),IRWS(I1),THETAS(1,1,ICELL),CDOSAT0(IE),&           ! LLY Lloyd
                  CDOSAT1(IE),IRMD,LMAX)                                          ! LLY Lloyd

               ! calculate contribution from free space

               CDOS2(I1) = CDOS2(I1) + CDOSAT1(IE)*WEZ(IE)                       ! LLY Lloyd
            enddo                                                                ! LLY Lloyd
            CDOS0(1:IELAST) = CDOS0(1:IELAST) + CDOSAT0(1:IELAST)                   ! LLY Lloyd
            CDOS1(1:IELAST) = CDOS1(1:IELAST) + CDOSAT1(1:IELAST)                   ! LLY Lloyd
         enddo                                                                   ! LLY Lloyd
         CDOS0(:) = -CDOS0(:) / PI                                               ! LLY Lloyd
         CDOS1(:) = -CDOS1(:) / PI                                               ! LLY Lloyd
         ! LLY Lloyd
         if(myrank==master) then
            open(701,FILE='freedos.dat',FORM='FORMATTED')                        ! LLY Lloyd
            do IE = 1,IELAST                                                     ! LLY Lloyd
               write(701,FMT='(10E16.8)') EZ(IE),CDOS0(IE),CDOS1(IE)             ! LLY Lloyd
            enddo                                                                ! LLY Lloyd
            close(701)                                                           ! LLY Lloyd
            open(701,FILE='singledos.dat',FORM='FORMATTED')                      ! LLY Lloyd
            do I1 = 1,NATYP                                                      ! LLY Lloyd
               write(701,FMT='(I5,10E16.8)') I1,CDOS2(I1)                        ! LLY Lloyd
            enddo                                                                ! LLY Lloyd
            close(701)                                                           ! LLY Lloyd
         end if ! myrank==master

         ! energy integration to get cdos_lly
         CDOS_LLY(1:ielast,1:NSPIN) = CZERO                                      ! LLY Lloyd
         if (.not.OPT('NEWSOSOL')) then
            if(t_lloyd%cdos_diff_lly_to_file) then
               if(myrank==master) then
                  open (701,FILE='cdosdiff_lly.dat',FORM='FORMATTED')               ! LLY Lloyd
                  do ISPIN = 1,NSPIN                                                ! LLY Lloyd
                     do IE = 1,IELAST                                               ! LLY Lloyd
                        read(701,FMT='(10E25.16)') EREAD,CDOS_LLY(IE,ISPIN)         ! LLY Lloyd
                     enddo                                                          ! LLY Lloyd
                  enddo                                                             ! LLY Lloyd
                  close(701)                                                        ! LLY Lloyd
               end if
#ifdef CPP_MPI
               call MPI_Bcast(CDOS_LLY, ielast*nspin, MPI_DOUBLE_COMPLEX, master, t_mpi_c_grid%myMPI_comm_at, ierr)
#endif
            else   !(t_lloyd%cdos_diff_lly_to_file)
#ifdef CPP_MPI
               ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
               ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
#else
               ie_start = 0
               ie_end = IELAST
#endif
               do ISPIN = 1,NSPIN                                                ! LLY Lloyd
                  do ie_num=1,ie_end
                     IE = ie_start+ie_num
                     CDOS_LLY(IE,ISPIN) = t_lloyd%cdos(ie_num,ispin)             ! LLY Lloyd
                  enddo  !ie_num                                                 ! LLY Lloyd
               enddo    !ispin
#ifdef CPP_MPI
               ! MPI gather cdos_lly on all processors
               ihelp = ielast*nspin  !IELAST*NSPIN
               allocate(work(ielast,nspin), stat=i_stat)
               call memocc(i_stat,product(shape(work))*kind(work),'work','main1c')
               work = (0.d0, 0.d0)
               call MPI_ALLREDUCE(cdos_lly(1:ielast,1:nspin),work,ihelp,   &
                  MPI_DOUBLE_COMPLEX,MPI_SUM,t_mpi_c_grid%myMPI_comm_at,ierr)
               call zcopy(ihelp,work,1,cdos_lly(1:ielast,1:nspin),1)

               i_all=-product(shape(work))*kind(work)
               deallocate(work, stat=i_stat)
               call memocc(i_stat,i_all,'work','main1c')
#endif
            end if   !(t_lloyd%cdos_diff_lly_to_file)
            ! LLY Lloyd
            ! Add free-space contribution cdos0                                  ! LLY Lloyd
            do ISPIN = 1,NSPIN                                                   ! LLY Lloyd
               CDOS_LLY(1:IELAST,ISPIN) = CDOS_LLY(1:IELAST,ISPIN)&
               + CDOS0(1:IELAST)                                                 ! LLY Lloyd
            enddo                                                                ! LLY Lloyd
            ! LLY Lloyd
            CHARGE_LLY(1:NSPIND) = 0.D0                                          ! LLY Lloyd
            do ISPIN = 1,NSPIN                                                   ! LLY Lloyd
               CSUM = CZERO                                                      ! LLY Lloyd
               do IE = 1,IELAST                                                  ! LLY Lloyd
                  CSUM = CSUM + CDOS_LLY(IE,ISPIN) * WEZ(IE)                     ! LLY Lloyd
               enddo                                                             ! LLY Lloyd
               CHARGE_LLY(ISPIN) = -aimag(CSUM) * PI / NSPINPOT                  ! LLY Lloyd
            enddo                                                                ! LLY Lloyd
            ! LLY Lloyd
            ! LLY Lloyd
            if(myrank==master) then
               open (701,FILE='cdos_lloyd.dat',FORM='FORMATTED')                 ! LLY Lloyd
               do ISPIN=1,NSPIN                                                  ! LLY Lloyd
                  do IE=1,IELAST                                                 ! LLY Lloyd
                     write(701,FMT='(10E16.8)') EZ(IE),CDOS_LLY(IE,ISPIN)        ! LLY Lloyd
                  enddo                                                          ! LLY Lloyd
               enddo                                                             ! LLY Lloyd
               close(701)                                                        ! LLY Lloyd
               write(*,*) 'Valence charge from Lloyds formula:',&                ! LLY Lloyd
                  (CHARGE_LLY(ISPIN),ISPIN=1,NSPIN)                              ! LLY Lloyd
               if(t_inc%i_write>0) then                                          ! LLY Lloyd
                  write(1337,*) 'Valence charge from Lloyds formula:',&          ! LLY Lloyd
                     (CHARGE_LLY(ISPIN),ISPIN=1,NSPIN)                           ! LLY Lloyd
               endif
            end if ! myrank==master

         else !NEWSOSOL

            if(t_lloyd%cdos_diff_lly_to_file) then
               open (701,FILE='cdosdiff_lly.dat',FORM='FORMATTED')               ! LLY
               do IE = 1,IELAST                                                  ! LLY
                  read(701,FMT='(10E25.16)') EREAD,CDOS_LLY(IE,1)                ! LLY
               enddo                                                             ! LLY
               close(701)                                                        ! LLY
            else  !(t_lloyd%cdos_diff_lly_to_file)
#ifdef CPP_MPI
               ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
               ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
#else
               ie_start = 0
               ie_end = IELAST
#endif
               do ie_num=1,ie_end
                  IE = ie_start+ie_num
                  CDOS_LLY(IE,1) = t_lloyd%cdos(ie_num,1)
               enddo
#ifdef CPP_MPI
               ! MPI gather cdos_lly on all processors
               ihelp = ielast*nspin  !IELAST*NSPIN
               allocate(work(ielast,nspin), stat=i_stat)
               call memocc(i_stat,product(shape(work))*kind(work),'work','main1c')
               work = (0.d0, 0.d0)
               CALL MPI_ALLREDUCE(cdos_lly(1:ielast,1:nspin),work,ihelp,   &
                  MPI_DOUBLE_COMPLEX,MPI_SUM,t_mpi_c_grid%myMPI_comm_at,ierr)
               call zcopy(ihelp,work,1,cdos_lly(1:ielast,1:nspin),1)
               i_all=-product(shape(work))*kind(work)
               deallocate(work,stat=i_stat)
               call memocc(i_stat,i_all,'work','main1c')
#endif
            end if  !(t_lloyd%cdos_diff_lly_to_file)

            ! Add free-space contribution cdos0                                        ! LLY
            CDOS_LLY(1:ielast,1)=CDOS_LLY(1:ielast,1)+2D0*CDOS0(1:ielast)              ! LLY
            CSUM=CZERO                                                                 ! LLY
            do IE = 1,IELAST                                                           ! LLY
               CSUM = CSUM + CDOS_LLY(IE,1) * WEZ(IE)                                  ! LLY
            enddo                                                                      ! LLY
            CHARGE_LLY(1) = -aimag(CSUM) * PI                                          ! LLY
            if(myrank==master) then
               open (701,FILE='cdos_lloyd.dat',FORM='FORMATTED')                       ! LLY
               do IE=1,IELAST                                                          ! LLY
                  write(701,FMT='(10E16.8)') EZ(IE),CDOS_LLY(IE,1)                     ! LLY
               enddo                                                                   ! LLY
               close(701)                                                              ! LLY
               write(*,*) 'Valence charge from Lloyds formula:',(CHARGE_LLY(1))        ! LLY
               if(t_inc%i_write>0) then                                                ! LLY
                  write(1337,*) 'Valence charge from Lloyds formula:',(CHARGE_LLY(1))  ! LLY
               endif                                                                   ! LLY
            end if ! myrank==master

         endif ! NEWSOSOL

      endif !LLY<>0                                                     ! LLY
      !-------------------------------------------------------------------------
      ! LLY Lloyd
      !-------------------------------------------------------------------------

      if (.not.OPT('NEWSOSOL')) then
         !----------------------------------------------------------------------
         ! NATYP
         !----------------------------------------------------------------------
#ifdef CPP_MPI
         ntot1 = t_inc%NATYP
         call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie,  &
            t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,        &
            master,ntot1,ntot_pT,ioff_pT,.true.)

         i1_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
         i1_end = ioff_pT(t_mpi_c_grid%myrank_ie)+ntot_pT(t_mpi_c_grid%myrank_ie)
         t_mpi_c_grid%ntot1 = ntot_pT(t_mpi_c_grid%myrank_ie)

         if(.not.(allocated(t_mpi_c_grid%ntot_pT1).and.  &
            allocated(t_mpi_c_grid%ioff_pT1))) then
            allocate(t_mpi_c_grid%ntot_pT1(0:t_mpi_c_grid%nranks_ie-1),stat=i_stat)
            call memocc(i_stat,product(shape(t_mpi_c_grid%ntot_pT1))*kind(t_mpi_c_grid%ntot_pT1),'t_mpi_c_grid%ntot_pT1','main1c')
            allocate(t_mpi_c_grid%ioff_pT1(0:t_mpi_c_grid%nranks_ie-1),stat=i_stat)
            call memocc(i_stat,product(shape(t_mpi_c_grid%ioff_pT1))*kind(t_mpi_c_grid%ioff_pT1),'t_mpi_c_grid%ioff_pT1','main1c')
         endif
         t_mpi_c_grid%ntot_pT1 = ntot_pT
         t_mpi_c_grid%ioff_pT1 = ioff_pT
#else
         i1_start = 1
         i1_end   = NATYP
#endif
         do I1=i1_start, i1_end
            !-------------------------------------------------------------------
            ! SPIN
            !-------------------------------------------------------------------
            IQ = IQAT(I1)
            do ISPIN = 1,NSPIN
               ICELL = NTCELL(I1)
               IPOT = (I1-1) * NSPINPOT + ISPIN
               IPOT1 = (I1-1) * NSPINPOT + 1
#ifdef CPP_TIMING
               call timing_start('main1c - rhoval')
#endif
               call RHOVAL(IHOST,LDORHOEF,ICST,INS,IELAST,NSRA,ISPIN,NSPIN,   &
                  NSPINPOT,I1,EZ,WEZ,DRDI(1,I1),RMESH(1,I1),                  &
                  VINS(IRMIND,1,KNOSPH*IPOT+(1-KNOSPH)),VISP(1,IPOT),ZAT(I1), &
                  IPAN(I1),IRCUT(0,I1),IRMIN(I1),THETAS(1,1,ICELL),           &
                  IFUNM1(1,ICELL),LMSP1(1,ICELL),RHO2N1(1,1,IPOT1),           &
                  RHO2N2(1,1,IPOT1),RHOORB(1,I1),DEN(0,1,1,IPOT),             &
                  DENLM(1,1,1,IPOT),MUORB(0,1,I1),ESPV(0,IPOT1),CLEB,LOFLM,   &
                  ICLEB,IEND,JEND,SOLVER,SOCSCL(1,KREL*I1+(1-KREL)),          &
                  CSCL(1,KREL*I1+(1-KREL)),VTREL(1,I1),BTREL(1,I1),           &
                  RMREL(1,I1),DRDIREL(1,I1),R2DRDIREL(1,I1),ZREL(I1),         &
                  JWSREL(I1),IRSHIFT(I1),LMOMVEC,QMTET(IQ),QMPHI(IQ),MVEVIL1, &
                  MVEVIL2,NMVECMAX,IDOLDAU,LOPT(I1),PHILDAU(1,I1),            &
                  WLDAU(1,1,1,I1),DENMATC(1,1,IPOT),NATYP,NQDOS,LMAX)
#ifdef CPP_TIMING
               call timing_pause('main1c - rhoval')
#endif
            end do
            !-------------------------------------------------------------------
            ! SPIN
            !-------------------------------------------------------------------
            IPOT1 = (I1-1)*NSPINPOT + 1
            !
            do LM = 1,LMPOTD
               do IR = 1,IRMD
                  RHO2NS(IR,LM,I1,1) = RHO2N1(IR,LM,IPOT1)
                  R2NEF(IR,LM,I1,1)  = RHO2N2(IR,LM,IPOT1)
               end do
            end do
            !
            do L = 0,LMAXD1
               DENEF = DENEF - 2.0D0 * CONC(I1)*   &
                  aimag(DEN(L,IELAST,1,IPOT1))/PI/DBLE(NSPINPOT)
               DENEFAT(I1) = DENEFAT(I1) - 2.0D0*  &
                  aimag(DEN(L,IELAST,1,IPOT1))/PI/DBLE(NSPINPOT)
            end do
            !
            if (NSPINPOT.EQ.2) then
               do LM = 1,LMPOTD
                  do IR = 1,IRMD
                     RHO2NS(IR,LM,I1,2) = RHO2N1(IR,LM,IPOT1+1)
                     R2NEF(IR,LM,I1,2)  = RHO2N2(IR,LM,IPOT1+1)
                  end do
               end do
               !
               do L = 0,LMAXD1
                  DENEF = DENEF - 2.0D0 * CONC(I1)*   &
                     aimag(DEN(L,IELAST,1,IPOT1+1))/PI/DBLE(NSPINPOT)
                  DENEFAT(I1) = DENEFAT(I1) - 2.0D0*  &
                     aimag(DEN(L,IELAST,1,IPOT1+1))/PI/DBLE(NSPINPOT)
               end do
            end if

            if ( TEST('RHOVALW ') ) then !Bauer
               open(unit=324234,file='out_rhoval')
               write(324234,*) '#IATOM',I1
               write(324234,'(50000F14.7)') RHO2NS(:,:,I1,1)
               if (NSPIN==2) write(324234,'(50000F14.7)') RHO2NS(:,:,I1,2)
            end if

            !-----------------------------------------------------------------------
            ! itermdir/kmrot <> 0
            !-----------------------------------------------------------------------
            if (LMOMVEC) then
               do IS = 1, NMVECMAX
                  do LM=1, 3
                     MVEVI(I1,LM,IS) = (0.0D0,0.0D0)
                     MVEVIEF(I1,LM,IS) = (0.0D0,0.0D0)
                     do L = 0, LMAX
                        MVEVIL(L,I1,LM,IS) = MVEVIL1(L,LM,IS)
                        MVEVI(I1,LM,IS) = MVEVI(I1,LM,IS)+MVEVIL1(L,LM,IS)
                        !
                        MVEVIEF(I1,LM,IS) = MVEVIEF(I1,LM,IS)+MVEVIL2(L,LM,IS)
                     end do
                  end do
               end do
            end if
         end do
         !----------------------------------------------------------------------
         ! NATYP
         !----------------------------------------------------------------------
#ifdef CPP_TIMING
         !     call timing_stop('main1c - rhoval')
         if(i1_end>=i1_start .and. t_inc%i_time>0) call timing_stop('main1c - rhoval')
#endif
         close (69) !gmat file
#ifndef CPP_MPI
         close (30) !close lmdos file if no mpi is used
         close (31) !close qdos file if no mpi is used, otherwise writeout in the following
#endif
         close (96) !close gflle file

#ifdef CPP_MPI
         !move writeout of qdos file here                                           ! qdos
         if (OPT('qdos    ')) then                                                  ! qdos
            ! first communicate den array to write out qdos files                   ! qdos
            IDIM = (LMAXD1+1)*IELAST*NQDOS*NPOTD                                     ! qdos
            allocate(workc(0:LMAXD1,IELAST,NQDOS,NPOTD))                             ! qdos
            call memocc(i_stat,product(shape(workc))*kind(workc),'workc','main1c')  ! qdos
            workc = (0.d0, 0.d0)                                                    ! qdos
            ! qdos
            CALL MPI_ALLREDUCE(DEN,workc,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM, &         ! qdos
               MPI_COMM_WORLD,IERR)                                                 ! qdos
            CALL ZCOPY(IDIM,WORKC,1,DEN,1)                                          ! qdos
            i_all=-product(shape(workc))*kind(workc)                                ! qdos
            deallocate(workc,stat=i_stat)                                           ! qdos
            call memocc(i_stat,i_all,'workc','main1c')                              ! qdos
            ! qdos
            if(myrank==master) then                                                 ! qdos
               ! qdos
               do I1=1,NATYP                                                        ! qdos
                  ! qdos
                  do ISPIN = 1,NSPIN                                                ! qdos
                     ! qdos
                     if(NATYP.ge.100) then                                          ! qdos
                        open(31,                                                 &  ! qdos
                        FILE="qdos."//char(48+I1/100)//char(48+mod(I1/10,10))//  &  ! qdos
                        char(48+mod(I1,10))//"."//char(48+ISPIN)//".dat")           ! qdos
                     else                                                           ! qdos
                        open(31,                                                 &  ! qdos
                        FILE="qdos."//char(48+I1/10)//char(48+mod(I1,10))//"."// &  ! qdos
                        char(48+ISPIN)//".dat")                                     ! qdos
                     end if                                                         ! qdos
                     call version_print_header(31)                                  ! qdos
                     write (31,'(7(A,3X))') '#   Re(E)','Im(E)',                 &  ! qdos
                        'k_x','k_y','k_z','DEN_tot','DEN_s,p,...'                   ! qdos
                     ! qdos
                     IPOT = (I1-1) * NSPINPOT + ISPIN                               ! qdos
                     ! qdos
                     do IE=1,IELAST                                                 ! qdos
                        do IQ=1,NQDOS                                               ! qdos
                           DENTOT = CMPLX(0.D0,0.D0, kind=dp)                               ! qdos
                           do L = 0,LMAXD1                                          ! qdos
                              DENTOT = DENTOT + DEN(L,IE,IQ,IPOT)                   ! qdos
                           enddo                                                    ! qdos
                           write(31,9000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ), &
                              QVEC(3,IQ),-aimag(DENTOT)/PI,             &           ! qdos
                              (-aimag(DEN(L,IE,IQ,IPOT))/PI,L=0,LMAXD1)             ! qdos
                        end do ! IQ=1,NQDOS                                         ! qdos
                     end do ! IE=1,IELAST                                           ! qdos
                     close(31)                                                      ! qdos
                     ! qdos
                     if(test('compqdos')) then                                            ! complex qdos
                        if(NATYP.ge.100) then                                             ! complex qdos
                           open(31,                                                 &     ! complex qdos
                           FILE="cqdos."//char(48+I1/100)//char(48+mod(I1/10,10))// &     ! complex qdos
                           char(48+mod(I1,10))//"."//char(48+ISPIN)//".dat")              ! complex qdos
                        else                                                              ! complex qdos
                           open(31,                                                    &  ! complex qdos
                           FILE="cqdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//   &  ! complex qdos
                           char(48+ISPIN)//".dat")                                        ! complex qdos
                        end if                                                            ! complex qdos
                        call version_print_header(31)                                     ! complex qdos
                        write (31,'(A)') '#   lmax, natyp, nspin, nqdos, ielast:'         ! complex qdos
                        write (31,'(5I9)') lmax, natyp, nspin, nqdos, ielast              ! complex qdos
                        write (31,'(7(A,3X))') '#   Re(E)','Im(E)', &                     ! complex qdos
                           'k_x','k_y','k_z','DEN_tot','DEN_s,p,...'                      ! complex qdos
                        ! complex qdos
                        IPOT = (I1-1) * NSPINPOT + ISPIN                                  ! complex qdos
                        ! complex qdos
                        do IE=1,IELAST                                                    ! complex qdos
                           do IQ=1,NQDOS                                                  ! complex qdos
                              DENTOT = CMPLX(0.D0,0.D0, kind=dp)                                  ! complex qdos
                              do L = 0,LMAXD1                                             ! complex qdos
                                 DEN(L,IE,IQ,IPOT) = -2.0d0/pi*DEN(L,IE,IQ,IPOT)          ! complex qdos
                                 DENTOT = DENTOT + DEN(L,IE,IQ,IPOT)                      ! complex qdos
                              enddo                                                       ! complex qdos
                              write(31,9002) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),    &           ! complex qdos
                                 QVEC(3,IQ),DENTOT,(DEN(L,IE,IQ,IPOT),L=0,LMAXD1)         ! complex qdos
                           end do ! IQ=1,NQDOS                                            ! complex qdos
                        end do ! IE=1,IELAST                                              ! complex qdos
                        close(31)                                                         ! complex qdos
                     end if                                                         ! qdos
                     ! qdos
                  end do !ISPIN=1,NSPIN                                             ! qdos
               enddo !I1                                                            ! qdos
            end if ! myrank_at==master                                              ! qdos
         endif ! OPT('qdos    ')                                                    ! qdos
         9000 FORMAT(5F10.6,40E16.8)                                                ! qdos
         9002 FORMAT(5F10.6,80E16.8)                                                      ! complex qdos

         !reset NQDOS number to avoid loo large communication which is not needed anyways for qdos run
         NQDOS=1

#ifdef CPP_TIMING
         call timing_start('main1c - communication')
#endif
         !MPI: reduce these arrays, so that master processor has the results (MPI_REDUCE)
#ifdef CPP_TIMING
         call timing_start('main1c - communication')
#endif
         call mympi_main1c_comm(IRMD,LMPOTD,NATYP,LMAX,LMAXD1,LMMAXD, &
            NPOTD,IELAST,MMAXD,IDOLDAU,NATYP,KREL,LMOMVEC,           &
            NMVECMAX,NQDOS,rho2ns,r2nef,espv,den,denlm,denmatc,     &
            denef,denefat,rhoorb,muorb,mvevi,mvevil,mvevief,        &
            t_mpi_c_grid%mympi_comm_ie)
         call mympi_main1c_comm(IRMD,LMPOTD,NATYP,LMAX,LMAXD1,LMMAXD, &
            NPOTD,IELAST,MMAXD,IDOLDAU,NATYP,KREL,LMOMVEC,           &
            NMVECMAX,NQDOS,rho2ns,r2nef,espv,den,denlm,denmatc,     &
            denef,denefat,rhoorb,muorb,mvevi,mvevil,mvevief,        &
            t_mpi_c_grid%mympi_comm_at)
#ifdef CPP_TIMING
         call timing_stop('main1c - communication')
#endif

         !     lmdos writeout
         if(myrank==master) then                                                 ! lm-dos
            !        IF (.not.OPT('qdos    ')) THEN                              ! lm-dos
            if (OPT('lmdos    ')) then                                           ! lm-dos
               do I1=1,NATYP                                                     ! lm-dos
                  do ISPIN = 1,NSPIN                                             ! lm-dos
                     IPOT = (I1-1) * NSPINPOT + ISPIN                            ! lm-dos
                     if(NATYP.ge.100) then                                       ! lm-dos
                        open(30,FILE="lmdos."//char(48+I1/100)//           &     ! lm-dos
                        char(48+mod(I1/10,10))//char(48+mod(I1,10))//"."// &     ! lm-dos
                        char(48+ISPIN)//".dat")                                  ! lm-dos
                     else                                                        ! lm-dos
                        open(30,FILE="lmdos."//char(48+I1/10)//            &     ! lm-dos
                        char(48+mod(I1,10))//"."//char(48+ISPIN)//".dat")        ! lm-dos
                     end if                                                      ! lm-dos
                     call version_print_header(30)
                     write (30,*) ' '                                            ! lm-dos
                     write (30,8600) '# ISPIN=',ISPIN,' I1=',I1                  ! lm-dos
                     8600         FORMAT (a8,I3,a4,I5)                           ! lm-dos
                     do IE=1,IELAST                                              ! lm-dos
                        write(30,9001) REAL(EZ(IE), kind=dp), &                          ! lm-dos
                           (-aimag(DENLM(LM,IE,1,IPOT))/PI,LM=1,LMMAXD)          ! lm-dos
                     end do ! IE                                                 ! lm-dos
                  end do ! ISPIN                                                 ! lm-dos
                  9001       FORMAT(30E12.4)                                     ! lm-dos
                  close(30)                                                      ! lm-dos
               end do !I1                                                        ! lm-dos
            endif  ! not qdos option                                             ! lm-dos
         end if ! myrank==master                                                 ! lm-dos
#endif

      else ! new spin-orbit solver

         ! nonco angles
         call read_angles(t_params,NATYP,THETA,PHI)

         ! interpolate potential
         if ( IDOLDAU.EQ.1 ) then
            call CINIT(MMAXD*MMAXD*4*NATYP,DENMATN(1,1,1,1,1))
         endif

         call INTERPOLATE_POTEN(LPOTD,IRMD,IRNSD,NATYP,IPAND,       &
            LMPOTD,NSPOTD,NTOTD,NTOTD*(NCHEB+1),NSPIN,RMESH, &
            IRMIN,IRWS,IRCUT,VINS,VISP,NPAN_LOG_AT,NPAN_EQ_AT,   &
            NPAN_TOT,RNEW,IPAN_INTERVALL,VINSNEW)


#ifdef CPP_MPI
         ntot1 = t_inc%NATYP

         if(t_mpi_c_grid%dims(1)>1) then
            nranks_local = t_mpi_c_grid%nranks_ie
            if(t_mpi_c_grid%nranks_ie>t_mpi_c_grid%dims(1)) then
               nranks_local = t_mpi_c_grid%dims(1)
            endif
         else
            nranks_local = 1
         endif
         call distribute_linear_on_tasks(nranks_local,      &
            t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,  &
            master,ntot1,ntot_pT,ioff_pT,.true.,.true.)
         if(t_mpi_c_grid%nranks_ie<=t_mpi_c_grid%dims(1)) then
            i1_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
            i1_end   = ioff_pT(t_mpi_c_grid%myrank_ie)+ntot_pT(t_mpi_c_grid%myrank_ie)
            t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie)
         else
            myrank_ie_tmp = t_mpi_c_grid%myrank_ie
            if(t_mpi_c_grid%myrank_ie>(t_mpi_c_grid%dims(1)-1))then
               myrank_ie_tmp = myrank_ie_tmp - t_mpi_c_grid%myrank_ie &
               + (t_mpi_c_grid%dims(1)-1)
            endif
            i1_start = ioff_pT(myrank_ie_tmp)+1
            i1_end   = ioff_pT(myrank_ie_tmp)+ntot_pT(myrank_ie_tmp)
            t_mpi_c_grid%ntot1  = ntot_pT(myrank_ie_tmp)
         endif

         if(.not.(allocated(t_mpi_c_grid%ntot_pT1).and.allocated(t_mpi_c_grid%ioff_pT1))) then
            allocate(t_mpi_c_grid%ntot_pT1(0:t_mpi_c_grid%nranks_ie-1),stat=i_stat)
            call memocc(i_stat,product(shape(t_mpi_c_grid%ntot_pT1))*kind(t_mpi_c_grid%ntot_pT1),'t_mpi_c_grid%ntot_pT1','main1c')
            allocate(t_mpi_c_grid%ioff_pT1(0:t_mpi_c_grid%nranks_ie-1),stat=i_stat)
            call memocc(i_stat,product(shape(t_mpi_c_grid%ioff_pT1))*kind(t_mpi_c_grid%ioff_pT1),'t_mpi_c_grid%ioff_pT1','main1c')
         endif
         t_mpi_c_grid%ntot_pT1 = ntot_pT
         t_mpi_c_grid%ioff_pT1 = ioff_pT
#else
         i1_start = 1
         i1_end   = NATYP
#endif

         do I1 = i1_start,i1_end

            ICELL = NTCELL(I1)
            IPOT = (I1-1) * NSPIN + 1

#ifdef CPP_TIMING
            call timing_start('main1c - rhovalnew')
#endif

            call RHOVALNEW( &
               LDORHOEF,IELAST,NSRA,NSPIN,LMAX,EZ,WEZ,     &
               ZAT(I1),SOCSCALE(I1),CLEB(1,1),ICLEB,IEND,               &
               IFUNM1(1,ICELL),LMSP1(1,ICELL),                          &
               NCHEB,NPAN_TOT(I1),NPAN_LOG_AT(I1),                         &
               NPAN_EQ_AT(I1),RMESH(1,I1),IRWS(I1),RPAN_INTERVALL(0,I1),   &
               IPAN_INTERVALL(0,I1),RNEW(1,I1),VINSNEW,                 &
               THETASNEW(1,1,ICELL),THETA(I1),PHI(I1),I1,IPOT,          &
               DEN1(0,1,1),ESPV1(0,1),RHO2M1,RHO2M2,MUORB(0,1,I1),      &
               angles_new(:,i1),                                        &
               IDOLDAU,LOPT(I1),WLDAU(1,1,1,I1),          &  ! LDAU
               DENMATN(1,1,1,1,I1),NATYP)                                  ! LDAU

#ifdef CPP_TIMING
            call timing_pause('main1c - rhovalnew')
#endif

            do L = 0,LMAXD1
               ESPV(L,IPOT)=ESPV1(L,1)
               ESPV(L,IPOT+NSPIN-1)=ESPV1(L,NSPIN)
               do IE=1,IELAST
                  DEN(L,IE,1,IPOT)=DEN1(L,IE,1)
                  DEN(L,IE,1,IPOT+NSPIN-1)=DEN1(L,IE,NSPIN)
               enddo
            enddo
            do ISPIN=1,NSPIN
               do LM = 1,LMPOTD
                  do IR = 1,IRMD
                     RHO2NS(IR,LM,I1,ISPIN) = RHO2M1(IR,LM,ISPIN)
                     R2NEF(IR,LM,I1,ISPIN)  = RHO2M2(IR,LM,ISPIN)
                  end do
               end do
            end do

            do L = 0,LMAXD1
               DENEF = DENEF - 2.0D0 * CONC(I1)*   &
                  aimag(DEN(L,IELAST,1,IPOT))/PI/DBLE(NSPIN)
               DENEFAT(I1) = DENEFAT(I1) - 2.0D0*  &
                  aimag(DEN(L,IELAST,1,IPOT))/PI/DBLE(NSPIN)
            end do
            do L = 0,LMAXD1
               DENEF = DENEF - 2.0D0 * CONC(I1)*   &
                  aimag(DEN(L,IELAST,1,IPOT+1))/PI/DBLE(NSPIN)
               DENEFAT(I1) = DENEFAT(I1) - 2.0D0*  &
                  aimag(DEN(L,IELAST,1,IPOT+1))/PI/DBLE(NSPIN)
            end do
            do ISPIN=1,NSPIN
               do L=0,LMAXD1
                  MUORB(L,3,I1)=MUORB(L,3,I1)+MUORB(L,ISPIN,I1)
               enddo
            enddo
            do ISPIN=1,3
               do L=0,LMAXD1
                  MUORB(LMAXD1+1,ISPIN,I1)=MUORB(LMAXD1+1,ISPIN,I1)+MUORB(L,ISPIN,I1)
               enddo
            enddo
         end do !I1

#ifdef CPP_MPI
         !reset NQDOS to 1 to avoid endless communication
         NQDOS = 1
         call mympi_main1c_comm_newsosol2(LMAXD1,LMMAXD,IELAST,NQDOS,    &
            NPOTD,NATYP,LMPOTD,IRMD,MMAXD,den, denlm, muorb, espv, r2nef, &
            rho2ns, denefat, denef,denmatn,angles_new,t_mpi_c_grid%mympi_comm_ie)
#endif

#ifdef CPP_TIMING
         !call timing_stop('main1c - rhovalnew')
         if(i1_end>=i1_start) call timing_stop('main1c - rhovalnew')
#endif

         close (69)

         if(myrank==master) then
            ! rewrite new theta and phi to nonco_angle_out.dat, nonco_angle.dat is the input
            if (.not.TEST('FIXMOM  ')) then
               open(UNIT=13,file='nonco_angle_out.dat',form='formatted')
               call version_print_header(13)
               do I1=1,NATYP
                  ! save to file in converted units (degrees)
                  write(13,*) angles_new(1,I1)/(2.0D0*PI)*360.0D0,   &
                     angles_new(2,I1)/(2.0D0*PI)*360.0D0
                  ! use internal units here
                  t_params%THETA(I1) = angles_new(1,I1)
                  t_params%PHI(I1)   = angles_new(2,I1)
               enddo
               close(13)

            endif ! .not.test('FIXMOM  ')
         end if !(myrank==master)

         close (29) !close lmdos file
         close (30) !close lmdos file
         close (31) !close qdos file
         close (32) !close qdos file
         close (91) !close gflle file

      endif !NEWSOSOL

#ifdef CPP_MPI
      if(myrank==master) then
#endif
#ifdef CPP_TIMING
         call timing_start('main1c - serial part')
#endif

         ! In case of Lloyds formula renormalize valence charge                  ! LLY Lloyd
         if (LLY.GT.0) then                                                      ! LLY Lloyd
            LMAXP1 = LMAX                                                        ! LLY Lloyd
            if (INS.ne.0) LMAXP1 = LMAX + 1                                      ! LLY Lloyd

            call RENORM_LLY(CDOS_LLY,IELAST,NSPIN,NATYP,DEN(:,:,1,:),LMAXP1,CONC,   &
               1,IELAST,WEZ,IRCUT,IPAN,EZ,ZAT,RHO2NS,R2NEF,DENEF,DENEFAT,ESPV)

         endif                                                                   ! LLY Lloyd

         !----------------------------------------------------------------------
         ! NATYP
         !----------------------------------------------------------------------
         CHRGSEMICORE = 0D0
         do I1 = 1,NATYP
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! l/m_s/atom-resolved charges
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            do ISPIN = 1,NSPINPOT
               IPOT = (I1-1)*NSPINPOT + ISPIN
               do L = 0,LMAXD1
                  CHARGE(L,I1,ISPIN) = 0.0D0
                  !
                  do IE = 1,IELAST
                     CHARGE(L,I1,ISPIN) = CHARGE(L,I1,ISPIN) + &
                        aimag(WEZ(IE)*DEN(L,IE,1,IPOT))/DBLE(NSPINPOT)
                     if ( IE.eq.IESEMICORE ) then
                        CHRGSEMICORE = CHRGSEMICORE + CONC(I1)*CHARGE(L,I1,ISPIN)
                     endif
                  end do
                  !
               end do
            end do
            EU(I1) = 0D0
            EDC(I1) = 0D0
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Orbital magnetic moments (array initialised to 0.0D0 in rhoval)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (KREL.eq.1) then
               do ISPIN = 1,3
                  do L = 0,LMAX + 1
                     MUORB(LMAXD1+1,ISPIN,I1) = MUORB(LMAXD1+1,ISPIN,I1) + MUORB(L,ISPIN,I1)
                  end do
               end do
            end if
         end do
         !----------------------------------------------------------------------
         ! NATYP
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         ! LDA+U
         !----------------------------------------------------------------------
         if ( IDOLDAU.eq.1 ) then
            !-------------------------------------------------------------------
            ! Save old LDA+U interaction matrix for mixing
            !-------------------------------------------------------------------
            call DCOPY(MMAXD*MMAXD*NSPIND*NATYP,WLDAU,1,WLDAUOLD,1)
            !-------------------------------------------------------------------
            ! Construct LDA+U interaction matrix for next iteration
            !-------------------------------------------------------------------
            if (.not.OPT('NEWSOSOL')) then
               call WMATLDAU(NTLDAU,ITLDAU,NSPINPOT,DENMATC,LOPT,UEFF,JEFF,   &
                  ULDAU,WLDAU,EU,EDC,MMAXD,NPOTD,NATYP,NSPIN,LMAX)
            else
               call WMATLDAUSOC(NTLDAU,ITLDAU,NSPINPOT,DENMATN,LOPT,UEFF,JEFF,&
                  ULDAU,WLDAU,EU,EDC,MMAXD,NATYP, NSPIN, LMAX)
            endif
            ! -> Mix old and new LDA+U interaction matrices
            call MIXLDAU(MMAXD,NSPIND,NATYP,NATYP,NSPIN,LOPT,WLDAUOLD,WLDAU)
            !-------------------------------------------------------------------
            ! Update variables-file
            !-------------------------------------------------------------------
            ITRUNLDAU = ITRUNLDAU + 1
            open (67,FILE='ldau.unformatted',FORM='unformatted')
            write (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
            close(67)
            !-------------------------------------------------------------------
            ! Write full lda+u information in ascii file ldaupot_new
            !-------------------------------------------------------------------
            call WRLDAUPOT(ITRUNLDAU,LOPT,UEFF,JEFF,EREFLDAU,NATYP,WLDAU,  &
               ULDAU,PHILDAU,IRMD,NATYP,NSPIND,MMAXD,IRWS)
         end if
         !----------------------------------------------------------------------
         ! LDA+U
         !----------------------------------------------------------------------

         call WRMOMS(KREL+KORBIT,NATYP,NSPINPOT,TEXTS,TEXTL,TEXTNS,CHARGE, &
            MUORB,LMAX,LMAXD1)
         !----------------------------------------------------------------------
         ! ITERMDIR
         !----------------------------------------------------------------------
         if ( ( KREL.eq.1 ) .and. LMOMVEC ) then
            do I1=1,NATYP
               IQ = IQAT(I1)
               call MVECGLOBAL(I1,IQ,NATYP,QMPHI(IQ),QMTET(IQ),MVEVI,MVEVIL,  &
                  MVEVIEF,NATYP,LMAX,NMVECMAX)
            end do
         end if
         !----------------------------------------------------------------------
         ! TERMDIR
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         ! TEST BRAHIM
         !----------------------------------------------------------------------
         if (NPOL.eq.0 .or.TEST('DOS     '))  then
            call WRLDOS(DEN,EZ,WEZ,LMAXD1,IEMXD,NPOTD,ITITLE,EFERMI,EMIN,EMAX,   &
               ALAT,TK,NACLS1,NSPINPOT,NATYP,CONC,IELAST,INTERVX,INTERVY,INTERVZ,&
               DOSTOT)
         endif

         !----------------------------------------------------------------------
         ! CORE STATES
         !----------------------------------------------------------------------
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! RHO_core is calculated only if also RHO_valence was
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (NPOL.ne.0) then
            if(t_inc%i_write>0) then
               write (1337,*)
               write (1337,'(78("#"))')
               write (1337,'(33X,A)') 'CORE  STATES'
               write (1337,'(78("#"))')
            endif
            do I1 = 1,NATYP
               do ISPIN = 1,NSPIN
                  IPOT = (I1-1) * NSPINPOT + ISPIN
                  IPOT1 = (I1-1) * NSPINPOT + 1
                  !
                  call RHOCORE(NSRA,ISPIN,NSPIN,I1,DRDI(1,I1),RMESH(1,I1), &
                     VISP(1,IPOT),A(I1),B(I1),ZAT(I1),IRCUT(0,I1),         &
                     RHOC(1,IPOT1),ECORE(1,IPOT),NCORE(IPOT),LCORE(1,IPOT),&
                     CSCL(1,KREL*I1+(1-KREL)),VTREL(1,I1),BTREL(1,I1),     &
                     RMREL(1,I1),DRDIREL(1,I1),R2DRDIREL(1,I1),ZREL(I1),   &
                     JWSREL(I1),IRSHIFT(I1),ECOREREL(1,IPOT1),NKCORE(1,I1),&
                     KAPCORE(1,IPOT1))
                  !
               end do
            end do
            if(t_inc%i_write>0) then
               write (1337,*)
               write (1337,'(78("#"))')
               write (1337,*)
            endif
         end if
         !----------------------------------------------------------------------
         ! CORE STATES
         !----------------------------------------------------------------------
         call save_density(t_params,RHO2NS,R2NEF,RHOC,DENEF,DENEFAT,ESPV,  &
            ECORE,IDOLDAU,LOPT,EU,EDC,CHRGSEMICORE,RHOORB,ECOREREL,NKCORE, &
            KAPCORE,KREL,NATYP,NPOTD,IRMD,LMPOTD,LMAXD1)

         if (TEST('den-asci')) then
            open (67,FILE='densitydn.ascii',FORM='formatted')
            do I1 = 1,NATYP
               do LM = 1,LMPOTD
                  do IR = 1,IRMD
                     write(67,FMT='(I6,2I5,2E25.16)') &
                        I1,LM,IR,RHO2NS(IR,LM,I1,1),RHO2NS(IR,LM,I1,2)
                  enddo
               enddo
            enddo
            close(67)
         endif

         !----------------------------------------------------------------------
         ! TERMDIR
         !----------------------------------------------------------------------
         if (ITERMVDIR) then
            t_params%MVEVI   = MVEVI
            t_params%MVEVIEF = MVEVIEF
         end if
         !----------------------------------------------------------------------
         ! TERMDIR
         !----------------------------------------------------------------------

         !
         if(t_inc%i_write>0) then
            write (1337,'(79("="),/,30X,"< KKR1c finished >",/,79("="),/)')
         endif

#ifdef CPP_TIMING
         call timing_stop('main1c - serial part')
#endif
#ifdef CPP_MPI
      end if !myrank==master
#endif


      if (lly/=0) then
        ! cleanup allocations
        deallocate(CDOS0, stat=i_stat)
        call memocc(i_stat,-product(shape(CDOS0))*kind(CDOS0),'CDOS0','main1c')
        deallocate(CDOS1, stat=i_stat)
        call memocc(i_stat,-product(shape(CDOS1))*kind(CDOS1),'CDOS1','main1c')
        deallocate(CDOS2, stat=i_stat)
        call memocc(i_stat,-product(shape(CDOS2))*kind(CDOS2),'CDOS2','main1c')
        deallocate(CDOSAT0, stat=i_stat)
        call memocc(i_stat,-product(shape(CDOSAT0))*kind(CDOSAT0),'CDOSAT0','main1c')
        deallocate(CDOSAT1, stat=i_stat)            
        call memocc(i_stat,-product(shape(CDOSAT1))*kind(CDOSAT1),'CDOSAT1','main1c')
        deallocate(CDOS_LLY, stat=i_stat)            
        call memocc(i_stat,-product(shape(CDOS_LLY))*kind(CDOS_LLY),'CDOS_LLY','main1c')
      end if

   end subroutine main1c

end module MOD_MAIN1C
