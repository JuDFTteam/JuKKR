module mod_types

implicit none


   type :: type_tgmatices

       ! logical switches to control if matrices are stored in memory or written to files
      logical :: tmat_to_file = .false.
      logical :: gmat_to_file = .false.
      logical :: gref_to_file = .false.
      
      integer :: Nelements = 4 ! 3 arrays in this type, for mpi bcast
      
      ! allocatable arrays for tmat, gmat and gref
      double complex, allocatable :: tmat(:,:,:)       ! dimensions=LMMAXD, LMMAXD, IREC; IREC= IE+IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1) ;IE=1,...,IELAST, ISPIN=1,...,NSPIN, I1=1,...,NATYP)
      double complex, allocatable :: gmat(:,:,:)       ! dimensions=LMMAXD, LMMAXD, IREC; IREC= IQDOS+NQDOS*(IE-1)+NQDOS*IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1) ;IE=1,...,IELAST, ISPIN=1,...,NSPIN, I1=1,...,NATYP)
      double complex, allocatable :: gref(:,:,:,:)       !GINP(NACLSD*LMGF0D,LMGF0D,NCLSD) IREC=IE=1,...,IELAST

   end type type_tgmatices
   

   type :: type_cpa

       ! logical switches to control if matrices are stored in memory or written to files
      logical :: dmatproj_to_file = .false.

      integer :: Nelements = 3 ! 2 array in this type, for mpi bcast
      
      ! allocatable arrays for tmat, gmat and gref
      double complex, allocatable :: dmatts(:,:,:,:)       ! dimensions=LMMAXD, LMMAXD, NATYP, IREC; IREC= IE+IELAST*(ISPIN-1)+; IE=1,...,IELAST, ISPIN=1,...,NSPIN)
      double complex, allocatable :: dtilts(:,:,:,:)       ! dimensions=LMMAXD, LMMAXD, NATYP, IREC; IREC= IE+IELAST*(ISPIN-1)+; IE=1,...,IELAST, ISPIN=1,...,NSPIN)
   end type type_cpa

   !data type for the derivatives of the t-matrix with respect to changing the non-collinear angles in directions {x,y,z}
   type :: type_dtmatJijDij
     
      integer :: Nelements = 3
      logical :: calculate = .false.
      double complex, allocatable :: dtmat_xyz(:,:,:,:) !dimensions= LMMAXD, LMMAXD, 3, IELAST;  3={x,y,z}

   end type type_dtmatJijDij


   type :: type_inc
   
      integer :: Nparams = 18   ! number of parameters in type_inc, excluding allocatable array KMESH
      integer :: LMMAXD  = -1
      integer :: NSPIN   = -1
      integer :: IELAST  = -1
      integer :: NQDOS   = -1
      integer :: NATYP   = -1
      integer :: LMGF0D  = -1
      integer :: NCLSD   = -1
      integer :: NACLSD  = -1
      integer :: i_iteration = -1 
      integer :: N_iteration = -1
      integer :: mit_bry = 1
      integer :: NSHELL0 = -1
      integer :: NKMESH = -1
      logical :: NEWSOSOL = .false.  ! use new solver for SOC
      logical :: deci_out = .false.  ! use deci_out case
      integer :: i_write = 0 ! switch to control if things are written out or not (verbosity levels 0,1,2)
      integer :: i_time  = 1 ! switch to control if timing files are written (verbosity levels 0,1,2)
      integer, allocatable :: KMESH(:)
         
   end type type_inc
   
   
   type :: type_mpi_cartesian_grid_info
   
      integer :: Nparams = 12
      integer :: dims(2) = (/ -1, -1 /)
      integer :: myMPI_comm_ie  = -1
      integer :: myMPI_comm_at  = -1
      integer :: myrank_ie  = -1
      integer :: myrank_at  = -1
      integer :: myrank_atcomm  = -1
      integer :: nranks_ie  = -1
      integer :: nranks_at  = -1
      integer :: nranks_atcomm  = -1
      integer :: ntot1    = -1
      integer, allocatable :: ntot_pT1(:)
      integer, allocatable :: ioff_pT1(:)
      integer :: ntot2    = -1
      integer, allocatable :: ntot_pT2(:)
      integer, allocatable :: ioff_pT2(:)
   
   end type type_mpi_cartesian_grid_info


   type :: type_lloyd

       ! logical switches to control if matrices are stored in memory or written to files
      logical :: dtmat_to_file = .false.          ! unit 691
      logical :: tralpha_to_file = .false.        ! unit 692
      logical :: cdos_diff_lly_to_file = .false.  ! unit 701
      logical :: dgref_to_file = .false.          ! unit 681
      logical :: g0tr_to_file = .false.           ! unit 682
      
      integer :: N1 = 6 ! 5 logicals and 5 arrays this type, for mpi bcast
      
      ! allocatable arrays
      double complex, allocatable :: dtmat(:,:,:)      ! DOUBLE COMPLEX TMAT0(LMMAXD,LMMAXD), IREC = ie_num + ie_end*(ISPIN-1) + ie_end*NSPIN* (I1-1)
      double complex, allocatable :: tralpha(:)        ! DOUBLE COMPLEX TRALPHA, IREC = ie_num + ie_end*(ISPIN-1) + ie_end*NSPIN* (I1-1)
      double complex, allocatable :: cdos(:,:)         ! DOUBLE COMPLEX CDOS_LLY(IEMXD,NSPIND), irec=IE, aalready in dim 1 of cdos!
      double complex, allocatable :: dgref(:,:,:,:)    ! DOUBLE COMPLEX; ALLOCATE ( DGINP(NACLSMAX*LMGF0D,LMGF0D,NCLS) ), IREC=IE 
      double complex, allocatable :: g0tr(:)           ! DOUBLE COMPLEX LLY_G0TR_IE, irec=ie

   end type type_lloyd


   type (type_inc), save :: t_inc
   type (type_tgmatices), save :: t_tgmat
   type (type_mpi_cartesian_grid_info), save :: t_mpi_c_grid
   type (type_lloyd), save :: t_lloyd
   type (type_dtmatJijDij), allocatable, save :: t_dtmatJij(:) !dimensions I1=1,...,NATYP 
   type (type_cpa), save :: t_cpa

contains

   subroutine init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
   
      use mod_mympi, only: myrank, master, nranks
   
      implicit none
   
      type(type_inc), intent(in) :: t_inc
      type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
      type(type_tgmatices), intent(inout) :: t_tgmat
      
      integer :: ierr, nspin
      
      nspin = t_inc%NSPIN
      if(t_inc%NEWSOSOL) nspin = 1
    
      if (.not. allocated(t_tgmat%tmat)) then
         if (.not. t_tgmat%tmat_to_file) then
            if(nranks.eq.1) then
               !allocate tmat(lmmax,lmmax,irec_max) for irec_max=ielast*nspin*natyp
               allocate(t_tgmat%tmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%IELAST*nspin*t_inc%NATYP), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%tmat'
            else
               !allocate tmat(lmmax,lmmax,irec_max) for irec_max=iemax_local*nspin*natyp
               allocate(t_tgmat%tmat(t_inc%LMMAXD,t_inc%LMMAXD,t_mpi_c_grid%ntot2*nspin*t_inc%NATYP), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%tmat for mpi'
            end if
         else
            allocate(t_tgmat%tmat(1,1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_tgmat%tmat'
         end if
      
         t_tgmat%tmat(:,:,:) = (0.d0, 0.d0)
      endif
      
      if (.not. allocated(t_tgmat%gmat)) then
         if (.not. t_tgmat%gmat_to_file) then
            if(nranks.eq.1) then
               !allocate gmat(lmmax,lmmax,irec_max) for irec_max=NQDOS*IELAST*NSPIN*NATYP
               allocate(t_tgmat%gmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NQDOS*t_inc%IELAST*nspin*t_inc%NSHELL0), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%gmat'
            else
               !allocate gmat(lmmax,lmmax,irec_max) for irec_max=NQDOS*IEMAX_local*NSPIN*NATYP
               allocate(t_tgmat%gmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NQDOS*t_mpi_c_grid%ntot2*nspin*t_inc%NSHELL0), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%gmat for mpi'
            end if
            
         else
            allocate(t_tgmat%gmat(1,1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_tgmat%gmat'
         end if
      
         t_tgmat%gmat(:,:,:) = (0.d0, 0.d0)
      endif
      
      if (.not. allocated(t_tgmat%gref)) then
         if (.not. t_tgmat%gref_to_file) then
            if(nranks.eq.1) then
               !allocate gref(NACLSD*LMGF0D,LMGF0D,NCLSD,irec_max) for irec_max=IELAST
               allocate(t_tgmat%gref(t_inc%NACLSD*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_inc%IELAST), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%gref'
            else
               !allocate gref(NACLSD*LMGF0D,LMGF0D,NCLSD,irec_max) for irec_max=IEMAX_local (=ntot2)
               allocate(t_tgmat%gref(t_inc%NACLSD*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_mpi_c_grid%ntot2), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%gref for mpi'
            end if
         else
            allocate(t_tgmat%gref(1,1,1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_tgmat%gref'
         end if
   
      t_tgmat%gref(:,:,:,:) = (0.d0, 0.d0)
      endif

   end subroutine init_tgmat
   

   subroutine init_t_cpa(t_inc,t_cpa,t_mpi_c_grid)
   
      use mod_mympi, only: myrank, master, nranks
   
      implicit none
   
      type(type_inc), intent(in) :: t_inc
      type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
      type(type_cpa), intent(inout) :: t_cpa
      
      integer :: ierr, nspin, nenergy
      
      nspin = t_inc%NSPIN
      if(t_inc%NEWSOSOL) nspin = 1

      if(nranks==1)then
        nenergy=t_inc%IELAST
      else
        nenergy=t_mpi_c_grid%ntot2
      end if
    
      if (.not. allocated(t_cpa%dmatts)) then
         if (.not. t_cpa%dmatproj_to_file) then
           !allocate tmat(lmmax,lmmax,NATYP,irec_max) for irec_max=nenergy*nspin
           allocate(t_cpa%dmatts(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NATYP,nenergy*nspin), STAT=ierr)
           if(ierr/=0) stop 'Problem allocating t_cpa%dmatts'
         else
            allocate(t_cpa%dtilts(1,1,1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_cpa%dmatts'
         end if
      
         t_cpa%dmatts(:,:,:,:) = (0.d0, 0.d0)
      endif

      if (.not. allocated(t_cpa%dtilts)) then
         if (.not. t_cpa%dmatproj_to_file) then
           !allocate tmat(lmmax,lmmax,NATYP,irec_max) for irec_max=nenergy*nspin
           allocate(t_cpa%dtilts(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NATYP,nenergy*nspin), STAT=ierr)
           if(ierr/=0) stop 'Problem allocating t_cpa%dtilts'
         else
            allocate(t_cpa%dtilts(1,1,1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_cpa%dtilts'
         end if
      
         t_cpa%dtilts(:,:,:,:) = (0.d0, 0.d0)
      endif

   end subroutine init_t_cpa 


   subroutine init_t_dtmatJij(t_inc,t_mpi_c_grid,t_dtmatJij)

      implicit none

      type(type_inc), intent(in) :: t_inc
      type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
      type(type_dtmatJijDij), intent(inout), allocatable :: t_dtmatJij(:)

      integer :: ierr
     
      if(.not.t_inc%NEWSOSOL) stop 'in init_t_dtmatJij: should only be called with NEWSOSOL'

      if (.not. allocated(t_dtmatJij)) then
        allocate(t_dtmatJij(t_inc%NATYP),STAT=ierr)
        if(ierr/=0) stop 'Problem allocating t_dtmatJij%dtmat(NATYP)'
      endif

   end subroutine init_t_dtmatJij

   
   subroutine init_t_dtmatJij_at(t_inc,t_mpi_c_grid,t_dtmatJij_at)

      use mod_mympi, only: myrank, master, nranks

      implicit none

      type(type_inc), intent(in) :: t_inc
      type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
      type(type_dtmatJijDij), intent(inout) :: t_dtmatJij_at

      integer :: ierr

      if(.not.t_inc%NEWSOSOL) stop 'in init_t_dtmatJij_single: should only be called with NEWSOSOL'

      if (.not. allocated(t_dtmatJij_at%dtmat_xyz) .and. t_dtmatJij_at%calculate) then
!        if (.not. t_dtmatJij%dtmat_to_file) then
            if(nranks.eq.1) then
               !allocate dtmat_xyz(lmmax,lmmax,3,irec) for irec_max=ielast
               allocate(t_dtmatJij_at%dtmat_xyz(t_inc%LMMAXD,t_inc%LMMAXD,3,t_inc%IELAST), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_dtmatJij%dtmat'
            else
               !allocate dtmat_xyz(lmmax,lmmax,3,irec) for irec_max=iemax_local
               allocate(t_dtmatJij_at%dtmat_xyz(t_inc%LMMAXD,t_inc%LMMAXD,3,t_mpi_c_grid%ntot2), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_dtmatJij%dtmat for mpi'
            end if
!        else
!           allocate(t_tgmat%tmat(1,1,1), STAT=ierr)
!           if(ierr/=0) stop 'Problem allocating dummy t_tgmat%tmat'
!        end if
      
         t_dtmatJij_at%dtmat_xyz(:,:,:,:) = (0.d0, 0.d0)
      endif

   end subroutine init_t_dtmatJij_at
   



#ifdef CPP_MPI
   subroutine bcast_t_inc_tgmat(t_inc,t_tgmat,t_cpa)
    !ruess: after myBcast_impcls from Pkkr_sidebranch2D_2014_12_16 by Bernd Zimmermann

    use mpi
    use mod_mympi,   only: myrank, nranks, master
    implicit none

    type(type_inc), intent(inout) :: t_inc
    type(type_tgmatices), intent(inout) :: t_tgmat
    type(type_cpa), intent(inout) :: t_cpa

    integer :: blocklen1(t_inc%Nparams), etype1(t_inc%Nparams), myMPItype1 ! for parameter from t_inc
    integer :: blocklen2(t_tgmat%Nelements), etype2(t_tgmat%Nelements), myMPItype2 ! for logicals in t_tgmat
    integer :: blocklen3(t_tgmat%Nelements), etype3(t_tgmat%Nelements), myMPItype3 ! for allocatable arrays in t_tgmat
    integer :: ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp1(t_inc%Nparams), disp2(t_tgmat%Nelements), disp3(t_cpa%Nelements), base

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !broadcast parameters from t_inc
    call MPI_Get_address(t_inc%Nparams,       disp1(1), ierr)
    call MPI_Get_address(t_inc%LMMAXD,        disp1(2), ierr)
    call MPI_Get_address(t_inc%NSPIN ,        disp1(3), ierr)
    call MPI_Get_address(t_inc%IELAST,        disp1(4), ierr)
    call MPI_Get_address(t_inc%NQDOS ,        disp1(5), ierr)
    call MPI_Get_address(t_inc%NATYP ,        disp1(6), ierr)
    call MPI_Get_address(t_inc%LMGF0D,        disp1(7), ierr)
    call MPI_Get_address(t_inc%NCLSD ,        disp1(8), ierr)
    call MPI_Get_address(t_inc%NACLSD,        disp1(9), ierr)
    call MPI_Get_address(t_inc%i_iteration,  disp1(10), ierr)
    call MPI_Get_address(t_inc%N_iteration,  disp1(11), ierr)
    call MPI_Get_address(t_inc%mit_bry,      disp1(12), ierr)
    call MPI_Get_address(t_inc%NSHELL0,      disp1(13), ierr)
    call MPI_Get_address(t_inc%NKMESH,       disp1(14), ierr)
    call MPI_Get_address(t_inc%NEWSOSOL,     disp1(15), ierr)
    call MPI_Get_address(t_inc%deci_out,     disp1(16), ierr)
    call MPI_Get_address(t_inc%i_write,      disp1(17), ierr)
    call MPI_Get_address(t_inc%i_time,       disp1(18), ierr)
    base  = disp1(1)
    disp1 = disp1 - base

    blocklen1(1:18)=1

    etype1(1:18) = MPI_INTEGER
    etype1(15:16) = MPI_LOGICAL

    call MPI_Type_create_struct(t_inc%Nparams, blocklen1, disp1, etype1, myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_t_inc'

    call MPI_Type_commit(myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_t_inc'

    call MPI_Bcast(t_inc%Nparams, 1, myMPItype1, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting t_inc'

    call MPI_Type_free(myMPItype1, ierr)
    
    !broadcast allocatable array kmesh(nkmesh)
    if(.not. allocated(t_inc%kmesh)) allocate(t_inc%kmesh(t_inc%nkmesh))
    call MPI_Bcast(t_inc%KMESH, t_inc%nkmesh, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting t_inc%kmesh'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !brodcast allocatable arrays from t_tgmat
    !first broadcast logocal switches
    call MPI_Get_address(t_tgmat%Nelements,      disp2(1), ierr)
    call MPI_Get_address(t_tgmat%tmat_to_file,   disp2(2), ierr)
    call MPI_Get_address(t_tgmat%gmat_to_file,   disp2(3), ierr)
    call MPI_Get_address(t_tgmat%gref_to_file,   disp2(4), ierr)
    
    base  = disp2(1)
    disp2 = disp2 - base

    blocklen2(1:4)=1

    etype2(1) = MPI_INTEGER
    etype2(2:4) = MPI_LOGICAL

    call MPI_Type_create_struct(t_tgmat%Nelements, blocklen2, disp2, etype2, myMPItype2, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_tgmat_logicals'

    call MPI_Type_commit(myMPItype2, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_tgmat_logicals'
    
    call MPI_Bcast(t_tgmat%Nelements, 1, myMPItype2, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting logicals from t_tgmat'

    call MPI_Type_free(myMPItype2, ierr)

    call MPI_Bcast(t_cpa%dmatproj_to_file,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting logicals from t_cpa'

   end subroutine bcast_t_inc_tgmat
#endif


#ifdef CPP_MPI
   subroutine bcast_t_lly_1(t_inc,t_lloyd)

    use mpi
    use mod_mympi,   only: myrank, nranks, master
    implicit none

    type(type_inc), intent(in) :: t_inc
    type(type_lloyd), intent(inout) :: t_lloyd

    integer :: blocklen1(t_lloyd%N1), etype1(t_lloyd%N1), myMPItype1 ! for parameter from t_lloyd
    integer :: ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp1(t_lloyd%N1), base

    call MPI_Get_address(t_lloyd%N1,                     disp1(1), ierr)
    call MPI_Get_address(t_lloyd%dtmat_to_file,          disp1(2), ierr)
    call MPI_Get_address(t_lloyd%tralpha_to_file,        disp1(3), ierr)
    call MPI_Get_address(t_lloyd%cdos_diff_lly_to_file,  disp1(4), ierr)
    call MPI_Get_address(t_lloyd%dgref_to_file,          disp1(5), ierr)
    call MPI_Get_address(t_lloyd%g0tr_to_file,           disp1(6), ierr)
    
    base  = disp1(1)
    disp1 = disp1 - base

    blocklen1(1:6)=1

    etype1(1) = MPI_INTEGER
    etype1(2:6) = MPI_LOGICAL

    call MPI_Type_create_struct(t_lloyd%N1, blocklen1, disp1, etype1, myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_tgmat_logicals'

    call MPI_Type_commit(myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_tgmat_logicals'
    
    call MPI_Bcast(t_lloyd%N1, 1, myMPItype1, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting logicals from t_tgmat'

    call MPI_Type_free(myMPItype1, ierr)

   end subroutine bcast_t_lly_1
#endif


   subroutine init_tlloyd(t_inc,t_lloyd,t_mpi_c_grid)
   
      use mod_mympi, only: myrank, master, nranks
   
      implicit none
   
      type(type_inc), intent(in) :: t_inc
      type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
      type(type_lloyd), intent(inout) :: t_lloyd
      
      integer :: ierr, nspin
      
      nspin = t_inc%NSPIN
      if(t_inc%NEWSOSOL) nspin = 1 ! t_inc%NSPIN !1
    
    
      if (.not. allocated(t_lloyd%dtmat)) then
         if (.not. t_lloyd%dtmat_to_file) then
            if(nranks.eq.1) then
               !allocate dtmat(lmmax,lmmax,irec_max) for irec_max=ielast*nspin*natyp
               allocate(t_lloyd%dtmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%IELAST*nspin*t_inc%NATYP), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%tmat'
            else
               !allocate dtmat(lmmax,lmmax,irec_max) for irec_max=iemax_local*nspin*natyp
               allocate(t_lloyd%dtmat(t_inc%LMMAXD,t_inc%LMMAXD,t_mpi_c_grid%ntot2*nspin*t_inc%NATYP), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%tmat for mpi'
            end if
         else
            allocate(t_lloyd%dtmat(1,1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_tgmat%tmat'
         end if
         t_lloyd%dtmat(:,:,:) = (0.d0, 0.d0)
      endif
      
      if (.not. allocated(t_lloyd%tralpha)) then
         if (.not. t_lloyd%tralpha_to_file) then
            if(nranks.eq.1) then
               allocate(t_lloyd%tralpha(t_inc%IELAST*nspin*t_inc%NATYP), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_lloyd%tralpha'
            else
               allocate(t_lloyd%tralpha(t_mpi_c_grid%ntot2*nspin*t_inc%NATYP), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_lloyd%tralpha for mpi'
            end if
         else
            allocate(t_lloyd%tralpha(1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_lloyd%tralpha'
         end if
         t_lloyd%tralpha(:) = (0.d0, 0.d0)
      endif
      
      if (.not. allocated(t_lloyd%cdos)) then
         if (.not. t_lloyd%cdos_diff_lly_to_file) then
            if(nranks.eq.1) then
               allocate(t_lloyd%cdos(t_inc%IELAST,t_inc%NSPIN), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_lloyd%cdos'
            else
               allocate(t_lloyd%cdos(t_inc%IELAST,t_inc%NSPIN), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_lloyd%cdos for mpi'
            end if
         else
            allocate(t_lloyd%cdos(1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_lloyd%cdos'
         end if
         t_lloyd%cdos(:,:) = (0.d0, 0.d0)
      endif

      if (.not. allocated(t_lloyd%dgref)) then
         if (.not. t_lloyd%dgref_to_file) then
            if(nranks.eq.1) then
               allocate(t_lloyd%dgref(t_inc%NACLSD*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_inc%IELAST), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_lloyd%dgref'
            else
               allocate(t_lloyd%dgref(t_inc%NACLSD*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_mpi_c_grid%ntot2), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_lloyd%dgref for mpi'
            end if
         else
            allocate(t_lloyd%dgref(1,1,1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_lloyd%dgref'
         end if
         t_lloyd%dgref(:,:,:,:) = (0.d0, 0.d0)
      endif
      
      if (.not. allocated(t_lloyd%g0tr)) then
         if (.not. t_lloyd%g0tr_to_file) then
            if(nranks.eq.1) then
               allocate(t_lloyd%g0tr(t_inc%IELAST), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_lloyd%g0tr'
            else
               allocate(t_lloyd%g0tr(t_mpi_c_grid%ntot2), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_lloyd%g0tr for mpi'
            end if
         else
            allocate(t_lloyd%g0tr(1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_lloyd%g0tr'
         end if
         t_lloyd%g0tr(:) = (0.d0, 0.d0)
      endif
      

   end subroutine init_tlloyd




#ifdef CPP_MPI
   subroutine save_t_mpi_c_grid(t_mpi_c_grid,subarr_dim, myMPI_comm_ie, myMPI_comm_at, myrank_ie, myrank_at, myrank_atcomm, nranks_ie, nranks_at, nranks_atcomm)
   
      use mpi
      implicit none
      type(type_mpi_cartesian_grid_info), intent(inout) :: t_mpi_c_grid
      integer, intent(in) :: subarr_dim(2)
      integer, intent(in) :: myMPI_comm_ie, myMPI_comm_at, myrank_ie, myrank_at, nranks_ie, nranks_at, nranks_atcomm, myrank_atcomm

      t_mpi_c_grid%dims = subarr_dim
      t_mpi_c_grid%myMPI_comm_ie  = myMPI_comm_ie
      t_mpi_c_grid%myMPI_comm_at  = myMPI_comm_at
      t_mpi_c_grid%myrank_ie  = myrank_ie
      t_mpi_c_grid%myrank_at  = myrank_at
      t_mpi_c_grid%myrank_atcomm  = myrank_atcomm
      t_mpi_c_grid%nranks_ie  = nranks_ie
      t_mpi_c_grid%nranks_at  = nranks_at
      t_mpi_c_grid%nranks_atcomm  = nranks_atcomm
   
   end subroutine save_t_mpi_c_grid
#endif


#ifdef CPP_MPI
   subroutine get_ntot_pT_ioff_pT_2D(t_mpi_c_grid,ntot_all,ioff_all)
   
   use mpi
!    use mod_mympi, only: myrank, nranks, master
   implicit none
   type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
   integer, intent(out) :: ntot_all(t_mpi_c_grid%nranks_ie*t_mpi_c_grid%nranks_at), ioff_all(t_mpi_c_grid%nranks_ie*t_mpi_c_grid%nranks_at)
   
   integer :: ntot_pT1(t_mpi_c_grid%nranks_ie), ioff_pT1(t_mpi_c_grid%nranks_ie)
   integer :: ntot_pT2(t_mpi_c_grid%nranks_at), ioff_pT2(t_mpi_c_grid%nranks_at)
   integer :: N1, N2, i1, i2, i3
   
   ntot_pT1 = t_mpi_c_grid%ntot_pT1
   ioff_pT1 = t_mpi_c_grid%ioff_pT1
   ntot_pT2 = t_mpi_c_grid%ntot_pT2
   ioff_pT2 = t_mpi_c_grid%ioff_pT2
   N1 = t_mpi_c_grid%nranks_ie
   N2 = t_mpi_c_grid%nranks_at
   
   
   do i1=1,N1
      do i2=1,N2
         i3 = i2+N2*(i1-1)
         ntot_all(i3) = ntot_pT1(i1)*ntot_pT2(i2)
         if (i3==1) then
            ioff_all(i3) = 0
         else
            ioff_all(i3) = ioff_all(i3-1)+ntot_all(i3-1)
         end if
      end do
   end do 
      
   end subroutine get_ntot_pT_ioff_pT_2D
#endif


#ifdef CPP_MPI
   subroutine gather_tmat(t_inc, t_tgmat, t_mpi_c_grid, ntot_pT, ioff_pT, mytot, mympi_comm, nranks)

    use mpi
    use mod_mympi,   only: myrank, master
    implicit none

    type(type_inc), intent(in) :: t_inc
    type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
    type(type_tgmatices), intent(inout) :: t_tgmat
    integer, intent(in) :: nranks
    integer, intent(in) :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1), mytot, mympi_comm

    integer :: ihelp
    integer :: nspin, recvcounts(nranks), displs(nranks)
    double complex, allocatable :: work(:,:,:)
    integer :: ierr,idim
    
    
    !Gather tmat so that all processors the full matrix for their part of the energy contour
    if(t_mpi_c_grid%nranks_ie>1) then
       nspin = t_inc%NSPIN
       if(t_inc%NEWSOSOL) nspin = 1
      
       ihelp      = t_inc%LMMAXD**2*t_mpi_c_grid%ntot2*nspin!*t_inc%NATYP/mytot
       recvcounts = ntot_pT*ihelp
       displs     = ioff_pT*ihelp
       
       allocate(work(t_inc%LMMAXD,t_inc%LMMAXD,t_mpi_c_grid%ntot2*nspin*t_inc%NATYP))
       call MPI_Allgatherv( t_tgmat%tmat, mytot*ihelp, MPI_DOUBLE_COMPLEX,         &
                          & work, recvcounts, displs, MPI_DOUBLE_COMPLEX, &
                          & mympi_comm, ierr )
       idim = t_inc%LMMAXD**2*t_mpi_c_grid%ntot2*nspin*t_inc%NATYP
       call zcopy(idim,work,1,t_tgmat%tmat,1)
       deallocate(work)
    end if
                          
   end subroutine gather_tmat
#endif


#ifdef CPP_MPI
   subroutine gather_gref(t_inc, t_tgmat, t_mpi_c_grid, ntot_pT, ioff_pT, mytot, mympi_comm, nranks)

    use mpi
    use mod_mympi,   only: myrank, master
    implicit none

    type(type_inc), intent(in) :: t_inc
    type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
    type(type_tgmatices), intent(inout) :: t_tgmat
    integer, intent(in) :: nranks
    integer, intent(in) :: ntot_pT(0:t_mpi_c_grid%nranks_ie-1), ioff_pT(0:t_mpi_c_grid%nranks_ie-1), mytot, mympi_comm

    integer :: ihelp
    integer :: ierr
        
    !Gather gref so that all processors have the full matrix for their part of the energy countour
    if(t_mpi_c_grid%dims(1)>1) then
       ihelp      = t_inc%NACLSD*t_inc%LMGF0D*t_inc%LMGF0D*t_inc%NCLSD
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_Bcast( t_tgmat%gref, mytot*ihelp, MPI_DOUBLE_COMPLEX, 0 , &
                          & mympi_comm, ierr )
    end if
    
   end subroutine gather_gref
#endif


#ifdef CPP_MPI
   subroutine gather_gmat(t_inc,t_tgmat,ntot_pT,ioff_pT,mytot)

    use mpi
    use mod_mympi,   only: myrank, nranks, master
    implicit none

    type(type_inc), intent(in) :: t_inc
    type(type_tgmatices), intent(inout) :: t_tgmat
    integer, intent(in) :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1), mytot

    integer :: ihelp
    integer :: recvcounts(0:nranks-1), displs(0:nranks-1)
    integer :: ierr

    !Gather Pkk' so that all processors have the full matrix
    ihelp      = t_inc%LMMAXD*t_inc%LMMAXD*t_inc%NQDOS!*t_inc%IELAST*t_inc%NSPIN*t_inc%NATYP
    if(t_mpi_c_grid%dims(1)>1) then
       recvcounts = ntot_pT*ihelp
       displs     = ioff_pT*ihelp
       call MPI_Allgatherv( t_tgmat%gmat, mytot*ihelp, MPI_DOUBLE_COMPLEX,         &
                          & t_tgmat%gmat, recvcounts, displs, MPI_DOUBLE_COMPLEX, &
                          & MPI_COMM_WORLD, ierr )
    end if

   end subroutine gather_gmat
#endif

end module
