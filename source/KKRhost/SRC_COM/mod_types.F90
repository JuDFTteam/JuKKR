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
   



   type :: type_inc
   
      integer :: Nparams = 13   ! number of parameters in type_inc
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
      logical :: NEWSOSOL = .false.
         
   end type type_inc
   
   
   type :: type_mpi_cartesian_grid_info
   
      integer :: Nparams = 10
      integer :: dims(2) = (/ -1, -1 /)
      integer :: myMPI_comm_grid = -1
      integer :: myMPI_comm_row  = -1
      integer :: myMPI_comm_col  = -1
      integer :: myrank_grid = -1
      integer :: myrank_row  = -1
      integer :: myrank_col  = -1
      integer :: nranks_row  = -1
      integer :: nranks_col  = -1
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


contains

   subroutine init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
   
      use mod_mympi, only: myrank, master, nranks
   
      implicit none
   
      type(type_inc), intent(in) :: t_inc
      type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
      type(type_tgmatices), intent(inout) :: t_tgmat
      
      integer :: ierr, nspin
      
      logical :: test= .false.
      
      nspin = t_inc%NSPIN
      if(t_inc%NEWSOSOL) nspin = 1
    
    if(test .and. myrank==master) then
   
      if (.not. allocated(t_tgmat%tmat)) then
         if (.not. t_tgmat%tmat_to_file) then
 
               !allocate tmat(lmmax,lmmax,irec_max) for irec_max=ielast*nspin*natyp
               allocate(t_tgmat%tmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%IELAST*nspin*t_inc%NATYP), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%tmat'
 
         else
            allocate(t_tgmat%tmat(1,1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_tgmat%tmat'
         end if
      
         t_tgmat%tmat(:,:,:) = (0.d0, 0.d0)
      endif
      
      if (.not. allocated(t_tgmat%gmat)) then
         if (.not. t_tgmat%gmat_to_file) then
 
               !allocate gmat(lmmax,lmmax,irec_max) for irec_max=NQDOS*IELAST*NSPIN*NATYP
               allocate(t_tgmat%gmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NQDOS*t_inc%IELAST*nspin*t_inc%NATYP), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%gmat'
            
         else
            allocate(t_tgmat%gmat(1,1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_tgmat%gmat'
         end if
      
         t_tgmat%gmat(:,:,:) = (0.d0, 0.d0)
      endif
      
      if (.not. allocated(t_tgmat%gref)) then
         if (.not. t_tgmat%gref_to_file) then
               !allocate gref(NACLSD*LMGF0D,LMGF0D,NCLSD,irec_max) for irec_max=IEMAX_local (=ntot2)
               allocate(t_tgmat%gref(t_inc%NACLSD*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_inc%IELAST), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%gref for mpi'
         else
            allocate(t_tgmat%gref(1,1,1,1), STAT=ierr)
            if(ierr/=0) stop 'Problem allocating dummy t_tgmat%gref'
         end if
   
      t_tgmat%gref(:,:,:,:) = (0.d0, 0.d0)
      endif
    
    else
      
      if (.not. allocated(t_tgmat%tmat)) then
         if (.not. t_tgmat%tmat_to_file) then
            if(nranks.eq.1) then
               !allocate tmat(lmmax,lmmax,irec_max) for irec_max=ielast*nspin*natyp
               allocate(t_tgmat%tmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%IELAST*t_inc%NSPIN*t_inc%NATYP), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%tmat'
            else
               !allocate tmat(lmmax,lmmax,irec_max) for irec_max=iemax_local*nspin*natyp
               allocate(t_tgmat%tmat(t_inc%LMMAXD,t_inc%LMMAXD,t_mpi_c_grid%ntot2*t_inc%NSPIN*t_inc%NATYP), STAT=ierr)
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
               allocate(t_tgmat%gmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NQDOS*t_inc%IELAST*t_inc%NSPIN*t_inc%NATYP), STAT=ierr)
               if(ierr/=0) stop 'Problem allocating t_tgmat%gmat'
            else
               !allocate gmat(lmmax,lmmax,irec_max) for irec_max=NQDOS*IEMAX_local*NSPIN*NATYP
               allocate(t_tgmat%gmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NQDOS*t_mpi_c_grid%ntot2*t_inc%NSPIN*t_inc%NATYP), STAT=ierr)
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
      
      end if

   end subroutine init_tgmat
   
   
   
#ifdef CPP_MPI
   subroutine bcast_t_inc_tgmat(t_inc,t_tgmat)
    !ruess: after myBcast_impcls from Pkkr_sidebranch2D_2014_12_16 by Bernd Zimmermann

    use mpi
    use mod_mympi,   only: myrank, nranks, master
    implicit none

    type(type_inc), intent(inout) :: t_inc
    type(type_tgmatices), intent(inout) :: t_tgmat

    integer :: blocklen1(t_inc%Nparams), etype1(t_inc%Nparams), myMPItype1 ! for parameter from t_inc
    integer :: blocklen2(t_tgmat%Nelements), etype2(t_tgmat%Nelements), myMPItype2 ! for logicals in t_tgmat
    integer :: blocklen3(t_tgmat%Nelements), etype3(t_tgmat%Nelements), myMPItype3 ! for allocatable arrays in t_tgmat
    integer :: ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp1(t_inc%Nparams), disp2(t_tgmat%Nelements), disp3(t_tgmat%Nelements), base

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
    call MPI_Get_address(t_inc%NEWSOSOL,     disp1(13), ierr)
    base  = disp1(1)
    disp1 = disp1 - base

    blocklen1(1:13)=1

    etype1(1:12) = MPI_INTEGER
    etype1(13) = MPI_LOGICAL

    call MPI_Type_create_struct(t_inc%Nparams, blocklen1, disp1, etype1, myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_impcls_1'

    call MPI_Type_commit(myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_impcls_1'

    call MPI_Bcast(t_inc%Nparams, 1, myMPItype1, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting impcls_1'

    call MPI_Type_free(myMPItype1, ierr)
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
               allocate(t_lloyd%cdos(t_mpi_c_grid%ntot2,t_inc%NSPIN), STAT=ierr)
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
   subroutine save_t_mpi_c_grid(t_mpi_c_grid,subarr_dim, myMPI_comm_grid, myMPI_comm_row,myMPI_comm_col, myrank_grid, myrank_row, myrank_col, nranks_row, nranks_col)
      use mpi
      implicit none
      type(type_mpi_cartesian_grid_info), intent(inout) :: t_mpi_c_grid
      integer, intent(in) :: subarr_dim(2)
      integer, intent(in) :: myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col,myrank_grid, myrank_row, myrank_col, nranks_row, nranks_col

      t_mpi_c_grid%dims = subarr_dim
      t_mpi_c_grid%myMPI_comm_grid = myMPI_comm_grid
      t_mpi_c_grid%myMPI_comm_row  = myMPI_comm_row
      t_mpi_c_grid%myMPI_comm_col  = myMPI_comm_col
      t_mpi_c_grid%myrank_grid = myrank_grid
      t_mpi_c_grid%myrank_row  = myrank_row
      t_mpi_c_grid%myrank_col  = myrank_col
      t_mpi_c_grid%nranks_row  = nranks_row
      t_mpi_c_grid%nranks_col  = nranks_col   
   
   end subroutine save_t_mpi_c_grid
#endif


#ifdef CPP_MPI
   subroutine get_ntot_pT_ioff_pT_2D(t_mpi_c_grid,ntot_all,ioff_all)
   
   use mpi
!    use mod_mympi, only: myrank, nranks, master
   implicit none
   type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
   integer, intent(out) :: ntot_all(t_mpi_c_grid%nranks_row*t_mpi_c_grid%nranks_col), ioff_all(t_mpi_c_grid%nranks_row*t_mpi_c_grid%nranks_col)
   
   integer :: ntot_pT1(t_mpi_c_grid%nranks_row), ioff_pT1(t_mpi_c_grid%nranks_row)
   integer :: ntot_pT2(t_mpi_c_grid%nranks_col), ioff_pT2(t_mpi_c_grid%nranks_col)
   integer :: N1, N2, i1, i2, i3
   
   ntot_pT1 = t_mpi_c_grid%ntot_pT1
   ioff_pT1 = t_mpi_c_grid%ioff_pT1
   ntot_pT2 = t_mpi_c_grid%ntot_pT2
   ioff_pT2 = t_mpi_c_grid%ioff_pT2
   N1 = t_mpi_c_grid%nranks_row
   N2 = t_mpi_c_grid%nranks_col
   
   
   do i1=1,N1
      do i2=1,N2
         i3 = i2+N2*(i1-1)
         ntot_all(i3) = ntot_pT1(i1)*ntot_pT2(i2)
!          if(ioff_pT1(i1).eq.0) ioff_pT1(i1) = 1
!          ioff_all(i3) = ioff_pT1(i1)*ntot_pT2(i2)+ioff_pT2(i2)
         if (i3==1) then
            ioff_all(i3) = 0
         else
            ioff_all(i3) = ioff_all(i3-1)+ntot_all(i3-1)
         end if
      end do
   end do 
   
!    write(*,*) 'in get_ntot_pT_ioff_pT_2D',   ntot_pT1,'ioff_pT1',ioff_pT1,'ntot_pT2',ntot_pT2,'ioff_pT2',ioff_pT2,'N1/2',N1,N2,'ntot_all',ntot_all,'ioff_all',ioff_all
   
   end subroutine get_ntot_pT_ioff_pT_2D
#endif


#ifdef CPP_MPI
   subroutine gather_tmat(t_inc, t_tgmat, t_mpi_c_grid, ntot_pT, ioff_pT, mytot, mympi_comm)

    use mpi
    use mod_mympi,   only: myrank, nranks, master
    implicit none

    type(type_inc), intent(in) :: t_inc
    type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
    type(type_tgmatices), intent(inout) :: t_tgmat
    integer, intent(in) :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1), mytot, mympi_comm

    integer :: ihelp
    integer :: recvcounts(0:nranks-1), displs(0:nranks-1)
    integer :: ierr
    
    logical :: test= .false.
    
    if(test) then
       write(*,*) myrank,'sends to master'
!        ihelp      = t_inc%LMMAXD*t_inc%LMMAXD*t_inc%NSPIN*t_inc%NATYP!*t_mpi_c_grid%ntot2 !*t_inc%IELAST*t_inc%NSPIN*t_inc%NATYP
!        if(t_inc%NEWSOSOL) ihelp = ihelp/t_inc%NSPIN
       ihelp      = t_inc%LMMAXD*t_inc%LMMAXD!*t_inc%NSPIN!*t_inc%NATYP!*t_mpi_c_grid%ntot2 !*t_inc%IELAST*t_inc%NSPIN*t_inc%NATYP
!        if(t_inc%NEWSOSOL) ihelp = ihelp/t_inc%NSPIN
       recvcounts = ntot_pT*ihelp
       displs     = ioff_pT*ihelp
       call MPI_gatherv( t_tgmat%tmat, mytot*ihelp, MPI_DOUBLE_COMPLEX,         &
                          & t_tgmat%tmat, recvcounts, displs, MPI_DOUBLE_COMPLEX, &
                          & master,mympi_comm, ierr )
                          
!                           int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
!                 void *recvbuf, const int *recvcounts, const int *displs,
!                 MPI_Datatype recvtype, int root, MPI_Comm comm)
    else
    !Gather tmat so that all processors the full matrix for their part of the energy contour
    if(t_mpi_c_grid%ntot1>1) then
       ihelp      = t_inc%LMMAXD*t_inc%LMMAXD!*t_inc%NSPIN!*t_inc%NATYP!*t_mpi_c_grid%ntot2 !*t_inc%IELAST*t_inc%NSPIN*t_inc%NATYP
!        if(t_inc%NEWSOSOL) ihelp = ihelp/t_inc%NSPIN
       recvcounts = ntot_pT*ihelp
       displs     = ioff_pT*ihelp
       call MPI_Allgatherv( t_tgmat%tmat, mytot*ihelp, MPI_DOUBLE_COMPLEX,         &
                          & t_tgmat%tmat, recvcounts, displs, MPI_DOUBLE_COMPLEX, &
                          & mympi_comm, ierr )
    end if
    end if !test
                          
   end subroutine gather_tmat
#endif


#ifdef CPP_MPI
   subroutine gather_gref(t_inc, t_tgmat, t_mpi_c_grid, ntot_pT, ioff_pT, mytot, mympi_comm)

    use mpi
    use mod_mympi,   only: myrank, nranks, master
    implicit none

    type(type_inc), intent(in) :: t_inc
    type(type_mpi_cartesian_grid_info), intent(in) :: t_mpi_c_grid
    type(type_tgmatices), intent(inout) :: t_tgmat
    integer, intent(in) :: ntot_pT(0:t_mpi_c_grid%nranks_col-1), ioff_pT(0:t_mpi_c_grid%nranks_col-1), mytot, mympi_comm

    integer :: ihelp
    integer :: recvcounts(0:t_mpi_c_grid%nranks_col-1), displs(0:t_mpi_c_grid%nranks_col-1)
!     integer :: recvcounts(0:nranks-1), displs(0:nranks-1)
    integer :: ierr
        
    logical :: test= .false.
    
    if(test) then
       write(*,*) 'gather_gref:',myrank,t_inc%NACLSD,t_inc%LMGF0D,t_inc%NCLSD,mytot,ntot_pT,ioff_pT,shape(t_tgmat%gref)
       ihelp      = t_inc%NACLSD*t_inc%LMGF0D*t_inc%LMGF0D*t_inc%NCLSD!*t_mpi_c_grid%ntot2!*t_inc%IELAST
       recvcounts = ntot_pT*ihelp
       displs     = ioff_pT*ihelp
       call MPI_gatherv( t_tgmat%gref, mytot*ihelp, MPI_DOUBLE_COMPLEX,         &
                       & t_tgmat%gref, recvcounts, displs, MPI_DOUBLE_COMPLEX,  &
                       & master, mympi_comm, ierr )    
    else
    !Gather gref so that all processors have the full matrix for their part of the energy countour
    if(t_mpi_c_grid%ntot1>1) then
       ihelp      = t_inc%NACLSD*t_inc%LMGF0D*t_inc%LMGF0D*t_inc%NCLSD!*t_mpi_c_grid%ntot2!*t_inc%IELAST
       recvcounts = ntot_pT*ihelp
       displs     = ioff_pT*ihelp
       call MPI_Allgatherv( t_tgmat%gref, mytot*ihelp, MPI_DOUBLE_COMPLEX,         &
                          & t_tgmat%gref, recvcounts, displs, MPI_DOUBLE_COMPLEX,  &
                          & mympi_comm, ierr )
    end if
    end if !test
    
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
    recvcounts = ntot_pT*ihelp
    displs     = ioff_pT*ihelp
    call MPI_Allgatherv( t_tgmat%gmat, mytot*ihelp, MPI_DOUBLE_COMPLEX,         &
                       & t_tgmat%gmat, recvcounts, displs, MPI_DOUBLE_COMPLEX, &
                       & MPI_COMM_WORLD, ierr )
   end subroutine gather_gmat
#endif

end module
