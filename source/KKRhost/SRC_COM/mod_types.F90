module mod_types

implicit none

      
   type :: type_main0
      
      integer :: i_iteration = -1 
      integer :: N_iteration = -1
      integer :: mit_bry = 1
      
   end type type_main0



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
   
      integer :: Nparams = 9   ! number of parameters in type_inc
      integer :: LMMAXD  = -1
      integer :: NSPIN   = -1
      integer :: IELAST  = -1
      integer :: NQDOS   = -1
      integer :: NATYP   = -1
      integer :: LMGF0D  = -1
      integer :: NCLSD   = -1
      integer :: NACLSD  = -1
         
   end type type_inc


!    type :: type_lloyd
! 
!       integer ::
!    
! 
!    end type type_lloyd


   type (type_main0), save :: type0
   type (type_inc), save :: t_inc
   type (type_tgmatices), save :: t_tgmat
!    type (type_lloyd), save :: t_lloyd


contains

   subroutine init_tgmat(t_inc,t_tgmat)
   
      implicit none
   
      type(type_inc), intent(in) :: t_inc
      type(type_tgmatices), intent(inout) :: t_tgmat
      
      integer :: ierr
   
      if (.not. t_tgmat%tmat_to_file) then
         !allocate tmat(lmmax,lmmax,irec_max) for irec_max=ielast*nspin*natyp
         allocate(t_tgmat%tmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%IELAST*t_inc%NSPIN*t_inc%NATYP), STAT=ierr)
         if(ierr/=0) stop 'Problem allocating t_tgmat%tmat'
      else
         allocate(t_tgmat%tmat(1,1,1), STAT=ierr)
         if(ierr/=0) stop 'Problem allocating dummy t_tgmat%tmat'
      end if
   
      t_tgmat%tmat(:,:,:) = (0.d0, 0.d0)
   
      if (.not. t_tgmat%gmat_to_file) then
         !allocate gmat(lmmax,lmmax,irec_max) for irec_max=NQDOS*IELAST*NSPIN*NATYP
         allocate(t_tgmat%gmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NQDOS*t_inc%IELAST*t_inc%NSPIN*t_inc%NATYP), STAT=ierr)
         if(ierr/=0) stop 'Problem allocating t_tgmat%gmat'
      else
         allocate(t_tgmat%gmat(1,1,1), STAT=ierr)
         if(ierr/=0) stop 'Problem allocating dummy t_tgmat%gmat'
      end if
   
      t_tgmat%gmat(:,:,:) = (0.d0, 0.d0)
   
      if (.not. t_tgmat%gref_to_file) then
         !allocate gref(NACLSD*LMGF0D,LMGF0D,NCLSD,irec_max) for irec_max=IELAST
         allocate(t_tgmat%gref(t_inc%NACLSD*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_inc%IELAST), STAT=ierr)
         if(ierr/=0) stop 'Problem allocating t_tgmat%gref'
      else
         allocate(t_tgmat%gref(1,1,1,1), STAT=ierr)
         if(ierr/=0) stop 'Problem allocating dummy t_tgmat%gref'
      end if
   
      t_tgmat%gref(:,:,:,:) = (0.d0, 0.d0)

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
    base  = disp1(1)
    disp1 = disp1 - base

    blocklen1(1:9)=1

    etype1(1:9) = MPI_INTEGER

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

    !then broadcast arrays
    if(myrank/=master) call init_tgmat(t_inc,t_tgmat)

    call MPI_Get_address(t_tgmat%Nelements,      disp3(1), ierr)
    call MPI_Get_address(t_tgmat%tmat,           disp3(2), ierr)
    call MPI_Get_address(t_tgmat%gmat,           disp3(3), ierr)
    call MPI_Get_address(t_tgmat%gref,           disp3(4), ierr)

    base  = disp3(1)
    disp3 = disp3 - base

    blocklen3(1)=1
    blocklen3(2)=size(t_tgmat%tmat)
    blocklen3(3)=size(t_tgmat%gmat)
    blocklen3(4)=size(t_tgmat%gref)

    etype3(1) = MPI_INTEGER
    etype3(2) = MPI_DOUBLE_COMPLEX
    etype3(3) = MPI_DOUBLE_COMPLEX
    etype3(4) = MPI_DOUBLE_COMPLEX

    call MPI_Type_create_struct(t_tgmat%Nelements, blocklen3, disp3, etype3, myMPItype3, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_tgmat_arrays'

    call MPI_Type_commit(myMPItype3, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_tgmat_arrays'

    call MPI_Bcast(t_tgmat%Nelements, 1, myMPItype3, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting tgmat_arrays'

    call MPI_Type_free(myMPItype3, ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end subroutine bcast_t_inc_tgmat
#endif

end module
