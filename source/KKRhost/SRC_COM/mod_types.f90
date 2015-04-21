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
      
      ! allocatable arrays for tmat, gmat and gref
      double complex, allocatable :: tmat(:,:,:)       ! dimensions=LMMAXD, LMMAXD, IREC; IREC= IE+IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1) ;IE=1,...,IELAST, ISPIN=1,...,NSPIN, I1=1,...,NATYP)
      double complex, allocatable :: gmat(:,:,:)       ! dimensions=LMMAXD, LMMAXD, IREC; IREC= IQDOS+NQDOS*(IE-1)+NQDOS*IELAST*(ISPIN-1)+IELAST*NSPIN*(I1-1) ;IE=1,...,IELAST, ISPIN=1,...,NSPIN, I1=1,...,NATYP)
      double complex, allocatable :: gref(:,:,:,:)       !GINP(NACLSD*LMGF0D,LMGF0D,NCLSD) IREC=IE=1,...,IELAST
   

   end type type_tgmatices
   



   type :: type_inc
   
      integer :: LMMAXD = -1
      integer :: NSPIN  = -1
      integer :: IELAST = -1
      integer :: NQDOS  = -1
      integer :: NATYP  = -1
      integer :: LMGF0D = -1
      integer :: NCLSD  = -1
      integer :: NACLSD = -1
         
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
   
         t_tgmat%tmat(:,:,:) = (0.d0, 0.d0)
      end if
   
      if (.not. t_tgmat%gmat_to_file) then
         !allocate gmat(lmmax,lmmax,irec_max) for irec_max=NQDOS*IELAST*NSPIN*NATYP
         allocate(t_tgmat%gmat(t_inc%LMMAXD,t_inc%LMMAXD,t_inc%NQDOS*t_inc%IELAST*t_inc%NSPIN*t_inc%NATYP), STAT=ierr)
         if(ierr/=0) stop 'Problem allocating t_tgmat%gmat'
   
         t_tgmat%gmat(:,:,:) = (0.d0, 0.d0)
      end if
   
      if (.not. t_tgmat%gref_to_file) then
         !allocate gref(NACLSD*LMGF0D,LMGF0D,NCLSD,irec_max) for irec_max=IELAST
         allocate(t_tgmat%gref(t_inc%NACLSD*t_inc%LMGF0D,t_inc%LMGF0D,t_inc%NCLSD,t_inc%IELAST), STAT=ierr)
         if(ierr/=0) stop 'Problem allocating t_tgmat%gref'
   
         t_tgmat%gref(:,:,:,:) = (0.d0, 0.d0)
      end if

   end subroutine


end module