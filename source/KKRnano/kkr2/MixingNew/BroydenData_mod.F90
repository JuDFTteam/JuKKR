! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module BroydenData_mod
  implicit none
  private
  public :: BroydenData, create, destroy
  public :: createBroydenData, destroyBroydenData ! deprecated
  

  type BroydenData
    double precision, allocatable :: sm1s(:)
    double precision, allocatable :: fm1s(:)
    double precision, allocatable :: ui2(:,:)
    double precision, allocatable :: vi2(:,:)
    double precision, allocatable :: wit(:)

    integer :: ntird
    integer :: itdbryd

    integer :: iteration_count
    integer :: imix
    double precision :: mixing
  endtype

  interface create
    module procedure createBroydenData
  endinterface
  
  interface destroy
    module procedure destroyBroydenData
  endinterface
  
  contains

  !-----------------------------------------------------------------------------
  !> Constructs a BroydenData object.
  !> @param[inout] self    The BroydenData object to construct.
  !> @param[in]    ntird
  !> @param[in]    itdbryd
  subroutine createBroydenData(self, ntird, itdbryd, imix, mixing)
    type(BroydenData), intent(inout) :: self
    integer, intent(in) :: ntird
    integer, intent(in) :: itdbryd
    integer, intent(in) :: imix
    double precision, intent(in) :: mixing

    integer :: memory_stat

    self%ntird = ntird
    self%itdbryd = itdbryd

    self%iteration_count = 1
    self%imix = imix
    self%mixing = mixing

    ALLOCATECHECK(self%sm1s(ntird))
    ALLOCATECHECK(self%fm1s(ntird))
    ALLOCATECHECK(self%ui2(ntird,2:itdbryd))
    ALLOCATECHECK(self%vi2(ntird,2:itdbryd))
    ALLOCATECHECK(self%wit(2:itdbryd))

    self%sm1s = 0.d0
    self%fm1s = 0.d0
    self%ui2  = 0.d0
    self%vi2  = 0.d0
    self%wit  = 0.d0
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a BroydenData object.
  !> @param[inout] self    The BroydenData object to destroy.
  subroutine destroyBroydenData(self)
    type(BroydenData), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%sm1s)
    DEALLOCATECHECK(self%fm1s)
    DEALLOCATECHECK(self%ui2)
    DEALLOCATECHECK(self%vi2)
    DEALLOCATECHECK(self%wit)
  endsubroutine ! destroy
  
endmodule ! BroydenData_mod
