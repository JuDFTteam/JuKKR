!> Datastructure that contains data for non-collinear magnetism

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)

module NonCollinearMagnetismData_mod
  implicit none
  private
  public :: NOCOData, create, destroy, load, store

  type NOCOData
    double precision, allocatable :: theta_noco(:)   !< non-collinear magnetism angle
    double precision, allocatable :: phi_noco(:)     !< non-collinear magnetism angle
    integer (kind=1), allocatable :: angle_fixed(:) !< keep angles fixed (1) or not (0)
  
  endtype ! NOCOData

  interface create
    module procedure createNOCOData
  endinterface
  
  interface destroy
    module procedure destroyNOCOData
  endinterface

  interface load
    module procedure readNOCOData
  endinterface
  
  interface store
    module procedure writeNOCOData
  endinterface
  
  contains

  !-----------------------------------------------------------------------------
  !> Constructs a NOCOData object. (implementation, don't call directly!)
  !> @param[in,out] self    The NOCOData object to construct.
  !> @param[in]    naez
  subroutine createNOCOData(self, naez)
    type(NOCOData), intent(inout) :: self
    integer, intent(in) :: naez
    
    integer :: memory_stat

    ALLOCATECHECK(self%theta_noco(naez))
    ALLOCATECHECK(self%phi_noco(naez))
    ALLOCATECHECK(self%angle_fixed(naez))
    
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a NOCOData object.
  !> @param[in,out] self    The NOCOData object to destroy.
  elemental subroutine destroyNOCOData(self)
    type(NOCOData), intent(inout) :: self

    integer :: ist ! ignore status
#define DEALLOCATECHECK(X) deallocate(X, stat=ist)

    DEALLOCATECHECK(self%theta_noco)
    DEALLOCATECHECK(self%phi_noco)
    DEALLOCATECHECK(self%angle_fixed)
  endsubroutine ! destroy

  !-----------------------------------------------------------------------------
  !> Writes NOCOData data to file.
  !> @param[in] self    NOCOData object to write.
  subroutine writeNOCOData(self, filename)
    type(NOCOData), intent(in) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 67

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) self%theta_noco, &
              self%phi_noco, &
              self%angle_fixed
    close(fu)

  endsubroutine ! store

  !-----------------------------------------------------------------------------
  !> Reads NOCOData data from file.
  !> @param[in,out] self    The NOCOData object to read.
  subroutine readNOCOData(self, filename)
    type(NOCOData), intent(inout) :: self
    character(len=*), intent(in) :: filename

    integer, parameter :: fu = 67

    open(fu, file=filename, form='unformatted', action='read')

    read(fu) self%theta_noco, &
              self%phi_noco, &
              self%angle_fixed
    close(fu)

  endsubroutine ! read

endmodule ! NonCollinearMagnetismData_mod
