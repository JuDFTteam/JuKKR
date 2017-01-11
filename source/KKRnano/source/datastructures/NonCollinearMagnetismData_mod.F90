!> Datastructure that contains data for non-collinear magnetism

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)

module NonCollinearMagnetismData_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  
  implicit none
  private
  public :: NOCOData, create, destroy, load, store, loadascii

  type NOCOData
    double precision, allocatable :: theta_noco(:)   !< non-collinear magnetism angle, WARNING: not synchronized between MPI threads
    double precision, allocatable :: phi_noco(:)     !< non-collinear magnetism angle, WARNING: not synchronized between MPI threads
    integer (kind=1), allocatable :: angle_fixed(:) !< keep angles fixed (1) or not (0), WARNING: not synchronized between MPI threads
  
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
  
  interface loadascii
    module procedure readNOCODataascii
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

  !-----------------------------------------------------------------------------
  !>    Reads nonco_angle.dat file.
  subroutine readNOCODataascii(theta, phi, angle_fixed, naez)
   
    double precision,  intent(out) :: theta (naez)
    double precision,  intent(out) :: phi (naez)
    integer (kind=1),  intent(out) :: angle_fixed (naez) ! (1): keep angle fixed, (0): relax angle
    integer, intent(in)            :: naez
  
    logical :: lread, lcheckangles
    integer :: I1 
    double precision, parameter :: PI=4.d0*datan(1.d0), eps=1d-5
  
    theta(:) = 0.D0
    phi(:) = 0.D0
    lread = .FALSE.
    LCHECKANGLES = .false.
  
    inquire(file='nonco_angle.dat',EXIST=lread)
    if (lread) then
            open(UNIT=10,FILE='nonco_angle.dat',FORM='FORMATTED')
            do I1 = 1,naez
               read(10,*) theta(I1),phi(I1),angle_fixed(I1)
              ! if((abs(theta(I1)).lt.(pi+eps)   .and. abs(theta(I1)).gt.eps) .or. &
              ! (abs(phi(I1)).lt.(2*pi+eps) .and. abs(phi(I1)).gt.eps)) then
              !   LCHECKANGLES = .true.
              ! endif
               if(LCHECKANGLES .eqv. .true.) then
                 die_here('angles in nonco_angle.dat are not given in a correct format')       
               endif
               theta(i1)=theta(I1)*(pi/180.0d0)
               phi(i1)  =phi(i1)*(pi/180.0d0)
            enddo
    else
         die_here('failed to read "nonco_angle.dat"!')
    endif
  
  endsubroutine
endmodule ! NonCollinearMagnetismData_mod
