! Author: Elias Rabel
module DebugCheckArrayZ_mod
  implicit none
  public

  type DebugCheckArrayZ
    private
    double complex, dimension(:), pointer :: array_data
    integer :: num_elements
    character(len=32) :: array_name 
  end type

  contains

  subroutine createDebugCheckArrayZ(debug_array, array_to_check, num_elements, array_name)
    implicit none

    type (DebugCheckArrayZ), intent(inout) :: debug_array
    double complex, dimension(num_elements), intent(in) :: array_to_check
    integer, intent(in) :: num_elements
    character(len=*), intent(in) :: array_name

    integer :: ii
    
    allocate(debug_array%array_data(num_elements))

    debug_array%num_elements = num_elements
    debug_array%array_name = array_name

    do ii = 1, num_elements
      debug_array%array_data(ii) = array_to_check(ii)  
    end do 

  end subroutine

  logical function testDebugCheckArrayZ(debug_array, array_to_check, fail_message)
    implicit none

    type (DebugCheckArrayZ), intent(in) :: debug_array
    double complex, dimension(*), intent(in) :: array_to_check ! accept any array
    character(len=*), intent(in), optional :: fail_message

    integer :: ii

    testDebugCheckArrayZ = .false.

    do ii = 1, debug_array%num_elements
      if(debug_array%array_data(ii) /= array_to_check(ii)) then
      write(*,*) "testDebugCheckArrayZ: Arrays do not match. Element ", ii
        if (present(fail_message)) then
          write(*,*) debug_array%array_name, fail_message
        else
          write(*,*) debug_array%array_name 
        end if
      return
      end if  
    end do 
     
    testDebugCheckArrayZ = .true.

  end function testDebugCheckArrayZ

  subroutine destroyDebugCheckArrayZ(debug_array)
    implicit none

    type (DebugCheckArrayZ), intent(inout) :: debug_array
    
    deallocate(debug_array%array_data)
  end subroutine
end module


!program TryDebugCheckArrayZ
!  use DebugCheckArrayZ_mod
!  implicit none
!
!  integer, parameter :: dimx = 10
!  integer, parameter :: dimy = 10
!
!  double complex, dimension(dimx, dimy) :: my_array
!
!  integer :: x, y
!  logical :: flag
!
!  type (DebugCheckArrayZ) :: db
!
!  do y = 1, dimy
!    do x = 1, dimx
!      my_array(x,y) = x * y
!    end do
!  end do
!
!  call createDebugCheckArrayZ(db, my_array, dimx*dimy, "my_array")
!
!  ! .. do something
!
!  write(*,*) testDebugCheckArrayZ(db, my_array)
!
!  ! .. do something bad
!
!  my_array(3,5) = -3
!
!  write(*,*) testDebugCheckArrayZ(db, my_array)
!
!  ! use optional fail_message
!  write(*,*) testDebugCheckArrayZ(db, my_array, fail_message="location: main")
!
!
!  call destroyDebugCheckArrayZ(db)
!end program
