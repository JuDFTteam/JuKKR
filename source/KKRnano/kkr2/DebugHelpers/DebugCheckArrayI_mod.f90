
module DebugCheckArrayI_mod
  implicit none

  type DebugCheckArrayI
    private
    integer, dimension(:), pointer :: array_data
    integer :: num_elements
    character(len=32) :: array_name 
  end type

  contains

  subroutine createDebugCheckArrayI(debug_array, array_to_check, num_elements, array_name)
    implicit none

    type (DebugCheckArrayI), intent(inout) :: debug_array
    integer, dimension(num_elements), intent(in) :: array_to_check
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

  logical function testDebugCheckArrayI(debug_array, array_to_check, num_elements, fail_message)
    implicit none

    type (DebugCheckArrayI), intent(in) :: debug_array
    integer, dimension(num_elements), intent(in) :: array_to_check
    integer, intent(in), optional :: num_elements
    character(len=*), intent(in), optional :: fail_message

    integer :: ii

    testDebugCheckArrayI = .false.

    if (present(num_elements)) then
      if (num_elements /= debug_array%num_elements) then
        write(*,*) "testDebugCheckArrayI: Wrong number of elements specified. ", debug_array%array_name
        if (present(fail_message)) then
          write(*,*) fail_message
        end if
        return
      end if
    end if

    do ii = 1, debug_array%num_elements
      if(debug_array%array_data(ii) /= array_to_check(ii)) then
      write(*,*) "testDebugCheckArrayI: Arrays do not match. Element ", ii
        if (present(fail_message)) then
          write(*,*) debug_array%array_name, " ", fail_message
        else
          write(*,*) debug_array%array_name
        end if
      return
      end if  
    end do 
     
    testDebugCheckArrayI = .true.

  end function testDebugCheckArrayI

  subroutine destroyDebugCheckArrayI(debug_array)
    implicit none

    type (DebugCheckArrayI), intent(inout) :: debug_array
    
    deallocate(debug_array%array_data)
  end subroutine
end module


!program TryDebugCheckArrayI
!  use DebugCheckArrayI_mod
!  implicit none
!
!  integer, parameter :: dimx = 10
!  integer, parameter :: dimy = 10
!
!  integer, dimension(dimx, dimy) :: my_array
!
!  integer :: x, y
!  logical :: flag
!
!  type (DebugCheckArrayI) :: db
!
!  do y = 1, dimy
!    do x = 1, dimx
!      my_array(x,y) = x * y
!    end do
!  end do
!
!  call createDebugCheckArrayI(db, my_array, dimx*dimy, "my_array")
!
!  ! .. do something
!
!  write(*,*) testDebugCheckArrayI(db, my_array, dimx*dimy)
!
!  ! .. do something bad
!
!  my_array(3,5) = -3
!
!  write(*,*) testDebugCheckArrayI(db, my_array, dimx*dimy)
!
!  ! size argument is optional for test
!  write(*,*) testDebugCheckArrayI(db, my_array)
!
!  ! use optional fail_message
!  write(*,*) testDebugCheckArrayI(db, my_array, fail_message="location: main")
!
!
!  call destroyDebugCheckArrayI(db)
!end program
