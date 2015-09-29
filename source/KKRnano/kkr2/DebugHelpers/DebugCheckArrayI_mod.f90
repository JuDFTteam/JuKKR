! Author: Elias Rabel
module DebugCheckArrayI_mod
  implicit none
  public

  type DebugCheckArrayI
    private
    integer, allocatable :: array_data(:)
    integer :: num_elements
    character(len=32) :: array_name 
  endtype

  contains

  subroutine createDebugCheckArrayI(debug_array, array_to_check, num_elements, array_name)
    type(DebugCheckArrayI), intent(inout) :: debug_array
    integer, intent(in) :: array_to_check(num_elements)
    integer, intent(in) :: num_elements
    character(len=*), intent(in) :: array_name

    integer :: ii
    
    allocate(debug_array%array_data(num_elements))

    debug_array%num_elements = num_elements
    debug_array%array_name = array_name

    do ii = 1, num_elements
      debug_array%array_data(ii) = array_to_check(ii)  
    enddo ! ii

  endsubroutine

  logical function testDebugCheckArrayI(debug_array, array_to_check, fail_message)
    type(DebugCheckArrayI), intent(in) :: debug_array
    integer, intent(in) :: array_to_check(*) ! accept any array
    character(len=*), intent(in), optional :: fail_message

    integer :: ii

    testDebugCheckArrayI = .false.

    do ii = 1, debug_array%num_elements
      if (debug_array%array_data(ii) /= array_to_check(ii)) then
        write(*,*) "testDebugCheckArrayI: Arrays do not match. Element ", ii
        if (present(fail_message)) then
          write(*,*) debug_array%array_name, fail_message
        else
          write(*,*) debug_array%array_name 
        endif
        return
      endif  
    enddo ! ii
     
    testDebugCheckArrayI = .true.

  endfunction

  subroutine destroyDebugCheckArrayI(debug_array)
    type(DebugCheckArrayI), intent(inout) :: debug_array
    
    deallocate(debug_array%array_data)
  endsubroutine
  
endmodule

!
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
!  type(DebugCheckArrayI) :: db
!
!  do y = 1, dimy
!    do x = 1, dimx
!      my_array(x,y) = x * y
!    enddo
!  enddo
!
!  call createDebugCheckArrayI(db, my_array, dimx*dimy, "my_array")
!
!  ! .. do something
!
!  write(*,*) testDebugCheckArrayI(db, my_array)
!
!  ! .. do something bad
!
!  my_array(3,5) = -3
!
!  write(*,*) testDebugCheckArrayI(db, my_array)
!
!  ! use optional fail_message
!  write(*,*) testDebugCheckArrayI(db, my_array, fail_message="location: main")
!
!
!  call destroyDebugCheckArrayI(db)
!endprogram
