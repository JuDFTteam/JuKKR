! Author: Elias Rabel
module DebugCheckArrayD_mod
  implicit none
  public

  type DebugCheckArrayD
    private
    double precision, allocatable :: array_data(:)
    integer :: num_elements
    character(len=32) :: array_name 
  endtype

  contains

  subroutine createDebugCheckArrayD(debug_array, array_to_check, num_elements, array_name)
    type(DebugCheckArrayD), intent(inout) :: debug_array
    double precision, intent(in) :: array_to_check(num_elements)
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

  logical function testDebugCheckArrayD(debug_array, array_to_check, fail_message)
    type(DebugCheckArrayD), intent(in) :: debug_array
    double precision, intent(in) :: array_to_check(*) ! accept any array
    character(len=*), intent(in), optional :: fail_message

    integer :: ii

    testDebugCheckArrayD = .false.

    do ii = 1, debug_array%num_elements
      if (debug_array%array_data(ii) /= array_to_check(ii)) then
        write(*,*) "testDebugCheckArrayD: Arrays do not match. Element ", ii
        if (present(fail_message)) then
          write(*,*) debug_array%array_name, fail_message
        else
          write(*,*) debug_array%array_name 
        endif
        return
      endif  
    enddo ! ii
     
    testDebugCheckArrayD = .true.

  endfunction

  subroutine destroyDebugCheckArrayD(debug_array)
    type(DebugCheckArrayD), intent(inout) :: debug_array
    
    deallocate(debug_array%array_data)
  endsubroutine
  
endmodule

!
!program TryDebugCheckArrayD
!  use DebugCheckArrayD_mod
!  implicit none
!
!  integer, parameter :: dimx = 10
!  integer, parameter :: dimy = 10
!
!  double precision, dimension(dimx, dimy) :: my_array
!
!  integer :: x, y
!  logical :: flag
!
!  type(DebugCheckArrayD) :: db
!
!  do y = 1, dimy
!    do x = 1, dimx
!      my_array(x,y) = x * y
!    enddo
!  enddo
!
!  call createDebugCheckArrayD(db, my_array, dimx*dimy, "my_array")
!
!  ! .. do something
!
!  write(*,*) testDebugCheckArrayD(db, my_array)
!
!  ! .. do something bad
!
!  my_array(3,5) = -3
!
!  write(*,*) testDebugCheckArrayD(db, my_array)
!
!  ! use optional fail_message
!  write(*,*) testDebugCheckArrayD(db, my_array, fail_message="location: main")
!
!
!  call destroyDebugCheckArrayD(db)
!endprogram
