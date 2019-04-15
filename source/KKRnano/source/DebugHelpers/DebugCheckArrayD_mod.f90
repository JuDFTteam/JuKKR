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

  subroutine createDebugCheckArrayD(self, array_to_check, num_elements, array_name)
    type(DebugCheckArrayD), intent(inout) :: self
    double precision, intent(in) :: array_to_check(num_elements)
    integer, intent(in) :: num_elements
    character(len=*), intent(in) :: array_name

    integer :: ii
    
    allocate(self%array_data(num_elements))

    self%num_elements = num_elements
    self%array_name = array_name

    do ii = 1, num_elements
      self%array_data(ii) = array_to_check(ii)  
    enddo ! ii

  endsubroutine

  logical function testDebugCheckArrayD(self, array_to_check, fail_message)
    type(DebugCheckArrayD), intent(in) :: self
    double precision, intent(in) :: array_to_check(*) ! accept any array
    character(len=*), intent(in), optional :: fail_message

    integer :: ii

    testDebugCheckArrayD = .false.

    do ii = 1, self%num_elements
      if (self%array_data(ii) /= array_to_check(ii)) then
        write(*,*) "testDebugCheckArrayD: Arrays do not match. Element ", ii
        if (present(fail_message)) then
          write(*,*) self%array_name, fail_message
        else
          write(*,*) self%array_name 
        endif
        return
      endif  
    enddo ! ii
     
    testDebugCheckArrayD = .true.

  endfunction

  elemental subroutine destroyDebugCheckArrayD(self)
    type(DebugCheckArrayD), intent(inout) :: self
    integer :: ist
    deallocate(self%array_data, stat=ist) ! ignore status
  endsubroutine ! destroy
  
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
