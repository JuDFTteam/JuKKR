! Author: Elias Rabel
module DebugCheckArrayZ_mod
  implicit none
  public

  type DebugCheckArrayZ
    private
    double complex, allocatable :: array_data(:)
    integer :: num_elements
    character(len=32) :: array_name 
  endtype

  contains

  subroutine createDebugCheckArrayZ(self, array_to_check, num_elements, array_name)
    type(DebugCheckArrayZ), intent(inout) :: self
    double complex, intent(in) :: array_to_check(num_elements)
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

  logical function testDebugCheckArrayZ(self, array_to_check, fail_message)
    type(DebugCheckArrayZ), intent(in) :: self
    double complex, intent(in) :: array_to_check(*) ! accept any array
    character(len=*), intent(in), optional :: fail_message

    integer :: ii

    testDebugCheckArrayZ = .false.

    do ii = 1, self%num_elements
      if (self%array_data(ii) /= array_to_check(ii)) then
        write(*,*) "testDebugCheckArrayZ: Arrays do not match. Element ", ii
        if (present(fail_message)) then
          write(*,*) self%array_name, fail_message
        else
          write(*,*) self%array_name 
        endif
        return
      endif  
    enddo ! ii
     
    testDebugCheckArrayZ = .true.

  endfunction

  elemental subroutine destroyDebugCheckArrayZ(self)
    type(DebugCheckArrayZ), intent(inout) :: self
    integer :: ist
    deallocate(self%array_data, stat=ist) ! ignore status
  endsubroutine ! destroy
  
endmodule


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
!  type(DebugCheckArrayZ) :: db
!
!  do y = 1, dimy
!    do x = 1, dimx
!      my_array(x,y) = x * y
!    enddo
!  enddo
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
!endprogram
