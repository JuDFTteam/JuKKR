module MyPrecond_mod
  use Precond_mod
  implicit none

  type, extends(Precond) :: MyPrecond
    private
    integer :: member = 5
    contains
    procedure :: apply => apply_MyPrecond
  end type

  contains

  subroutine apply_MyPrecond(self)
    class(MyPrecond) :: self
    write(*,*) "MyPrecond applied."
  end subroutine
 
end module
