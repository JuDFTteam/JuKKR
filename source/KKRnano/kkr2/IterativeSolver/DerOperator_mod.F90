module DerOperator_mod
  use Operator_mod
  implicit none

  type, extends(Operator) :: DerOperator
    contains
    procedure :: apply => apply_DerOperator
  end type

  contains

  subroutine apply_DerOperator(self)
    class(DerOperator) :: self
    write(*,*) "DerOperator applied."
  end subroutine
 
end module
