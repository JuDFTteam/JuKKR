module DerOperator_mod
  use Operator_mod
  implicit none

  type, extends(Operator) :: DerOperator
    contains
    procedure :: apply => apply_DerOperator
  end type

  contains

  !----------------------------------------------------------------------------
  !> Applies Operator on mat_X and returns result in mat_AX
  subroutine apply_DerOperator(self, mat_X, mat_AX)
    class(DerOperator) :: self
    double complex, intent(in)  :: mat_X
    double complex, intent(out) :: mat_AX
    write(*,*) "DerOperator applied."
  end subroutine
 
end module
