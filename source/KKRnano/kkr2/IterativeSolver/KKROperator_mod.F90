module KKROperator_mod
  use OperatorT_mod
  implicit none

  type, extends(OperatorT) :: KKROperator
    contains
    procedure :: apply => apply_KKROperator
  end type

  contains

  !----------------------------------------------------------------------------
  !> Applies Operator on mat_X and returns result in mat_AX
  subroutine apply_KKROperator(self, mat_X, mat_AX)
    class(KKROperator) :: self
    double complex, intent(in)  :: mat_X
    double complex, intent(out) :: mat_AX
    write(*,*) "DerOperator applied."
  end subroutine
 
end module
