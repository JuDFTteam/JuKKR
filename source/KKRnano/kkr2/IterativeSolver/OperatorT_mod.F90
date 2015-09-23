!> Simple operator class.
!> An operator can be applied (method 'apply') to a matrix mat_X and yields
!> another matrix mat_AX
!>
!> @author Elias Rabel

module OperatorT_mod
  implicit none
  private
  public :: OperatorT

  type, abstract :: OperatorT
    contains
    procedure(apply), deferred :: apply
  endtype

  interface
    !----------------------------------------------------------------------------
    !> Applies Operator on mat_X and returns result in mat_AX
    subroutine apply(self, mat_X, mat_AX)
      import OperatorT
      class(OperatorT) :: self
      double complex, intent(in)  :: mat_X(:,:)
      double complex, intent(out) :: mat_AX(:,:)
    endsubroutine
  endinterface

endmodule ! OperatorT_mod
