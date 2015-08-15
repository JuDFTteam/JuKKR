!> Module to apply KKR coefficient matrix.
!>
!> The setup of the KKR coefficient matrix is rather complicated.
!> Therefore the needed data is stored in a 'MultScatData' struct.
!>
!> *) One has to get a reference (pointer) to the MultScatData workspace by
!>    using 'get_ms_workspace' and set up the workspace properly (routine kkrmat01)
!> *) Then one can apply the KKR coefficient matrix on any dense matrix using
!>    'apply'

module KKROperator_mod
  use MultScatData_mod, only: MultScatData
  use OperatorT_mod, only: OperatorT
  implicit none
  private
  public :: KKROperator, create, destroy
  public :: create_KKROperator, destroy_KKROperator ! deprecated
  public :: get_ms_workspace

  !> Represents the operator/matrix (1 - \Delta T G_ref).
  type, extends(OperatorT) :: KKROperator
    private
    type(MultScatData), pointer :: ms => null()
    contains
    procedure :: create => create_KKROperator
    procedure :: apply  => apply_KKROperator
    procedure :: get_ms_workspace
    procedure :: destroy => destroy_KKROperator
  end type

  interface create
    module procedure create_KKROperator
  endinterface
  
  interface destroy
    module procedure destroy_KKROperator
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  subroutine create_KKROperator(self)
    class(KKROperator) :: self

    allocate(self%ms)
  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroy_KKROperator(self)
    class(KKROperator) :: self

    deallocate(self%ms)
  end subroutine

  !----------------------------------------------------------------------------
  !> Applies Operator on mat_X and returns result in mat_AX.
  subroutine apply_KKROperator(self, mat_X, mat_AX)
    use vbrmv_mat_mod, only: multiply_vbr
    class(KKROperator) :: self
    double complex, intent(in)  :: mat_X(:,:)
    double complex, intent(out) :: mat_AX(:,:)

    double complex, parameter :: CZERO = (0.0D0,0.0D0)
    mat_AX = CZERO

    ! perform sparse VBR matrix * dense matrix
    call multiply_vbr(self%ms%GLLH, mat_X, mat_AX, self%ms%sparse)
  end subroutine

  !----------------------------------------------------------------------------
  !> Return a reference to the multiple scattering workspace.
  function get_ms_workspace(self) result(ms)
    class(KKROperator) :: self
    type(MultScatData), pointer :: ms

    ms => self%ms
  end function

end module
