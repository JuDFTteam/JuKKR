#include "../DebugHelpers/test_macros.h"

!> Define the block-circulant preconditioning matrix

! SEVERE restriction of preconditioner code:
! The number of non-zero blocks must be the same in each row
! this is not the case for all crystal structures (e.g. perovskite)
! or amorphous structures

module BCPOperator_mod
  use OperatorT_mod
  use SolverOptions_mod
  use ClusterInfo_mod

  implicit none

  !> Represents the Block-Circulant preconditioning matrix
  type, extends(OperatorT) :: BCPOperator
    ! private
    double complex, allocatable :: GLLHBLCK(:,:)
    type (SolverOptions) :: solver_opts
    type (ClusterInfo), pointer :: cluster_info
    integer :: lmmaxd
    contains
    procedure :: create => create_BCPOperator
    procedure :: calc   => calc_BCPOperator
    procedure :: apply  => apply_BCPOperator
    procedure :: destroy => destroy_BCPOperator
  end type

  contains

  !----------------------------------------------------------------------------
  !> Setup of BCP preconditioner
  subroutine create_BCPOperator(self, solver_opts, cluster_info, lmmaxd)
    class(BCPOperator) :: self
    type(SolverOptions), intent(in) :: solver_opts
    type (ClusterInfo), target  :: cluster_info
    integer, intent(in) :: lmmaxd

    integer naezd
    integer blocks_per_row

    naezd = cluster_info%naez_trc

    CHECKASSERT(naezd == solver_opts%NATBLD*solver_opts%XDIM * solver_opts%YDIM*solver_opts%ZDIM)

    allocate(self%GLLHBLCK(solver_opts%NATBLD*LMMAXD, naezd*LMMAXD))

    blocks_per_row = cluster_info%numn0_trc(1)

    ! SEVERE restriction of preconditioner code:
    ! The number of non-zero blocks must be the same in each row
    ! this is not the case for all crystal structures (e.g. perovskite)
    CHECKASSERT(all(cluster_info%numn0_trc == blocks_per_row ))

    self%solver_opts = solver_opts
    self%cluster_info => cluster_info
    self%lmmaxd = lmmaxd

  end subroutine

  !----------------------------------------------------------------------------
  !> Calculation of BCP preconditioner.
  !>
  !> uses BCPWUPPER
  subroutine calc_BCPOperator(self, GLLH)
    class(BCPOperator) :: self
    double complex, intent(in) :: GLLH(:)

    integer naezd
    integer blocks_per_row

    if (self%solver_opts%bcp /= 1) return ! no preconditioning selected

    naezd = self%cluster_info%naez_trc
    blocks_per_row = self%cluster_info%numn0_trc(1)

    call BCPWUPPER(GLLH, self%GLLHBLCK,NAEZD,self%cluster_info%NUMN0_trc,self%cluster_info%INDN0_trc, &
                   self%lmmaxd, self%solver_opts%natbld, self%solver_opts%xdim, &
                   self%solver_opts%ydim, self%solver_opts%zdim, &
                   blocks_per_row)

  end subroutine

  !----------------------------------------------------------------------------
  !> Applies Preconditioner/Operator on mat_X and returns result in mat_AX.
  subroutine apply_BCPOperator(self, mat_X, mat_AX)
    class(BCPOperator) :: self
    double complex, intent(in)  :: mat_X(:,:)
    double complex, intent(out) :: mat_AX(:,:)

    integer :: num_columns, naez
    integer :: natbld, xdim, ydim, zdim

    natbld = self%solver_opts%natbld
    xdim = self%solver_opts%xdim
    ydim = self%solver_opts%ydim
    zdim = self%solver_opts%zdim

    num_columns = size(mat_X, 2)
    naez = self%cluster_info%naez_trc

    ASSERT(naez == natbld*xdim*ydim*zdim)

    call APPBLCKCIRC  (mat_X, mat_AX, self%GLLHBLCK, &
                       naez, self%lmmaxd, &
                       natbld, xdim, ydim, zdim, num_columns)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroy_BCPOperator(self)
    class(BCPOperator) :: self

    if (allocated(self%GLLHBLCK)) deallocate(self%GLLHBLCK)
    nullify(self%cluster_info)
  end subroutine

end module
