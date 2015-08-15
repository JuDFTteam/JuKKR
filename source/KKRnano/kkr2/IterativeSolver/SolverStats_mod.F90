!> Module to collect statistics of iterative solver
module SolverStats_mod
  implicit none
  private
  public :: SolverStats, sum_stats, reset_stats

  type SolverStats
    integer :: iterations = 0
    integer :: sum_iterations = 0
    integer :: max_iterations = 0
    double precision :: max_residual = 0.0d0
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  !> Add statistics of an iterative solver run to total statistics.
  subroutine sum_stats(self, total)
    type (SolverStats), intent(in) :: self
    type (SolverStats), intent(inout) :: total

    total%sum_iterations = total%sum_iterations + self%iterations
    total%max_residual = max(total%max_residual, self%max_residual)
    total%max_iterations = max(total%max_iterations, self%iterations)
  end subroutine

  !----------------------------------------------------------------------------
  !> Reset stats to zero.
  subroutine reset_stats(self)
    type (SolverStats), intent(inout) :: self
    self%iterations = 0
    self%sum_iterations = 0
    self%max_iterations = 0
    self%max_residual = 0.0d0
  end subroutine
  
end module SolverStats_mod
