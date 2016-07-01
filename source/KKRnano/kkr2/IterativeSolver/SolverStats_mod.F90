!> Module to collect statistics of iterative solver
module SolverStats_mod
  implicit none
  private
  public :: SolverStats, reduce, reset, represent

  type SolverStats
    integer :: max_iterations = 0
    integer :: sum_iterations(0:2) = 0
    double precision :: max_residual = 0.d0
  endtype
  
  interface reduce
    module procedure reduce_stats
  endinterface
  
  interface reset
    module procedure reset_stats
  endinterface

  interface represent
    module procedure repr_stats
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  !> Add statistics of an iterative solver run to total statistics.
  subroutine reduce_stats(self, iterations, residual)
    type(SolverStats), intent(inout) :: self
    integer, intent(in) :: iterations
    double precision, intent(in) :: residual
    
    self%sum_iterations = self%sum_iterations + max(1, iterations)**[0, 1, 2]
    self%max_residual = max(residual, self%max_residual)
    self%max_iterations = max(iterations, self%max_iterations)
  endsubroutine ! sum

  !----------------------------------------------------------------------------
  !> Reset stats to zero.
  subroutine reset_stats(self)
    type(SolverStats), intent(inout) :: self
    self%sum_iterations = 0
    self%max_iterations = 0
    self%max_residual = 0.d0
  endsubroutine ! reset
  
  character(len=64) function repr_stats(self) result(string)
    type(SolverStats), intent(in) :: self
    double precision :: avg_var(2)
    string = ''; if (self%sum_iterations(0) < 1) return
    avg_var(1:2) = self%sum_iterations(1:2) / dble(self%sum_iterations(0))
    avg_var(2) = sqrt(avg_var(2) - avg_var(1)**2) ! Gaussian variance
    write(unit=string, fmt='(2(a,f0.1),9(a,i0))') 'stats ',avg_var(1),' +/- ',avg_var(2),' its, max ',self%max_iterations
  endfunction ! represent
  
endmodule ! SolverStats_mod
