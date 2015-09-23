module SolverOptions_mod
  implicit none
  private
  public :: SolverOptions

  !> type that provides additional options to the solver
  type SolverOptions
    ! for block-circulant preconditioner
    integer :: bcp = 0
    integer :: xdim = 1
    integer :: ydim = 1
    integer :: zdim = 1
    integer :: natbld = 1
  endtype
  
endmodule ! SolverOptions_mod
