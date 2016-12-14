!> Iterative solver.
!>

module IterativeSolver_mod
  use KKROperator_mod, only: KKROperator
  use BCPOperator_mod, only: BCPOperator
  use SolverStats_mod, only: SolverStats
  implicit none
  private
  
  public :: IterativeSolver, create, solve, destroy, free_memory
  
  type :: IterativeSolver
    type(KKROperator), pointer :: op => null()
    type(BCPOperator), pointer :: precond => null()

    double complex, allocatable :: vecs(:,:,:,:) ! dim(lmsa,lmsd,nnzb,3:9) workspace

    logical :: use_precond = .false.
    logical :: initial_zero = .true.
    double precision :: qmrbound = 1.d-9
    type(SolverStats) :: stats
  endtype

  interface create
    module procedure create_solver
  endinterface
  
  interface solve
    module procedure solve_with_solver
  endinterface

  interface free_memory
    module procedure free_memory_solver
  endinterface

  interface destroy
    module procedure destroy_solver
  endinterface

  contains

  subroutine create_solver(self, qmrbound, op, precond)
    type(IterativeSolver), intent(inout) :: self
    double precision, intent(in) :: qmrbound
    type(KKROperator), intent(in), target :: op
    type(BCPOperator), intent(in), target, optional :: precond
    ! local
    type(BCPOperator), target, save :: precond_null
    self%qmrbound = qmrbound
    self%op => op
    self%precond => precond_null
    self%use_precond = .false.
    if (present(precond)) then
      self%precond => precond
      self%use_precond = .true.
    endif ! preconditioning
  endsubroutine ! init

  !----------------------------------------------------------------------------
  !> Solves problem for right hand side mat_B, solution in mat_X.
  !>
  !> The workspace is allocated on demand and stays allocated.
  !> Deallocate with call IterativeSolver%destroy
  subroutine solve_with_solver(self, timer)
    use tfQMR_mod, only: solve
    use TimerMpi_mod, only: TimerMpi
    use SolverStats_mod, only: reduce
    type(IterativeSolver), intent(inout) :: self
    type(TimerMpi), intent(inout) :: timer
    
    integer :: ndim(4), iterations_needed, ist, nRHSs
    double precision :: largest_residual
    integer(kind=8) :: nFlops

    if (.not. associated(self%op)) stop "IterativeSolver error: No matrix/operator set."
    
    ! adopt the dimensions of mat_X
    ndim(1:3) = shape(self%op%mat_X)

    ndim(4) = 7 ; if(self%use_precond) ndim(4) = ndim(4) + 1 ! need one more for preconditioning
    if (any(shape(self%vecs) /= ndim)) then
      ! resize
      deallocate(self%vecs, stat=ist) ! ignore status
      allocate(self%vecs(ndim(1),ndim(2),ndim(3),ndim(4)), stat=ist)
      if (ist /= 0) then
        write(*, '(9(a,f0.3))') "IterativeSolver error: (Re-)Allocation of workspace failed! requested ", &
         product(ndim(1:3)/1024.)*ndim(4)*16.," GiByte" ! 16 Byte per dp-complex
        stop
      endif ! failed
    endif ! needs resize

#ifdef  TRANSPOSE_TO_ROW_MAJOR
    nRHSs = ndim(1)
#else
    nRHSs = ndim(2)
#endif

    nFlops = 0
    ! call tfQMR solver
    call solve(self%op, self%op%mat_X, self%op%mat_B, self%qmrbound, nRHSs, self%op%bsr_X%nCols, &
               self%initial_zero, self%precond, self%use_precond, &
               self%vecs, timer, iterations_needed, largest_residual, nFlops)

    call reduce(self%stats, iterations_needed, largest_residual, nFlops)
  endsubroutine ! solve

  !----------------------------------------------------------------------------
  !> The workspace is deallocated and needs to be allocated on next use
  integer(kind=8) function free_memory_solver(self) result(nBytes)
    type(IterativeSolver), intent(inout) :: self
    integer :: ist
    nBytes = 16 * size(self%vecs)
    deallocate(self%vecs, stat=ist) ! ignore status
  endfunction ! free_memory

  !----------------------------------------------------------------------------
  elemental subroutine destroy_solver(self)
    type(IterativeSolver), intent(inout) :: self
    
    integer :: ist
    deallocate(self%vecs, stat=ist) ! ignore status
    nullify(self%precond)
    nullify(self%op)
  endsubroutine ! destroy
  
endmodule ! IterativeSolver_mod

