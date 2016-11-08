!> Iterative solver.
!>

module IterativeSolver_mod
  use KKROperator_mod, only: KKROperator
  use BCPOperator_mod, only: BCPOperator
  use SolverStats_mod, only: SolverStats
  implicit none
  private
  
  public :: IterativeSolver, create, solve, destroy
  
  type :: IterativeSolver
    type(KKROperator), pointer :: op => null()
    type(BCPOperator), pointer :: precond => null()

    double complex, allocatable :: vecs(:,:,:,:) ! workspace

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
  subroutine solve_with_solver(self)
    use TFQMR_mod, only: solve
    use SolverStats_mod, only: reduce
    type(IterativeSolver), intent(inout) :: self

    integer :: nrow, ncol, iterations_needed, nvecs, ist, nnzb
    double precision :: largest_residual
    integer(kind=8) :: nFlops

    if (.not. associated(self%op)) stop "IterativeSolver error: No matrix/operator set."
    
    ! adopt the dimensions of mat_X
#define mat_X self%op%mat_X    
#define mat_B self%op%mat_B
    nrow = size(mat_X, 1)
    ncol = size(mat_X, 2)
    nnzb = size(mat_X, 3)

    nvecs = 7 ; if(self%use_precond) nvecs = nvecs + 1 ! need one more for preconditioning
    if (any(shape(self%vecs) /= [nrow,ncol,nnzb,nvecs])) then
      ! resize
      deallocate(self%vecs, stat=ist) ! ignore status
      allocate(self%vecs(nrow,ncol,nnzb,nvecs), stat=ist)
      if (ist /= 0) then
        write(*,*) "IterativeSolver error: (Re-)Allocation of workspace failed! requested ", &
         (nrow/1024.)*(ncol/1024.)*(nnzb/1024.)*(nvecs*16.)," GiByte" ! 16 Byte per dp-complex
        stop
      endif ! failed
    endif ! needs resize

#ifdef useBSR
#define nRHSs self%op%bsr_X%nCols
#else
#define nRHSs 1
#endif

    nFlops = 0
    ! call TFQMR solver
    call solve(self%op, mat_X, mat_B, self%qmrbound, ncol, nRHSs, &
               self%initial_zero, self%precond, self%use_precond, &
               self%vecs, iterations_needed, largest_residual, nFlops)

    call reduce(self%stats, iterations_needed, largest_residual, nFlops)
  endsubroutine ! solve

  !----------------------------------------------------------------------------
  elemental subroutine destroy_solver(self)
    type(IterativeSolver), intent(inout) :: self
    
    integer :: ist ! ignore status
    deallocate(self%vecs, stat=ist)
    nullify(self%precond)
    nullify(self%op)
  endsubroutine ! destroy
  
endmodule ! IterativeSolver_mod

