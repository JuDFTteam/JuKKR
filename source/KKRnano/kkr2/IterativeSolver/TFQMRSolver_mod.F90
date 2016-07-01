!> TFQMR solver.
!>
!> *) Set coefficent matrix with 'init'.
!> *) Set preconditioner with 'init_precond'
!> *) Use 'set_qmrbound' and
!> 'set_initial_zero' (true: start with zero vector, false: start with vector passed to solve)
!> *) use 'solve' to solve for right hand side mat_B (memory is allocated at first use)
!> *) call 'destroy' after last use
!>
!> @author Elias Rabel

module TFQMRSolver_mod
  use Solver_mod, only: Solver
  use OperatorT_mod, only: OperatorT
  use SolverStats_mod, only: SolverStats
  implicit none
  private
  
  public :: TFQMRSolver, solve, destroy, init
  
  type, extends(Solver) :: TFQMRSolver
    class(OperatorT), pointer :: op => null()
    class(OperatorT), pointer :: precond => null()

    double complex, allocatable :: vecs(:,:,:)

    type(SolverStats) :: stats
    logical :: use_precond = .false.
    logical :: initial_zero = .true.
    double precision :: qmrbound = 1.d-9

!     contains
!       procedure :: solve => solve_with_solver
!       procedure :: init => init_solver
!       procedure :: init_precond => init_precond_solver
!       procedure :: set_qmrbound
!       procedure :: set_initial_zero
!       procedure :: get_stats
!       procedure :: reset_stats
!       procedure :: represent_stats
  endtype
  
  interface init
    module procedure init_solver
  endinterface
  
  interface solve
    module procedure solve_with_solver
  endinterface

  interface destroy
    module procedure destroy_solver
  endinterface

  contains

  subroutine init_solver(self, qmrbound, op, precond)
    type(TFQMRSolver), intent(inout) :: self
    double precision, intent(in) :: qmrbound
    class(OperatorT), target :: op
    class(OperatorT), target, optional :: precond
    self%qmrbound = qmrbound
    self%op => op
    self%precond => op ! also init the pointer to precond with op, may be overwritten later
    self%use_precond = .false.
    if (present(precond)) then
      self%precond => precond
      self%use_precond = .true.
    endif
  endsubroutine ! init

  !----------------------------------------------------------------------------
  !> Solves problem for right hand side mat_B, solution in mat_X.
  !>
  !> The workspace is allocated on demand and stays allocated.
  !> Deallocate with call TFQMRSolver%destroy
  subroutine solve_with_solver(self, mat_X, mat_B)
    use mminvmod_oop_mod, only: mminvmod_oop
    use SolverStats_mod, only: reduce
    type(TFQMRSolver) :: self
    double complex, intent(inout) :: mat_X(:,:)
    double complex, intent(in)    :: mat_B(:,:)

    integer :: nrow, ncol, iterations_needed, nvecs, ist
    double precision :: largest_residual

    nrow = size(mat_B, 1)
    ncol = size(mat_B, 2)

    nvecs = 7; if(self%use_precond) nvecs = nvecs+1 ! need only 7 without preconditioning
    if (size(self%vecs, 2) /= nvecs) then
      deallocate(self%vecs, stat=ist) ! ignore status
      
      allocate(self%vecs(nrow,ncol,nvecs), stat=ist)
      
      if (ist /= 0) then
        write(*,*) "TFQMRSolver error: Allocation of workspace failed! requested ",(nrow/1024.)*(ncol/1024.)*(nvecs/64.)," GiByte" ! 16 Byte per dp-complex
        stop
      endif
      
    else
      if (size(self%vecs, 1) /= nrow .or. size(self%vecs, 2) /= ncol) then
        ! when problem size has changed, one should destroy the solver and create a new one
        write(*,*) "TFQMRSolver error: Problem size has changed. Cannot reuse solver."
        stop
      endif
    endif

    if (.not. associated(self%op)) then
      write(*,*) "TFQMRSolver error: No matrix/operator set."
      stop
    endif

    call mminvmod_oop(self%op, mat_X, mat_B, self%qmrbound, ncol, nrow, &
                      .true., self%precond, self%use_precond, &
                      self%vecs, iterations_needed, largest_residual)

    call reduce(self%stats, iterations_needed, largest_residual)

  endsubroutine ! solve
  
  !----------------------------------------------------------------------------
  elemental subroutine destroy_solver(self)
    type(TFQMRSolver), intent(inout) :: self
    
    integer :: ist ! ignore status

    deallocate(self%vecs, stat=ist)
    nullify(self%op)
    nullify(self%precond)
  endsubroutine ! destroy

! 
!   subroutine init_precond_solver(self
!     type(TFQMRSolver) :: self
!     class(OperatorT), target :: precond
!   endsubroutine ! init with preconditioner
  
!   !----------------------------------------------------------------------------
!   subroutine set_qmrbound(self, qmrbound)
!     type(TFQMRSolver) :: self
!     double precision, intent(in) :: qmrbound
!     self%qmrbound = qmrbound
!   endsubroutine
! 
!   !----------------------------------------------------------------------------
!   subroutine set_initial_zero(self, initial_zero)
!     type(TFQMRSolver) :: self
!     logical, intent(in) :: initial_zero
!     self%initial_zero = initial_zero
!   endsubroutine
! 
!   !----------------------------------------------------------------------------
!   !> Get statistics of last solver run.
!   type(SolverStats) function get_stats(self)
!     type(TFQMRSolver) :: self
! 
!     get_stats = self%stats
!   endfunction
! 
!   !----------------------------------------------------------------------------
!   !> Reset accumulated statistics.
!   subroutine reset_stats(self)
!     use SolverStats_mod, only: reset
!     type(TFQMRSolver) :: self
! 
!     call reset(self%stats)
!   endsubroutine
! 
!   !----------------------------------------------------------------------------
!   character(len=64) function represent_stats(self) result(string)
!     use SolverStats_mod, only: represent
!     type(TFQMRSolver) :: self
!     string = represent(self%stats)
!   endfunction
  
  
endmodule ! TFQMRSolver_mod

