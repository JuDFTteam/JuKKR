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
  
  public :: TFQMRSolver
  
  type, extends(Solver) :: TFQMRSolver
    PRIVATE
    class(OperatorT), pointer :: op => null()
    class(OperatorT), pointer :: precond => null()

    double complex, allocatable :: vecs(:,:,:)
    double complex, allocatable :: temp(:,:)

    type(SolverStats) :: stats
    type(SolverStats) :: total_stats
    logical :: use_precond = .false.
    logical :: initial_zero = .true.
    double precision :: qmrbound = 1d-6

    contains
      procedure :: init => init_solver
      procedure :: init_precond => init_precond_solver
      procedure :: solve => solve_with_solver
      procedure :: destroy => destroy_solver
      procedure :: set_qmrbound
      procedure :: set_initial_zero
      procedure :: get_stats
      procedure :: get_total_stats
      procedure :: reset_total_stats
  end type

  contains

  subroutine init_solver(self, op)
    class(TFQMRSolver) :: self
    class(OperatorT), target :: op
    self%op => op
    self%precond => op ! also init the pointer to precond with op, may be overwritten later
    self%use_precond = .false.
  end subroutine

  subroutine init_precond_solver(self, precond)
    class(TFQMRSolver) :: self
    class(OperatorT), target :: precond
    self%precond => precond
    self%use_precond = .true.
  end subroutine

  !----------------------------------------------------------------------------
  !> Solves problem for right hand side mat_B, solution in mat_X.
  !>
  !> The workspace is allocated on demand and stays allocated.
  !> Deallocate with call TFQMRSolver%destroy
  subroutine solve_with_solver(self, mat_X, mat_B)
    use mminvmod_oop_mod, only: mminvmod_oop
    use SolverStats_mod, only: sum_stats
    class(TFQMRSolver) :: self
    double complex, intent(inout) :: mat_X(:,:)
    double complex, intent(inout) :: mat_B(:,:)

    integer nlen
    integer num_columns

    nlen = size(mat_B, 1)
    num_columns = size(mat_B, 2)

    if (.not. allocated(self%vecs)) then
      allocate(self%vecs(nlen, num_columns, 7))
      allocate(self%temp(nlen, num_columns))
    else
      if (size(self%vecs, 1) /= nlen .or. size(self%vecs, 2) /= num_columns) then
        ! when problem size has changed, one should destroy the solver and create a new one
        write(*,*) "TFQMRSolver error: Problem size has changed. Cannot reuse solver."
        STOP
      endif
    endif

    if (.not. associated(self%op)) then
      write(*,*) "TFQMRSolver error: No matrix/operator set."
      STOP
    endif

    call mminvmod_oop(self%op, mat_X, mat_B, self%qmrbound, num_columns, NLEN, &
                      .true., self%stats, self%precond, self%use_precond, &
                      self%vecs, self%temp)

    call sum_stats(self%stats, self%total_stats)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroy_solver(self)
    class(TFQMRSolver) :: self

    if (allocated(self%vecs)) deallocate(self%vecs)
    if (allocated(self%temp)) deallocate(self%temp)

    nullify(self%op)
    nullify(self%precond)
  end subroutine

  !----------------------------------------------------------------------------
  subroutine set_qmrbound(self, qmrbound)
    class(TFQMRSolver) :: self
    double precision, intent(in) :: qmrbound
    self%qmrbound = qmrbound
  end subroutine

  !----------------------------------------------------------------------------
  subroutine set_initial_zero(self, initial_zero)
    class(TFQMRSolver) :: self
    logical, intent(in) :: initial_zero
    self%initial_zero = initial_zero
  end subroutine

  !----------------------------------------------------------------------------
  !> Get statistics of last solver run.
  function get_stats(self)
    class(TFQMRSolver) :: self
    type(SolverStats) :: get_stats

    get_stats = self%stats
  end function

  !----------------------------------------------------------------------------
  !> Get accumulated statistics of all solver runs.
  function get_total_stats(self)
    class(TFQMRSolver) :: self
    type(SolverStats) :: get_total_stats

    get_total_stats = self%total_stats
  end function

  !----------------------------------------------------------------------------
  !> Reset accumulated statistics.
  subroutine reset_total_stats(self)
    use SolverStats_mod, only: reset_stats
    class(TFQMRSolver) :: self

    call reset_stats(self%total_stats)
  end subroutine
  
end module

