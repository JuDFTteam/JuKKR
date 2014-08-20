module TFQMRSolver_mod
  use Solver_mod
  use OperatorT_mod
  use SolverStats_mod
  implicit none

  type TFQMRSolver
    PRIVATE
    class(OperatorT), pointer :: op => null()
    class(OperatorT), pointer :: precond => null()

    double complex, allocatable :: vecs(:,:,:)
    double complex, allocatable :: temp(:,:)

    type(SolverStats) :: stats
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
  end type

  contains

  subroutine init_solver(self, op)
    class(TFQMRSolver) :: self
    class(OperatorT), target :: op
    self%op => op
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
    use mminvmod_oop_mod
    class(TFQMRSolver) :: self
    double complex, intent(inout) :: mat_X(:,:)
    double complex, intent(inout)    :: mat_B(:,:)

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

    call mminvmod_oop(self%op, mat_X, mat_B, 1d-6, num_columns, NLEN, &
                      .true., self%stats, self%precond, self%use_precond, &
                      self%vecs, self%temp)

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
end module

