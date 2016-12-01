!> Multiple Scattering problem: Solving the Dyson equation for all k-points,
!> using *real-space* G_ref and the t-matrices
!>
!> Output: Brillouin-zone integrated diagonal elements of structural Green's
!> function

#include "../DebugHelpers/logging_macros.h"
#include "../DebugHelpers/test_array_log.h"
#include "../DebugHelpers/test_macros.h"

#define SPLIT_REFERENCE_FOURIER_COM

module KKRmat_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  use arraytest2_mod, only: !import no name here, just mention it for the module dependency
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  
  implicit none
  private
  public :: MultipleScattering

  double complex, parameter :: zero=(0.d0, 0.d0), cone=(1.d0, 0.d0)
  

  contains

  !------------------------------------------------------------------------------
  !> Solves multiple scattering problem for every k-point.
  !>
  !> Returns diagonal k-integrated part of Green's function in GS.
  subroutine MultipleScattering(solver, op, preconditioner, kpoints, nkpoints, kpointweight, GS, tmatLL, alat, nsymat, RR, &
                          Ginp, global_atom_id, communicator, xTable, iguess_data, ienergy, ispin, &
                          mssq, tr_alph, lly_grdt, volbz, global_atom_idx_lly, Lly, & ! LLY
                          solver_type, kpoint_timer, kernel_timer)
    !   performs k-space integration,
    !   determines scattering path operator (g(k,e)-t**-1)**-1 and
    !   Greens function of the real system -> GS(*,*,*),
    USE_LOGGING_MOD
    USE_ARRAYLOG_MOD
    use InitialGuess_mod, only: InitialGuess
    use jij_calc_mod, only: global_jij_data, kkrjij
    use SolverStats_mod, only: reset, GiFlops
    use IterativeSolver_mod, only: IterativeSolver
    use BCPOperator_mod, only: BCPOperator
    use KKROperator_mod, only: KKROperator
    use TimerMpi_mod, only: TimerMpi, startTimer, stopTimer
    use fillKKRMatrix_mod, only: getGreenDiag ! retrieve result
    use MPI, only: MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD!, MPI_Allreduce 
    
    use ExchangeTable_mod, only: ExchangeTable
    use two_sided_commZ_mod, only: reference_sys_com

    type(IterativeSolver), intent(inout) :: solver
    type(KKROperator), intent(inout) :: op
    type(BCPOperator), intent(inout) :: preconditioner
    
    integer, intent(in) :: nkpoints !< number of k-points
    double precision, intent(in) :: kpoints(:,:) !< list of k-points dim(3,nkpoints)
    double precision, intent(in) :: kpointweight(:) !< k-point weights dim(nkpoints)

    double complex, intent(out) :: GS(:,:,:) ! dim(lmsd,lmsd,num_local_atoms) result
    double complex, intent(in) :: tmatLL(:,:,:,0:) ! dim(lmsd,lmsd,num_trunc_atoms,0:Lly)
    double precision, intent(in) :: alat
    integer, intent(in) :: nsymat ! needed only for Jij-calculation and Lloyds formula
    double precision, intent(in) :: RR(:,0:)
    double complex, intent(inout) :: Ginp(:,:,0:,:,:) ! dim(lmmaxd,lmmaxd,0:Lly,naclsd,num_local_atoms), contains dG_ref/dE if Lly>0
    integer, intent(in) :: global_atom_id(:) ! dim(num_trunc_atoms)
    integer, intent(in) :: communicator
    type(ExchangeTable), intent(in) :: xTable
    type(InitialGuess), intent(inout) :: iguess_data
    integer, intent(in) :: ienergy, ispin

    !LLY
    double complex, intent(in)   :: mssq(:,:,:) !< inverted T-matrix
!   double complex, intent(in)   :: dtde(:,:,:) !< dT/dE
    double complex, intent(in)   :: tr_alph(:) 
    double complex, intent(out)  :: lly_grdt
    double precision, intent(in) :: volbz
    integer, intent(in)          :: global_atom_idx_lly
    integer, intent(in)          :: Lly

    integer, intent(in) :: solver_type
    type(TimerMpi), intent(inout) :: kpoint_timer, kernel_timer

    ! locals
    double complex, allocatable :: G_diag(:,:) ! dim(lmsd,lmsd)
    double complex :: bztr2, trace ! LLY
    integer :: num_local_atoms, num_trunc_atoms, ikpoint, ila, ierr, ist, lmsd, nd(5)

#ifdef SPLIT_REFERENCE_FOURIER_COM
    ! needs more memory
    double complex, allocatable :: Gref_buffer(:,:,:,:,:) ! split_reference_fourier_com uses more memory but calls the communication routine only 1x per energy point
#endif

    ! array dimensions
    num_trunc_atoms = op%cluster%naez_trc
    lmsd = size(GS, 2)
    num_local_atoms = size(op%atom_indices)

    ! WARNING: Symmetry assumptions might have been used that are
    ! not valid in cases of non-local potential (e.g. for Spin-Orbit coupling)
    ! ---> use sit
    !      G(n,n',L,L')(-k) = G(n',n,L',L)(k)

    GS = zero ! init result with zero
    bztr2 = zero ! init zero

    TESTARRAYLOG(3, Ginp)

    if (global_jij_data%do_jij_calculation) global_jij_data%GSXIJ = zero

    call reset(solver%stats)

#ifdef SPLIT_REFERENCE_FOURIER_COM

!     ! get the required reference Green functions from the other MPI processes
!     call referenceFourier_com_part1(Gref_buffer, num_trunc_atoms, Ginp, global_atom_id, communicator)

    nd = shape(Ginp) ! ToDo maybe also here, we can move the Lly dimension to be 3rd instead of two different arrays
    allocate(Gref_buffer(nd(1),nd(2),0:nd(3)-1,nd(4),num_trunc_atoms), stat=ist) ! ToDo: move 0:Lly to be the 3rd dimension
    if (ist /= 0) die_here("failed to allocate Gref_buffer with"+(product(nd(1:4))*.5**26*num_trunc_atoms)+"GiByte for reference_sys_com") 

    ! get the required reference Green functions from the other MPI processes
    call reference_sys_com(xTable, product(nd(1:4)), Ginp, Gref_buffer)
#endif

    allocate(G_diag(lmsd,lmsd))
    !==============================================================================
    do ikpoint = 1, nkpoints ! K-POINT-LOOP
    !==============================================================================
      call startTimer(kpoint_timer)

      WRITELOG(4, *) "k-point ", ikpoint

      ! Get the scattering path operator for k-point kpoints(:,ikpoint)
      ! output: op%mat_X
      call kloopbody(solver, op, preconditioner, kpoints(1:3,ikpoint), tmatLL, alat, RR, &
#ifdef SPLIT_REFERENCE_FOURIER_COM
                     Gref_buffer, &
#else
                     Ginp, global_atom_id, communicator, &
#endif
                     iguess_data, ienergy, ispin, ikpoint, &
                     mssq, bztr2, kpointweight, global_atom_idx_lly, Lly, & !LLY
                     solver_type, kernel_timer)

      do ila = 1, num_local_atoms
        call getGreenDiag(G_diag, op%mat_X, op%bsr_X, op%atom_indices(ila), ila) ! extract solution

        ! ----------- Integrate Scattering Path operator over k-points --> GS -----
        ! Note: here k-integration only in irreducible wedge
        GS(:,:,ila) = GS(:,:,ila) + kpointweight(ikpoint)*G_diag(:,:) 
        ! -------------------------------------------------------------------------
      enddo ! ila
      
      ! ToDo: use mat_X to calculate Jij -- needs to be updated
      if (global_jij_data%do_jij_calculation) then
        ! communicate off-diagonal elements and multiply with exp-factor
        ila = op%atom_indices(1) ! convert to default integer kind
        call KKRJIJ(kpoints(1:3,ikpoint), kpointweight(ikpoint), nsymat, num_trunc_atoms, ila, &
                    global_jij_data%NXIJ, global_jij_data%IXCP,global_jij_data%ZKRXIJ, &
                    op%mat_X(:,:,1), global_jij_data%GSXIJ, communicator, lmsd, global_jij_data%nxijd)
                    stop __LINE__ ! invalid argument is passed, data layout of mat_X has changed
      endif ! jij

      call stopTimer(kpoint_timer)
    !==============================================================================
    enddo ! ikpoint = 1, nkpoints
    !==============================================================================
    
#ifdef SPLIT_REFERENCE_FOURIER_COM
    deallocate(Gref_buffer, stat=ist) ! ignore status
#endif
    deallocate(G_diag, stat=ist) ! ignore status

    if (Lly == 1) then   
      bztr2 = bztr2*nsymat/volbz + tr_alph(1)
      trace = zero
      call MPI_Allreduce(bztr2, trace, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, communicator, ierr)
      lly_grdt = trace
    endif ! Lly == 1

    do ila = 1, num_local_atoms
      TESTARRAYLOG(3, GS(:,:,ila))
    enddo ! ila

    WRITELOG(3, *) "Max. TFQMR residual for this E-point: ", solver%stats%max_residual
    WRITELOG(3, *) "Max. num iterations for this E-point: ", solver%stats%max_iterations
    WRITELOG(3, *) "Sum of iterations for this E-point:   ", solver%stats%sum_iterations
    WRITELOG(2, *) "useful Floating point operations:     ", GiFlops," GiFlop"

  endsubroutine ! MultipleScattering


  !------------------------------------------------------------------------------
  !> Calculate scattering path operator for 'kpoint'.
  !>
  !> Input are the \Delta T and the realspace G_ref (Ginp).
  !> Solution is stored in op%mat_X.
  !> Scattering path operator is calculated for atoms given in
  !> op%atom_indices(:)
  subroutine kloopbody(iterative_solver, op, preconditioner, kpoint, tmatLL, alat, RR, Ginp, &
#ifndef SPLIT_REFERENCE_FOURIER_COM
                       global_atom_id, communicator, &
#endif
                       iguess_data, ienergy, ispin, ikpoint, &
                       mssq, bztr2, volcub, global_atom_idx_lly, Lly, & !LLY
                       solver_type, kernel_timer)

    use fillKKRMatrix_mod, only: buildKKRCoeffMatrix, buildRightHandSide
    use fillKKRMatrix_mod, only: convertBSRToFullMatrix! LLY , convertFullMatrixToBSR!, solveFull
    use fillKKRMatrix_mod, only: dump
    use IterativeSolver_mod, only: IterativeSolver, solve
    use DirectSolver_mod, only: DirectSolver, solve
    use SparseMatrixDescription_mod, only: dump
    use InitialGuess_mod, only: InitialGuess, load, store
    use KKROperator_mod, only: KKROperator
    use BCPOperator_mod, only: BCPOperator, calc
    use TimerMpi_mod, only: TimerMpi
    
    USE_ARRAYLOG_MOD
    USE_LOGGING_MOD

    type(DirectSolver), save :: direct_solver
    type(IterativeSolver), intent(inout) :: iterative_solver
    type(KKROperator), intent(inout) :: op
    type(BCPOperator), intent(inout) :: preconditioner
    double precision, intent(in) :: kpoint(3)
    double complex, intent(in) :: tmatLL(:,:,:,0:) !> tmatLL(lmsd,lmsd,num_trunc_atoms,0:Lly)
    double complex, intent(in) :: Ginp(:,:,0:,:,:) !> Ginp(lmmaxd,lmmaxd,0:Lly,naclsd,N) where N=num_local_atoms or num_trunc_atoms if defined(SPLIT_REFERENCE_FOURIER_COM)
    double precision, intent(in) :: alat
    double precision, intent(in)  :: RR(:,0:)
#ifndef SPLIT_REFERENCE_FOURIER_COM    
    integer, intent(in) :: global_atom_id(:) ! becomes redundant with SPLIT_REFERENCE_FOURIER_COM
    integer, intent(in) :: communicator      ! becomes redundant with SPLIT_REFERENCE_FOURIER_COM
#else     
#define referenceFourier_com(A,B,C,D,E,F,G,H,I) referenceFourier_part2(A,B,C,D,E,F,G)
#endif
    type(InitialGuess), intent(inout) :: iguess_data
    integer, intent(in) :: ienergy, ispin, ikpoint
    ! LLY
    double complex, intent(in)         :: mssq(:,:,:) !< inverted T-matrix
    double complex, intent(out)        :: bztr2
    double precision, intent(in)       :: volcub (:)
    integer, intent(in)                :: global_atom_idx_lly !< includes the global index of local atom so that atom-specific entries in global arrays can be accessed, e.g. dtde, tmatLL
    integer, intent(in)                :: Lly             !< LLY=1/0, turns Lloyd's formula on/off
    integer, intent(in) :: solver_type
    type(TimerMpi), intent(inout) :: kernel_timer
    
    double complex, allocatable :: dPdE_local(:,:), gllke_x(:,:), dgde(:,:), gllke_x_t(:,:), dgde_t(:,:), gllke_x2(:,:), dgde2(:,:) ! LLY
    double complex :: tracek, gtdPdE ! LLY

    integer :: num_trunc_atoms, nacls, alm, lmmaxd, ist, matrix_index, lm1, lm2, lm3, il1
    double complex :: cfctorinv

    cfctorinv = (cone*8.d0*atan(1.d0))/alat
    
    lmmaxd = op%lmmaxd
    num_trunc_atoms = size(tmatLL, 3) ! number of atoms in the unified truncation zones
    nacls = op%cluster%naclsd


    !=======================================================================
    ! ---> fourier transformation
    !
    !     added by h.hoehler 3.7.2002
    !                                                     n   0          n
    !     define fourier transform as g mu mu'= ( sum_n g mu mu' exp(-iKR )
    !                                   L  L'             L   L'
    !
    !                                             n   0           n
    !                                 +   sum_n g mu'mu exp(-iK(-R ))) *0.5
    !                                             L'  L
    !
    !     this operation has to be done to satisfy e.g. the point symmetry!
    !     application of fourier transformation is just an approximation
    !     for the tb system, since the transl. invariance is not satisfied.
    !
    ! The same calculation as with Lloyds formula is done all over again ???
    ! - NO! eikrm and eikrp are SWAPPED in call to DLKE0 !!!!


    ! if the following macro is defined, don't use MPI RMA locks
    ! not using locks does not scale well
  
    call referenceFourier_com(op%mat_A, op%bsr_A, kpoint, alat, RR, op%cluster, Ginp, global_atom_id, communicator)

    TESTARRAYLOG(3, op%mat_A)

    ! ToDo: merge the referenceFourier_part2 with buildKKRCoeffMatrix

    if (Lly == 1) then ! LLY
      ! Allocate additional arrays for Lloyd's formula
      alm = lmmaxd*num_trunc_atoms

      allocate(gllke_x(lmmaxd*num_trunc_atoms,lmmaxd*nacls))
      allocate(dgde(lmmaxd*num_trunc_atoms,lmmaxd*nacls))
      allocate(gllke_x_t(lmmaxd*nacls,lmmaxd*num_trunc_atoms))
      allocate(dgde_t(lmmaxd*nacls,lmmaxd*num_trunc_atoms))
      allocate(gllke_x2(lmmaxd*num_trunc_atoms,lmmaxd))
      allocate(dgde2(lmmaxd*num_trunc_atoms,lmmaxd))
      allocate(dPdE_local(lmmaxd*num_trunc_atoms,lmmaxd))
    
!     call referenceFourier_com(op%mat_A(:,:,:,Lly), op%bsr_A, kpoint, alat, RR, op%cluster, dGinp, global_atom_id, communicator)

#undef     referenceFourier_com

      TESTARRAYLOG(3, op%mat_A(:,:,:,Lly))

      call convertBSRToFullMatrix(gllke_x, op%bsr_A, op%mat_A(:,:,:,0))
      call convertBSRToFullMatrix(   dgde, op%bsr_A, op%mat_A(:,:,:,Lly))

      !--------------------------------------------------------
      ! dP(E,k)   dG(E,k)                   dT(E)
      ! ------- = ------- * T(E) + G(E,k) * -----
      !   dE        dE                       dE
  
      matrix_index = (global_atom_idx_lly - 1)*lmmaxd

      gllke_x_t = transpose(gllke_x)
      dgde_t = transpose(dgde)

      gllke_x2 = gllke_x_t(:,matrix_index+ 1:lmmaxd +matrix_index)
      dgde2    = dgde_t   (:,matrix_index+ 1:lmmaxd +matrix_index)

      dPdE_local = zero ! init

      call zgemm('n','n',alm,lmmaxd,lmmaxd,cone, dgde2,alm, tmatLL(1,1,global_atom_idx_lly,0),lmmaxd,zero, dPdE_local,alm)

      call zgemm('n','n',alm,lmmaxd,lmmaxd,cfctorinv, gllke_x2,alm, tmatLL(1,1,global_atom_idx_lly,Lly),lmmaxd,cone, dPdE_local,alm)
      !--------------------------------------------------------

    endif ! LLY
    
    ! ToDo: merge the referenceFourier_part2 with buildKKRCoeffMatrix

    !----------------------------------------------------------------------------
    call buildKKRCoeffMatrix(op%mat_A(:,:,:,0), tmatLL(:,:,:,0), op%bsr_A)
    !----------------------------------------------------------------------------

    TESTARRAYLOG(3, op%mat_A(:,:,:,0))

    ! ==> now mat_A holds (Delta_t * G_ref - 1)

    ! Now solve the linear matrix equation A*X = B (B is also a matrix),
    ! where A = (Delta_t*G_ref - 1) (inverse of scattering path operator)
    ! and B = Delta_t

    !===================================================================
    ! 3) solve linear set of equations by iterative TFQMR scheme
    !    solve (\Delta t * G_ref - 1) X = \Delta t
    !    the solution X is the scattering path operator

    ! ToDo: move buildRightHandSide out of the k-loop (independent of k-points)
    call buildRightHandSide(op%mat_B, op%bsr_B, op%atom_indices, tmatLL(:,:,:,0)) ! construct RHS with t-matrices
!   call buildRightHandSide(op%mat_B, op%bsr_B, op%atom_indices) ! construct RHS as unity matrices

    call calc(preconditioner, op%mat_A(:,:,:,0)) ! calculate preconditioner from sparse matrix data ! should be BROKEN due to variable block row format ! ToDo: check

    selectcase(solver_type)
    case (4) ! direct solution with LAPACK, should only be used for small systems

      call solve(direct_solver, op, op%mat_X)
      ! warning, the memory of the direct solver can only be freed here as this is a save variable of this procedure

    case (0, 3) ! iterative solver
      if (solver_type == 0) warn(6, "solver_type ="+solver_type+"is deprecated, please use 3")

      ! only iterative solvers need an initial guess
      if (iguess_data%prec > 0) then
        iterative_solver%initial_zero = .false.
        call load(iguess_data, op%mat_X, ik=ikpoint, is=ispin, ie=ienergy)
      else
        iterative_solver%initial_zero = .true.
      endif

      call solve(iterative_solver, kernel_timer) ! use iterative solver

#ifdef DEBUG_dump_matrix
      call dump(op%bsr_A, "matrix_descriptor.dat") ! SparseMatrixDescription
      call dump(op%mat_A,  "bin.matrix", formatted=.false.)
      call dump(op%mat_A,  "matrix_form.dat", formatted=.true.)
      call dump(op%mat_X, "bin.solution", formatted=.false.)
      call dump(op%mat_X, "solution_form.dat", formatted=.true.)
      call dump(op%mat_B, "bin.rhs", formatted=.false.)
      call dump(op%mat_B, "rhs_form.dat", formatted=.true.)
#endif

      ! store the initial guess
      call store(iguess_data, op%mat_X, ik=ikpoint, is=ispin, ie=ienergy)

    case default
      warn(6, "No solver selected! Problem is not solved, solver_type ="+solver_type)
    endselect ! solver_type
    
    TESTARRAYLOG(3, op%mat_B)
    TESTARRAYLOG(3, op%mat_X)
    ! RESULT: mat_X

    if (Lly == 1) then ! LLY
      !--------------------------------------------------------
      !                /  -1    dM  \
      ! calculate  Tr  | M   * ---- | 
      !                \        dE  /

      tracek = zero

      do lm2 = 1, lmmaxd
        do lm1 = 1, lmmaxd
          gtdPdE = zero
          do il1 = 1, num_trunc_atoms
            do lm3 = 1, lmmaxd
              gtdPdE = gtdPdE + op%mat_X(lm3,lm2,il1)*dPdE_local(lmmaxd*(il1 - 1) + lm3,lm1)
            enddo ! lm3
          enddo ! il1
          tracek = tracek + mssq(lm1,lm2,1)*gtdPdE
        enddo ! lm1
      enddo ! lm2

      bztr2 = bztr2 + tracek*volcub(ikpoint)
      !--------------------------------------------------------
      deallocate(gllke_x, dgde, gllke_x_t, dgde_t, gllke_x2, dgde2, dPdE_local, stat=ist)
    endif ! LLY

  endsubroutine ! kloopbody



  !------------------------------------------------------------------------------
  !> See H. Hoehler
  !=======================================================================
  ! ---> fourier transformation
  !
  !     added by h.hoehler 3.7.2002
  !                                                     n   0          n
  !     define fourier transform as g mu mu'= ( sum_n g mu mu' exp(-iKR )
  !                                   L  L'             L   L'
  !
  !                                             n   0           n
  !                                 +   sum_n g mu'mu exp(-iK(-R ))) *0.5
  !                                             L'  L
  !
  !     this operation has to be done to satisfy e.g. the point symmetry!
  !     application of fourier transformation is just an approximation
  !     for the tb system, since the transl. invariance is not satisfied.
  !
  ! The same calculation as with lloyds formula is done all over again ???
  ! - NO! eikrm and eikrp are SWAPPED in call to DLKE0 !!!!
  !
  !> Alternative implementation of 'referenceFourier_com' to prevent a bug that occured on
  !> BGQ regarding MPI RMA locks.
  !>
  !> Uses fence calls instead of locks.
  !> Might not perform and scale as well as referenceFourier_com
  subroutine referenceFourier_com(Grefk, sparse, kpoint, alat, RR, cluster, Ginp, global_atom_id, communicator)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    use ChunkIndex_mod, only: getRankAndLocalIndex
    use one_sided_commZ_mod, only: exposeBufferZ, copyChunksNoSyncZ, hideBufferZ
#ifdef NO_LOCKS_MPI
    use one_sided_commZ_mod, only: fenceZ
#endif
    use ClusterInfo_mod, only: ClusterInfo
    include 'mpif.h'
    double complex, intent(out) :: Grefk(:,:,:,:) ! (lmsa,lmsd,nnzb,0:Lly)
    type(SparseMatrixDescription), intent(in) :: sparse
    double precision, intent(in) :: kpoint(3)
    double precision, intent(in) :: alat
    double precision, intent(in) :: RR(:,0:)
    type(ClusterInfo), intent(in) :: cluster
    double complex, intent(in) :: Ginp(:,:,:,:,:) ! (lmmaxd,lmmaxd,0:Lly,naclsd,num_local_atoms)
    integer, intent(in) :: global_atom_id(:) !> mapping trunc. index -> global atom index
    integer, intent(in) :: communicator

    ! locals
    integer :: site_index, num_trunc_atoms, naclsd, nLly, lmmaxd, ist
    integer :: num_local_atoms, atom_requested
    double complex, allocatable :: Gref_buffer(:,:,:,:), eikRR(:,:)
    integer(kind=4) :: chunk_inds(2,1)
    integer :: win, nranks, ierr
    integer :: naez_max
    
    num_local_atoms = size(Ginp, 5)
    
    if (num_local_atoms == 1) then
      ! this version using point-to-point MPI communication is so far only implemented for 1 atom per process
      call referenceFourier_mpi(Grefk, sparse, kpoint, alat, RR, cluster, Ginp, global_atom_id, communicator)
      return
    endif

    num_trunc_atoms = cluster%naez_trc
    ASSERT(num_trunc_atoms == size(global_atom_id))
    lmmaxd = size(Ginp, 1)
    ASSERT(lmmaxd == size(Ginp, 2))
    nLly = size(Ginp, 3)
    naclsd = size(Ginp, 4)

    allocate(eikRR(0:1,0:size(RR, 2)-1))
    call Bloch_factors(alat, RR, kpoint, eikRR)
    
    ! Note: some MPI implementations might need the use of MPI_Alloc_mem
    allocate(Gref_buffer(lmmaxd,lmmaxd,0:nLly-1,naclsd))

    call MPI_Comm_size(communicator, nranks, ierr)

#ifdef NO_LOCKS_MPI
    ! get maximum number of atoms of all truncation zones
    call MPI_Allreduce(num_trunc_atoms, naez_max, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)
#else
    naez_max = num_trunc_atoms
#endif

    ! share Ginp with all other processes in 'communicator'
    call exposeBufferZ(win, Ginp, lmmaxd*lmmaxd*nLly*naclsd*num_local_atoms, lmmaxd*lmmaxd*nLly*naclsd, communicator)

    Grefk = zero ! init

    ! loop up to naez_max to ensure that each rank does the same amount of fence calls
    do site_index = 1, naez_max

#ifdef NO_LOCKS_MPI
      call fenceZ(win)
      if (site_index <= num_trunc_atoms) then
#endif

        ! get Ginp(:,:,:)[global_atom_id(site_index)]

        atom_requested = global_atom_id(site_index)
        chunk_inds(:,1) = getRankAndLocalIndex(atom_requested, num_local_atoms*nranks, nranks)

#ifdef NO_LOCKS_MPI
        call copyChunksNoSyncZ(Gref_buffer, win, chunk_inds, lmmaxd*lmmaxd*nLly*naclsd)
      endif ! site_index <= num_trunc_atoms
      
      call fenceZ(win) ! ensures that data has arrived in Gref_buffer
      
      if (site_index <= num_trunc_atoms) then
#else

#ifdef IDENTICAL_REF
      Gref_buffer(:,:,:,:) = Ginp(:,:,:,:,1) ! use this if all Grefs are the same
#else
      call MPI_Win_Lock(MPI_LOCK_SHARED, chunk_inds(1,1), 0, win, ierr)
      CHECKASSERT(ierr == 0)

      call copyChunksNoSyncZ(Gref_buffer, win, chunk_inds, lmmaxd*lmmaxd*nLly*naclsd)

      call MPI_Win_Unlock(chunk_inds(1,1), win, ierr)
      CHECKASSERT(ierr == 0)
#endif

      if (.true.) then
#endif
        call dlke0_smat(Grefk, site_index, sparse, eikRR, cluster, Gref_buffer)

      endif ! site_index in bounds
      
    enddo ! site_index

    call hideBufferZ(win)

    deallocate(Gref_buffer, eikRR, stat=ist) ! ignore status

  endsubroutine ! referenceFourier_com

  
  subroutine referenceFourier_com_part1(Gref_buffer, num_trunc_atoms, Ginp, global_atom_id, communicator)
    !! distributes the Ginp
    use ChunkIndex_mod, only: getRankAndLocalIndex
    use one_sided_commZ_mod, only: exposeBufferZ, copyChunksNoSyncZ, hideBufferZ
#ifdef NO_LOCKS_MPI
    use one_sided_commZ_mod, only: fenceZ
#endif
    include 'mpif.h'
    double complex, allocatable, intent(out) :: Gref_buffer(:,:,:,:,:) !> dim(lmmaxd,lmmaxd,0:Lly,naclsd,num_trunc_atoms)
    integer, intent(in) :: num_trunc_atoms ! ToDo: remove from interface
    double complex, intent(in) :: Ginp(:,:,:,:,:) !> dim(lmmaxd,lmmaxd,0:Lly,naclsd,num_local_atoms)
    integer, intent(in) :: global_atom_id(:) !> dim(num_trunc_atoms) mapping trunc. index -> global atom index
    integer, intent(in) :: communicator

    ! locals
    integer :: site_index, naclsd, lmmaxd, nLly, ist
    integer :: num_local_atoms, atom_requested
    integer(kind=4) :: chunk_inds(2,1)
    integer :: win, nranks, ierr
    integer :: naez_max
    
    ASSERT( num_trunc_atoms == size(global_atom_id) )
    lmmaxd          = size(Ginp, 1)
    ASSERT( lmmaxd == size(Ginp, 2) )
    nLly            = size(Ginp, 3)
    naclsd          = size(Ginp, 4)
    num_local_atoms = size(Ginp, 5)

    ! Note: some MPI implementations might need the use of MPI_Alloc_mem
    if (any(shape(Gref_buffer) /= [lmmaxd,lmmaxd,nLly,naclsd,num_trunc_atoms])) then
      deallocate(Gref_buffer, stat=ist)
      allocate(Gref_buffer(lmmaxd,lmmaxd,0:nLly-1,naclsd,num_trunc_atoms), stat=ist)
      ! Note: some MPI implementations might need the use of MPI_Alloc_mem
      if (ist /= 0) die_here("failed to allocate Gref_buffer in referenceFourier_com_part1 with"+(lmmaxd*lmmaxd*nLly*naclsd*.5**26*num_trunc_atoms)+"GiByte") 
    endif ! shape matches

    if (num_local_atoms == 1) then
      ! this version using point-to-point MPI communication is so far only implemented for 1 atom per process
!     call referenceFourier_mpi_part1(Gref_buffer, num_trunc_atoms, Ginp, global_atom_id, communicator)
!     return
      warn(6, "referenceFourier_mpi_part1 is not called")
    endif
    
    call MPI_Comm_size(communicator, nranks, ierr)

#ifdef NO_LOCKS_MPI
    ! get maximum number of atoms of all truncation zones
    call MPI_Allreduce(num_trunc_atoms, naez_max, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)
#else
    naez_max = num_trunc_atoms
#endif

    ! share Ginp with all other processes in 'communicator'
    call exposeBufferZ(win, Ginp, lmmaxd*lmmaxd*nLly*naclsd*num_local_atoms, lmmaxd*lmmaxd*nLly*naclsd, communicator)

    ! loop up to naez_max to ensure that each rank does the same amount of fence calls
    do site_index = 1, naez_max

#ifdef NO_LOCKS_MPI
      call fenceZ(win)
      if (site_index <= num_trunc_atoms) then
#endif

        atom_requested = global_atom_id(site_index)
        chunk_inds(:,1) = getRankAndLocalIndex(atom_requested, num_local_atoms*nranks, nranks)

#ifdef NO_LOCKS_MPI
        call copyChunksNoSyncZ(Gref_buffer(:,:,:,site_index), win, chunk_inds, lmmaxd*lmmaxd*nLly*naclsd)
      endif ! site_index <= num_trunc_atoms
      
      call fenceZ(win) ! ensures that data has arrived in Gref_buffer
#else

#ifdef IDENTICAL_REF
      Gref_buffer(:,:,:,:,site_index) = Ginp(:,:,:,:,1) ! use this if all Grefs are the same
#else
      call MPI_Win_Lock(MPI_LOCK_SHARED, chunk_inds(1,1), 0, win, ierr)
      CHECKASSERT(ierr == 0)

      call copyChunksNoSyncZ(Gref_buffer(:,:,:,:,site_index), win, chunk_inds, lmmaxd*lmmaxd*nLly*naclsd)

      call MPI_Win_Unlock(chunk_inds(1,1), win, ierr)
      CHECKASSERT(ierr == 0)
#endif

#endif
    enddo ! site_index

    call hideBufferZ(win)

  endsubroutine ! referenceFourier_com_part1



  subroutine referenceFourier_mpi(Grefk, sparse, kpoint, alat, RR, cluster, Ginp, global_atom_id, comm)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    use ClusterInfo_mod, only: ClusterInfo
    double complex, intent(out) :: Grefk(:,:,:,:) !> dim(lmsa,lmsd,nnzb,0:Lly)
    type(SparseMatrixDescription), intent(in) :: sparse
    double precision, intent(in) :: kpoint(3)
    double precision, intent(in) :: alat
    double precision, intent(in) :: RR(:,0:)
    type(ClusterInfo), intent(in) :: cluster
    double complex, intent(in) :: Ginp(:,:,:,:,:) !> dim(lmmaxd,lmmaxd,0:Lly,naclsd,num_local_atoms)
    integer, intent(in) :: global_atom_id(:) !> mapping trunc. index -> atom index
    integer, intent(in) :: comm

    ! locals
    integer :: num_local_atoms, num_trunc_atoms, naclsd, lmmaxd, nLly, ist
    double complex, allocatable :: Gref_buffer(:,:,:,:,:) !> dim(lmmaxd,lmmaxd,0:Lly,naclsd,num_trunc_atoms)
    integer, parameter :: TAGMOD = 2**15
    include 'mpif.h'

    num_trunc_atoms = cluster%naez_trc
    ASSERT( num_trunc_atoms == size(global_atom_id) )
#ifdef IDENTICAL_REF
    num_trunc_atoms = 1
#endif
    lmmaxd          = size(Ginp, 1)
    ASSERT( lmmaxd == size(Ginp, 2) )
    nLly            = size(Ginp, 3)
    naclsd          = size(Ginp, 4)
    num_local_atoms = size(Ginp, 5)

    ASSERT( num_local_atoms == 1 ) ! only 1 atom per MPI process

    ! Note: some MPI implementations might need the use of MPI_Alloc_mem
    allocate(Gref_buffer(lmmaxd,lmmaxd,0:nLly-1,naclsd,num_trunc_atoms))
    
    call referenceFourier_mpi_part1(Gref_buffer, Ginp, global_atom_id, comm)

    call referenceFourier_part2(Grefk, sparse, kpoint, alat, RR, cluster, Gref_buffer)

    deallocate(Gref_buffer, stat=ist) ! ignore status
  endsubroutine ! referenceFourier_mpi

  subroutine referenceFourier_part2(Grefk, sparse, kpoint, alat, RR, cluster, Gref_buffer)
    !! this operation will be performed for every k-point and does not include MPI communication
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    use ClusterInfo_mod, only: ClusterInfo
    double complex, intent(out) :: Grefk(:,:,:,:) !> dim(lmsa,lmsd,nnzb,0:Lly)
    type(SparseMatrixDescription), intent(in) :: sparse
    double precision, intent(in) :: kpoint(3)
    double precision, intent(in) :: alat
    double precision, intent(in) :: RR(:,0:)
    type(ClusterInfo), intent(in) :: cluster
    double complex, intent(in) :: Gref_buffer(:,:,:,:,:) !> dim(lmmaxd,lmmaxd,0:Lly,naclsd,num_trunc_atoms)

    ! locals
    integer :: site_index, ist
    double complex, allocatable :: eikRR(:,:) !> dim(0:1,0:nr) Bloch phase factors
    
    ASSERT ( size(Gref_buffer, 2) == size(Gref_buffer, 1) )
    ASSERT ( size(Gref_buffer, 4) >= cluster%naclsd )

    allocate(eikRR(0:1,0:size(RR, 2)-1))
    call Bloch_factors(alat, RR, kpoint, eikRR)

    Grefk = zero ! init

    ! loop up to naez_max to ensure that each rank does the same amount of fence calls
    do site_index = 1, cluster%naez_trc
       
#ifdef IDENTICAL_REF
#define SITE_INDEX 1
#endif

      call dlke0_smat(Grefk, site_index, sparse, eikRR, cluster, Gref_buffer(:,:,:,:,SITE_INDEX)) !! site_index must be ALL_CAPS
      
#ifdef IDENTICAL_REF
#undef SITE_INDEX
#endif

    enddo ! site_index

    deallocate(eikRR, stat=ist) ! ignore status

  endsubroutine ! referenceFourier_part2

  

  subroutine referenceFourier_mpi_part1(Gref_buffer, Ginp, global_atom_id, comm)
    use ChunkIndex_mod, only: getRankAndLocalIndex
    double complex, intent(out) :: Gref_buffer(:,:,:,:,:) !> dim(lmmaxd,lmmaxd,0:Lly,naclsd,num_trunc_atoms)
    double complex, intent(in)  :: Ginp(:,:,:,:,:)        !< dim(lmmaxd,lmmaxd,0:Lly,naclsd,num_local_atoms)
    integer, intent(in) :: global_atom_id(:) !> mapping trunc. index -> atom index
    integer, intent(in) :: comm

#ifdef IDENTICAL_REF
    Gref_buffer(:,:,:,:,1) = Ginp(:,:,:,:,1)
#else

    ! locals
    integer(kind=4) :: chunk_inds(2,1)
    integer :: site_index, naclsd, nLly, lmmaxd, ist
    integer :: num_local_atoms, num_trunc_atoms, atom_requested
    integer :: rank, tag, myrank, nranks, ierr, ncount
    integer, allocatable :: reqs(:,:), stats(:,:,:)
    integer, parameter :: TAGMOD = 2**15
    include 'mpif.h'

    lmmaxd          = size(Ginp, 1)
    ASSERT( lmmaxd == size(Ginp, 2) )
    nLly            = size(Ginp, 3)
    naclsd          = size(Ginp, 4)
    num_local_atoms = size(Ginp, 5)
    num_trunc_atoms = size(Gref_buffer, 5)
    ASSERT( num_trunc_atoms == size(global_atom_id) )

    ASSERT( num_local_atoms == 1 ) ! only 1 atom per MPI process

    call MPI_Comm_size(comm, nranks, ierr)
    call MPI_Comm_rank(comm, myrank, ierr)
    
    allocate(reqs(2,num_trunc_atoms), stats(MPI_STATUS_SIZE,2,num_trunc_atoms), stat=ist)
    reqs(:,:) = MPI_REQUEST_NULL

    ncount = lmmaxd*lmmaxd*nLly*naclsd

    ! ===============================================================================
    ! part 1: exchange the Ginp arrays between the MPI processes 
    ! ===============================================================================
    
    ! loop over the sites sending the information
    do site_index = 1, num_trunc_atoms
      atom_requested = global_atom_id(site_index) ! get the global atom id

!     rank = (atom_requested - 1)/num_local_atoms ! block distribution of atoms to ranks
      chunk_inds(:,1) = getRankAndLocalIndex(atom_requested, num_local_atoms*nranks, nranks)
      rank = chunk_inds(1,1)
      assert( chunk_inds(2,1) == 1 ) ! since there may only be one local atom, its local index must be one

      if (rank /= myrank) then

        tag = modulo(myrank, TAGMOD)
        call MPI_Isend(Ginp(:,:,:,:,1),                 ncount, MPI_DOUBLE_COMPLEX, rank, tag, comm, reqs(1,site_index), ierr)

        tag = modulo(atom_requested - 1, TAGMOD)
        call MPI_Irecv(Gref_buffer(:,:,:,:,site_index), ncount, MPI_DOUBLE_COMPLEX, rank, tag, comm, reqs(2,site_index), ierr)

      else

        reqs(:,site_index) = MPI_REQUEST_NULL
        Gref_buffer(:,:,:,:,site_index) = Ginp(:,:,:,:,1) ! copy locally

      endif ! distant rank

    enddo ! site_index

    call MPI_Waitall(2*num_trunc_atoms, reqs, stats, ierr) ! wait until all sends and all receives have finished
#endif
  endsubroutine ! referenceFourier_mpi_part1

  
  subroutine Bloch_factors(alat, RR, kpoint, eikRR)
  use Constants_mod, only: pi
    ! ----------------------------------------------------------------------
    !     Fourier transformation of the cluster Greens function
    !     Prepares the calculation (calculates Fourier factors) for dlke0
    ! ----------------------------------------------------------------------
    double precision, intent(in)  :: alat !< lattice constant
    double precision, intent(in)  :: kpoint(1:3) !< k-point (vector in the Brillouin zone)
    double precision, intent(in)  :: RR(1:,0:) !< dim(1:3,0:nr-1) periodic image vectors
    double complex,   intent(out) :: eikRR(0:,0:) !< dim(0=m:1=p,0:nr-1)

    ! .. locals
    double complex, parameter     :: ci=(0.d0, 1.d0)
    double precision    :: convpuh, tpi
    double complex      :: tt, exparg
    integer             :: ezoa, nr

    tpi = 2.d0*pi
    convpuh = alat/tpi * 0.5d0 ! the factor 0.5 comes in because we anti-hermitize the A operator
    nr   =  size(RR, 2) ! number of periodic image vectors
    assert( size(eikRR, 2) == nr )

    do ezoa = 0, nr-1
       
      ! Here we do       --                  nn'
      !                  \                   ii'          ii'
      !                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
      !                  --          n'  n   LL'          LL'
      !                  n'
      ! Be careful about the minus sign included here. RR is not
      ! symmetric around each atom. The minus comes from the fact that
      ! the repulsive potential GF is calculated for 0n and not n0!
  
      tt = -ci*tpi*dot_product(kpoint(1:3), RR(1:3,ezoa)) ! purely imaginary number

      ! convert to p.u. and multiply with 1/2 (done above)
      exparg = exp(tt)
      eikRR(0,ezoa) = conjg(exparg) * convpuh ! we can re-use exparg here instead of exp(-tt) since tt is purely imaginary
      eikRR(1,ezoa) =       exparg  * convpuh
    enddo ! ezoa

  endsubroutine ! Bloch_factors

  
  subroutine dlke0_smat(smat, isa, sparse, eikRR, c, Ginp)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    use ClusterInfo_mod, only: ClusterInfo
    ! assume a block sparse row matrix description
    double complex, intent(inout) :: smat(:,:,:,0:) !< dim(lmsd,lmsd,nnzb,0:Lly)
    integer, intent(in) :: isa !> local_atom_index of the source atom
    type(SparseMatrixDescription), intent(in) :: sparse
    double complex, intent(in) :: eikRR(0:,0:) ! (0=m:1=p,0:nr-1)
    double complex, intent(in) :: Ginp(:,:,0:,:) !< dim(lmmaxd,lmmaxd,0:Lly,nacls+)
    type(ClusterInfo), intent(in) :: c ! cluster
!!! members of ClusterInfo c now:
!     integer,        :: nacls(:)   !< dim(naez_trc) number of atoms in the cluster around each site
!     integer(kind=2) :: atom(:,:)  !< dim(maxval(nacls)+,naez_trc)
!     integer,        :: numn0(:)   !< dim(naez_trc)
!     integer(kind=2) :: indn0(:,:) !< dim(maxval(numn0)+,naez_trc)
!     integer(kind=2) :: ezoa(:,:)  !< dim(maxval(nacls)+,naez_trc) ! index of periodic image into array lattice_vectors%rr or array eikRR

    ! .. locals
    integer :: ita, iacls, in0, jCol, Aind, iLly
   

    ! symmetrization of the reference Green function with simultaneous Fourier transformation (applying Bloch factors)
    
    do iLly = 0, min(ubound(smat, 4), ubound(Ginp, 3))
      
      do iacls = 1, c%nacls(isa) ! loop over all atoms in the reference cluster around source atom isa
        ita =  c%atom(iacls,isa) ! local_atom_index of the target atom iacls in the cluster around source atom isa

        if (ita > 0) then ! target atoms that do not exists inside the truncation zone are not treated

          do in0 = 1,  c%numn0(isa) ! loop over the set of inequivalent atoms in the reference cluster around source atom isa
            jCol = c%indn0(in0,isa) ! local_atom_index of the inequivalent target atom
            if (jCol == ita) then ! see which one of the inequivalent atoms is hit (should only be true once)

              assert( in0 < sparse%RowStart(isa + 1) )
              Aind = sparse%RowStart(isa) - 1 + in0
              assert( jCol == sparse%ColIndex(Aind) )
              smat(:,:,Aind,iLly) = smat(:,:,Aind,iLly) + eikRR(0,c%ezoa(iacls,isa)) * transpose(Ginp(:,:,iLly,iacls))

            endif ! ita == jCol
          enddo ! in0

          do in0 = 1,  c%numn0(ita) ! loop over the set of inequivalent atoms in the reference cluster around target atom ita
            jCol = c%indn0(in0,ita) ! local_atom_index of the inequivalent target atom
            if (jCol == isa) then ! see which one of the inequivalent atoms is hit (should only be true once)

              assert( in0 < sparse%RowStart(ita + 1) )
              Aind = sparse%RowStart(ita) - 1 + in0
              assert( jCol == sparse%ColIndex(Aind) )
              smat(:,:,Aind,iLly) = smat(:,:,Aind,iLly) + eikRR(1,c%ezoa(iacls,isa)) * Ginp(:,:,iLly,iacls)

            endif ! isa == jCol
          enddo ! in0

        endif ! ita > 0

      enddo ! iacls

    enddo ! iLly

  endsubroutine ! dlke0_smat

endmodule ! KKRmat_mod
