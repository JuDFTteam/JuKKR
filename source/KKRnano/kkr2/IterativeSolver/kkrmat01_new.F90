!> Multiple Scattering problem: Solving the Dyson equation for all k-points,
!> using *real-space* G_ref and the t-matrices
!>
!> Output: Brillouin-zone integrated diagonal elements of structural Green's
!> function

#include "../DebugHelpers/logging_macros.h"
#include "../DebugHelpers/test_array_log.h"
#include "../DebugHelpers/test_macros.h"

module kkrmat_new_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  use arraytest2_mod, only: !import no name here, just mention it for the module dependency
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: kkrmat01_new

  double complex, allocatable :: full(:,:)
  double complex, parameter :: zero=(0.d0, 0.d0), cone=(1.d0, 0.d0)

  contains

  !------------------------------------------------------------------------------
  !> Solves multiple scattering problem for every k-point.
  !>
  !> Returns diagonal k-integrated part of Green's function in GS.
  subroutine kkrmat01_new(solver, kkr_op, preconditioner, Bzkp, NofKs, volcub, GS, tmatLL, alat, nsymat, RR, &
                          Ginp, lmmaxd, trunc2atom_index, communicator, iguess_data)
    !   performs k-space integration,
    !   determines scattering path operator (g(k,e)-t**-1)**-1 and
    !   Greens function of the real system -> GS(*,*,*,*),
    USE_LOGGING_MOD
    USE_ARRAYLOG_MOD
    use InitialGuess_mod, only: InitialGuess, iguess_set_k_ind
    use jij_calc_mod, only: global_jij_data, kkrjij
    use SolverStats_mod, only: SolverStats, reset
    use TFQMRSolver_mod, only: TFQMRSolver
    use BCPOperator_mod, only: BCPOperator
    use KKROperator_mod, only: KKROperator

    class(TFQMRSolver), intent(inout) :: solver
    class(KKROperator), intent(inout) :: kkr_op
    class(BCPOperator), intent(inout) :: preconditioner
    
    double precision, intent(in) :: Bzkp(:,:) !< list of k-points
    integer, intent(in) :: NofKs !< number of k-points
    double precision, intent(in) :: volcub(:) !< k-point weights

    double complex, intent(out) ::  GS(:,:,:,:) ! (lmmaxd,lmmaxd,nsymat,num_local_atoms)
    double complex, intent(in) :: tmatLL(:,:,:) ! (lmmaxd,lmmaxd,naez)
    double precision, intent(in) :: alat
    integer, intent(in) :: nsymat
    double precision, intent(in) :: RR(:,0:)
    double complex, intent(inout) :: Ginp(:,:,:,:) ! (lmmaxd,lmmaxd,naclsd,nclsd)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: trunc2atom_index(:)
    integer, intent(in) :: communicator
    type(InitialGuess), intent(inout) :: iguess_data

    ! locals
    double complex, allocatable :: G_diag(:,:,:)
    integer :: site_lm_size, ilocal, num_local_atoms, naclsd, naez, k_point_index
    type(SolverStats) :: stats

#define ms kkr_op%ms
#define cluster_info ms%cluster_info

    ! array dimensions
    naez = cluster_info%naez_trc
    naclsd = cluster_info%naclsd

    site_lm_size = naez*LMMAXD

    num_local_atoms = size(ms%atom_indices)

    !-----------------------------------------------------------------------
    ! Allocate arrays
    !-----------------------------------------------------------------------

    allocate(G_diag(lmmaxd,lmmaxd,num_local_atoms))

    ! WARNING: Symmetry assumptions might have been used that are
    ! not valid in cases of non-local potential (e.g. for Spin-Orbit coupling)
    ! ---> use sit
    !      G(n,n',L,L')(-k) = G(n',n,L',L)(k)

    GS = zero ! init zero

    TESTARRAYLOG(3, Ginp)

    if (global_jij_data%do_jij_calculation) global_jij_data%GSXIJ = zero

    call solver%reset_stats()

    !==============================================================================
    do k_point_index = 1, NofKs ! K-POINT-LOOP
    !==============================================================================

      WRITELOG(4, *) "k-point ", k_point_index

      ! select right slot for storing initial guess
      call iguess_set_k_ind(iguess_data, k_point_index)

      ! Get the scattering path operator for k-point Bzkp(:, k_point_index)
      ! output: ms%mat_X

      call kloopbody(solver, kkr_op, preconditioner, Bzkp(:,k_point_index), tmatLL, Ginp, alat, RR, trunc2atom_index, communicator, iguess_data)

      call getGreenDiag(G_diag, ms%mat_X, ms%atom_indices, ms%sparse%kvstr)

      ! TODO: use mat_X to calculate Jij

      ! ----------- Integrate Scattering Path operator over k-points --> GS -----
      ! Note: here k-integration only in irreducible wedge
      call greenKSummation(G_diag, GS, volcub(k_point_index), num_local_atoms, nsymat, lmmaxd)
      ! -------------------------------------------------------------------------

      if (global_jij_data%do_jij_calculation) then
        ! communicate off-diagonal elements and multiply with exp-factor
        call KKRJIJ(Bzkp(:,k_point_index), volcub(k_point_index), nsymat, naez, ms%atom_indices(1), &
                    global_jij_data%NXIJ, global_jij_data%IXCP,global_jij_data%ZKRXIJ, &
                    ms%mat_X, global_jij_data%GSXIJ, communicator, lmmaxd, global_jij_data%nxijd)
      endif ! jij

      do ilocal = 1, num_local_atoms
        TESTARRAYLOG(3, GS(:,:,:,ilocal))
      enddo ! ilocal

    !==============================================================================
    enddo ! k_point_index = 1, NofKs
    !==============================================================================

    
    deallocate(G_diag, stat=ilocal) ! Cleanup

    stats = solver%get_stats()

    WRITELOG(3, *) "Max. TFQMR residual for this E-point: ", stats%max_residual
    WRITELOG(3, *) "Max. num iterations for this E-point: ", stats%max_iterations
    WRITELOG(3, *) "Sum of iterations for this E-point:   ", stats%sum_iterations

#undef cluster_info
#undef ms
  endsubroutine ! kkrmat01_new



  !------------------------------------------------------------------------------
  !> Copy the diagonal elements G_{LL'}^{nn'} of the Green's-function,
  !> dependent on (k,E) into matrix G_diag
  subroutine getGreenDiag(G_diag, mat_X, atom_indices, kvstr)
    double complex, intent(out) :: G_diag(:,:,:) ! dim lmmaxd*lmmaxd*num_local_atoms
    double complex, intent(in) :: mat_X(:,:)
    integer, intent(in) :: atom_indices(:)
    integer, intent(in) :: kvstr(:)

    integer :: atom_index, lmmax1, start, ii !< local atom index

    !                                      nn
    !         Copy the diagonal elements G_LL' of the Green's-function,
    !         dependent on (k,E) into matrix G_diag
    !         (n = n' = atom_index)

    ASSERT(size(atom_indices) == size(G_diag, 3))

    G_diag = zero

    do ii = 1, size(atom_indices)

      atom_index = atom_indices(ii)

      start  = kvstr(atom_index) - 1
      lmmax1 = kvstr(atom_index+1) - kvstr(atom_index)

      ASSERT(lmmax1 == size(G_diag, 1))
      ASSERT(lmmax1 == size(G_diag, 2))

      G_diag(1:lmmax1,1:lmmax1,ii) = mat_X(start+ 1:lmmax1 +start,(ii-1)*lmmax1+1:ii*lmmax1)

    enddo ! ii

  endsubroutine ! getGreenDiag


  !------------------------------------------------------------------------------
  !> Summation of Green's function over k-points. Has to be called for every k-point
  !> TODO: it would be better to do the k-space-symmetry treatment separately ???
  !> this routine creates nsymat copies of the same solution
  !> set GS to 0 before first call
  !> in: gllke1
  !> inout: GS (set to 0 before first call)
  subroutine greenKSummation(G_diag, GS, k_point_weight, natoms, nsymat, lmmaxd)
    double complex, intent(in) :: G_diag(lmmaxd,lmmaxd,natoms)
    double complex, intent(inout) :: GS(lmmaxd,lmmaxd,nsymat,natoms)
    double precision, intent(in) :: k_point_weight
    integer, intent(in) :: lmmaxd, natoms, nsymat

    integer :: isym, ila !< local atom index

    do ila = 1, natoms
      ! perform the Brillouin-zone integration for diagonal block entries of the Green's function
      do isym = 1, nsymat
        GS(:,:,isym,ila) = GS(:,:,isym,ila) + k_point_weight*G_diag(:,:,ila)
      enddo ! isym
    enddo ! ila

  endsubroutine ! greenKSummation

  !------------------------------------------------------------------------------
  !> Calculate scattering path operator for 'kpoint'.
  !>
  !> Input are the \Delta T and the realspace G_ref (Ginp).
  !> Solution is stored in ms%mat_X.
  !> Scattering path operator is calculated for atoms given in
  !> ms%atom_indices(:)
  subroutine kloopbody(solver, kkr_op, preconditioner, kpoint, tmatLL, Ginp, alat, RR, trunc2atom_index, communicator, iguess_data)
    use fillKKRMatrix_mod, only: buildKKRCoeffMatrix, buildRightHandSide, solveFull, convertToFullMatrix
    use fillKKRMatrix_mod, only: dump
    use TFQMRSolver_mod, only: TFQMRSolver, solve
    use SparseMatrixDescription_mod, only: dump
    use InitialGuess_mod, only: InitialGuess, load, store
    use TEST_lcutoff_mod, only: cutoffmode, DEBUG_dump_matrix
    use KKROperator_mod, only: KKROperator
    use BCPOperator_mod, only: BCPOperator, calc

    USE_ARRAYLOG_MOD
    USE_LOGGING_MOD

    class(TFQMRSolver), intent(inout) :: solver
    class(KKROperator), intent(inout) :: kkr_op
    class(BCPOperator), intent(inout) :: preconditioner
    double precision, intent(in) :: kpoint(3)
    double complex, intent(in) :: tmatLL(:,:,:)
    double complex, intent(in) :: Ginp(:,:,:,:) !> Ginp(lmmaxd,lmmaxd,naclsd,nclsd) independent of the kpoint, ToDo: try not to communicate it again for every kpoint
    double precision, intent(in) :: alat
    double precision, intent(in)  :: RR(:,0:)
    integer, intent(in) :: trunc2atom_index(:)
    integer, intent(in) :: communicator
    type(InitialGuess), intent(inout) :: iguess_data

    integer :: naez, lmmaxd, n, ist
    logical :: initial_zero

#define ms kkr_op%ms
#define cluster ms%cluster_info
    lmmaxd = ms%lmmaxd
    naez = ms%naez

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

    !TODO: use an alternative referenceFourier to work around a bug on BGQ
    !call referenceFourier_com(ms%GLLh, ms%sparse, kpoint, alat, &

  ! if the following macro is defined, don't use MPI RMA locks
  ! not using locks does not scale well
  
  
    call referenceFourier_com(ms%GLLh, ms%sparse, kpoint, alat, &
             cluster%nacls_trc, cluster%atom_trc,  cluster%numn0_trc, cluster%indn0_trc, &
             rr, cluster%ezoa_trc, Ginp, trunc2atom_index, communicator)

    TESTARRAYLOG(3, ms%GLLh)

    !----------------------------------------------------------------------------
    call buildKKRCoeffMatrix(ms%GLLh, tmatLL, ms%lmmaxd, naez, ms%sparse)
    !----------------------------------------------------------------------------

    TESTARRAYLOG(3, ms%GLLh)

    ! ==> now GLLh holds (1 - Delta_t * G_ref)

    ! Now solve the linear matrix equation A*X = b (b is also a matrix),
    ! where A = (1 - Delta_t*G_ref) (inverse of scattering path operator)
    ! and b = Delta_t

    !===================================================================
    ! 3) solve linear set of equations by iterative TFQMR scheme
    !    solve (1 - \Delta t * G_ref) X = \Delta t
    !    the solution X is the scattering path operator

    ! call buildRightHandSide(ms%mat_B, lmmaxd, ms%atom_indices, ms%sparse%kvstr, tmatLL=tmatLL) ! construct RHS with t-matrices
    call buildRightHandSide(ms%mat_B, lmmaxd, ms%atom_indices, ms%sparse%kvstr) ! construct RHS as negative unity

    initial_zero = .true.
    if (iguess_data%iguess == 1) then
      initial_zero = .false.
      call load(iguess_data, ms%mat_X)
    endif

    call solver%set_initial_zero(initial_zero)

    call calc(preconditioner, ms%GLLh) ! calculate preconditioner from sparse matrix data

    if (cutoffmode == 3 .or. cutoffmode == 0) then

      call solve(solver, ms%mat_X, ms%mat_B) ! use iterative solver

      if (DEBUG_dump_matrix) then
        call dump(ms%sparse, "matrix_descriptor.dat") ! SparseMatrixDescription
        call dump(ms%GLLh,  "matrix.unf", formatted=.false.)
        call dump(ms%GLLh,  "matrix_form.dat", formatted=.true.)
        call dump(ms%mat_X, "solution.unf", formatted=.false.)
        call dump(ms%mat_X, "solution_form.dat", formatted=.true.)
        call dump(ms%mat_B, "rhs.unf", formatted=.false.)
        call dump(ms%mat_B, "rhs_form.dat", formatted=.true.)
      endif ! DEBUG_dump_matrix
      
    endif ! cutoffmode in {0,3}

    TESTARRAYLOG(3, ms%mat_B)

    ! ALTERNATIVE: direct solution with LAPACK
    if (cutoffmode == 4) then
      n = size(ms%mat_B,1)
      if (any(shape(full) /= [n,n])) then
        deallocate(full, stat=ist)
        allocate(full(n,n), stat=ist)
        if (ist /= 0) die_here("failed to allocate dense matrix with"+(n*n*.5**16)+"MiByte!")
      endif
      call convertToFullMatrix(ms%GLLh, ms%sparse%ia, ms%sparse%ja, ms%sparse%ka, ms%sparse%kvstr, ms%sparse%kvstr, full)
      TESTARRAYLOG(3, full)
      call solveFull(full, ms%mat_B, ms%mat_X)
    endif ! cutoffmode == 4

    ! store the initial guess in previously selected slot (selected with 'iguess_set_k_ind')
    call store(iguess_data, ms%mat_X)

    TESTARRAYLOG(3, ms%mat_X)
    
    ! RESULT: mat_X
#undef cluster
#undef ms
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
  subroutine referenceFourier_com(GLLh, sparse, kpoint, alat, nacls, atom, numn0, &
                indn0, rr, ezoa, Ginp, trunc2atom_index, communicator)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    use ChunkIndex_mod, only: getRankAndLocalIndex
    use one_sided_commZ_mod, only: exposeBufferZ, copyChunksNoSyncZ, hideBufferZ
#ifdef NO_LOCKS_MPI
    use one_sided_commZ_mod, only: fenceZ
#endif
    include 'mpif.h'
    double complex, intent(out) :: GLLh(:)
    type(SparseMatrixDescription), intent(in) :: sparse
    double precision, intent(in) :: kpoint(3)
    double precision, intent(in) :: alat
    integer, intent(in) :: nacls(:)
    integer, intent(in) :: atom(:,:)
    integer, intent(in) :: numn0(:)
    integer, intent(in) :: indn0(:,:)
    double precision, intent(in) :: rr(:,0:)
    integer, intent(in) :: ezoa(:,:)
    double complex, intent(in) :: Ginp(:,:,:,:)
    integer, intent(in) :: trunc2atom_index(:) !> mapping trunc. index -> global atom index
    integer, intent(in) :: communicator

    ! locals
    integer :: site_index, naez, naclsd, lmmaxd, ist
    integer :: num_local_atoms, atom_requested
    double complex, allocatable :: Gref_buffer(:,:,:), eikrm(:), eikrp(:) ! dim: naclsd
    integer(kind=4) :: chunk_inds(2,1)
    integer :: win, nranks, ierr
    integer :: naez_max
    
    num_local_atoms = size(Ginp, 4)
    
    if (num_local_atoms == 1) then
      ! this version using point-to-point MPI communication is so far only implemented for 1 atom per process
      call referenceFourier_mpi(GLLh, sparse, kpoint, alat, nacls, atom, numn0, &
                indn0, rr, ezoa, Ginp, trunc2atom_index, communicator)
      return
    endif

    naez = size(nacls)
    ASSERT(naez == size(trunc2atom_index))
    lmmaxd = size(Ginp, 1)
    ASSERT(lmmaxd == size(Ginp, 2))
    naclsd = size(Ginp, 3)
    
    allocate(eikrm(naclsd), eikrp(naclsd))

    ! Note: some MPI implementations might need the use of MPI_Alloc_mem
    allocate(Gref_buffer(lmmaxd,lmmaxd,naclsd))

    call MPI_Comm_size(communicator, nranks, ierr)

#ifdef NO_LOCKS_MPI
    ! get maximum number of atoms of all truncation zones
    call MPI_Allreduce(naez, naez_max, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)
#else
    naez_max = naez
#endif

    ! share Ginp with all other processes in 'communicator'
    call exposeBufferZ(win, Ginp, lmmaxd*lmmaxd*naclsd*num_local_atoms, lmmaxd*lmmaxd*naclsd, communicator)

    GLLh = zero ! init

    ! loop up to naez_max to ensure that each rank does the same amount of fence calls
    do site_index = 1, naez_max

#ifdef NO_LOCKS_MPI
      call fenceZ(win)
      if (site_index <= naez) then
#endif
        call dlke1(alat, nacls(site_index), rr, ezoa(:,site_index), kpoint, eikrm, eikrp)

        ! get Ginp(:,:,:)[trunc2atom_index(site_index)]

        atom_requested = trunc2atom_index(site_index)
!         chunk_inds(1)%owner = getOwner(atom_requested, num_local_atoms*nranks, nranks)
!         chunk_inds(1)%local_ind = getLocalInd(atom_requested, num_local_atoms*nranks, nranks)
        chunk_inds(:,1) = getRankAndLocalIndex(atom_requested, num_local_atoms*nranks, nranks)

#ifdef NO_LOCKS_MPI
        call copyChunksNoSyncZ(Gref_buffer, win, chunk_inds, lmmaxd*lmmaxd*naclsd)
      endif ! site_index <= naez
      
      call fenceZ(win) ! ensures that data has arrived in Gref_buffer
      
      if (site_index <= naez) then
#else

#ifdef IDENTICAL_REF
      Gref_buffer(:,:,:) = Ginp(:,:,:,1) ! use this if all Grefs are the same
#else
!       call MPI_Win_Lock(MPI_LOCK_SHARED, chunk_inds(1)%owner, 0, win, ierr)
      call MPI_Win_Lock(MPI_LOCK_SHARED, chunk_inds(1,1), 0, win, ierr)
      CHECKASSERT(ierr == 0)

      call copyChunksNoSyncZ(Gref_buffer, win, chunk_inds, lmmaxd*lmmaxd*naclsd)

!       call MPI_Win_Unlock(chunk_inds(1)%owner, win, ierr)
      call MPI_Win_Unlock(chunk_inds(1,1), win, ierr)
      CHECKASSERT(ierr == 0)
#endif

      if (.true.) then
#endif
        call dlke0_smat(site_index, GLLh, sparse%ia, sparse%ka, sparse%kvstr, eikrm, eikrp, &
                        nacls(site_index), atom(:,site_index), numn0, indn0, Gref_buffer, naez, lmmaxd)
      endif ! site_index in bounds
      
    enddo ! site_index

    call hideBufferZ(win)

    deallocate(Gref_buffer, eikrm, eikrp, stat=ist)

  endsubroutine ! referenceFourier_com


  
  
  
  subroutine referenceFourier_mpi(GLLh, sparse, kpoint, alat, nacls, atom, numn0, &
                indn0, rr, ezoa, Ginp, trunc2atom_index, comm)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    include 'mpif.h'
    double complex, intent(out) :: GLLh(:)
    type(SparseMatrixDescription), intent(in) :: sparse
    double precision, intent(in) :: kpoint(3)
    double precision, intent(in) :: alat
    integer, intent(in) :: nacls(:)
    integer, intent(in) :: atom(:,:)
    integer, intent(in) :: numn0(:)
    integer, intent(in) :: indn0(:,:)
    double precision, intent(in) :: rr(:,0:)
    integer, intent(in) :: ezoa(:,:)
    double complex, intent(in) :: Ginp(:,:,:,:)
    integer, intent(in) :: trunc2atom_index(:) !> mapping trunc. index -> atom index
    integer, intent(in) :: comm

    ! locals
    integer :: site_index, naez, naclsd, lmmaxd, ist
    integer :: num_local_atoms, atom_requested
    double complex, allocatable :: Gref_buffer(:,:,:,:), eikrm(:), eikrp(:) ! dim: naclsd
    integer :: rank, tag, myrank, nranks, ierr, ncount
    integer, allocatable :: reqs(:,:), stats(:,:,:)

    naez = size(nacls)
    ASSERT(naez == size(trunc2atom_index))
    lmmaxd = size(Ginp, 1)
    ASSERT(lmmaxd == size(Ginp, 2))
    naclsd = size(Ginp, 3)
    num_local_atoms = size(Ginp, 4)
    
    ASSERT(num_local_atoms == 1) ! only 1 atom per MPI process
    
    allocate(eikrm(naclsd), eikrp(naclsd))
    
    call MPI_Comm_size(comm, nranks, ierr)
    call MPI_Comm_rank(comm, myrank, ierr)

    GLLh = zero ! init
    
#ifndef IDENTICAL_REF
    ! Note: some MPI implementations might need the use of MPI_Alloc_mem
    allocate(Gref_buffer(lmmaxd,lmmaxd,naclsd,naez))
    allocate(reqs(2,naez), stats(MPI_STATUS_SIZE,2,naez))
    reqs(:,:) = MPI_REQUEST_NULL

    ncount = lmmaxd*lmmaxd*naclsd

    ! loop up to naez sending the information
    do site_index = 1, naez
      atom_requested = trunc2atom_index(site_index) ! get the global atom id

      rank = (atom_requested - 1)/num_local_atoms ! block distribution of atoms to ranks

      if (rank /= myrank) then

        tag  = modulo(myrank, 2**15)
        call MPI_Isend(Ginp(:,:,:,1),                 ncount, MPI_DOUBLE_COMPLEX, rank, tag, comm, reqs(1,site_index), ierr)

        tag = modulo(atom_requested - 1, 2**15)
        call MPI_Irecv(Gref_buffer(:,:,:,site_index), ncount, MPI_DOUBLE_COMPLEX, rank, tag, comm, reqs(2,site_index), ierr)

      else
        reqs(:,site_index) = MPI_REQUEST_NULL
        Gref_buffer(:,:,:,site_index) = Ginp(:,:,:,1) ! copy locally
      endif ! distant rank

    enddo ! site_index

    call MPI_Waitall(2*naez, reqs, stats, ierr) ! wait until all sends and all receives have finished
#endif
    do site_index = 1, naez
    
      call dlke1(alat, nacls(site_index), rr, ezoa(:,site_index), kpoint, eikrm, eikrp)

      call dlke0_smat(site_index, GLLh, sparse%ia, sparse%ka, sparse%kvstr, eikrm, eikrp, &
                        nacls(site_index), atom(:,site_index), numn0, indn0, &
#ifndef IDENTICAL_REF
                        Gref_buffer(:,:,:,site_index), &
#else
                        Ginp(:,:,:,1), &
#endif
                        naez, lmmaxd)

    enddo ! site_index
    
    deallocate(Gref_buffer, eikrm, eikrp, stat=ist)

  endsubroutine ! referenceFourier_mpi
  
  
  
  
  !     >> Input parameters
  !>    @param     alat    lattice constant a
  !>    @param     nacls   number of atoms in cluster
  !>    @param     RR      array of real space vectors
  !>    @param     EZOA
  !>    @param     Bzkp
  !>    @param     nrd     There are nrd+1 real space vectors in RR
  !>    @param     naclsd. maximal number of atoms in cluster

  !     << Output parameters
  !>    @param     eikrm   Fourier exponential factor with minus sign
  !>    @param     eikrp   Fourier exponential factor with plus sign
  subroutine dlke1(alat, nacls, rr, ezoa, Bzkp, eikrm, eikrp)
  use Constants_mod, only: pi
    ! ----------------------------------------------------------------------
    !     Fourier transformation of the cluster Greens function
    !     Prepares the calculation (calculates Fourier factors) for dlke0
    ! ----------------------------------------------------------------------
    double precision, intent(in) :: alat
    integer, intent(in) :: nacls !< number of vectors in the cluster
    integer, intent(in) :: ezoa(1:) !< index list of ...
    double precision, intent(in) :: Bzkp(1:3) !< k-point (vector in the Brillouin zone)
    double precision, intent(in) :: rr(1:,0:) !< dim(1:3,0:nrd) real space cluster vectors
    double complex, intent(out) :: eikrp(:), eikrm(:)
    
    double complex, parameter :: ci=(0.d0,1.d0)
    double precision :: convpuh, tpi
    double complex :: tt, exparg
    integer :: iacls

    tpi = 2.d0*pi
    convpuh = alat/tpi * 0.5d0

    do iacls = 1, nacls
       
  !     Here we do   --                  nn'
  !                  \                   ii'          ii'
  !                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
  !                  --          n'  n   LL'          LL'
  !                  n'
  !  Be careful about the minus sign included here. RR is not
  !  symmetric around each atom. The minus comes from the fact that
  !  the repulsive potential GF is calculated for 0n and not n0!                   
  
      tt = -ci*tpi*dot_product(Bzkp(1:3), rr(1:3,ezoa(iacls))) ! purely imaginary number

  !  convert to p.u. and multiply with 1/2 (done above)
      exparg = exp(tt)
      eikrp(iacls) =       exparg  * convpuh
      eikrm(iacls) = conjg(exparg) * convpuh ! we can re-use exparg here instead of exp(-tt) since tt is purely imaginary
    enddo ! iacls
  
  endsubroutine ! dlke1
  
  
  subroutine dlke0_smat(ind, smat, ia, ka, kvstr, eikrm, eikrp, nacls, atom, numn0, indn0, Ginp, naez, lmmaxd)
    integer, intent(in) :: ind !> site_index
    double complex, intent(inout) :: smat(:)
    integer, intent(in) :: ia(:)
    integer, intent(in) :: ka(:)
    integer, intent(in) :: kvstr(:)
    double complex, intent(in) :: eikrm(nacls), eikrp(nacls)
    integer, intent(in) :: nacls, atom(nacls) 
    integer, intent(in) :: numn0(naez)
    integer, intent(in) :: indn0(naez,nacls)
    integer, intent(in) :: naez, lmmaxd
    
    double complex, intent(in) :: Ginp(lmmaxd,lmmaxd,nacls)
    
    integer :: jat, lm1, lm2, iacls, ni, jnd, lmmax1, lmmax2, is

    do iacls = 1, nacls
      jat = atom(iacls)
      if (jat < 1) cycle

      do ni = 1, numn0(ind)
        jnd = indn0(ind,ni)
        if (jat == jnd) then

          lmmax1 = kvstr(ind+1) - kvstr(ind)
          lmmax2 = kvstr(jnd+1) - kvstr(jnd)

          do lm2 = 1, lmmax2
            do lm1 = 1, lmmax1

              is = ka(ia(ind) + ni-1) + lmmax1*(lm2-1) + (lm1-1)

              smat(is) = smat(is) + eikrm(iacls) * Ginp(lm2,lm1,iacls)

            enddo ! lm1
          enddo ! lm2

        endif ! jat == jnd
      enddo ! ni

      do ni = 1, numn0(jat)
        jnd = indn0(jat,ni)
        if (ind == jnd) then

          lmmax1 = kvstr(jat+1) - kvstr(jat)
          lmmax2 = kvstr(jnd+1) - kvstr(jnd)

          do lm2 = 1, lmmax2
            do lm1 = 1, lmmax1

              is = ka(ia(jat) + ni-1) + lmmax1*(lm2-1) + (lm1-1)

              smat(is) = smat(is) + eikrp(iacls) * Ginp(lm1,lm2,iacls)

            enddo ! lm1
          enddo ! lm2

        endif ! ind == jnd
      enddo ! ni

    enddo ! iacls

  endsubroutine ! dlke0_smat

endmodule ! kkrmat_new_mod
