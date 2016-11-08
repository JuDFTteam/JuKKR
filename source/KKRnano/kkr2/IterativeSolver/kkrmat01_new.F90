!> Multiple Scattering problem: Solving the Dyson equation for all k-points,
!> using *real-space* G_ref and the t-matrices
!>
!> Output: Brillouin-zone integrated diagonal elements of structural Green's
!> function

#include "../DebugHelpers/logging_macros.h"
#include "../DebugHelpers/test_array_log.h"
#include "../DebugHelpers/test_macros.h"


#define SPLIT_REFERENCE_FOURIER_COM

module kkrmat_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  use arraytest2_mod, only: !import no name here, just mention it for the module dependency
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: kkrmat01

  double complex, allocatable :: full_A(:,:), full_X(:,:)
  double complex, parameter :: zero=(0.d0, 0.d0), cone=(1.d0, 0.d0)

  contains

  !------------------------------------------------------------------------------
  !> Solves multiple scattering problem for every k-point.
  !>
  !> Returns diagonal k-integrated part of Green's function in GS.
  subroutine kkrmat01(solver, op, preconditioner, kpoints, nkpoints, kpointweight, GS, tmatLL, alat, nsymat, RR, &
                          Ginp, global_atom_id, communicator, iguess_data, ienergy, ispin, &
                          mssq, dGinp, dtde, tr_alph, lly_grdt, volbz, global_atom_idx_lly, Lly, solver_type, kpoint_timer) ! LLY
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
    use MPI, only: MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD!, MPI_Allreduce 

    type(IterativeSolver), intent(inout) :: solver
    type(KKROperator), intent(inout) :: op
    type(BCPOperator), intent(inout) :: preconditioner
    
    integer, intent(in) :: nkpoints !< number of k-points
    double precision, intent(in) :: kpoints(:,:) !< list of k-points dim(3,nkpoints)
    double precision, intent(in) :: kpointweight(:) !< k-point weights dim(nkpoints)

    double complex, intent(out) :: GS(:,:,:) ! dim(lmmaxd,lmmaxd,num_local_atoms) result
    double complex, intent(in) :: tmatLL(:,:,:) ! dim(lmmaxd,lmmaxd,naez_trc)
    double precision, intent(in) :: alat
    integer, intent(in) :: nsymat ! needed only for Jij-calculation and Lloyds formula
    double precision, intent(in) :: RR(:,0:)
    double complex, intent(inout) :: Ginp(:,:,:,:) ! dim(lmmaxd,lmmaxd,naclsd,num_local_atoms)
    integer, intent(in) :: global_atom_id(:)
    integer, intent(in) :: communicator
    type(InitialGuess), intent(inout) :: iguess_data
    integer, intent(in) :: ienergy, ispin

    !LLY
    double complex, intent(in)   :: mssq (:,:,:)    !< inverted T-matrix
    double complex, intent(inout) :: dGinp(:,:,:,:)  !< dG_ref/dE dim(lmmaxd,lmmaxd,naclsd,num_local_atoms)
    double complex, intent(in)   :: dtde(:,:,:)     !< dT/dE
    double complex, intent(in)   :: tr_alph(:) 
    double complex, intent(out)  :: lly_grdt
    double precision, intent(in) :: volbz
    integer, intent(in)          :: global_atom_idx_lly
    integer, intent(in)          :: Lly

    integer, intent(in) :: solver_type
    type(TimerMpi), intent(inout) :: kpoint_timer

    ! locals
    double complex, allocatable :: G_diag(:,:) ! dim(lmmaxd,lmmaxd)
    double complex :: bztr2, trace ! LLY
    integer :: num_local_atoms, naez, ikpoint, ila, ierr, lmmaxd, ist!, naclsd

#ifdef SPLIT_REFERENCE_FOURIER_COM
    double complex, allocatable :: Gref_buffer(:,:,:,:) ! split_reference_fourier_com uses more memory but calls the communication routine only 1x per energy point
    double complex, allocatable :: dGref_buffer(:,:,:,:) ! LLY
#endif

    ! array dimensions
    naez = op%cluster_info%naez_trc
!   naclsd = op%cluster_info%naclsd
    lmmaxd = size(GS, 2)
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
    ! get the required reference Green functions from the other MPI processes
    call referenceFourier_com_part1(Gref_buffer, naez, Ginp, global_atom_id, communicator)
#define Ginp Gref_buffer
    if (Lly == 1) then ! LLY
      call referenceFourier_com_part1(dGref_buffer, naez, dGinp, global_atom_id, communicator)
#define dGinp dGref_buffer
    endif
#endif
    allocate(G_diag(lmmaxd,lmmaxd))
    !==============================================================================
    do ikpoint = 1, nkpoints ! K-POINT-LOOP
    !==============================================================================
      call startTimer(kpoint_timer)

      WRITELOG(4, *) "k-point ", ikpoint

      ! Get the scattering path operator for k-point kpoints(:,ikpoint)
      ! output: op%mat_X
      call kloopbody(solver, op, preconditioner, kpoints(1:3,ikpoint), tmatLL, Ginp, &
                     alat, RR, global_atom_id, communicator, iguess_data, ienergy, ispin, &
                     mssq, dtde, dGinp, bztr2, kpointweight, ikpoint, &
                     global_atom_idx_lly, Lly, & !LLY
                     solver_type)

      do ila = 1, num_local_atoms
        call getGreenDiag(G_diag, op%mat_X, op%atom_indices(ila), ila, op%bsr_X) ! extract solution

        ! ----------- Integrate Scattering Path operator over k-points --> GS -----
        ! Note: here k-integration only in irreducible wedge
        GS(:,:,ila) = GS(:,:,ila) + kpointweight(ikpoint)*G_diag(:,:) 
        ! -------------------------------------------------------------------------
      enddo ! ila
      
      ! TODO: use mat_X to calculate Jij
      if (global_jij_data%do_jij_calculation) then
        ! communicate off-diagonal elements and multiply with exp-factor
        ila = op%atom_indices(1) ! convert to default integer kind
        call KKRJIJ(kpoints(1:3,ikpoint), kpointweight(ikpoint), nsymat, naez, ila, &
                    global_jij_data%NXIJ, global_jij_data%IXCP,global_jij_data%ZKRXIJ, &
                    op%mat_X(:,:,1), global_jij_data%GSXIJ, communicator, lmmaxd, global_jij_data%nxijd)
                    stop __LINE__ ! invalid argument is passed, data layout of mat_X has changed
      endif ! jij

      call stopTimer(kpoint_timer)
    !==============================================================================
    enddo ! ikpoint = 1, nkpoints
    !==============================================================================

    if (Lly == 1) then   
      bztr2 = bztr2*nsymat/volbz + tr_alph(1)
      trace = zero
      call MPI_Allreduce(bztr2, trace, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
      lly_grdt = trace
    endif ! Lly == 1

#ifdef SPLIT_REFERENCE_FOURIER_COM
#undef Ginp
    deallocate(Gref_buffer, stat=ist) ! ignore status
#undef dGinp
    deallocate(dGref_buffer, stat=ist) ! LLY, ignore status
#endif

    do ila = 1, num_local_atoms
      TESTARRAYLOG(3, GS(:,:,ila))
    enddo ! ila
    
    WRITELOG(3, *) "Max. TFQMR residual for this E-point: ", solver%stats%max_residual
    WRITELOG(3, *) "Max. num iterations for this E-point: ", solver%stats%max_iterations
    WRITELOG(3, *) "Sum of iterations for this E-point:   ", solver%stats%sum_iterations
    WRITELOG(2, *) "useful Floating point operations:     ", GiFlops," GiFlop"
    deallocate(G_diag, stat=ist) ! ignore status

  endsubroutine ! kkrmat01



  !------------------------------------------------------------------------------
  !> Copy the diagonal elements G_{LL'}^{nn'} of the Green's-function,
  !> dependent on (k,E) into matrix G_diag
  subroutine getGreenDiag(G_diag, mat_X, atom_index, iRHS, bsr_X)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription, exists
    double complex, intent(out) :: G_diag(:,:) ! dim(lmsd,lmsd)
    double complex, intent(in) :: mat_X(:,:,:) !> dim(lmsa,lmsd*nRHS,nrows) with useBSR: dim(lmsa,lmsd,nnzb)
    integer(kind=2), intent(in) :: atom_index
    integer, intent(in) :: iRHS
    type(SparseMatrixDescription), intent(in) :: bsr_X

    integer :: lmsd, Xind, row
    !                                      nn
    !         Copy the diagonal elements G_LL' of the Green's-function,
    !         dependent on (k,E) into matrix G_diag
    !         (n = n' = atom_index)
    lmsd = size(mat_X, 1)
    ASSERT( lmsd == size(G_diag, 1))
    ASSERT( lmsd == size(G_diag, 2))
    G_diag = ZERO
#ifndef useBSR
    ASSERT( 1 <= iRHS .and. iRHS <= size(mat_X, 2)/lmsd )
    G_diag(:,:) = mat_X(:lmsd,(iRHS - 1)*lmsd + 1:iRHS*lmsd,atom_index)
#else
    row = atom_index
    Xind = exists(bsr_X, row, col=iRHS)
    if (Xind < 1) die_here("diagonal element not contained in X, row="-row-", col="-iRHS)
    G_diag(:,:) = mat_X(:lmsd,:,Xind)
#endif
  endsubroutine ! getGreenDiag

  !------------------------------------------------------------------------------
  !> Calculate scattering path operator for 'kpoint'.
  !>
  !> Input are the \Delta T and the realspace G_ref (Ginp).
  !> Solution is stored in op%mat_X.
  !> Scattering path operator is calculated for atoms given in
  !> op%atom_indices(:)
  subroutine kloopbody(solver, op, preconditioner, kpoint, tmatLL, Ginp, alat, RR, global_atom_id, communicator, iguess_data, ienergy, ispin, &
                       mssq, dtde, dGinp, bztr2, volcub, ikpoint, &
                       global_atom_idx_lly, Lly, solver_type) !LLY

    use fillKKRMatrix_mod, only: buildKKRCoeffMatrix, buildRightHandSide, solveFull, convertBSRToFullMatrix, convertFullMatrixToBSR
    use fillKKRMatrix_mod, only: dump
    use IterativeSolver_mod, only: IterativeSolver, solve
    use SparseMatrixDescription_mod, only: dump
    use InitialGuess_mod, only: InitialGuess, load, store
    use KKROperator_mod, only: KKROperator
    use BCPOperator_mod, only: BCPOperator, calc

    USE_ARRAYLOG_MOD
    USE_LOGGING_MOD

    type(IterativeSolver), intent(inout) :: solver
    type(KKROperator), intent(inout) :: op
    type(BCPOperator), intent(inout) :: preconditioner
    double precision, intent(in) :: kpoint(3)
    double complex, intent(in) :: tmatLL(:,:,:) !> tmatLL(lmsd,lmsd,naez_trc)
    double complex, intent(in) :: Ginp(:,:,:,:) !> Ginp(lmmaxd,lmmaxd,naclsd,N) where N=num_local_atoms or naez_trc if defined(SPLIT_REFERENCE_FOURIER_COM)
    double precision, intent(in) :: alat
    double precision, intent(in)  :: RR(:,0:)
    integer, intent(in) :: global_atom_id(:) ! becomes redundant with SPLIT_REFERENCE_FOURIER_COM
    integer, intent(in) :: communicator      ! becomes redundant with SPLIT_REFERENCE_FOURIER_COM
    type(InitialGuess), intent(inout) :: iguess_data
    integer, intent(in) :: ienergy, ispin
    ! LLY
    double complex, intent(in)         :: mssq(:,:,:)    !< inverted T-matrix
    double complex, intent(in)         :: dtde(:,:,:)     !< energy derivative of T-matrix
    double complex, intent(in)         :: dGinp(:,:,:,:)  !< dG_ref/dE dim(lmmaxd,lmmaxd,naclsd,N), compare N from above
    double complex, intent(out)        :: bztr2
    double precision, intent(in)       :: volcub (:)
    integer, intent(in)                :: ikpoint 
    integer, intent(in)                :: global_atom_idx_lly !< includes the global index of local atom so that atom-specific entries in global arrays can be accessed, e.g. dtde, tmatll
    integer, intent(in)                :: Lly             !< LLY=1/0, turns Lloyd's formula on/off
    integer, intent(in) :: solver_type

    double complex, allocatable :: dPdE_local(:,:), gllke_x(:,:), dgde(:,:), gllke_x_t(:,:), dgde_t(:,:), gllke_x2(:,:), dgde2(:,:) ! LLY
    double complex :: tracek, gtdPdE ! LLY

    integer :: naez, nacls, alm, lmmaxd, nRHSs, ist, matrix_index, lm1, lm2, lm3, il1
    integer :: i, n, nB, Bd
    double complex :: cfctorinv

    cfctorinv = (cone*8.d0*atan(1.d0))/alat
    
#define cluster op%cluster_info
    lmmaxd = op%lmmaxd
    naez = size(tmatLL, 3) ! number of atoms in the unified truncation zones
    nacls = cluster%naclsd


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
  
#ifndef SPLIT_REFERENCE_FOURIER_COM
    call referenceFourier_com(op%mat_A(:,:,:,0), op%bsr_A, kpoint, alat, &
             cluster%nacls_trc, cluster%atom_trc,  cluster%numn0_trc, cluster%indn0_trc, &
             RR, cluster%ezoa_trc, Ginp, global_atom_id, communicator)
#else
    call referenceFourier_part2(op%mat_A(:,:,:,0), op%bsr_A, kpoint, alat, &
             cluster%nacls_trc, cluster%atom_trc,  cluster%numn0_trc, cluster%indn0_trc, &
             RR, cluster%ezoa_trc, Ginp)
#endif

    TESTARRAYLOG(3, op%mat_A)

    ! ToDo: merge the referenceFourier_part2 with buildKKRCoeffMatrix

    if (Lly == 1) then ! LLY
      ! Allocate additional arrays for Lloyd's formula
      alm = lmmaxd*naez

      allocate(gllke_x(lmmaxd*naez,lmmaxd*nacls))
      allocate(dgde(lmmaxd*naez,lmmaxd*nacls))
      allocate(gllke_x_t(lmmaxd*nacls,lmmaxd*naez))
      allocate(dgde_t(lmmaxd*nacls,lmmaxd*naez))
      allocate(gllke_x2(lmmaxd*naez,lmmaxd))
      allocate(dgde2(lmmaxd*naez,lmmaxd))
      allocate(dPdE_local(lmmaxd*naez,lmmaxd))
    
#ifndef SPLIT_REFERENCE_FOURIER_COM
      call referenceFourier_com(op%mat_A(:,:,:,Lly), op%bsr_A, kpoint, alat, &
             cluster%nacls_trc, cluster%atom_trc,  cluster%numn0_trc, cluster%indn0_trc, &
             RR, cluster%ezoa_trc, dGinp, global_atom_id, communicator)
#else
      call referenceFourier_part2(op%mat_A(:,:,:,Lly), op%bsr_A, kpoint, alat, &
             cluster%nacls_trc, cluster%atom_trc,  cluster%numn0_trc, cluster%indn0_trc, &
             RR, cluster%ezoa_trc, dGinp)
#endif

      TESTARRAYLOG(3, op%mat_A(:,:,:,Lly))

      call convertBSRToFullMatrix(op%mat_A(:,:,:,0),   op%bsr_A, gllke_x)
      call convertBSRToFullMatrix(op%mat_A(:,:,:,Lly), op%bsr_A, dgde) 

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

      call zgemm('n','n',alm,lmmaxd,lmmaxd,cone, dgde2,alm, tmatll(1,1,global_atom_idx_lly),lmmaxd,zero, dPdE_local,alm)

      call zgemm('n','n',alm,lmmaxd,lmmaxd,cfctorinv, gllke_x2,alm, dtde(:,:,global_atom_idx_lly),lmmaxd,cone, dPdE_local,alm)
      !--------------------------------------------------------
 
    endif ! LLY
    
    ! TODO: merge the referenceFourier_part2 with buildKKRCoeffMatrix

    !----------------------------------------------------------------------------
    call buildKKRCoeffMatrix(op%mat_A(:,:,:,0), tmatLL, op%bsr_A)
    !----------------------------------------------------------------------------

    TESTARRAYLOG(3, op%mat_A(:,:,:,0))

    ! ==> now GLLh holds (1 - Delta_t * G_ref)

    ! Now solve the linear matrix equation A*X = b (b is also a matrix),
    ! where A = (1 - Delta_t*G_ref) (inverse of scattering path operator)
    ! and b = Delta_t

    !===================================================================
    ! 3) solve linear set of equations by iterative TFQMR scheme
    !    solve (1 - \Delta t * G_ref) X = \Delta t
    !    the solution X is the scattering path operator

    ! ToDo: use bsr_B instead of bsr_X
    call buildRightHandSide(op%mat_B, op%bsr_X, op%atom_indices, tmatLL=tmatLL) ! construct RHS with t-matrices
!   call buildRightHandSide(op%mat_B, op%bsr_X, op%atom_indices) ! construct RHS as unity

    if (iguess_data%prec == 1) then
      solver%initial_zero = .false.
      call load(iguess_data, op%mat_X, ik=ikpoint, is=ispin, ie=ienergy)
    else
      solver%initial_zero = .true.
    endif

    call calc(preconditioner, op%mat_A(:,:,:,0)) ! calculate preconditioner from sparse matrix data ! should be BROKEN due to variable block row format ! TODO: check

    selectcase(solver_type)
    case (4) ! direct solution with LAPACK, should only be used for small systems

#ifndef useBSR
      nB = size(op%mat_X, 3)
      n  = size(op%mat_A, 2)*nB
      Bd = size(op%mat_A, 1)
      nRHSs = size(op%mat_X, 2)
      ASSERT( nRHSs == size(op%mat_B, 2) )
      ASSERT( nB    == size(op%mat_B, 3) )
#else
      Bd = size(op%mat_A, 2)
      n =     Bd*op%bsr_A%nRows
      nRHSs = Bd*op%bsr_X%nCols
#endif

      if (any(shape(full_A) /= [n,n])) then
        deallocate(full_A, full_X, stat=ist) ! ignore status
        allocate(full_A(n,n), full_X(n,nRHSs), stat=ist)
        if (ist /= 0) die_here("failed to allocate dense matrix with"+(n*.5**26*n)+"GiByte!")
      endif
      call convertBSRToFullMatrix(op%mat_A(:,:,:,0), op%bsr_A, full_A)
      TESTARRAYLOG(3, full_A)

      ! convert op%mat_B to full_B
#ifndef useBSR
      do i = 1, nB
        full_X(Bd*(i - 1) + 1:Bd*i,:) = op%mat_B(:,:,i)
      enddo ! i
#else
      call convertBSRToFullMatrix(op%mat_B, op%bsr_X, full_X)
#endif

      ist = solveFull(full_A, full_X) ! on entry, full_X contains mat_B, compute the direct solution using LAPACK
      if (ist /= 0) die_here("failed to directly invert a matrix of dim"+n+"with"+nRHSs+"right hand sides!")
      
      ! convert back full_X to op%mat_X
#ifndef useBSR
      do i = 1, nB
        op%mat_X(:,:,i) = full_X(Bd*(i - 1) + 1:Bd*i,:) 
      enddo ! i
#else
      call convertFullMatrixToBSR(op%mat_X, op%bsr_X, full_X)
#endif

    case (0, 3) ! iterative solver
      if(solver_type == 0) warn(6, "solver_type ="+solver_type+"is deprecated, please use 3")

      call solve(solver)!, op%mat_X, op%mat_B) ! use iterative solver

#ifdef DEBUG_dump_matrix
        call dump(op%bsr_A, "matrix_descriptor.dat") ! SparseMatrixDescription
        call dump(op%mat_A,  "bin.matrix", formatted=.false.)
        call dump(op%mat_A,  "matrix_form.dat", formatted=.true.)
        call dump(op%mat_X, "bin.solution", formatted=.false.)
        call dump(op%mat_X, "solution_form.dat", formatted=.true.)
        call dump(op%mat_B, "bin.rhs", formatted=.false.)
        call dump(op%mat_B, "rhs_form.dat", formatted=.true.)
#endif
    case default
      warn(6, "No solver selected! Problem is not solved, solver_type ="+solver_type)
    endselect ! solver_type
    
    TESTARRAYLOG(3, op%mat_B)

    ! store the initial guess in previously selected slot (selected with 'iguess_set_k_ind')
    call store(iguess_data, op%mat_X, ik=ikpoint, is=ispin, ie=ienergy)

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
          do il1 = 1, naez
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

#undef cluster
  endsubroutine ! kloopbody


#ifndef SPLIT_REFERENCE_FOURIER_COM

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
                indn0, rr, ezoa, Ginp, global_atom_id, communicator)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    use ChunkIndex_mod, only: getRankAndLocalIndex
    use one_sided_commZ_mod, only: exposeBufferZ, copyChunksNoSyncZ, hideBufferZ
#ifdef NO_LOCKS_MPI
    use one_sided_commZ_mod, only: fenceZ
#endif
    include 'mpif.h'
    double complex, intent(out) :: GLLh(:,:,:)
    type(SparseMatrixDescription), intent(in) :: sparse
    double precision, intent(in) :: kpoint(3)
    double precision, intent(in) :: alat
    integer, intent(in) :: nacls(:) ! dim(naez_trc)
    integer(kind=2), intent(in) :: atom(:,:)
    integer, intent(in) :: numn0(:)
    integer(kind=2), intent(in) :: indn0(:,:)
    double precision, intent(in) :: rr(:,0:)
    integer(kind=2), intent(in) :: ezoa(:,:)
    double complex, intent(in) :: Ginp(:,:,:,:)
    integer, intent(in) :: global_atom_id(:) !> mapping trunc. index -> global atom index
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
                indn0, rr, ezoa, Ginp, global_atom_id, communicator)
      return
    endif

    naez = size(nacls)
    ASSERT(naez == size(global_atom_id))
    lmmaxd = size(Ginp, 1)
    ASSERT(lmmaxd == size(Ginp, 2))
    naclsd = size(Ginp, 3)
    
    allocate(eikrm(naclsd), eikrp(naclsd))

    ! Note: some MPI implementations might need the use of MPI_Alloc_mem
    allocate(Gref_buffer(lmmaxd,lmmaxd,naclsd))
!    allocate(dGref_buffer(lmmaxd,lmmaxd,naclsd))

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

        ! get Ginp(:,:,:)[global_atom_id(site_index)]

        atom_requested = global_atom_id(site_index)
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
        call dlke1(alat, nacls(site_index), rr, ezoa(:,site_index), kpoint, eikrm, eikrp)
        
        call dlke0_smat(GLLh, site_index, sparse, eikrm, eikrp, nacls(site_index), atom(:,site_index), numn0, indn0, Gref_buffer)

      endif ! site_index in bounds
      
    enddo ! site_index

    call hideBufferZ(win)

    deallocate(Gref_buffer, eikrm, eikrp, stat=ist)

  endsubroutine ! referenceFourier_com

#endif
  
  
#ifdef SPLIT_REFERENCE_FOURIER_COM
  
  subroutine referenceFourier_com_part1(Gref_buffer, naez, Ginp, global_atom_id, communicator)
    !! distributes the Ginp
    use ChunkIndex_mod, only: getRankAndLocalIndex
    use one_sided_commZ_mod, only: exposeBufferZ, copyChunksNoSyncZ, hideBufferZ
#ifdef NO_LOCKS_MPI
    use one_sided_commZ_mod, only: fenceZ
#endif
    include 'mpif.h'
    double complex, allocatable, intent(out) :: Gref_buffer(:,:,:,:)
    integer, intent(in) :: naez ! = size(nacls)
    double complex, intent(in) :: Ginp(:,:,:,:)
    integer, intent(in) :: global_atom_id(:) !> mapping trunc. index -> global atom index
    integer, intent(in) :: communicator

    ! locals
    integer :: site_index, naclsd, lmmaxd, ist
    integer :: num_local_atoms, atom_requested
    integer(kind=4) :: chunk_inds(2,1)
    integer :: win, nranks, ierr
    integer :: naez_max
    
    num_local_atoms = size(Ginp, 4)
    ASSERT(naez == size(global_atom_id))
    lmmaxd = size(Ginp, 1)
    ASSERT(lmmaxd == size(Ginp, 2))
    naclsd = size(Ginp, 3)

    ! Note: some MPI implementations might need the use of MPI_Alloc_mem
    if (any(shape(Gref_buffer) /= [lmmaxd,lmmaxd,naclsd,naez])) then
      deallocate(Gref_buffer, stat=ist)
      allocate(Gref_buffer(lmmaxd,lmmaxd,naclsd,naez), stat=ist)
      ! Note: some MPI implementations might need the use of MPI_Alloc_mem
      if (ist /= 0) stop 'failed to allocate Gref_buffer in referenceFourier_com_part1!'
    endif ! shape matches
   
    if (num_local_atoms == 1) then
      ! this version using point-to-point MPI communication is so far only implemented for 1 atom per process
      call referenceFourier_mpi_part1(Gref_buffer, naez, Ginp, global_atom_id, communicator)
      return
    endif
    
    call MPI_Comm_size(communicator, nranks, ierr)

#ifdef NO_LOCKS_MPI
    ! get maximum number of atoms of all truncation zones
    call MPI_Allreduce(naez, naez_max, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)
#else
    naez_max = naez
#endif

    ! share Ginp with all other processes in 'communicator'
    call exposeBufferZ(win, Ginp, lmmaxd*lmmaxd*naclsd*num_local_atoms, lmmaxd*lmmaxd*naclsd, communicator)

    ! loop up to naez_max to ensure that each rank does the same amount of fence calls
    do site_index = 1, naez_max

#ifdef NO_LOCKS_MPI
      call fenceZ(win)
      if (site_index <= naez) then
#endif

        atom_requested = global_atom_id(site_index)
        chunk_inds(:,1) = getRankAndLocalIndex(atom_requested, num_local_atoms*nranks, nranks)

#ifdef NO_LOCKS_MPI
        call copyChunksNoSyncZ(Gref_buffer(:,:,:,site_index), win, chunk_inds, lmmaxd*lmmaxd*naclsd)
      endif ! site_index <= naez
      
      call fenceZ(win) ! ensures that data has arrived in Gref_buffer
#else

#ifdef IDENTICAL_REF
      Gref_buffer(:,:,:,site_index) = Ginp(:,:,:,1) ! use this if all Grefs are the same
#else
      call MPI_Win_Lock(MPI_LOCK_SHARED, chunk_inds(1,1), 0, win, ierr)
      CHECKASSERT(ierr == 0)

      call copyChunksNoSyncZ(Gref_buffer(:,:,:,site_index), win, chunk_inds, lmmaxd*lmmaxd*naclsd)

      call MPI_Win_Unlock(chunk_inds(1,1), win, ierr)
      CHECKASSERT(ierr == 0)
#endif

#endif
    enddo ! site_index

    call hideBufferZ(win)

  endsubroutine ! referenceFourier_com_part1

  subroutine referenceFourier_part2(GLLh, sparse, kpoint, alat, nacls, atom, numn0, indn0, rr, ezoa, Gref_buffer)
    !! this operation will be performed for every k-point and does not include MPI communication
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    double complex, intent(out) :: GLLh(:,:,:)
    type(SparseMatrixDescription), intent(in) :: sparse
    double precision, intent(in) :: kpoint(3)
    double precision, intent(in) :: alat
    integer, intent(in) :: nacls(:)
    integer(kind=2), intent(in) :: atom(:,:)
    integer, intent(in) :: numn0(:)
    integer(kind=2), intent(in) :: indn0(:,:)
    double precision, intent(in) :: rr(:,0:)
    integer(kind=2), intent(in) :: ezoa(:,:)
    double complex, intent(in) :: Gref_buffer(:,:,:,:)

    ! locals
    integer :: site_index, naez, lmmaxd, naclsd, ist
    double complex, allocatable :: eikrm(:), eikrp(:)

    naez = size(nacls)
    lmmaxd = size(Gref_buffer, 1)
    ASSERT ( size(Gref_buffer, 2) == lmmaxd )

    naclsd = maxval(nacls)
    ASSERT ( size(Gref_buffer, 3) >= naclsd )

    allocate(eikrm(naclsd), eikrp(naclsd))

    GLLh = zero ! init

    ! loop up to naez_max to ensure that each rank does the same amount of fence calls
    do site_index = 1, naez

      call dlke1(alat, nacls(site_index), rr, ezoa(:,site_index), kpoint, eikrm, eikrp)

      call dlke0_smat(GLLh, site_index, sparse, eikrm, eikrp, nacls(site_index), atom(:,site_index), numn0, indn0, Gref_buffer(:,:,:,site_index))

    enddo ! site_index

    deallocate(eikrm, eikrp, stat=ist)

  endsubroutine ! referenceFourier_part2
  
#endif  
  
  
  
  
  
#ifndef SPLIT_REFERENCE_FOURIER_COM
  
  subroutine referenceFourier_mpi(GLLh, sparse, kpoint, alat, nacls, atom, numn0, &
                indn0, rr, ezoa, Ginp, global_atom_id, comm)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    use ChunkIndex_mod, only: getRankAndLocalIndex
    double complex, intent(out) :: GLLh(:,:,:)
    type(SparseMatrixDescription), intent(in) :: sparse
    double precision, intent(in) :: kpoint(3)
    double precision, intent(in) :: alat
    integer, intent(in) :: nacls(:)
    integer(kind=2), intent(in) :: atom(:,:)
    integer, intent(in) :: numn0(:)
    integer(kind=2), intent(in) :: indn0(:,:)
    double precision, intent(in) :: rr(:,0:)
    integer(kind=2), intent(in) :: ezoa(:,:)
    double complex, intent(in) :: Ginp(:,:,:,:)
    integer, intent(in) :: global_atom_id(:) !> mapping trunc. index -> atom index
    integer, intent(in) :: comm

    ! locals
    integer(kind=4) :: chunk_inds(2,1)
    integer :: site_index, naez, naclsd, lmmaxd, ist
    integer :: num_local_atoms, atom_requested
    double complex, allocatable :: Gref_buffer(:,:,:,:), eikrm(:), eikrp(:) ! dim: naclsd
    integer :: rank, tag, myrank, nranks, ierr, ncount
    integer, allocatable :: reqs(:,:), stats(:,:,:)
    integer, parameter :: TAGMOD = 2**15
    include 'mpif.h'

    naez = size(nacls)
    ASSERT(naez == size(global_atom_id))
    lmmaxd = size(Ginp, 1)
    ASSERT(lmmaxd == size(Ginp, 2))
    naclsd = size(Ginp, 3)
    num_local_atoms = size(Ginp, 4)
    
    ASSERT(num_local_atoms == 1) ! only 1 atom per MPI process
    
    
    call MPI_Comm_size(comm, nranks, ierr)
    call MPI_Comm_rank(comm, myrank, ierr)

#ifndef IDENTICAL_REF
    ! Note: some MPI implementations might need the use of MPI_Alloc_mem
    allocate(Gref_buffer(lmmaxd,lmmaxd,naclsd,naez))
    
    allocate(reqs(2,naez), stats(MPI_STATUS_SIZE,2,naez))
    reqs(:,:) = MPI_REQUEST_NULL

    ncount = lmmaxd*lmmaxd*naclsd

    ! ===============================================================================
    ! part 1: exchange the Ginp arrays between the MPI processes 
    ! ===============================================================================
    
    ! loop over the sites sending the information
    do site_index = 1, naez
      atom_requested = global_atom_id(site_index) ! get the global atom id

!     rank = (atom_requested - 1)/num_local_atoms ! block distribution of atoms to ranks
      chunk_inds(:,1) = getRankAndLocalIndex(atom_requested, num_local_atoms*nranks, nranks)
      rank = chunk_inds(1,1)
      assert( chunk_inds(2,1) == 1 ) ! since there may only be one local atom, its local index must be one

      if (rank /= myrank) then

        tag = modulo(myrank, TAGMOD)
        call MPI_Isend(Ginp(:,:,:,1),                 ncount, MPI_DOUBLE_COMPLEX, rank, tag, comm, reqs(1,site_index), ierr)

        tag = modulo(atom_requested - 1, TAGMOD)
        call MPI_Irecv(Gref_buffer(:,:,:,site_index), ncount, MPI_DOUBLE_COMPLEX, rank, tag, comm, reqs(2,site_index), ierr)

      else
      
        reqs(:,site_index) = MPI_REQUEST_NULL
        Gref_buffer(:,:,:,site_index) = Ginp(:,:,:,1) ! copy locally
        
      endif ! distant rank

    enddo ! site_index

    call MPI_Waitall(2*naez, reqs, stats, ierr) ! wait until all sends and all receives have finished
#else
    allocate(Gref_buffer(lmmaxd,lmmaxd,naclsd,1))
    Gref_buffer(:,:,:,1) = Ginp(:,:,:,1)
#define SITE_INDEX 1
#endif

    ! ===============================================================================
    ! part 2: perform the fourier transformation and prepare the entries of GLLh
    ! ===============================================================================
    GLLh = zero ! init
    allocate(eikrm(naclsd), eikrp(naclsd))

    do site_index = 1, naez
    
      call dlke1(alat, nacls(site_index), rr, ezoa(:,site_index), kpoint, eikrm, eikrp)

      call dlke0_smat(GLLh, site_index, sparse, eikrm, eikrp, nacls(site_index), atom(:,site_index), numn0, indn0, Gref_buffer(:,:,:,SITE_INDEX))

    enddo ! site_index
    
    deallocate(Gref_buffer, eikrm, eikrp, stat=ist)

  endsubroutine ! referenceFourier_mpi

#endif

#ifdef SPLIT_REFERENCE_FOURIER_COM

  subroutine referenceFourier_mpi_part1(Gref_buffer, naez, Ginp, global_atom_id, comm)
    use ChunkIndex_mod, only: getRankAndLocalIndex
    double complex, intent(out) :: Gref_buffer(:,:,:,:) ! (lmmaxd,lmmaxd,naclsd,naez)
    integer, intent(in) :: naez
    double complex, intent(in) :: Ginp(:,:,:,:)
    integer, intent(in) :: global_atom_id(:) !> mapping trunc. index -> atom index
    integer, intent(in) :: comm

    ! locals
    integer(kind=4) :: chunk_inds(2,1)
    integer :: site_index, naclsd, lmmaxd, ist
    integer :: num_local_atoms, atom_requested
    integer :: rank, tag, myrank, nranks, ierr, ncount
    integer, allocatable :: reqs(:,:), stats(:,:,:)
    integer, parameter :: TAGMOD = 2**15
    include 'mpif.h'

    ASSERT(naez == size(global_atom_id))
    lmmaxd = size(Ginp, 1)
    ASSERT(lmmaxd == size(Ginp, 2))
    naclsd = size(Ginp, 3)
    num_local_atoms = size(Ginp, 4)
    
    ASSERT(num_local_atoms == 1) ! only 1 atom per MPI process
    
    call MPI_Comm_size(comm, nranks, ierr)
    call MPI_Comm_rank(comm, myrank, ierr)

#ifdef IDENTICAL_REF
    Gref_buffer(:,:,:,1) = Ginp(:,:,:,1)
#else
    allocate(reqs(2,naez), stats(MPI_STATUS_SIZE,2,naez), stat=ist)
    reqs(:,:) = MPI_REQUEST_NULL

    ncount = lmmaxd*lmmaxd*naclsd

    ! ===============================================================================
    ! part 1: exchange the Ginp arrays between the MPI processes 
    ! ===============================================================================
    
    ! loop over the sites sending the information
    do site_index = 1, naez
      atom_requested = global_atom_id(site_index) ! get the global atom id

!     rank = (atom_requested - 1)/num_local_atoms ! block distribution of atoms to ranks
      chunk_inds(:,1) = getRankAndLocalIndex(atom_requested, num_local_atoms*nranks, nranks)
      rank = chunk_inds(1,1)
      assert( chunk_inds(2,1) == 1 ) ! since there may only be one local atom, its local index must be one

      if (rank /= myrank) then

        tag = modulo(myrank, TAGMOD)
        call MPI_Isend(Ginp(:,:,:,1),                 ncount, MPI_DOUBLE_COMPLEX, rank, tag, comm, reqs(1,site_index), ierr)

        tag = modulo(atom_requested - 1, TAGMOD)
        call MPI_Irecv(Gref_buffer(:,:,:,site_index), ncount, MPI_DOUBLE_COMPLEX, rank, tag, comm, reqs(2,site_index), ierr)

      else

        reqs(:,site_index) = MPI_REQUEST_NULL
        Gref_buffer(:,:,:,site_index) = Ginp(:,:,:,1) ! copy locally

      endif ! distant rank

    enddo ! site_index

    call MPI_Waitall(2*naez, reqs, stats, ierr) ! wait until all sends and all receives have finished
#endif
  endsubroutine ! referenceFourier_mpi_part1

#endif


  !     >> Input parameters
  !>    @param     alat    lattice constant a
  !>    @param     nacls   number of atoms in cluster
  !>    @param     RR      periodic image vectors
  !>    @param     ezoa
  !>    @param     kpoint
  !>    @param     naclsd. maximal number of atoms in cluster

  !     << Output parameters
  !>    @param     eikrm   Fourier exponential factor with minus sign
  !>    @param     eikrp   Fourier exponential factor with plus sign
  subroutine dlke1(alat, nacls, RR, ezoa, kpoint, eikrm, eikrp)
  use Constants_mod, only: pi
    ! ----------------------------------------------------------------------
    !     Fourier transformation of the cluster Greens function
    !     Prepares the calculation (calculates Fourier factors) for dlke0
    ! ----------------------------------------------------------------------
    double precision, intent(in) :: alat
    integer, intent(in) :: nacls !< number of vectors in the cluster
    integer(kind=2), intent(in) :: ezoa(1:) !< index list of periodic image vector
    double precision, intent(in) :: kpoint(1:3) !< k-point (vector in the Brillouin zone)
    double precision, intent(in) :: RR(1:,0:) !< dim(1:3,0:) periodic image vectors
    double complex, intent(out) :: eikrp(:), eikrm(:) !< dim(nacls)

    double complex, parameter :: ci=(0.d0, 1.d0)
    double precision :: convpuh, tpi
    double complex :: tt, exparg
    integer :: iacls

    tpi = 2.d0*pi
    convpuh = alat/tpi * 0.5d0

    do iacls = 1, nacls
       
      ! Here we do       --                  nn'
      !                  \                   ii'          ii'
      !                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
      !                  --          n'  n   LL'          LL'
      !                  n'
      ! Be careful about the minus sign included here. RR is not
      ! symmetric around each atom. The minus comes from the fact that
      ! the repulsive potential GF is calculated for 0n and not n0!

      if (ezoa(iacls) == 0) then
        ! the periodic image vector is (0,0,0)
        ASSERT( all(RR(1:3,ezoa(iacls)) == 0.d0) )
        eikrp(iacls) = dcmplx(convpuh, 0.d0)
        eikrm(iacls) = dcmplx(convpuh, 0.d0)
      else
  
        tt = -ci*tpi*dot_product(kpoint(1:3), RR(1:3,ezoa(iacls))) ! purely imaginary number

        ! convert to p.u. and multiply with 1/2 (done above)
        exparg = exp(tt)
        eikrp(iacls) =       exparg  * convpuh
        eikrm(iacls) = conjg(exparg) * convpuh ! we can re-use exparg here instead of exp(-tt) since tt is purely imaginary
      endif
    enddo ! iacls

  endsubroutine ! dlke1
  
  
  subroutine dlke0_smat(smat, iat, sparse, eikrm, eikrp, nacls, atom, numn0, indn0, Ginp)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    ! assume a variable block row sparse matrix description
    double complex, intent(inout) :: smat(:,:,:)
    integer, intent(in) :: iat !> local_atom_idx of the source atom
    type(SparseMatrixDescription), intent(in) :: sparse
    double complex, intent(in) :: eikrm(nacls), eikrp(nacls) ! todo: many of these phase factors are real
    integer, intent(in) :: nacls !< number of atoms in the cluster around site iat == nacls(iat)
    integer(kind=2), intent(in) :: atom(:) !< dim(nacls) == atom(:,iat)
    integer, intent(in) :: numn0(:) !< dim(naez)
    integer(kind=2), intent(in) :: indn0(:,:) !< dims(nacls+,naez)
    double complex, intent(in) :: Ginp(:,:,:) !< dims(lmmaxd,lmmaxd,nacls)

    integer :: jat, iacls, ni, ind, ist
    double complex, allocatable :: GinT(:,:) ! (size(Ginp, 2),size(Ginp, 1)) !< dims(lmmaxd,lmmaxd)
    
    allocate(GinT(size(Ginp, 2),size(Ginp, 1)))

    do iacls = 1, nacls ! loop over all atoms in the reference cluster around iat
      jat = atom(iacls) ! local_atom_idx of the target atom
      ASSERT( jat > 0 ) 

      do ni = 1, numn0(iat) ! loop over the set of inequivalent atoms in the reference cluster around iat
        ind = indn0(ni,iat)
        if (ind == jat) then ! see which one of the inequivalent atoms is hit (should only be true once)

          GinT = transpose(Ginp(:,:,iacls))
          call modify_smat(sparse, iat, ni, eikrm(iacls), GinT, smat)

        endif ! jat == ind
      enddo ! ni

      do ni = 1, numn0(jat) ! loop over the set of inequivalent atoms in the reference cluster around jat
        ind = indn0(ni,jat)
        if (ind == iat) then ! see which one of the inequivalent atoms is hit (should only be true once)
 
          call modify_smat(sparse, jat, ni, eikrp(iacls), Ginp(:,:,iacls), smat)

        endif ! iat == ind
      enddo ! ni

    enddo ! iacls
    
    deallocate(GinT, stat=ist) ! ignore status
  endsubroutine ! dlke0_smat

  
  subroutine modify_smat(sparse, ind, ni, eikr, Gin, smat)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    ! assume a variable block row sparse matrix description
    type(SparseMatrixDescription), intent(in) :: sparse
    integer, intent(in) :: ind !> source site_index
    integer, intent(in) :: ni  !> ???
    double complex, intent(in) :: eikr ! phase factor
    double complex, intent(in) :: Gin(:,:) !< dims(lmmaxd,lmmaxd)
    double complex, intent(inout) :: smat(:,:,:) !< dim(lmmaxd,lmmaxd,nnzb)

    integer :: Aind

    Aind = sparse%RowStart(ind) + ni - 1
    smat(:,:,Aind) = smat(:,:,Aind) + eikr * Gin(:,:)

  endsubroutine ! modify_smat
  
endmodule ! kkrmat_mod
