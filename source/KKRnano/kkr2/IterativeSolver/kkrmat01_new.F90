!> Multiple Scattering problem: Solving the Dyson equation for all k-points,
!> using *real-space* G_ref and the t-matrices
!>
!> Output: Brillouin-zone integrated diagonal elements of structural Green's
!> function

#include "../DebugHelpers/logging_macros.h"
#include "../DebugHelpers/test_array_log.h"
#include "../DebugHelpers/test_macros.h"


!#define SPLIT_REFERENCE_FOURIER_COM

module kkrmat_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  use arraytest2_mod, only: !import no name here, just mention it for the module dependency
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: kkrmat01

  double complex, allocatable :: full(:,:)
  double complex, parameter :: zero=(0.d0, 0.d0), cone=(1.d0, 0.d0)

  contains

  !------------------------------------------------------------------------------
  !> Solves multiple scattering problem for every k-point.
  !>
  !> Returns diagonal k-integrated part of Green's function in GS.
  subroutine kkrmat01(solver, op, preconditioner, kpoints, nkpoints, kpointweight, GS, tmatLL, alat, nsymat, RR, &
                          Ginp, lmmaxd, global_atom_id, communicator, iguess_data, &
                          mssq, dginp, dtde, tr_alph, lly_grdt, volcub, volbz, global_atom_idx_lly, lly) !LLY
    !   performs k-space integration,
    !   determines scattering path operator (g(k,e)-t**-1)**-1 and
    !   Greens function of the real system -> GS(*,*,*),
    USE_LOGGING_MOD
    USE_ARRAYLOG_MOD
    use InitialGuess_mod, only: InitialGuess, iguess_set_k_ind
    use jij_calc_mod, only: global_jij_data, kkrjij
    use SolverStats_mod, only: SolverStats, reset
    use IterativeSolver_mod, only: IterativeSolver
    use BCPOperator_mod, only: BCPOperator
    use KKROperator_mod, only: KKROperator
    use mpi

    type(IterativeSolver), intent(inout) :: solver
    type(KKROperator), intent(inout) :: op
    type(BCPOperator), intent(inout) :: preconditioner
    
    integer, intent(in) :: nkpoints !< number of k-points
    double precision, intent(in) :: kpoints(:,:) !< list of k-points dim(3,nkpoints)
    double precision, intent(in) :: kpointweight(:) !< k-point weights dim(nkpoints)

    double complex, intent(out) :: GS(:,:,:) ! (lmmaxd,lmmaxd,num_local_atoms)
    double complex, intent(in) :: tmatLL(:,:,:) ! (lmmaxd,lmmaxd,naez)
    double precision, intent(in) :: alat
    integer, intent(in) :: nsymat ! needed only for Jij-calculation
    double precision, intent(in) :: RR(:,0:)
    double complex, intent(inout) :: Ginp(:,:,:,:) ! (lmmaxd,lmmaxd,naclsd,nclsd)
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: global_atom_id(:)
    integer, intent(in) :: communicator
    type(InitialGuess), intent(inout) :: iguess_data

    !LLY
    double complex, intent(in)   :: mssq (:,:,:)    !< inverted T-matrix
    double complex, intent(inout) :: dginp(:,:,:,:)  !< dG_ref/dE,  dim: lmmaxd, lmmaxd, naclsd, nclsd
    double complex, intent(in)   :: dtde(:,:,:)     !< dT/dE
    double complex, intent(in)   :: tr_alph(:) 
    double complex, intent(out)  :: lly_grdt
    double precision, intent(in) :: volcub (:)
    double precision, intent(in) :: volbz
    integer, intent(in)          :: global_atom_idx_lly
    integer, intent(in)          :: lly

    ! locals
    double complex :: G_diag(LMMAXD,LMMAXD)
    double complex :: bztr2, trace ! LLY
    integer :: site_lm_size, num_local_atoms, naclsd, naez, ikpoint, ila, ierr
    type(SolverStats) :: stats
#ifdef SPLIT_REFERENCE_FOURIER_COM
    double complex, allocatable :: Gref_buffer(:,:,:,:) ! split_reference_fourier_com uses more memory but calls the communication routine only 1x per energy point
!    double complex, allocatable :: DGref_buffer(:,:,:,:) ! LLY
#endif

#define cluster op%cluster_info

    ! array dimensions
    naez = cluster%naez_trc
    naclsd = cluster%naclsd

    site_lm_size = naez*LMMAXD

    num_local_atoms = size(op%atom_indices)

    ! WARNING: Symmetry assumptions might have been used that are
    ! not valid in cases of non-local potential (e.g. for Spin-Orbit coupling)
    ! ---> use sit
    !      G(n,n',L,L')(-k) = G(n',n,L',L)(k)

    GS = zero ! init zero
    bztr2 = zero ! init zero

    TESTARRAYLOG(3, Ginp)

    if (global_jij_data%do_jij_calculation) global_jij_data%GSXIJ = zero

    call reset(solver%stats)

#ifdef SPLIT_REFERENCE_FOURIER_COM
    ! get the required reference Green functions from the other MPI processes
    call referenceFourier_com_part1(Gref_buffer, naez, Ginp, global_atom_id, communicator)
! #define Ginp Gref_buffer
!     if (lly == 1) then ! LLY
!       call referenceFourier_com_part1(DGref_buffer, naez, DGinp, global_atom_id, communicator)
!     endif ! LLY
! #define dginp DGref_buffer
#endif
   
    !==============================================================================
    do ikpoint = 1, nkpoints ! K-POINT-LOOP
    !==============================================================================

      WRITELOG(4, *) "k-point ", ikpoint

      ! select right slot for storing initial guess
      call iguess_set_k_ind(iguess_data, ikpoint)

      ! Get the scattering path operator for k-point kpoints(:,ikpoint)
      ! output: op%mat_X
      call kloopbody(solver, op, preconditioner, kpoints(1:3,ikpoint), tmatLL, Ginp, &
                     alat, RR, global_atom_id, communicator, iguess_data, &
                     mssq, dtde, dginp, bztr2, volcub, ikpoint, &
                     global_atom_idx_lly ,lly) !LLY

      do ila = 1, num_local_atoms
        call getGreenDiag(G_diag, op%mat_X, op%atom_indices(ila), op%sparse%kvstr, ila) ! extract solution

        ! ----------- Integrate Scattering Path operator over k-points --> GS -----
        ! Note: here k-integration only in irreducible wedge
        GS(:,:,ila) = GS(:,:,ila) + kpointweight(ikpoint)*G_diag(:,:) 
        ! -------------------------------------------------------------------------
      enddo ! ila
      
      ! TODO: use mat_X to calculate Jij
      if (global_jij_data%do_jij_calculation) then
        ! communicate off-diagonal elements and multiply with exp-factor
        call KKRJIJ(kpoints(1:3,ikpoint), kpointweight(ikpoint), nsymat, naez, op%atom_indices(1), &
                    global_jij_data%NXIJ, global_jij_data%IXCP,global_jij_data%ZKRXIJ, &
                    op%mat_X, global_jij_data%GSXIJ, communicator, lmmaxd, global_jij_data%nxijd)
      endif ! jij

    !==============================================================================
    enddo ! ikpoint = 1, nkpoints
    !==============================================================================

    !--------------------- LLY ----------------------------------------------------
    if (lly == 1) then   
       bztr2 = bztr2*nsymat/volbz + tr_alph(1)
       trace = zero
       CALL MPI_ALLREDUCE(bztr2,trace,1, &
                         MPI_DOUBLE_COMPLEX,MPI_SUM, &
                         MPI_COMM_WORLD,ierr)
       lly_grdt = trace
    endif    
    !------------------------------------------------------------------------------

#ifdef SPLIT_REFERENCE_FOURIER_COM
#undef Ginp
#undef dginp
    deallocate(Gref_buffer, stat=ila) ! ignore status
!    deallocate(DGref_buffer, stat=ila) ! LLY, ignore status
#endif    
    
    do ila = 1, num_local_atoms
      TESTARRAYLOG(3, GS(:,:,ila))
    enddo ! ila
    
    WRITELOG(3, *) "Max. TFQMR residual for this E-point: ", solver%stats%max_residual
    WRITELOG(3, *) "Max. num iterations for this E-point: ", solver%stats%max_iterations
    WRITELOG(3, *) "Sum of iterations for this E-point:   ", solver%stats%sum_iterations

#undef cluster
#undef ms
  endsubroutine ! kkrmat01



  !------------------------------------------------------------------------------
  !> Copy the diagonal elements G_{LL'}^{nn'} of the Green's-function,
  !> dependent on (k,E) into matrix G_diag
  subroutine getGreenDiag(G_diag, mat_X, atom_index, kvstr, local_atom_index)
    double complex, intent(out) :: G_diag(:,:) ! dim(lmmaxd,lmmaxd)
    double complex, intent(in) :: mat_X(:,:)
    integer, intent(in) :: atom_index 
    integer, intent(in) :: kvstr(:)
    integer, intent(in) :: local_atom_index 

    integer :: lmmax1, start, ila !< local atom index

    !                                      nn
    !         Copy the diagonal elements G_LL' of the Green's-function,
    !         dependent on (k,E) into matrix G_diag
    !         (n = n' = atom_index)

    ila = local_atom_index
    G_diag = zero

    start  = kvstr(atom_index) - 1
    lmmax1 = kvstr(atom_index+1) - kvstr(atom_index)

    ASSERT(lmmax1 == size(G_diag, 1))
    ASSERT(lmmax1 == size(G_diag, 2))

    G_diag(1:lmmax1,1:lmmax1) = mat_X(start+ 1:lmmax1 +start,(ila-1)*lmmax1+1:ila*lmmax1)
    
  endsubroutine ! getGreenDiag

  !------------------------------------------------------------------------------
  !> Calculate scattering path operator for 'kpoint'.
  !>
  !> Input are the \Delta T and the realspace G_ref (Ginp).
  !> Solution is stored in op%mat_X.
  !> Scattering path operator is calculated for atoms given in
  !> op%atom_indices(:)
  subroutine kloopbody(solver, op, preconditioner, kpoint, tmatLL, Ginp, alat, RR, global_atom_id, communicator, iguess_data, &
                       mssq, dtde, dginp, bztr2, volcub, ikpoint, &
                       global_atom_idx_lly, lly) !LLY

    use fillKKRMatrix_mod, only: buildKKRCoeffMatrix, buildRightHandSide, solveFull, convertToFullMatrix
    use fillKKRMatrix_mod, only: dump
    use IterativeSolver_mod, only: IterativeSolver, solve
    use SparseMatrixDescription_mod, only: dump
    use InitialGuess_mod, only: InitialGuess, load, store
    use TEST_lcutoff_mod, only: cutoffmode, DEBUG_dump_matrix
    use KKROperator_mod, only: KKROperator
    use BCPOperator_mod, only: BCPOperator, calc

    USE_ARRAYLOG_MOD
    USE_LOGGING_MOD

    type(IterativeSolver), intent(inout) :: solver
    type(KKROperator), intent(inout) :: op
    type(BCPOperator), intent(inout) :: preconditioner
    double precision, intent(in) :: kpoint(3)
    double complex, intent(in) :: tmatLL(:,:,:)
    double complex, intent(in) :: Ginp(:,:,:,:) !> Ginp(lmmaxd,lmmaxd,naclsd,nclsd) independent of the kpoint, ToDo: try not to communicate it again for every kpoint
    double precision, intent(in) :: alat
    double precision, intent(in)  :: RR(:,0:)
    integer, intent(in) :: global_atom_id(:) ! becomes redundant with SPLIT_REFERENCE_FOURIER_COM
    integer, intent(in) :: communicator      ! becomes redundant with SPLIT_REFERENCE_FOURIER_COM
    type(InitialGuess), intent(inout) :: iguess_data
    ! LLY
    double complex, intent(in)         :: mssq(:,:,:)    !< inverted T-matrix
    double complex, intent(in)         :: dtde(:,:,:)     !< energy derivative of T-matrix
    double complex, intent(in)         :: dginp(:,:,:,:)  !< dG_ref/dE dim: lmmaxd, lmmaxd, naclsd, nclsd 
    double complex, intent(out)        :: bztr2
    double precision, intent(in)       :: volcub (:)
    integer, intent(in)                :: ikpoint 
    integer, intent(in)                :: global_atom_idx_lly !< includes the global index of local atom so that atom-specific entries in global arrays can be accessed, e.g. dtde, tmatll
    integer, intent(in)                :: lly             !< LLY=1/0, turns Lloyd's formula on/off

    ! Local LLY
    double complex, allocatable :: dpde_local(:,:)
    double complex, allocatable :: gllke_x(:,:)
    double complex, allocatable :: dgde(:,:)
    double complex, allocatable :: gllke_x_t(:,:)
    double complex, allocatable :: dgde_t(:,:)
    double complex, allocatable :: gllke_x2(:,:)
    double complex, allocatable :: dgde2(:,:)
    double complex :: tracek
    double complex :: gtdpde
    
    integer :: naez, nacls, alm, lmmaxd, n, ist, matrix_index, lm1, lm2, il1
    logical :: initial_zero
    double complex :: cfctorinv

    cfctorinv = (cone*8.d0*atan(1.d0))/alat
    
#define cluster op%cluster_info
    lmmaxd = op%lmmaxd
    naez = op%naez
    nacls = cluster%naclsd
    alm = naez*lmmaxd

    ! Allocate additional arrays for Lloyd's formula    
    if (.not. allocated(gllke_x)) then
      allocate(gllke_x(naez*lmmaxd, nacls*lmmaxd))
    end if
    if (.not. allocated(dgde)) then
      allocate(dgde(naez*lmmaxd, nacls*lmmaxd))
    end if
    if (.not. allocated(gllke_x_t)) then
      allocate(gllke_x_t(nacls*lmmaxd, naez*lmmaxd))
    end if
    if (.not. allocated(dgde_t)) then
      allocate(dgde_t(nacls*lmmaxd, naez*lmmaxd))
    end if
    if (.not. allocated(gllke_x2)) then
      allocate(gllke_x2(naez*lmmaxd,lmmaxd))
    end if
    if (.not. allocated(dgde2)) then
      allocate(dgde2(naez*lmmaxd,lmmaxd))
    end if
    if (.not. allocated(dpde_local)) then
      allocate(dpde_local(naez*lmmaxd,lmmaxd))
    end if
    
    

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
    call referenceFourier_com(op%GLLh, op%sparse, kpoint, alat, &
             cluster%nacls_trc, cluster%atom_trc,  cluster%numn0_trc, cluster%indn0_trc, &
             RR, cluster%ezoa_trc, Ginp, global_atom_id, communicator)
#else
    call referenceFourier_part2(op%GLLh, op%sparse, kpoint, alat, &
             cluster%nacls_trc, cluster%atom_trc,  cluster%numn0_trc, cluster%indn0_trc, &
             RR, cluster%ezoa_trc, Ginp)
#endif

    TESTARRAYLOG(3, op%GLLh)

   ! TODO: merge the referenceFourier_part2 with buildKKRCoeffMatrix

    if (lly == 1) then ! LLY
#ifndef SPLIT_REFERENCE_FOURIER_COM
      call referenceFourier_com(op%DGLLh, op%sparse, kpoint, alat, &
             cluster%nacls_trc, cluster%atom_trc,  cluster%numn0_trc, cluster%indn0_trc, &
             RR, cluster%ezoa_trc, DGinp, global_atom_id, communicator)
#else
      call referenceFourier_part2(op%DGLLh, op%sparse, kpoint, alat, &
             cluster%nacls_trc, cluster%atom_trc,  cluster%numn0_trc, cluster%indn0_trc, &
             RR, cluster%ezoa_trc, DGinp)
#endif

      TESTARRAYLOG(3, op%DGLLh)

      call convertToFullMatrix(op%GLLH, op%sparse%ia, op%sparse%ja, op%sparse%ka, &
                           op%sparse%kvstr, op%sparse%kvstr, GLLKE_X)
      call convertToFullMatrix(op%DGLLH, op%sparse%ia, op%sparse%ja, op%sparse%ka, &
                           op%sparse%kvstr, op%sparse%kvstr, DGDE) 

      !--------------------------------------------------------
      ! dP(E,k)   dG(E,k)                   dT(E)
      ! ------- = ------- * T(E) + G(E,k) * -----
      !   dE        dE                       dE
  
      matrix_index = (global_atom_idx_lly-1)*lmmaxd+1

      gllke_x_t = transpose(gllke_x)
      dgde_t = transpose(dgde)

      gllke_x2 = gllke_x_t(:, matrix_index:matrix_index+lmmaxd)
      dgde2 = dgde_t(:, matrix_index:matrix_index+lmmaxd)

      call cinit(naez*lmmaxd*lmmaxd,dpde_local)

      call zgemm('n','n',alm,lmmaxd,lmmaxd,cone,&
                  dgde2,alm,&
                  tmatll(1,1,global_atom_idx_lly),lmmaxd,zero,&
                  dpde_local,alm)

      call zgemm('n','n',alm,lmmaxd,lmmaxd,cfctorinv,&
                  gllke_x2,alm,&
                  dtde(:,:,global_atom_idx_lly),lmmaxd,cone,dpde_local,alm)
      !--------------------------------------------------------
 
    endif ! LLY
    
    ! TODO: merge the referenceFourier_part2 with buildKKRCoeffMatrix

    !----------------------------------------------------------------------------
    call buildKKRCoeffMatrix(op%GLLh, tmatLL, op%lmmaxd, naez, op%sparse)
    !----------------------------------------------------------------------------

    TESTARRAYLOG(3, op%GLLh)

    ! ==> now GLLh holds (1 - Delta_t * G_ref)

    ! Now solve the linear matrix equation A*X = b (b is also a matrix),
    ! where A = (1 - Delta_t*G_ref) (inverse of scattering path operator)
    ! and b = Delta_t

    !===================================================================
    ! 3) solve linear set of equations by iterative TFQMR scheme
    !    solve (1 - \Delta t * G_ref) X = \Delta t
    !    the solution X is the scattering path operator

    call buildRightHandSide(op%mat_B, lmmaxd, op%atom_indices, op%sparse%kvstr, tmatLL=tmatLL) ! construct RHS with t-matrices
    ! call buildRightHandSide(op%mat_B, lmmaxd, op%atom_indices, op%sparse%kvstr) ! construct RHS as unity

    if (iguess_data%iguess == 1) then
      solver%initial_zero = .false.
      call load(iguess_data, op%mat_X)
    else
      solver%initial_zero = .true.
    endif

    call calc(preconditioner, op%GLLh) ! calculate preconditioner from sparse matrix data ! should be BROKEN due to variable block row format ! TODO: check

    if (cutoffmode == 3 .or. cutoffmode == 0) then

      call solve(solver, op%mat_X, op%mat_B) ! use iterative solver

      if (DEBUG_dump_matrix) then
        call dump(op%sparse, "matrix_descriptor.dat") ! SparseMatrixDescription
        call dump(op%GLLh,  "bin.matrix", formatted=.false.)
        call dump(op%GLLh,  "matrix_form.dat", formatted=.true.)
        call dump(op%mat_X, "bin.solution", formatted=.false.)
        call dump(op%mat_X, "solution_form.dat", formatted=.true.)
        call dump(op%mat_B, "bin.rhs", formatted=.false.)
        call dump(op%mat_B, "rhs_form.dat", formatted=.true.)
      endif ! DEBUG_dump_matrix
      
    endif ! cutoffmode in {0,3}

    TESTARRAYLOG(3, op%mat_B)

    ! ALTERNATIVE: direct solution with LAPACK
    if (cutoffmode == 4) then
      n = size(op%mat_B, 1)
      if (any(shape(full) /= [n,n])) then
        deallocate(full, stat=ist) ! ignore status
        allocate(full(n,n), stat=ist)
        if (ist /= 0) die_here("failed to allocate dense matrix with"+(n*.5**26*n)+"GiByte!")
      endif
      call convertToFullMatrix(op%GLLh, op%sparse%ia, op%sparse%ja, op%sparse%ka, op%sparse%kvstr, op%sparse%kvstr, full)
      TESTARRAYLOG(3, full)
      call solveFull(full, op%mat_B, op%mat_X)
    endif ! cutoffmode == 4

    ! store the initial guess in previously selected slot (selected with 'iguess_set_k_ind')
    call store(iguess_data, op%mat_X)

    TESTARRAYLOG(3, op%mat_X)
    
    ! RESULT: mat_X


    if (lly == 1) then ! LLY
      !--------------------------------------------------------
      !                /  -1    dM  \
      ! calculate  Tr  | M   * ---- | 
      !                \        dE  /
    
      tracek=zero

      do lm1 = 1, lmmaxd
        do lm2 = 1, lmmaxd
          gtdpde = zero
          do il1 = 1, naez*lmmaxd
            gtdpde = gtdpde + op%mat_x(il1,lm2)*dpde_local(il1,lm1)
          enddo
          tracek = tracek + mssq(lm1,lm2,1)*gtdpde
        enddo
      enddo

      bztr2 = bztr2 + tracek*volcub(ikpoint)
      !--------------------------------------------------------
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
!    allocate(DGref_buffer(lmmaxd,lmmaxd,naclsd))

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
    if (any(shape(Gref_buffer) /= [lmmaxd,lmmaxd,naclsd,naez])) deallocate(Gref_buffer, stat=ist)
    
    ! Note: some MPI implementations might need the use of MPI_Alloc_mem
    allocate(Gref_buffer(lmmaxd,lmmaxd,naclsd,naez), stat=ist)
    if (ist /= 0) stop 'failed to allocate Gref_buffer in referenceFourier_com_part1!'
    
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
    integer, intent(in) :: global_atom_id(:) !> mapping trunc. index -> atom index
    integer, intent(in) :: comm

    ! locals
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

      rank = (atom_requested - 1)/num_local_atoms ! block distribution of atoms to ranks

      if (rank /= myrank) then

        tag  = modulo(myrank, TAGMOD)
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
    double complex, intent(out) :: Gref_buffer(:,:,:,:) ! (lmmaxd,lmmaxd,naclsd,naez)
    integer, intent(in) :: naez
    double complex, intent(in) :: Ginp(:,:,:,:)
    integer, intent(in) :: global_atom_id(:) !> mapping trunc. index -> atom index
    integer, intent(in) :: comm

    ! locals
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

      rank = (atom_requested - 1)/num_local_atoms ! block distribution of atoms to ranks

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
    integer, intent(in) :: ezoa(1:) !< index list of periodic image vector
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
  
  
  subroutine dlke0_smat(smat, ind, sparse, eikrm, eikrp, nacls, atom, numn0, indn0, Ginp)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
  ! assume a variable block row sparse matrix description
    double complex, intent(inout) :: smat(:)
    integer, intent(in) :: ind !> site_index
    type(SparseMatrixDescription), intent(in) :: sparse
    double complex, intent(in) :: eikrm(nacls), eikrp(nacls) ! todo: many of these phase factors are real
    integer, intent(in) :: nacls !< number of atoms in the cluster around site ind
    integer, intent(in) :: atom(:) !< dim(nacls) == atom(:,ind)
    integer, intent(in) :: numn0(:) !< dim(naez)
    integer, intent(in) :: indn0(:,:) !< dims(naez,nacls)
    double complex, intent(in) :: Ginp(:,:,:) !< dims(lmmaxd,lmmaxd,nacls)

    integer :: jat, iacls, ni, jnd, ist, gint_iacls
    double complex, allocatable :: GinT(:,:) ! (size(Ginp, 2),size(Ginp, 1)) !< dims(lmmaxd,lmmaxd)

    gint_iacls = -1

    do iacls = 1, nacls
      jat = atom(iacls)
      if (jat < 1) cycle

      do ni = 1, numn0(ind)
        jnd = indn0(ind,ni)
        if (jat == jnd) then

          if (gint_iacls /= iacls) then
            if (gint_iacls == -1) allocate(GinT(size(Ginp, 2),size(Ginp, 1)))
            GinT = transpose(Ginp(:,:,iacls))
            gint_iacls = iacls
          endif

          ist = modify_smat(sparse, ind, jnd, ni, eikrm(iacls), GinT(:,:), smat)

        endif ! jat == jnd
      enddo ! ni

      do ni = 1, numn0(jat)
        jnd = indn0(jat,ni)
        if (ind == jnd) then

          ist = modify_smat(sparse, jat, jnd, ni, eikrp(iacls), Ginp(:,:,iacls), smat)

        endif ! ind == jnd
      enddo ! ni

    enddo ! iacls
    
    deallocate(GinT, stat=ist) ! ignore status
  endsubroutine ! dlke0_smat
  
  
  integer function modify_smat(sparse, ind, jnd, ni, eikr, Gin, smat) result(nOps)
    use SparseMatrixDescription_mod, only: SparseMatrixDescription
    ! assume a variable block row sparse matrix description
    type(SparseMatrixDescription), intent(in) :: sparse
    integer, intent(in) :: ind, jnd, ni !> source site_index, target site_index, ni
    double complex, intent(in) :: eikr ! phase factor
    double complex, intent(in) :: Gin(:,:) !< dims(lmmaxd,lmmaxd)
    double complex, intent(inout) :: smat(:)

    integer :: lm2, lmmax1, lmmax2, is0

    lmmax1 = sparse%kvstr(ind+1) - sparse%kvstr(ind)
    lmmax2 = sparse%kvstr(jnd+1) - sparse%kvstr(jnd)
    
    do lm2 = 1, lmmax2
      is0 = sparse%ka(sparse%ia(ind) + ni-1) + lmmax1*(lm2-1) - 1
      
      smat(is0 + 1:lmmax1 + is0) = smat(is0 + 1:lmmax1 + is0) + eikr * Gin(1:lmmax1,lm2)

    enddo ! lm2

    nOps = lmmax1*lmmax2
  endfunction ! modify_smat
  
endmodule ! kkrmat_mod
