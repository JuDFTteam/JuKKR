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
  implicit none
  private
  public :: KKRMAT01_new

  double complex, allocatable :: full(:,:)
  double complex, parameter :: zero=(0.d0,0.d0), cone=(1.d0,0.d0)

  CONTAINS

!------------------------------------------------------------------------------
!> Solves multiple scattering problem for every k-point.
!>
!> Returns diagonal k-integrated part of Green's function in GS.
subroutine KKRMAT01_new(solv, kkr_op, precond, BZKP,NOFKS,GS,VOLCUB, TMATLL, ALAT,NSYMAT,RR, &
                        GINP, lmmaxd, trunc2atom_index, communicator, iguess_data)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD
  use InitialGuess_mod, only: InitialGuess, iguess_set_k_ind
  use jij_calc_mod, only: global_jij_data, kkrjij
  use SolverStats_mod, only: SolverStats, reset_stats
  use TFQMRSolver_mod, only: TFQMRSolver
  use BCPOperator_mod, only: BCPOperator
  use KKROperator_mod, only: KKROperator, get_ms_workspace
  use MultScatData_mod, only: MultScatData
  use ClusterInfo_mod, only: ClusterInfo
  use one_sided_commZ_mod, only: ChunkIndex
  ! ************************************************************************
  !   performs k-space integration,
  !   determines scattering path operator (g(k,e)-t**-1)**-1 and
  !   Greens function of the real system -> GS(*,*,*,*),
  ! ------------------------------------------------------------------------

  class (TFQMRSolver), intent(inout) :: solv
  class (KKROperator), intent(inout) :: kkr_op
  class (BCPOperator), intent(inout) :: precond
  integer, intent(in) :: lmmaxd
  integer, intent(in) :: trunc2atom_index(:)
  integer, intent(in) :: communicator
  type (InitialGuess), intent(inout) :: iguess_data
  double precision, intent(in) :: ALAT
  integer, intent(in) :: NOFKS
  integer, intent(in) :: NSYMAT
  double complex, intent(inout) :: TMATLL(:,:,:) ! (lmmaxd,lmmaxd,naez)
  doublecomplex, intent(inout) :: GINP(:,:,:,:) ! (lmmaxd,lmmaxd,naclsd,nclsd)
  double complex, intent(out) ::  GS(:,:,:,:) ! (lmmaxd,lmmaxd,NSYMAXD,num_local_atoms)
  double precision, intent(in) :: BZKP(:,:)
  double precision, intent(in) :: VOLCUB(:)
  double precision, intent(in) :: RR(:,0:) 

  type (MultScatData), pointer :: ms
  type (ClusterInfo), pointer :: cluster_info
  integer :: k_point_index
  double complex, allocatable :: G_diag(:,:,:)
  integer::  site_lm_size
  integer :: iat
  integer :: num_local_atoms
  integer :: naclsd
  integer :: naez
  type (SolverStats) :: total_stats

  ! array dimensions

  ms => get_ms_workspace(kkr_op)
  cluster_info => ms%cluster_info

  naez = cluster_info%naez_trc
  naclsd = cluster_info%naclsd

  site_lm_size = NAEZ*LMMAXD

  num_local_atoms = size(ms%atom_indices)

  !-----------------------------------------------------------------------
  ! Allocate arrays
  !-----------------------------------------------------------------------

  allocate(G_diag(lmmaxd,lmmaxd,num_local_atoms))

  ! WARNING: Symmetry assumptions might have been used that are
  ! not valid in cases of non-local potential (e.g. for Spin-Orbit coupling)
  ! ---> use sit
  !      G(n,n',L,L')(-k) = G(n',n,L',L)(k)

  GS = zero

  TESTARRAYLOG(3, GINP)

  if (global_jij_data%do_jij_calculation) then
    global_jij_data%GSXIJ = zero
  end if

  call reset_stats(total_stats)

!==============================================================================
  do k_point_index = 1, NOFKS                       ! K-POINT-LOOP
!==============================================================================

    WRITELOG(4, *) "k-point ", k_point_index

    ! select right slot for storing initial guess
    call iguess_set_k_ind(iguess_data, k_point_index)

    ! Get the scattering path operator for k-point BZKP(:, k_point_index)
    ! output: ms%mat_X

    call kloopbody(solv, kkr_op, precond, BZKP(:,k_point_index), TMATLL, GINP, ALAT, RR, trunc2atom_index, communicator, iguess_data)

    call getGreenDiag(G_diag, ms%mat_X, ms%atom_indices, ms%sparse%kvstr)

    ! TODO: use mat_X to calculate Jij

    ! ----------- Integrate Scattering Path operator over k-points --> GS -----
    ! Note: here k-integration only in irreducible wedge
    call greenKSummation(G_diag, GS, VOLCUB(k_point_index), ms%atom_indices, NSYMAT, lmmaxd)
    ! -------------------------------------------------------------------------

    if (global_jij_data%do_jij_calculation) then
      !communicate off-diagonal elements and multiply with exp-factor
      call KKRJIJ( BZKP(:,k_point_index),VOLCUB(k_point_index), NSYMAT,NAEZ,ms%atom_indices(1), global_jij_data%NXIJ, global_jij_data%IXCP,global_jij_data%ZKRXIJ, &
      ms%mat_X, global_jij_data%GSXIJ, communicator, lmmaxd, global_jij_data%nxijd)
    endif

    do iat = 1, size(ms%atom_indices)
      TESTARRAYLOG(3, GS(:,:,:,iat))
    enddo ! iat

!==============================================================================
  enddo ! KPT = 1,NOFKS
!==============================================================================

  ! Cleanup
  deallocate(G_diag)

  total_stats = solv%get_total_stats()

  WRITELOG(3, *) "Max. TFQMR residual for this E-point: ", total_stats%max_residual
  WRITELOG(3, *) "Max. num iterations for this E-point: ", total_stats%max_iterations
  WRITELOG(3, *) "Sum of iterations for this E-point:   ", total_stats%sum_iterations

  call solv%reset_total_stats()

endsubroutine KKRMAT01_new

  
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
! - NO! EIKRM and EIKRP are SWAPPED in call to DLKE0 !!!!
subroutine referenceFourier_com(GLLH, sparse, kpoint, alat, nacls, atom, numn0, &
              indn0, rr, ezoa, GINP, EIKRM, EIKRP, trunc2atom_index, communicator)
  use dlke0_smat_mod, only: DLKE0_smat
  use SparseMatrixDescription_mod, only: SparseMatrixDescription
  use one_sided_commZ_mod, only: ChunkIndex, getOwner, getLocalInd, exposeBufferZ, copyChunksNoSyncZ, hideBufferZ

  include 'mpif.h'
  double complex, intent(inout) :: GLLH(:)
  type(SparseMatrixDescription), intent(in) :: sparse
  double precision, intent(in) :: kpoint(3)
  double precision, intent(in) :: alat
  integer, intent(in) :: nacls(:)
  integer, intent(in) :: atom(:,:)
  integer, intent(in) :: numn0(:)
  integer, intent(in) :: indn0(:,:)

  double precision, intent(in) :: rr(:,:)
  integer, intent(in) :: ezoa(:,:)
  double complex, intent(inout) :: GINP(:,:,:,:)

  ! work arrays
  double complex, intent(inout) :: EIKRM(:)   ! dim: naclsd
  double complex, intent(inout) :: EIKRP(:)

  !> mapping trunc. index -> atom index
  integer, intent(in) :: trunc2atom_index(:)
  integer, intent(in) :: communicator

  ! local
  integer site_index
  integer naez
  integer nrd
  integer naclsd
  integer lmmaxd
  integer num_local_atoms
  integer atom_requested
  double complex, allocatable :: Gref_buffer(:,:,:)
  type(ChunkIndex) :: chunk_inds(1)
  integer :: win
  integer :: nranks
  integer :: ierr

  naez = size(nacls)
  nrd = size(rr, 2) - 1  ! because rr has dim (0:nrd)
  lmmaxd = size(GINP,1)
  naclsd = size(GINP, 3)
  num_local_atoms = size(GINP, 4)

  ! checks
  ASSERT(lmmaxd == size(GINP,2))
  ASSERT(naclsd == size(eikrm))
  ASSERT(naclsd == size(eikrp))
  ASSERT(naez == size(trunc2atom_index))

  ! Note: some MPI implementations might need
  ! the use of MPI_Alloc_mem (for GINP)
  allocate(Gref_buffer(lmmaxd, lmmaxd, naclsd))

  call MPI_Comm_size(communicator, nranks, ierr)

  ! share GINP with all other processes in 'communicator'
  call exposeBufferZ(win, GINP, lmmaxd*lmmaxd*naclsd*num_local_atoms, lmmaxd*lmmaxd*naclsd, communicator)

  GLLH = zero
  do site_index = 1, NAEZ

    call DLKE1(ALAT,NACLS(site_index),RR,EZOA(:,site_index), kpoint,EIKRM,EIKRP, nrd, naclsd)

    ! get GINP(:,:,:)[trunc2atom_index(site_index)]

    atom_requested = trunc2atom_index(site_index)
    chunk_inds(1)%owner = getOwner(atom_requested, num_local_atoms * nranks, nranks)
    chunk_inds(1)%local_ind = getLocalInd(atom_requested, num_local_atoms * nranks, nranks)

#ifndef IDENTICAL_REF
    call MPI_Win_Lock(MPI_LOCK_SHARED, chunk_inds(1)%owner, 0, win, ierr)
    CHECKASSERT(ierr == 0)

    call copyChunksNoSyncZ(Gref_buffer, win, chunk_inds, lmmaxd*lmmaxd*naclsd)

    call MPI_Win_Unlock(chunk_inds(1)%owner, win, ierr)
    CHECKASSERT(ierr == 0)
#else
    Gref_buffer(:,:,:) = GINP(:,:,:,1) ! use this if all Grefs are the same
#endif

    call DLKE0_smat(site_index,GLLH,sparse%ia,sparse%ka,sparse%kvstr,EIKRM,EIKRP, &
                    NACLS(site_index), ATOM(:,site_index),NUMN0,INDN0, Gref_buffer, naez, lmmaxd, naclsd)
  enddo ! site_index

  call hideBufferZ(win)

  deallocate(Gref_buffer)

endsubroutine referenceFourier_com

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

    start = kvstr(atom_index)-1
    lmmax1 = kvstr(atom_index+1) - kvstr(atom_index)

    ASSERT(lmmax1 == size(G_diag, 1))
    ASSERT(lmmax1 == size(G_diag, 2))

    G_diag(1:lmmax1,:,ii) = mat_X(start+1:start+lmmax1,((ii-1)*lmmax1+1):(ii*lmmax1))

  enddo ! ii

endsubroutine getGreenDiag


!------------------------------------------------------------------------------
!> Summation of Green's function over k-points. Has to be called for every k-point
!> TODO: it would be better to do the k-space-symmetry treatment separately ???
!> this routine creates nsymat copies of the same solution
!> set gs to 0 before first call
!> in: gllke1
!> inout: gs (set to 0 before first call)
subroutine greenKSummation(g_diag, gs, k_point_weight, atom_indices, nsymat, lmmaxd)
  integer, parameter :: nsymaxd = 48

  integer, intent(in) :: lmmaxd
  integer, intent(in) :: atom_indices(:) ! not used except for its array length information

  double complex, intent(in) :: g_diag(lmmaxd,lmmaxd,size(atom_indices))
  double complex, intent(inout) :: gs(lmmaxd,lmmaxd,nsymaxd,size(atom_indices))

  integer, intent(in) :: nsymat
  double precision, intent(in) :: k_point_weight

  ! -------- local ------------------

  integer :: isym
  integer :: ii !< local atom index
! integer :: iat

  do ii = 1, size(atom_indices)
!   iat = atom_indices(ii)

    ! perform the k-space integration for diagonal element of green's function of atom iat

    do isym = 1, nsymat
      gs(:,:,isym,ii) = gs(:,:,isym,ii) + k_point_weight * g_diag(:,:,ii)
    enddo ! isym
  enddo ! ii
  
endsubroutine ! greenKSummation

!------------------------------------------------------------------------------
!> Calculate scattering path operator for 'kpoint'.
!>
!> Input are the \Delta T and the realspace G_ref (GINP).
!> Solution is stored in ms%mat_X.
!> Scattering path operator is calculated for atoms given in
!> ms%atom_indices(:)
subroutine kloopbody(solv, kkr_op, precond, kpoint, TMATLL, GINP, ALAT, RR, trunc2atom_index, communicator, iguess_data)
  use fillKKRMatrix_mod, only: buildKKRCoeffMatrix, buildRightHandSide, solveFull, convertToFullMatrix
  use fillKKRMatrix_mod, only: dumpDenseMatrix, dumpDenseMatrixFormatted, dumpSparseMatrixData, dumpSparseMatrixDataFormatted
  use TFQMRSolver_mod, only: TFQMRSolver
  use SparseMatrixDescription_mod, only: dumpSparseMatrixDescription
  use InitialGuess_mod, only: InitialGuess, iguess_load, iguess_save
  use TEST_lcutoff_mod, only: cutoffmode, DEBUG_dump_matrix
  use KKROperator_mod, only: KKROperator
  use BCPOperator_mod, only: BCPOperator
  use MultScatData_mod, only: MultScatData
  use ClusterInfo_mod, only: ClusterInfo

  USE_ARRAYLOG_MOD
  USE_LOGGING_MOD

  class(TFQMRSolver), intent(inout) :: solv
  class(KKROperator), intent(inout) :: kkr_op
  class(BCPOperator), intent(inout) :: precond
  integer, intent(in) :: trunc2atom_index(:)
  integer, intent(in) :: communicator
  type(InitialGuess), intent(inout) :: iguess_data
  double precision, intent(in) :: ALAT
  double precision, intent(in) :: kpoint(3)
  doublecomplex, intent(inout) :: GINP(:,:,:,:) ! todo: check intent ! dim: lmmaxd, lmmaxd, naclsd, nclsd
  double precision, intent(in) :: RR(:,0:)
  doublecomplex, intent(inout) :: TMATLL(:,:,:) ! todo: check intent

  integer :: NAEZ
  logical :: initial_zero
  type (MultScatData), pointer :: ms
  type (ClusterInfo), pointer :: cluster_info
  integer :: lmmaxd

  ms => kkr_op%get_ms_workspace()

  lmmaxd = ms%lmmaxd
  naez = ms%naez
  cluster_info => ms%cluster_info

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
  ! - NO! EIKRM and EIKRP are SWAPPED in call to DLKE0 !!!!

  !TODO: use an alternative referenceFourier to work around a bug on BGQ
  !call referenceFourier_com(ms%GLLH, ms%sparse, kpoint, alat, &

! if the following macro is defined, don't use MPI RMA locks
! not using locks does not scale well

#ifdef NO_LOCKS_MPI
#define REF_FOURIER referenceFourier_com_fenced
#else
#define REF_FOURIER referenceFourier_com
#endif

  call REF_FOURIER (ms%GLLH, ms%sparse, kpoint, alat, cluster_info%nacls_trc, cluster_info%atom_trc,  cluster_info%numn0_trc, cluster_info%indn0_trc, &
                    rr, cluster_info%ezoa_trc, GINP, ms%EIKRM, ms%EIKRP, trunc2atom_index, communicator)

  TESTARRAYLOG(3, ms%GLLH)

  !----------------------------------------------------------------------------
  call buildKKRCoeffMatrix(ms%GLLH, TMATLL, ms%lmmaxd, naez, ms%sparse)
  !----------------------------------------------------------------------------

  TESTARRAYLOG(3, ms%GLLH)

  ! ==> now GLLH holds (1 - Delta_t * G_ref)

  ! Now solve the linear matrix equation A*X = b (b is also a matrix),
  ! where A = (1 - Delta_t*G_ref) (inverse of scattering path operator)
  ! and b = Delta_t

  !===================================================================
  ! 3) solve linear set of equations by iterative TFQMR scheme
  !    solve (1 - \Delta t * G_ref) X = \Delta t
  !    the solution X is the scattering path operator

  call buildRightHandSide(ms%mat_B, TMATLL, lmmaxd, ms%atom_indices, ms%sparse%kvstr)

  initial_zero = .true.
  if (iguess_data%iguess == 1) then
    initial_zero = .false.
    call iguess_load(iguess_data, ms%mat_X)
  end if

  call solv%set_initial_zero(initial_zero)

  call precond%calc(ms%GLLH)  ! calculate preconditioner from sparse matrix data

  if (cutoffmode == 3 .or. cutoffmode == 0) then

    call solv%solve(ms%mat_X, ms%mat_B)

    if (DEBUG_dump_matrix) then
      call dumpSparseMatrixDescription(ms%sparse, "matrix_desc.dat")
      call dumpSparseMatrixData(ms%GLLH, "matrix.unf")
      call dumpSparseMatrixDataFormatted(ms%GLLH, "matrix_form.dat")
      call dumpDenseMatrix(ms%mat_X, "solution.unf")
      call dumpDenseMatrixFormatted(ms%mat_X, "solution_form.dat")
      call dumpDenseMatrix(ms%mat_B, "rhs.unf")
      call dumpDenseMatrixFormatted(ms%mat_B, "rhs_form.dat")
    endif ! DEBUG_dump_matrix
  endif

  ! store the initial guess in previously selected slot
  ! (selected with 'iguess_set_k_ind')
  call iguess_save(iguess_data, ms%mat_X)

  TESTARRAYLOG(3, ms%mat_B)

  ! ALTERNATIVE: direct solution with LAPACK
  if (cutoffmode == 4) then
    if (.not. allocated(full)) then
      allocate(full(size(ms%mat_B,1), size(ms%mat_B,1)))
    end if
    call convertToFullMatrix(ms%GLLH, ms%sparse%ia, ms%sparse%ja, ms%sparse%ka, ms%sparse%kvstr, ms%sparse%kvstr, full)
    TESTARRAYLOG(3, full)
    call solveFull(full, ms%mat_B)
    ms%mat_X = ms%mat_B
  endif

  TESTARRAYLOG(3, ms%mat_X)
  ! RESULT: mat_X

endsubroutine kloopbody


!------------------------------------------------------------------------------
!> Alternative implementation of 'referenceFourier_com' to prevent a bug that occured on
!> BGQ regarding MPI RMA locks.
!>
!> Uses fence calls instead of locks.
!> Might not perform and scale as well as referenceFourier_com
subroutine referenceFourier_com_fenced(GLLH, sparse, kpoint, alat, nacls, atom, numn0, indn0, rr, ezoa, GINP, EIKRM, EIKRP, trunc2atom_index, communicator)
  use dlke0_smat_mod, only: DLKE0_smat
  use SparseMatrixDescription_mod, only: SparseMatrixDescription
  use one_sided_commZ_mod, only: ChunkIndex, getOwner, getLocalInd, exposeBufferZ, fenceZ, copyChunksNoSyncZ, hideBufferZ

  include 'mpif.h'
  double complex, intent(inout) :: GLLH(:)
  type(SparseMatrixDescription), intent(in) :: sparse
  double precision, intent(in) :: kpoint(3)
  double precision, intent(in) :: alat
  integer, intent(in) :: nacls(:)
  integer, intent(in) :: atom(:,:)
  integer, intent(in) :: numn0(:)
  integer, intent(in) :: indn0(:,:)

  double precision, intent(in) :: rr(:,:)
  integer, intent(in) :: ezoa(:,:)
  double complex, intent(inout) :: GINP(:,:,:,:)

  ! work arrays
  double complex, intent(inout) :: EIKRM(:)   ! dim: naclsd
  double complex, intent(inout) :: EIKRP(:)

  !> mapping trunc. index -> atom index
  integer, intent(in) :: trunc2atom_index(:)
  integer, intent(in) :: communicator

  ! local
  integer site_index
  integer naez
  integer nrd
  integer naclsd
  integer lmmaxd
  integer num_local_atoms
  integer atom_requested
  double complex, allocatable :: Gref_buffer(:,:,:)
  type (ChunkIndex) :: chunk_inds(1)
  integer :: win
  integer :: nranks
  integer :: ierr
  integer :: naez_max

  naez = size(nacls)
  nrd = size(rr, 2) - 1  ! because rr has dim (0:nrd)
  lmmaxd = size(GINP,1)
  naclsd = size(GINP, 3)
  num_local_atoms = size(GINP, 4)

  ! checks
  ASSERT(lmmaxd == size(GINP,2))
  ASSERT(naclsd == size(eikrm))
  ASSERT(naclsd == size(eikrp))
  ASSERT(naez == size(trunc2atom_index))

  ! Note: some MPI implementations might need
  ! the use of MPI_Alloc_mem
  allocate(Gref_buffer(lmmaxd, lmmaxd, naclsd))

  call MPI_Comm_size(communicator, nranks, ierr)

  ! get maximum number of atoms of all truncation zones
  call MPI_Allreduce(naez, naez_max, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)

  ! share GINP with all other processes in 'communicator'
  call exposeBufferZ(win, GINP, lmmaxd*lmmaxd*naclsd*num_local_atoms, lmmaxd*lmmaxd*naclsd, communicator)

  GLLH = zero

  ! loop up to naez_max to ensure that each rank does the same amount of fence calls
  do site_index = 1, naez_max

    call fenceZ(win)

    if (site_index <= naez) then
      call DLKE1(ALAT,NACLS(site_index),RR,EZOA(:,site_index), kpoint,EIKRM,EIKRP, nrd, naclsd)

      ! get GINP(:,:,:)[trunc2atom_index(site_index)]

      atom_requested = trunc2atom_index(site_index)
      chunk_inds(1)%owner = getOwner(atom_requested, num_local_atoms * nranks, nranks)
      chunk_inds(1)%local_ind = getLocalInd(atom_requested, num_local_atoms * nranks, nranks)

      call copyChunksNoSyncZ(Gref_buffer, win, chunk_inds, lmmaxd*lmmaxd*naclsd)
      !!!Gref_buffer(:,:,:) = GINP(:,:,:,1) ! use this if all Grefs are the same
    end if

    call fenceZ(win) ! ensures that data has arrived in Gref_buffer

    if (site_index <= naez) then
      call DLKE0_smat(site_index,GLLH,sparse%ia,sparse%ka,sparse%kvstr,EIKRM,EIKRP, NACLS(site_index), ATOM(:,site_index),NUMN0,INDN0, Gref_buffer, naez, lmmaxd, naclsd)
    endif ! si < N
  enddo ! site_index

  call hideBufferZ(win)

  deallocate(Gref_buffer)

endsubroutine referenceFourier_com_fenced

endmodule kkrmat_new_mod
