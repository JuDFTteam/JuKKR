!> Multiple Scattering problem: Solving the Dyson equation for all k-points,
!> using *real-space* G_ref and the t-matrices
!>
!> Output: Brillouin-zone integrated diagonal elements of structural Green's
!> function

#include "../DebugHelpers/logging_macros.h"
#include "../DebugHelpers/test_array_log.h"
#include "../DebugHelpers/test_macros.h"

module kkrmat_new_mod

use SparseMatrixDescription_mod
use ClusterInfo_mod

implicit none

double complex, allocatable, dimension(:, :), save :: full

type MultScatData
  type (SparseMatrixDescription) :: sparse

  double complex, dimension(:,:), allocatable :: mat_B
  double complex, dimension(:,:), allocatable :: mat_X
  double complex, allocatable :: EIKRP(:)
  double complex, allocatable :: EIKRM(:)
  double complex, allocatable, dimension(:) :: GLLH

  integer, allocatable :: atom_indices(:)
  integer :: lmmaxd
  integer :: naez
  type (ClusterInfo), pointer :: cluster_info

end type

CONTAINS

!------------------------------------------------------------------------------
!> Create workspace for multiple scattering calculation.
subroutine createMultScatData(ms, cluster_info, lmmaxd, atom_indices)
  use TEST_lcutoff_mod, only: lmarray
  use fillKKRMatrix_mod, only: getKKRMatrixStructure
  implicit none
  type (MultScatData), intent(inout) :: ms
  type (ClusterInfo), target  :: cluster_info
  integer, intent(in) :: lmmaxd
  integer, intent(in) :: atom_indices(:)

  integer :: sum_cluster
  integer :: naez
  integer :: naclsd

  ms%cluster_info => cluster_info

  sum_cluster = sum(cluster_info%numn0_trc)
  naez = size(cluster_info%indn0_trc, 1)
  naclsd = size(cluster_info%indn0_trc, 2)
  ms%lmmaxd = lmmaxd
  ms%naez = naez

  allocate(ms%atom_indices, source=atom_indices)

  call createSparseMatrixDescription(ms%sparse, naez, sum_cluster)

  call getKKRMatrixStructure(lmarray, cluster_info%numn0_trc, &
                             cluster_info%indn0_trc, ms%sparse)

  allocate(ms%mat_B(ms%sparse%kvstr(naez+1)-1,LMMAXD * size(atom_indices)))
  allocate(ms%mat_X(ms%sparse%kvstr(naez+1)-1,LMMAXD * size(atom_indices)))

  ! allocate memory for sparse matrix
  allocate(ms%GLLH(getNNZ(ms%sparse)))

  allocate(ms%eikrm(naclsd))
  allocate(ms%eikrp(naclsd))
end subroutine

!------------------------------------------------------------------------------
!> Destroys workspace for multiple scattering calculation.
subroutine destroyMultScatData(ms)
  implicit none
  type (MultScatData), intent(inout) :: ms

  deallocate(ms%eikrp)
  deallocate(ms%eikrm)
  deallocate(ms%GLLH)
  deallocate(ms%mat_X)
  deallocate(ms%mat_B)

  call destroySparseMatrixDescription(ms%sparse)

  deallocate(ms%atom_indices)
end subroutine

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
                                indn0, rr, ezoa, GINP, EIKRM, EIKRP, &
                                trunc2atom_index, communicator)
  use dlke0_smat_mod
  use SparseMatrixDescription_mod
  use one_sided_commZ_mod
  implicit none

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
  double complex, parameter :: CZERO= ( 0.0D0,0.0D0)
  type (ChunkIndex) :: chunk_inds(1)
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
  ! the use of MPI_Alloc_mem
  allocate(Gref_buffer(lmmaxd, lmmaxd, naclsd))

  call MPI_Comm_size(communicator, nranks, ierr)

  ! share GINP with all other processes in 'communicator'
  call exposeBufferZ(win, GINP, lmmaxd*lmmaxd*naclsd*num_local_atoms, &
                     lmmaxd*lmmaxd*naclsd, communicator)

  GLLH = CZERO
  do site_index = 1,NAEZ

    call DLKE1(ALAT,NACLS(site_index),RR,EZOA(:,site_index), &
               kpoint,EIKRM,EIKRP, &
               nrd, naclsd)

    ! get GINP(:,:,:)[trunc2atom_index(site_index)]

    atom_requested = trunc2atom_index(site_index)
    chunk_inds(1)%owner = getOwner(atom_requested, num_local_atoms * nranks, nranks)
    chunk_inds(1)%local_ind = getLocalInd(atom_requested, num_local_atoms * nranks, nranks)

    call MPI_Win_Lock(MPI_LOCK_SHARED, chunk_inds(1)%owner, 0, win, ierr)
    CHECKASSERT(ierr == 0)

    call copyChunksNoSyncZ(Gref_buffer, win, chunk_inds, lmmaxd*lmmaxd*naclsd)
    !!!Gref_buffer(:,:,:) = GINP(:,:,:,1) ! use this if all Grefs are the same

    call MPI_Win_Unlock(chunk_inds(1)%owner, win, ierr)
    CHECKASSERT(ierr == 0)

    call DLKE0_smat(site_index,GLLH,sparse%ia,sparse%ka,sparse%kvstr,EIKRM,EIKRP, &
                    NACLS(site_index), ATOM(:,site_index),NUMN0,INDN0, &
                    Gref_buffer, &
                    naez, lmmaxd, naclsd)
  end do

  call hideBufferZ(win)

end subroutine

!------------------------------------------------------------------------------
!> Copy the diagonal elements G_{LL'}^{nn'} of the Green's-function,
!> dependent on (k,E) into matrix G_diag
subroutine getGreenDiag(G_diag, mat_X, atom_indices, kvstr)
  implicit none

  double complex, intent(out) :: G_diag (:,:,:) ! dim lmmaxd*lmmaxd*num_local_atoms
  double complex, intent(in) :: mat_X (:,:)
  integer, intent(in) :: atom_indices(:)
  integer, intent(in) :: kvstr(:)

  double complex, parameter :: CZERO =(0.0D0,0.0D0)
  integer :: atom_index
  integer :: lm1
  integer :: lmmax1
  integer :: start
  integer :: ii !< local atom index

  !                                      nn
  !         Copy the diagonal elements G_LL' of the Green's-function,
  !         dependent on (k,E) into matrix G_diag
  !         (n = n' = atom_index)

  ASSERT(size(atom_indices) == size(G_diag, 3))

  G_diag = CZERO

  do ii = 1, size(atom_indices)

    atom_index = atom_indices(ii)

    start = kvstr(atom_index) - 1
    lmmax1 = kvstr(atom_index + 1) - kvstr(atom_index)

    ASSERT(lmmax1 == size(G_diag, 1))
    ASSERT(lmmax1 == size(G_diag, 2))

    do lm1 = 1, lmmax1
       G_diag(lm1, :, ii) = mat_X(start + lm1, ((ii - 1) * lmmax1 + 1) : (ii * lmmax1))
    end do

  end do

end subroutine


!------------------------------------------------------------------------------
!> Summation of Green's function over k-points. Has to be called for every k-point
!> TODO: it would be better to do the k-space-symmetry treatment separately ???
!> This routine creates NSYMAT copies of the same solution
!> Set GS to 0 before first call
!> in: GLLKE1
!> inout: GS (set to 0 before first call)
subroutine greenKSummation(G_diag, GS, k_point_weight, &
                           atom_indices, NSYMAT, lmmaxd)
  implicit none
  integer, parameter :: NSYMAXD = 48

  integer, intent(in) :: lmmaxd
  integer, intent(in) :: atom_indices(:)

  double complex :: G_diag(lmmaxd,lmmaxd,size(atom_indices))
  double complex :: GS(lmmaxd,lmmaxd,NSYMAXD,size(atom_indices))

  integer :: NSYMAT
  double precision :: k_point_weight

  ! -------- local ------------------

  integer :: LM1
  integer :: LM2
  integer :: ISYM
  integer :: IAT
  integer :: ii !< local atom index

  do ii = 1, size(atom_indices)

    iat = atom_indices(ii)

      !         Perform the k-space integration for diagonal element of
      !         Green's function of atom IAT

      do ISYM = 1,NSYMAT
        do LM2=1,LMMAXD
          do LM1=1,LMMAXD
            GS(LM1,LM2,ISYM, ii) = GS(LM1,LM2,ISYM, ii) + k_point_weight * G_diag(LM1,LM2,ii)
          end do
        end do
      end do        ! ISYM = 1,NSYMAT

    end do !ii
end subroutine

!------------------------------------------------------------------------------
!> Calculate scattering path operator for 'kpoint'.
!>
!> Input are the \Delta T and the realspace G_ref (GINP).
!> Solution is stored in ms%mat_X.
!> Scattering path operator is calculated for atoms given in
!> ms%atom_indices(:)
subroutine kloopbody(ms, kpoint, &
                     TMATLL, GINP, ALAT, &
                     RR, QMRBOUND, &
                     trunc2atom_index, communicator, iguess_data, &
                     solver_opts, stats)

  use fillKKRMatrix_mod
  use mminvmod_mod
  use dlke0_smat_mod
  use SparseMatrixDescription_mod
  use InitialGuess_mod
  use TEST_lcutoff_mod, only: cutoffmode, DEBUG_dump_matrix
  use SolverOptions_mod
  use SolverStats_mod

  USE_ARRAYLOG_MOD
  USE_LOGGING_MOD
  implicit none

  type (MultScatData), intent(inout) :: ms

  integer, intent(in), dimension(:) :: trunc2atom_index
  integer, intent(in) :: communicator
  type(InitialGuess), intent(inout) :: iguess_data
  type(SolverOptions), intent(in) :: solver_opts
  type (SolverStats), intent(inout) :: stats

  double precision :: ALAT
  double precision :: kpoint(3)
  doublecomplex :: GINP(:,:,:,:) ! dim: lmmaxd, lmmaxd, naclsd, nclsd

  double precision :: QMRBOUND
  double precision :: RR(:,0:)
  doublecomplex :: TMATLL(:,:,:)

  !-------- local ---------
  double complex, parameter :: CONE = ( 1.0D0,0.0D0)
  double complex, parameter :: CZERO= ( 0.0D0,0.0D0)

  integer :: NAEZ
  logical :: initial_zero
  type (ClusterInfo), pointer :: cluster_info


  integer :: lmmaxd

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
   call referenceFourier_com_fenced(ms%GLLH, ms%sparse, kpoint, alat, &
                            cluster_info%nacls_trc, cluster_info%atom_trc, &
                            cluster_info%numn0_trc, cluster_info%indn0_trc, &
                            rr, cluster_info%ezoa_trc, GINP, ms%EIKRM, ms%EIKRP, &
                            trunc2atom_index, communicator)

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

  if (cutoffmode == 3) then
    call MMINVMOD_new(ms%GLLH, ms%sparse, ms%mat_X, ms%mat_B, &
                      QMRBOUND, size(ms%mat_B, 2), size(ms%mat_B, 1), initial_zero, stats)

    if (DEBUG_dump_matrix) then
      call dumpSparseMatrixDescription(ms%sparse, "matrix_desc.dat")
      call dumpSparseMatrixData(ms%GLLH, "matrix.unf")
      call dumpSparseMatrixDataFormatted(ms%GLLH, "matrix_form.dat")
      call dumpDenseMatrix(ms%mat_X, "solution.unf")
      call dumpDenseMatrixFormatted(ms%mat_X, "solution_form.dat")
      call dumpDenseMatrix(ms%mat_B, "rhs.unf")
      call dumpDenseMatrixFormatted(ms%mat_B, "rhs_form.dat")
    end if
  endif

  if (cutoffmode == 0) then
    call bcp_solver(ms%GLLH, ms%mat_X, ms%mat_B, qmrbound, cluster_info, solver_opts, ms%sparse, initial_zero)
  end if

  ! store the initial guess in previously selected slot
  ! (selected with 'iguess_set_k_ind')
  call iguess_save(iguess_data, ms%mat_X)

  TESTARRAYLOG(3, ms%mat_B)

  ! ALTERNATIVE: direct solution with LAPACK
  if (cutoffmode == 4) then
    if (.not. allocated(full)) then
      allocate(full(size(ms%mat_B,1), size(ms%mat_B,1)))
    end if
    call convertToFullMatrix(ms%GLLH, ms%sparse%ia, ms%sparse%ja, ms%sparse%ka, &
                                   ms%sparse%kvstr, ms%sparse%kvstr, full)
    TESTARRAYLOG(3, full)
    call solveFull(full, ms%mat_B)
    ms%mat_X = ms%mat_B
  endif

  TESTARRAYLOG(3, ms%mat_X)
  ! RESULT: mat_X

end subroutine

!------------------------------------------------------------------------------
!> Solves multiple scattering problem for every k-point.
!>
!> Returns diagonal k-integrated part of Green's function in GS.
subroutine KKRMAT01_new(BZKP,NOFKS,GS,VOLCUB, &
TMATLL, ALAT,NSYMAT,RR, &
GINP, atom_indices, QMRBOUND, &
lmmaxd, trunc2atom_index, communicator, iguess_data, cluster_info, &
solver_opts)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD
  use InitialGuess_mod
  use jij_calc_mod, only: global_jij_data, kkrjij
  use SolverOptions_mod
  use SolverStats_mod
  implicit none

  !     .. parameters ..
  double complex, parameter :: CZERO= ( 0.0D0,0.0D0)

  ! ************************************************************************
  !   performs k-space integration,
  !   determines scattering path operator (g(k,e)-t**-1)**-1 and
  !   Greens function of the real system -> GS(*,*,*,*),
  ! ------------------------------------------------------------------------

  integer, intent(in) :: lmmaxd
  integer, dimension(:), intent(in) :: atom_indices !< indices of atoms treated at once
  integer, intent(in) :: trunc2atom_index(:)
  integer, intent(in) :: communicator
  type (InitialGuess), intent(inout) :: iguess_data
  type (ClusterInfo), target         :: cluster_info ! in
  type (SolverOptions), intent(in)   :: solver_opts

  double precision:: ALAT

  integer::NOFKS
  integer::NSYMAT

  !double complex :: TMATLL(lmmaxd,lmmaxd,naez)
  double complex :: TMATLL(:,:,:)

  doublecomplex :: GINP(:,:,:,:) ! dim: lmmaxd, lmmaxd, naclsd, nclsd
  !double complex :: GS   (lmmaxd,lmmaxd,NSYMAXD,num_local_atoms)
  double complex ::  GS(:,:,:,:)

  double precision::BZKP(:,:)
  double precision::VOLCUB(:)
  double precision::RR(:,0:)

  double precision::QMRBOUND

! ------- local ----------

  type (MultScatData) :: ms
  integer::k_point_index

  double complex, allocatable, dimension(:,:,:) ::G_diag

  integer::        site_lm_size

  integer :: iat
  integer :: num_local_atoms
  integer :: naclsd
  integer :: naez

  type (SolverStats) :: stats, total_stats

  ! array dimensions

  naez = cluster_info%naez_trc
  naclsd = cluster_info%naclsd

  site_lm_size = NAEZ*LMMAXD

  num_local_atoms = size(atom_indices)

  !-----------------------------------------------------------------------
  ! Allocate arrays
  !-----------------------------------------------------------------------

  allocate(G_diag(lmmaxd,lmmaxd,num_local_atoms))

  ! WARNING: Symmetry assumptions might have been used that are
  ! not valid in cases of non-local potential (e.g. for Spin-Orbit coupling)
  ! ---> use sit
  !      G(n,n',L,L')(-k) = G(n',n,L',L)(k)

  GS = CZERO

  TESTARRAYLOG(3, GINP)

  call createMultScatData(ms, cluster_info, lmmaxd, atom_indices)

  if (global_jij_data%do_jij_calculation) then
    global_jij_data%GSXIJ = CZERO
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

    call kloopbody(ms, BZKP(:, k_point_index), &
                   TMATLL, GINP, ALAT, &
                   RR, QMRBOUND, &
                   trunc2atom_index, communicator, iguess_data, &
                   solver_opts, stats)

    call sum_stats(stats, total_stats)

    call getGreenDiag(G_diag, ms%mat_X, ms%atom_indices, ms%sparse%kvstr)

    ! TODO: use mat_X to calculate Jij

    ! ----------- Integrate Scattering Path operator over k-points --> GS -----
    ! Note: here k-integration only in irreducible wedge
    call greenKSummation(G_diag, &
                         GS, VOLCUB(k_point_index), &
                         atom_indices, NSYMAT, lmmaxd)
    ! -------------------------------------------------------------------------

    if (global_jij_data%do_jij_calculation) then
      !communicate off-diagonal elements and multiply with exp-factor
      call KKRJIJ( BZKP(:,k_point_index),VOLCUB(k_point_index), &
      NSYMAT,NAEZ,atom_indices(1), &
      global_jij_data%NXIJ, global_jij_data%IXCP,global_jij_data%ZKRXIJ, &
      ms%mat_X, &
      global_jij_data%GSXIJ, &
      communicator, &
      lmmaxd, global_jij_data%nxijd)
    end if

    do iat = 1, size(atom_indices)
      TESTARRAYLOG(3, GS(:,:,:,iat))
    end do

!==============================================================================
  end do ! KPT = 1,NOFKS
!==============================================================================

  ! Cleanup
  deallocate(G_diag)

  call destroyMultScatData(ms)

  WRITELOG(3, *) "Max. TFQMR residual for this E-point: ", total_stats%max_residual
  WRITELOG(3, *) "Max. num iterations for this E-point: ", total_stats%max_iterations
  WRITELOG(3, *) "Sum of iterations for this E-point:   ", total_stats%sum_iterations

end subroutine KKRMAT01_new

!------------------------------------------------------------------------------
!> Wrapper for solver with BCP preconditioner. (cutoffmode=0)
subroutine bcp_solver(GLLH, mat_X, mat_B, qmrbound, cluster_info, solver_opts, sparse, initial_zero)

  USE_ARRAYLOG_MOD
  USE_LOGGING_MOD

  use SolverOptions_mod
  use ClusterInfo_mod
  use SparseMatrixDescription_mod
  use mminvmod_bcp_mod
  use mminvmod_mod

  use SolverStats_mod
  implicit none

  double complex GLLH(:)
  double complex mat_X(:,:)
  double complex mat_B(:,:)
  double precision, intent(in) :: qmrbound
  type (ClusterInfo), intent(in)  :: cluster_info
  type (SolverOptions), intent(in) :: solver_opts
  type(SparseMatrixDescription) :: sparse
  logical :: initial_zero

  type(SolverStats) :: stats

  integer naezd
  integer lmmaxd
  integer nlen
  integer num_columns
  integer blocks_per_row
  double complex, allocatable :: GLLHBLCK(:,:)
  double complex, allocatable :: temp(:,:)

  naezd = cluster_info%naez_trc
  lmmaxd = sparse%kvstr(2) - sparse%kvstr(1)
  num_columns = size(mat_X, 2)

  nlen = naezd*lmmaxd

  CHECKASSERT(nlen == size(mat_X, 1))
  CHECKASSERT(nlen == size(mat_B, 1))
  CHECKASSERT(num_columns == size(mat_B, 2))

  allocate(temp(nlen, num_columns))

  if (solver_opts%BCP == 1) then
    CHECKASSERT(naezd == solver_opts%NATBLD*solver_opts%XDIM * &
                         solver_opts%YDIM*solver_opts%ZDIM)

    allocate(GLLHBLCK(solver_opts%NATBLD*LMMAXD, naezd*LMMAXD))

    blocks_per_row = cluster_info%numn0_trc(1)

    ! SEVERE restriction of preconditioner code:
    ! The number of non-zero blocks must be the same in each row
    ! this is not the case for all crystal structures (e.g. perovskite)
    CHECKASSERT(all(cluster_info%numn0_trc == blocks_per_row ))

    ! setup preconditioning matrix
    call BCPWUPPER(GLLH,GLLHBLCK,NAEZD,cluster_info%NUMN0_trc,cluster_info%INDN0_trc, &
                   lmmaxd, solver_opts%natbld, solver_opts%xdim, solver_opts%ydim, solver_opts%zdim, &
                   blocks_per_row)

  endif

  call MMINVMOD_bcp(GLLH, sparse, mat_X, mat_B, qmrbound, num_columns, nlen, initial_zero, &
                          solver_opts%bcp, lmmaxd, GLLHBLCK, cluster_info%numn0_trc, &
                          cluster_info%indn0_trc, temp, solver_opts%natbld, &
                          solver_opts%xdim,solver_opts%ydim, solver_opts%zdim)

end subroutine

!------------------------------------------------------------------------------
!> Alternative implementation of 'referenceFourier_com' to prevent a bug that occured on
!> BGQ regarding MPI RMA locks.
!>
!> Uses fence calls instead of locks.
!> Might not perform and scale as well as referenceFourier_com
subroutine referenceFourier_com_fenced(GLLH, sparse, kpoint, alat, nacls, atom, numn0, &
                                indn0, rr, ezoa, GINP, EIKRM, EIKRP, &
                                trunc2atom_index, communicator)
  use dlke0_smat_mod
  use SparseMatrixDescription_mod
  use one_sided_commZ_mod
  implicit none

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
  double complex, parameter :: CZERO= ( 0.0D0,0.0D0)
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
  call exposeBufferZ(win, GINP, lmmaxd*lmmaxd*naclsd*num_local_atoms, &
                     lmmaxd*lmmaxd*naclsd, communicator)

  GLLH = CZERO

  ! loop up to naez_max to ensure that each rank does the same amount of fence calls
  do site_index = 1, naez_max

    call fenceZ(win)

    if (site_index <= naez) then
      call DLKE1(ALAT,NACLS(site_index),RR,EZOA(:,site_index), &
                 kpoint,EIKRM,EIKRP, &
                 nrd, naclsd)

      ! get GINP(:,:,:)[trunc2atom_index(site_index)]

      atom_requested = trunc2atom_index(site_index)
      chunk_inds(1)%owner = getOwner(atom_requested, num_local_atoms * nranks, nranks)
      chunk_inds(1)%local_ind = getLocalInd(atom_requested, num_local_atoms * nranks, nranks)

      call copyChunksNoSyncZ(Gref_buffer, win, chunk_inds, lmmaxd*lmmaxd*naclsd)
      !!!Gref_buffer(:,:,:) = GINP(:,:,:,1) ! use this if all Grefs are the same
    end if

    call fenceZ(win) ! ensures that data has arrived in Gref_buffer

    if (site_index <= naez) then
      call DLKE0_smat(site_index,GLLH,sparse%ia,sparse%ka,sparse%kvstr,EIKRM,EIKRP, &
                      NACLS(site_index), ATOM(:,site_index),NUMN0,INDN0, &
                      Gref_buffer, &
                      naez, lmmaxd, naclsd)
    end if
  end do

  call hideBufferZ(win)

end subroutine

end module
