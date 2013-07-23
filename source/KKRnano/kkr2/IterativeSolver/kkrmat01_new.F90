! cleaned version of kkrmat

#include "../DebugHelpers/logging_macros.h"
#include "../DebugHelpers/test_array_log.h"
#include "../DebugHelpers/test_macros.h"

module kkrmat_new_mod

double complex, allocatable, dimension(:, :), save :: full
!IBM* ALIGN(32, full)

CONTAINS

! WARNING: Symmetry assumptions might have been used that are
! not valid in cases of non-local potential (e.g. for Spin-Orbit coupling)

subroutine KKRMAT01_new(BZKP,NOFKS,GS,VOLCUB, &
TMATLL, &
ALAT,NSYMAT,NAEZ,NACLS,RR,EZOA,ATOM, &
GINP, &
NUMN0,INDN0,atom_indices, &
QMRBOUND, &
lmmaxd, naclsd, trunc2atom_index, communicator)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD
  implicit none

  !     .. parameters ..
  integer, parameter :: NSYMAXD = 48
  double complex, parameter :: CZERO= ( 0.0D0,0.0D0)

  ! ************************************************************************
  !   performs k-space integration,
  !   determines scattering path operator (g(k,e)-t**-1)**-1 and
  !   Greens function of the real system -> GS(*,*,*,*),

  !   NEW VERSION 10.99
  !   up -> left , down -> right, for decimation
  ! ------------------------------------------------------------------------

  integer, intent(in) :: lmmaxd
  integer, intent(in) :: naclsd  !< max. number of atoms in reference cluster
  integer, dimension(:), intent(in) :: atom_indices !< indices of atoms treated at once
  integer, intent(in) :: trunc2atom_index(:)
  integer, intent(in) :: communicator

  !     .. SCALAR ARGUMENTS ..
  double precision:: ALAT

  integer::NAEZ
  integer::NOFKS
  integer::NSYMAT

  double complex :: TMATLL(lmmaxd,lmmaxd,naez)

  doublecomplex :: GINP(:,:,:,:) ! dim: lmmaxd, lmmaxd, naclsd, nclsd
  !double complex :: GS   (lmmaxd,lmmaxd,NSYMAXD,num_local_atoms)
  double complex ::  GS(:,:,:,:)

  double precision::BZKP(:,:)
  double precision::VOLCUB(:)
  double precision::RR(:,0:)

  integer:: NUMN0(:)
  integer:: INDN0(:,:)
  integer:: ATOM(:,:) ! dim naclsd, naez?
  integer:: EZOA(:,:) ! dim naclsd, naez?
  integer:: NACLS(:)

  double precision::QMRBOUND

! ------- local ----------

  double complex :: EIKRP(naclsd)
  double complex :: EIKRM(naclsd)
  integer::k_point_index

  double complex, allocatable, dimension(:,:,:) ::G_diag
  double complex, allocatable, dimension(:) ::GLLH

!IBM* ALIGN(32, GLLH)

  integer::        site_lm_size

  integer :: memory_stat
  logical :: memory_fail

  integer :: iat
  integer :: num_local_atoms

  ! array dimensions

  site_lm_size = NAEZ*LMMAXD

  num_local_atoms = size(atom_indices)

  !-----------------------------------------------------------------------
  ! Allocate arrays
  !-----------------------------------------------------------------------
  memory_stat = 0
  memory_fail = .false.

  allocate(G_diag(lmmaxd,lmmaxd,num_local_atoms), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  if (memory_fail .eqv. .true.) then
    write(*,*) "KKRMAT01: FATAL Error, failure to allocate memory."
    write(*,*) "       Probably out of memory."
    stop
  end if


  ! WARNING: Symmetry assumptions might have been used that are
  ! not valid in cases of non-local potential (e.g. for Spin-Orbit coupling)
  ! ---> use sit
  !      G(n,n',L,L')(-k) = G(n',n,L',L)(k)

  GS = CZERO

  TESTARRAYLOG(3, GINP)

!==============================================================================
  do k_point_index = 1, NOFKS                       ! K-POINT-LOOP
!==============================================================================

    WRITELOG(4, *) "k-point ", k_point_index

    ! Get the scattering path operator for k-point BZKP(:, k_point_index)
    ! output: GLLKE1, NOITER
    ! inout: PRSC
    ! inout (temporary arrays): GLLH, GLLHBLCK
    call kloopbody( G_diag, BZKP(:, k_point_index), TMATLL, GINP, ALAT, &
                   NAEZ, ATOM, EZOA, RR, INDN0, &
                   NUMN0, EIKRM, EIKRP, GLLH, &
                   atom_indices, QMRBOUND, NACLS, lmmaxd, trunc2atom_index, &
                   communicator)

    ! ----------- Integrate Scattering Path operator over k-points --> GS -----
    ! Note: here k-integration only in irreducible wedge
    call greenKSummation(G_diag, &
                         GS, VOLCUB(k_point_index), &
                         atom_indices, NSYMAT, lmmaxd)
    ! -------------------------------------------------------------------------

    do iat = 1, size(atom_indices)
      TESTARRAYLOG(3, GS(:,:,:,iat))
    end do

!==============================================================================
  end do ! KPT = 1,NOFKS
!==============================================================================

  ! ----------------------------------------------------------------
  ! Deallocate arrays
  ! ----------------------------------------------------------------

  deallocate(G_diag)
  if (allocated(GLLH)) deallocate(GLLH)

end subroutine KKRMAT01_new

!------------------------------------------------------------------------------
!> Calculate scattering path operator for 'kpoint'
subroutine kloopbody( G_diag, kpoint, &
                      TMATLL, GINP, ALAT, &
                      NAEZ, ATOM, EZOA, RR, INDN0, &
                      NUMN0, EIKRM, EIKRP, GLLH, &
                      atom_indices, QMRBOUND, NACLS, &
                      lmmaxd, trunc2atom_index, communicator)

  !use initialGuess_store_mod
  use fillKKRMatrix_mod
  use mminvmod_mod
  use dlke0_smat_mod
  use SparseMatrixDescription_mod
  use TEST_lcutoff_mod, only: lmarray, cutoffmode, DEBUG_dump_matrix !TODO: remove

  USE_ARRAYLOG_MOD
  USE_LOGGING_MOD
  implicit none

  integer, intent(in), dimension(:) :: atom_indices !< indices of local atoms
  integer, intent(in), dimension(:) :: trunc2atom_index
  integer, intent(in) :: communicator
  integer :: NAEZ
  double precision :: ALAT
  integer :: ATOM(:,:)         ! dim: naclsd, *
  double precision :: kpoint(3)
  double complex :: EIKRM(:)   ! dim: naclsd
  double complex :: EIKRP(:)
  integer :: EZOA(:,:) ! dim naclsd,*
  doublecomplex :: GINP(:,:,:,:) ! dim: lmmaxd, lmmaxd, naclsd, nclsd
  double complex, allocatable :: GLLH(:)

  double complex :: G_diag(:,:,:)

  integer :: INDN0(:,:)
  integer, intent(in) :: lmmaxd
  integer :: NACLS(:)
  integer :: NUMN0(:)
  double precision :: QMRBOUND
  double precision :: RR(:,0:)
  doublecomplex :: TMATLL(:,:,:)

  !-------- local ---------
  double complex, parameter :: CONE = ( 1.0D0,0.0D0)
  double complex, parameter :: CZERO= ( 0.0D0,0.0D0)

  type (SparseMatrixDescription) :: sparse

  double complex, dimension(:,:), allocatable :: mat_B
  double complex, dimension(:,:), allocatable :: mat_X

  integer :: sum_cluster
  logical :: initial_zero

  integer :: lm1, lm2

  sum_cluster = sum(numn0)

  call createSparseMatrixDescription(sparse, naez, sum_cluster)

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

  call getKKRMatrixStructure(lmarray, numn0, indn0, sparse)

  allocate(mat_B(sparse%kvstr(naez+1)-1,LMMAXD * size(atom_indices)))
  allocate(mat_X(sparse%kvstr(naez+1)-1,LMMAXD * size(atom_indices)))

  if (.not. allocated(GLLH)) then
    allocate(GLLH(getNNZ(sparse)))
  endif

  call referenceFourier_com(GLLH, sparse, kpoint, alat, nacls, atom, numn0, &
                            indn0, rr, ezoa, GINP, EIKRM, EIKRP, &
                            trunc2atom_index, communicator)

!  do site_index = 1, naclsd ! this was just a test
!  do lm2 = 1, lmmaxd
!    do lm1 = 1, lm2
!      if (abs(GINP(lm1, lm2, site_index, 1) - GINP(lm2, lm1, site_index, 1)) > 1e-4) then
!        write(*,*) lm1, lm2, site_index, abs(GINP(lm1, lm2, site_index, 1) - GINP(lm2, lm1, site_index, 1))
!        !STOP
!      end if
!    end do
!  end do
!  end do

  TESTARRAYLOG(3, GLLH)

  !----------------------------------------------------------------------------
  call buildKKRCoeffMatrix(GLLH, TMATLL, lmmaxd, naez, sparse)
  !----------------------------------------------------------------------------

  TESTARRAYLOG(3, GLLH)

  ! ==> now GLLH holds (1 - Delta_t * G_ref)

  ! Now solve the linear matrix equation A*X = b (b is also a matrix),
  ! where A = (1 - Delta_t*G_ref) (inverse of scattering path operator)
  ! and b = Delta_t

  !===================================================================
  ! 3) solve linear set of equations by iterative TFQMR scheme
  !    solve (1 - \Delta t * G_ref) X = \Delta t
  !    the solution X is the scattering path operator

  call buildRightHandSide(mat_B, TMATLL, lmmaxd, atom_indices, sparse%kvstr)

  initial_zero = .true.

  if (cutoffmode == 3) then
    call MMINVMOD_new(GLLH, sparse, mat_X, mat_B, &
                      QMRBOUND, size(mat_B, 2), size(mat_B, 1), initial_zero)

    if (DEBUG_dump_matrix) then
      call dumpSparseMatrixDescription(sparse, "matrix_desc.dat")
      call dumpSparseMatrixData(GLLH, "matrix.unf")
      call dumpSparseMatrixDataFormatted(GLLH, "matrix_form.dat")
      call dumpDenseMatrix(mat_X, "solution.unf")
      call dumpDenseMatrixFormatted(mat_X, "solution_form.dat")
      call dumpDenseMatrix(mat_B, "rhs.unf")
      call dumpDenseMatrixFormatted(mat_B, "rhs_form.dat")
    end if

  end if

  TESTARRAYLOG(4, mat_B)

  ! solve full matrix equation
  if (cutoffmode == 4) then
    if (.not. allocated(full)) then
      allocate(full(size(mat_B,1), size(mat_B,1)))
    end if
    call convertToFullMatrix(GLLH, sparse%ia, sparse%ja, sparse%ka, &
                                   sparse%kvstr, sparse%kvstr, full)
    TESTARRAYLOG(3, full)
    call solveFull(full, mat_B)
    mat_X = mat_B
  endif

  TESTARRAYLOG(4, mat_X)
  !call toOldSolutionFormat(GLLKE1, mat_X, lmmaxd, sparse%kvstr)

  call getGreenDiag(G_diag, mat_X, atom_indices, sparse%kvstr)

  ! solved. Result in G_diag

  !TESTARRAYLOG(3, GLLKE1)

  call destroySparseMatrixDescription(sparse)

  deallocate(mat_B)
  deallocate(mat_X)

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
  use one_sided_commZ_mod, only: copyFromZ_com
  implicit none

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
  integer atom_requested(1)
  double complex, allocatable :: Gref_buffer(:,:,:)
  double complex, parameter :: CZERO= ( 0.0D0,0.0D0)

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

  allocate(Gref_buffer(lmmaxd, lmmaxd, naclsd))

  GLLH = CZERO
  do site_index = 1,NAEZ

    call DLKE1(ALAT,NACLS(site_index),RR,EZOA(:,site_index), &
               kpoint,EIKRM,EIKRP, &
               nrd, naclsd)

    ! get GINP(:,:,:)[trunc2atom_index(site_index)]

    atom_requested(1) = trunc2atom_index(site_index)
    call copyFromZ_com(Gref_buffer, GINP, atom_requested, &
                       lmmaxd*lmmaxd*naclsd, num_local_atoms, communicator)
    !!!Gref_buffer(:,:,:) = GINP(:,:,:) ! use this if all Grefs are the same

    call DLKE0_smat(site_index,GLLH,sparse%ia,sparse%ka,sparse%kvstr,EIKRM,EIKRP, &
                    NACLS(site_index), ATOM(:,site_index),NUMN0,INDN0, &
                    Gref_buffer, &
                    naez, lmmaxd, naclsd)
  end do

end subroutine

!------------------------------------------------------------------------------
!>
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

end module
