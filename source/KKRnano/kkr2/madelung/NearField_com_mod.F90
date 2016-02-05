!> Calculation of near-field electrostatic corrections:
!> MPI communication routines.
!>
!> @author Elias Rabel

#include "../DebugHelpers/logging_macros.h"

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

module NearField_com_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  USE_LOGGING_MOD
  use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_SIGNALING_NAN
  implicit none
  private
  public :: LocalCellInfo!, create, destroy
  public :: NearFieldCorrection!, create, destroy
!   public :: createNearFieldCorrection, destroyNearFieldCorrection ! deprecated
!   public :: createLocalCellInfo, destroyLocalCellInfo ! deprecated
  public :: calc_nf_correction
  
  type LocalCellInfo
    double precision, allocatable :: charge_moments(:)
    double precision, allocatable :: v_intra(:,:)
    double precision, allocatable :: radial_points(:)
    integer :: critical_index !< apply correction only starting at 'critical_index'
    integer, allocatable :: near_cell_indices(:)
    !> vectors pointing from near cell center to local cell
    double precision, allocatable :: near_cell_dist_vec(:,:)

    contains
    procedure :: create => createLocalCellInfo
    procedure :: destroy => destroyLocalCellInfo

  endtype
  
  type NearFieldCorrection
    double precision, allocatable :: delta_potential(:,:)

    contains
    procedure :: create => createNearFieldCorrection
    procedure :: destroy => destroyNearFieldCorrection
  endtype
  
!   interface create
!     module procedure createNearFieldCorrection, createLocalCellInfo
!   endinterface
!   
  interface destroy
    module procedure destroyNearFieldCorrection, destroyLocalCellInfo
  endinterface
  
  contains
  
  !----------------------------------------------------------------------------
  subroutine calc_nf_correction(nf_correction, local_cells, gaunt, communicator)
    use NearField_kkr_mod, only: IntracellPotential
    use MadelungCalculator_mod, only: MadelungClebschData
!     use ChunkIndex_mod, only: ChunkIndex, getChunkIndex
    use ChunkIndex_mod, only: getRankAndLocalIndex
    use one_sided_commD_mod, only: exposeBufferD, hideBufferD, copyChunksNoSyncD
  
    type(NearFieldCorrection), intent(inout) :: nf_correction(:)
    type(LocalCellInfo), intent(in) :: local_cells(:)
    type(MadelungClebschData), intent(in) :: gaunt
    integer, intent(in) :: communicator
    
    include 'mpif.h'
    
    type(IntracellPotential) :: intra_pot
!    type(ChunkIndex) :: chunk(1)
    integer(kind=4) :: chunk(2,1)
    integer :: npoints, lmpotd
    integer :: num_local_atoms
    integer :: max_npoints
    integer :: ierr, ii
    integer :: irmd
    integer :: nranks
    integer :: ilocal, icell
    double precision, allocatable :: send_buffer(:,:,:)
    double precision, allocatable :: recv_buffer(:,:)
    double precision :: nan
    
    integer :: win
    
    nan = ieee_value(nan, IEEE_SIGNALING_NAN)
    
    num_local_atoms = size(local_cells)
    
    lmpotd = size(local_cells(1)%v_intra, 2)
    
    npoints = 0
    do ii = 1, num_local_atoms
      npoints = max(npoints, size(local_cells(ii)%radial_points))
    enddo ! ii

    call MPI_Comm_size(communicator, nranks, ierr)
    call MPI_Allreduce(npoints, max_npoints, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)

    allocate(send_buffer(max_npoints+1, lmpotd+1, num_local_atoms))
    allocate(recv_buffer(max_npoints+1, lmpotd+1))
    
#ifndef NDEBUG
    send_buffer = nan
    recv_buffer = nan
#endif
    
    do ii = 1, num_local_atoms
      irmd = size(local_cells(ii)%radial_points)
      send_buffer(1:irmd, 1:lmpotd, ii) = local_cells(ii)%v_intra
      send_buffer(max_npoints+1, 1:lmpotd, ii) = local_cells(ii)%charge_moments
      send_buffer(1:irmd, lmpotd+1, ii) = local_cells(ii)%radial_points
      
      ! we also have to store the actual number of mesh points for each local atom
      ! small hack: convert integer number to double to simplify communication
      send_buffer(max_npoints+1, lmpotd+1, ii) = dble(irmd) 
    enddo ! ii

    call exposeBufferD(win, send_buffer, &
                       (max_npoints+1)*(lmpotd+1)*num_local_atoms, &
                       (max_npoints+1)*(lmpotd+1), communicator)

    do ilocal = 1, num_local_atoms
      nf_correction(ilocal)%delta_potential = 0.0d0
      do icell = 1, size(local_cells(ilocal)%near_cell_indices)
      
!       chunk(1) = getChunkIndex(local_cells(ilocal)%near_cell_indices(icell), nranks*num_local_atoms, nranks)
        chunk(:,1) = getRankAndLocalIndex(local_cells(ilocal)%near_cell_indices(icell), nranks*num_local_atoms, nranks)
        
!         call MPI_Win_Lock(MPI_LOCK_SHARED, chunk(1)%owner, 0, win, ierr)
        call MPI_Win_Lock(MPI_LOCK_SHARED, chunk(1,1), 0, win, ierr)
        CHECKASSERT(ierr == 0)

        call copyChunksNoSyncD(recv_buffer, win, chunk, (max_npoints+1)*(lmpotd+1))
        
!         call MPI_Win_Unlock(chunk(1)%owner, win, ierr)
        call MPI_Win_Unlock(chunk(1,1), win, ierr)
        CHECKASSERT(ierr == 0)

        irmd = int(recv_buffer(max_npoints+1, lmpotd+1) + 0.1d0)  ! convert back to integer

        WRITELOG(2,*) "cell, irmd(icell) ", local_cells(ilocal)%near_cell_indices(icell), irmd

        call intra_pot%create(lmpotd, irmd)
        intra_pot%v_intra_values = recv_buffer(1:irmd, 1:lmpotd)
        intra_pot%charge_moments = recv_buffer(max_npoints+1, 1:lmpotd)
        intra_pot%radial_points = recv_buffer(1:irmd, lmpotd+1)

        call intra_pot%init() ! setup intra-cell pot. interpolation

        call add_potential_correction(nf_correction(ilocal)%delta_potential, intra_pot, &
                                       local_cells(ilocal)%radial_points, &
                                       local_cells(ilocal)%near_cell_dist_vec(:,icell), gaunt, &
                                       local_cells(ilocal)%critical_index)
        
        call intra_pot%destroy()
      enddo ! icell

      WRITELOG(2,*) "Near field corrections for local atom ", ilocal
      WRITELOG(2,*) "Near field - delta potential (point/norm):"
      do ii = 1, size(nf_correction(ilocal)%delta_potential, 1)
        WRITELOG(2,*) ii, local_cells(ilocal)%radial_points(ii), &
                      sqrt(dot_product(nf_correction(ilocal)%delta_potential(ii,:), &
                                       nf_correction(ilocal)%delta_potential(ii,:)))
      enddo ! ii

      WRITELOG(2,*) "Near field - delta potential (LM/norm):"
      do ii = 1, size(nf_correction(ilocal)%delta_potential, 2)
        WRITELOG(2,*) ii, sqrt(dot_product(nf_correction(ilocal)%delta_potential(:,ii), &
                               nf_correction(ilocal)%delta_potential(:,ii)))
      enddo ! ii

    enddo ! ilocal

    call hideBufferD(win)
  endsubroutine
  
  !> Set delta_potential to 0.0d0 before first call!
  !>
  !> Apply correction from 'critical_index' to size(radial_points)
  !> Set critical_index to 1 to calculate the near field correction everywhere
  subroutine add_potential_correction(delta_potential, intra_pot, radial_points, dist_vec, gaunt, critical_index)
    use NearField_kkr_mod, only: IntracellPotential
    use MadelungCalculator_mod, only: MadelungClebschData
    use NearField_mod, only: calc_near_field, calc_wrong_contribution_coeff
    
    double precision, intent(inout) :: delta_potential(:,:)
    type(IntracellPotential), intent(inout) :: intra_pot
    double precision, intent(in) :: radial_points(:)
    double precision, intent(in) :: dist_vec(3)
    type(MadelungClebschData), intent(in) :: gaunt
    integer, intent(in) :: critical_index
    
    double precision, allocatable :: coeffs(:)
    integer :: ii, lm, L, M
    
    allocate(coeffs(size(delta_potential, 2)))
    
    call calc_wrong_contribution_coeff(coeffs, dist_vec, intra_pot%charge_moments, gaunt)
    
    ! substract wrong potential
    L = 0
    M = 0
    do lm = 1, size(coeffs)
      do ii = critical_index, size(radial_points)
        delta_potential(ii, lm) = delta_potential(ii, lm) - coeffs(lm) * (-radial_points(ii))**L
      enddo ! ii
      M = M + 1
      if (M > L) then
        L = L+1
        M = -L
      endif
    enddo ! lm
    
    ! add correct contribution
    do ii = critical_index, size(radial_points)
      call calc_near_field(coeffs, radial_points(ii), dist_vec, intra_pot)
      delta_potential(ii, :) = delta_potential(ii, :) + coeffs
    enddo ! ii

  endsubroutine ! add

  !----------------------------------------------------------------------------
  subroutine createLocalCellInfo(self, irmd, lmpotd, critical_index)
    class (LocalCellInfo), intent(inout) :: self
    integer, intent(in) :: irmd
    integer, intent(in) :: lmpotd
    integer, intent(in) :: critical_index

    allocate(self%charge_moments(lmpotd))
    allocate(self%v_intra(irmd, lmpotd))
    allocate(self%radial_points(irmd))
    self%critical_index = critical_index

  endsubroutine ! create

  !----------------------------------------------------------------------------
  elemental subroutine destroyLocalCellInfo(self)
    class (LocalCellInfo), intent(inout) :: self
    integer :: ist
    deallocate(self%charge_moments, stat=ist)
    deallocate(self%v_intra, stat=ist)
    deallocate(self%radial_points, stat=ist)

    if (allocated(self%near_cell_indices)) deallocate(self%near_cell_indices, stat=ist)
    if (allocated(self%near_cell_dist_vec)) deallocate(self%near_cell_dist_vec, stat=ist)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  subroutine createNearFieldCorrection(self, irmd, lmpotd)
    class (NearFieldCorrection), intent(inout) :: self
    integer, intent(in) :: irmd
    integer, intent(in) :: lmpotd

    allocate(self%delta_potential(irmd, lmpotd))
  endsubroutine ! create

  !----------------------------------------------------------------------------
  elemental subroutine destroyNearFieldCorrection(self)
    class (NearFieldCorrection), intent(inout) :: self
    integer :: ist
    deallocate(self%delta_potential, stat=ist)
  endsubroutine ! destroy
  
endmodule
