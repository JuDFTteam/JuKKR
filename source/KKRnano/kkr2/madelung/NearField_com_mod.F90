!> Calculation of near-field electrostatic corrections:
!> MPI communication routines.
!>
!> @author Elias Rabel

#include "../DebugHelpers/logging_macros.h"
#include "macros.h"

module NearField_com_mod
  use Logging_mod, only:    !import no name here, just mention it for the module dependency 
  USE_LOGGING_MOD
  implicit none
  private
  public :: LocalCellInfo, NearFieldCorrection
  public :: create, calculate, destroy

  type LocalCellInfo
    double precision, allocatable :: charge_moments(:)
    double precision, allocatable :: v_intra(:,:)
    double precision, allocatable :: radial_points(:)
    integer :: critical_index !< apply correction only starting at 'critical_index'
    integer, allocatable :: near_cell_indices(:)
    !> vectors pointing from near cell center to local cell
    double precision, allocatable :: near_cell_dist_vec(:,:)
  endtype
  
  type NearFieldCorrection
    double precision, allocatable :: delta_potential(:,:)
  endtype

  interface create
    module procedure createNearFieldCorrection, createLocalCellInfo
  endinterface

  interface calculate
    module procedure calc_nf_correction
  endinterface
   
  interface destroy
    module procedure destroyNearFieldCorrection, destroyLocalCellInfo
  endinterface
  
  contains
  
  !----------------------------------------------------------------------------
  subroutine calc_nf_correction(nfc, lci, gaunt, communicator)
    use NearField_kkr_mod, only: IntracellPotential, create, init, destroy
    use MadelungCalculator_mod, only: MadelungClebschData
    use ChunkIndex_mod, only: getRankAndLocalIndex
    use one_sided_commD_mod, only: exposeBufferD, hideBufferD, copyChunksNoSyncD
#ifndef __GFORTRAN__    
    use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_SIGNALING_NAN
#endif
    type(NearFieldCorrection), intent(inout) :: nfc(:)
    type(LocalCellInfo), intent(in) :: lci(:)
    type(MadelungClebschData), intent(in) :: gaunt
    integer, intent(in) :: communicator
    
    include 'mpif.h'
    
    type(IntracellPotential) :: intra_pot
    integer(kind=4) :: chunk(2,1)
    integer :: num_local_atoms, ila, icell
    integer :: max_npoints, npoints, lmpotd
    integer :: ierr, i1, i2
    integer :: irmd
    integer :: nranks
    double precision, allocatable :: send_buffer(:,:,:), recv_buffer(:,:)
    integer :: win
#ifndef __GFORTRAN__    
    double precision :: nan
    nan = ieee_value(nan, IEEE_SIGNALING_NAN)
#else
    double precision, parameter :: nan = -99.d9
#endif
    
    num_local_atoms = size(lci)
    
    lmpotd = size(lci(1)%v_intra, 2)
    
    npoints = 0
    do ila = 1, num_local_atoms
      npoints = max(npoints, size(lci(ila)%radial_points))
    enddo ! ila

    call MPI_Comm_size(communicator, nranks, ierr)
    call MPI_Allreduce(npoints, max_npoints, 1, MPI_INTEGER, MPI_MAX, communicator, ierr)

    allocate(send_buffer(0:max_npoints,0:lmpotd,num_local_atoms))
    allocate(recv_buffer(0:max_npoints,0:lmpotd))
    
#ifndef NDEBUG
    send_buffer = nan
    recv_buffer = nan
#endif

    do ila = 1, num_local_atoms
      irmd = size(lci(ila)%radial_points)
      ! begin packing
      !--------------------------------------------------------------------------------------------
      send_buffer(1:irmd,1:lmpotd,ila) = lci(ila)%v_intra(:,:)
      send_buffer(0,1:lmpotd,ila) = lci(ila)%charge_moments(:)
      send_buffer(1:irmd,0,ila) = lci(ila)%radial_points(:)
      
      ! we also have to store the actual number of mesh points for each local atom
      ! small hack: convert integer number to double to simplify communication
      send_buffer(0,0,ila) = dble(irmd) 
      !--------------------------------------------------------------------------------------------
      ! end packing
    enddo ! ila

    call exposeBufferD(win, send_buffer, (1+max_npoints)*(1+lmpotd)*num_local_atoms, &
                                         (1+max_npoints)*(1+lmpotd), communicator)

    do ila = 1, num_local_atoms
      nfc(ila)%delta_potential = 0.d0
      do icell = 1, size(lci(ila)%near_cell_indices)

        chunk(:,1) = getRankAndLocalIndex(lci(ila)%near_cell_indices(icell), nranks*num_local_atoms, nranks)

        call MPI_Win_Lock(MPI_LOCK_SHARED, chunk(1,1), 0, win, ierr)
        CHECKASSERT(ierr == 0)

        call copyChunksNoSyncD(recv_buffer, win, chunk, (1+max_npoints)*(1+lmpotd))

        call MPI_Win_Unlock(chunk(1,1), win, ierr)
        CHECKASSERT(ierr == 0)

        ! begin unpacking
        !--------------------------------------------------------------------------------------------
        irmd = int(recv_buffer(0,0) + 0.1d0) ! convert back to integer

        WRITELOG(2,*) "cell, irmd(icell) ", lci(ila)%near_cell_indices(icell), irmd

        call create(intra_pot, lmpotd, irmd)
        intra_pot%v_intra_values(:,:) = recv_buffer(1:irmd,1:lmpotd)
        intra_pot%charge_moments(:) = recv_buffer(0,1:lmpotd)
        intra_pot%radial_points(:) = recv_buffer(1:irmd,0)
        !--------------------------------------------------------------------------------------------
        ! end unpacking

        call init(intra_pot) ! setup intra-cell pot. interpolation

        call add_potential_correction(nfc(ila)%delta_potential, intra_pot, &
                                      lci(ila)%radial_points, &
                                      lci(ila)%near_cell_dist_vec(:,icell), gaunt, &
                                      lci(ila)%critical_index)
        
        call destroy(intra_pot)
      enddo ! icell

      WRITELOG(2,*) "Near field corrections for local atom ", ila
      WRITELOG(2,*) "Near field - delta potential (point/norm):"
      do i1 = 1, size(nfc(ila)%delta_potential, 1)
        WRITELOG(2,*) i1, lci(ila)%radial_points(i1), &
                      sqrt(dot_product(nfc(ila)%delta_potential(i1,:), &
                                       nfc(ila)%delta_potential(i1,:)))
      enddo ! i1

      WRITELOG(2,*) "Near field - delta potential (LM/norm):"
      do i2 = 1, size(nfc(ila)%delta_potential, 2)
        WRITELOG(2,*) i2, sqrt(dot_product(nfc(ila)%delta_potential(:,i2), nfc(ila)%delta_potential(:,i2)))
      enddo ! i2

    enddo ! ila

    call hideBufferD(win)
  endsubroutine ! calc
  
  !> Set delta_potential to 0.d0 before first call!
  !>
  !> Apply correction from 'critical_index' to size(radial_points)
  !> Set critical_index to 1 to calculate the near field correction everywhere
  subroutine add_potential_correction(delta_potential, intra_pot, radial_points, dist_vec, gaunt, critical_index)
    use NearField_kkr_mod, only: IntracellPotential
    use MadelungCalculator_mod, only: MadelungClebschData
    use NearField_mod, only: calculate
    
    double precision, intent(inout) :: delta_potential(:,:)
    type(IntracellPotential), intent(inout) :: intra_pot
    double precision, intent(in) :: radial_points(:)
    double precision, intent(in) :: dist_vec(3)
    type(MadelungClebschData), intent(in) :: gaunt
    integer, intent(in) :: critical_index
    
    double precision, allocatable :: coeffs(:)
    integer :: ii, lm, ell, emm
    
    allocate(coeffs(size(delta_potential, 2)))
    
    call calculate(coeffs, dist_vec, intra_pot%charge_moments, gaunt) ! calc_wrong_contribution_coeff
    
    ! substract wrong potential
    ell = 0
    emm = 0
    do lm = 1, size(coeffs)
      do ii = critical_index, size(radial_points)
        delta_potential(ii,lm) = delta_potential(ii,lm) - coeffs(lm) * (-radial_points(ii))**ell
      enddo ! ii
      emm = emm + 1
      if (emm > ell) then
        ell = ell + 1
        emm = -ell
      endif
    enddo ! lm

    ! add correct contribution
    do ii = critical_index, size(radial_points)
      call calculate(coeffs, radial_points(ii), dist_vec, intra_pot) ! calc_near_field
      delta_potential(ii,:) = delta_potential(ii,:) + coeffs
    enddo ! ii

  endsubroutine ! add

  !----------------------------------------------------------------------------
  subroutine createLocalCellInfo(self, irmd, lmpotd, critical_index)
    type(LocalCellInfo), intent(inout) :: self
    integer, intent(in) :: irmd
    integer, intent(in) :: lmpotd
    integer, intent(in) :: critical_index

    allocate(self%charge_moments(lmpotd))
    allocate(self%v_intra(irmd,lmpotd))
    allocate(self%radial_points(irmd))
    self%critical_index = critical_index
  endsubroutine ! create

  !----------------------------------------------------------------------------
  elemental subroutine destroyLocalCellInfo(self)
    type(LocalCellInfo), intent(inout) :: self
    
    integer :: ist
    deallocate(self%charge_moments, self%v_intra, self%radial_points, stat=ist)
    if (allocated(self%near_cell_indices)) deallocate(self%near_cell_indices, stat=ist)
    if (allocated(self%near_cell_dist_vec)) deallocate(self%near_cell_dist_vec, stat=ist)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  subroutine createNearFieldCorrection(self, irmd, lmpotd)
    type(NearFieldCorrection), intent(inout) :: self
    integer, intent(in) :: irmd
    integer, intent(in) :: lmpotd

    allocate(self%delta_potential(irmd,lmpotd))
  endsubroutine ! create

  !----------------------------------------------------------------------------
  elemental subroutine destroyNearFieldCorrection(self)
    type(NearFieldCorrection), intent(inout) :: self
    integer :: ist
    deallocate(self%delta_potential, stat=ist)
  endsubroutine ! destroy

endmodule ! NearField_com_mod
