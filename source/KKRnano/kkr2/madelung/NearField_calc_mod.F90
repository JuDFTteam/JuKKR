!> This module implements the near field corrections for KKRnano.
!>
!> @author Elias Rabel

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

! TODO: There are a lot of optimisation possibilities in this near field calculation:
! *) Do not look for near cells at each iteration - especially since it is O(N**2)
! *) Do not recalculate spherical harmonics at Lebedev-points each time
! *) Rotate distance vector between cells parallel to z-axis, use Wigner matrices
!    to rotate back

module NearField_calc_mod
  implicit none
  private
  public :: add_near_field_corr

  contains

  subroutine add_near_field_corr(calc_data, arrays, alat, mpi_comm)
    use Main2Arrays_mod, only: Main2Arrays
    use DensityResults_mod, only: DensityResults
    use CalculationData_mod, only: CalculationData, getAtomData, getDensities
    use BasisAtom_mod, only: BasisAtom
    use RadialMeshData_mod, only: RadialMeshData
    use NearField_com_mod, only: LocalCellInfo, NearFieldCorrection, calculate, create, destroy
    
    type(CalculationData), intent(inout) :: calc_data
    type(Main2Arrays), intent(in) :: arrays
    double precision, intent(in) :: alat
    integer, intent(in) :: mpi_comm

    type(LocalCellInfo), allocatable :: local_cell(:)
    type(NearFieldCorrection), allocatable :: nf_correction(:)
    integer :: num_local_atoms
    integer :: ilocal
    integer :: ispin
    integer :: atom_id
    type(BasisAtom), pointer :: atomdata
    type(RadialMeshData), pointer :: mesh
    type(DensityResults), pointer :: densities

    num_local_atoms = calc_data%num_local_atoms
    allocate(local_cell(num_local_atoms))
    allocate(nf_correction(num_local_atoms))

    do ilocal = 1, num_local_atoms
      atomdata => getAtomData(calc_data, ilocal)
      mesh => atomdata%mesh_ptr
      densities => getDensities(calc_data, ilocal)
      atom_id = calc_data%atom_ids(ilocal)

      call create(nf_correction(ilocal), mesh%irmd, atomdata%potential%lmpot)

      ! setup information on local cells
      ! calculate near-field corrections for each radial point
      call create(local_cell(ilocal), mesh%irmd, atomdata%potential%lmpot, 1)

      local_cell(ilocal)%charge_moments = densities%cmom + densities%cminst
      local_cell(ilocal)%v_intra = atomdata%potential%vons(:,:,1)
      local_cell(ilocal)%radial_points = mesh%r

      call find_near_cells(local_cell(ilocal)%near_cell_indices, &
                           local_cell(ilocal)%near_cell_dist_vec, &
                           arrays%rbasis, arrays%bravais, &
                           atom_id, mesh%rws / alat)

      ! VERY IMPORTANT: convert distance vectors to units of alat!
      local_cell(ilocal)%near_cell_dist_vec = local_cell(ilocal)%near_cell_dist_vec * alat

    enddo ! ilocal

    call calculate(nf_correction, local_cell, calc_data%madelung_calc%clebsch, mpi_comm) ! calc_nf_correction

    do ilocal = 1, num_local_atoms
      atomdata => getAtomData(calc_data, ilocal)

      ! nf correction has to be added to each spin channel
      do ispin = 1, atomdata%potential%nspin
        atomdata%potential%vons(:,:,ispin) = atomdata%potential%vons(:,:,ispin) &
                                           + nf_correction(ilocal)%delta_potential
      enddo ! ispin
    enddo ! ilocal

    call destroy(nf_correction)
    call destroy(local_cell)

  endsubroutine


  !----------------------------------------------------------------------------
  !> Find near cells taking into account periodic boundary conditions.
  !> Returns atom-indices of near cells and distance vectors (pointing to center)
  !>
  !> Returns allocated near_inds, dist_vecs
  !> Caution with units!!!
  !> Criterion for near cells: distance < 2*(Radius bounding sphere)
  !> NOTE: This is only approximately valid - but should be enough for
  !> realistic lattice structures
  !> Note: O(N**2) scaling!
  subroutine find_near_cells(near_inds, dist_vecs, rbasis, bravais, center_ind, radius_bounding)
    integer, allocatable, intent(out) :: near_inds(:)
    double precision, allocatable, intent(out) :: dist_vecs(:,:)

    double precision, intent(in) :: rbasis(:,:)
    double precision, intent(in) :: bravais(3,3)
    integer, intent(in) :: center_ind
    double precision, intent(in) :: radius_bounding

    integer :: ii, num, nx, ny, nz, count_near
    double precision :: center(3), vec(3), four_rws_squared, dist_sq
    integer, parameter :: MAX_NEAR = 64 ! there should be never more than 64 near cells
    double precision :: near_inds_temp(MAX_NEAR), dist_vecs_temp(3,MAX_NEAR)

    double precision, parameter :: TOL = 1.d-6

    CHECKASSERT( .not. allocated(near_inds) )
    CHECKASSERT( .not. allocated(dist_vecs) )

    center = rbasis(:, center_ind)
    four_rws_squared = 4.d0*radius_bounding**2

    count_near = 0
    num = size(rbasis, 2)
    do ii = 1, num

      do nx = -1, 1
        do ny = -1, 1
          do nz = -1, 1
            vec = distance_vec(rbasis(:,ii), center, bravais, nx, ny, nz)
            dist_sq = vec(1)**2 + vec(2)**2 + vec(3)**2
            ! do check with squared distances
            if (dist_sq < four_rws_squared .and. dist_sq > TOL) then
              ! it is important to skip the center!
              count_near = count_near + 1

              CHECKASSERT (count_near <= MAX_NEAR )
              near_inds_temp(count_near) = ii
              dist_vecs_temp(:,count_near) = vec
            endif

          enddo ! nz
        enddo ! ny
      enddo ! nx

    enddo ! ii

    allocate(near_inds(count_near), dist_vecs(3,count_near))
    near_inds = near_inds_temp(1:count_near)
    dist_vecs = dist_vecs_temp(:,1:count_near)
  endsubroutine ! find_near_cells

  !----------------------------------------------------------------------------
  !> Calculate (point2 - point1) taking into account addition of a lattice vector
  !> as specified by integer indices nx, ny, nz
  function distance_vec(point1, point2, bravais, nx, ny, nz) result(dv)
    double precision :: dv(3)
    
    double precision, intent(in) :: point1(3), point2(3), bravais(3,3)
    integer, intent(in) :: nx, ny, nz

    dv(1:3) = point2(1:3) - point1(1:3) + nx*bravais(:,1) + ny*bravais(:,2) + nz*bravais(:,3)

  endfunction ! distance_vec

endmodule ! NearField_calc_mod
