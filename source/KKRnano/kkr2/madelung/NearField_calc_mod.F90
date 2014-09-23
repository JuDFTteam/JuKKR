!> This module implements the near field corrections for KKRnano.
!>
!> @author Elias Rabel

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

! TODO: There are a lot of optimisation possibilities in this near field calculation:
! *) Do not look for near cells at each iteration - especially since it is O(N**2)
! *) Do not calculate near field correction in non-critical region
! *) Do not recalculate spherical harmonics at Lebedev-points each time
! *) Rotate distance vector between cells parallel to z-axis, use Wigner matrices
!    to rotate back

module NearField_calc_mod
  use NearField_com_mod
  use CalculationData_mod
  use KKRnanoParallel_mod
  use BasisAtom_mod
  use RadialMeshData_mod
  implicit none

  contains

  subroutine add_near_field_corr(calc_data, arrays, alat, my_mpi)
    implicit none
    type(CalculationData), intent(inout) :: calc_data
    type(Main2Arrays), intent(in) :: arrays
    double precision, intent(in) :: alat
    type(KKRnanoParallel), intent(in) :: my_mpi

    type (LocalCellInfo), allocatable :: local_cell(:)
    type (NearFieldCorrection), allocatable :: nf_correction(:)
    integer :: num_local_atoms
    integer :: ilocal
    integer :: ispin
    integer :: atom_ind
    type(BasisAtom), pointer :: atomdata
    type(RadialMeshData), pointer :: mesh
    type(DensityResults), pointer :: densities
    type(MadelungCalculator), pointer :: madelung_calc

    num_local_atoms  = getNumLocalAtoms(calc_data)
    allocate(local_cell(num_local_atoms))
    allocate(nf_correction(num_local_atoms))

    madelung_calc => getMadelungCalculator(calc_data)

    do ilocal = 1, num_local_atoms
      atomdata => getAtomData(calc_data, ilocal)
      mesh => atomdata%mesh_ptr
      densities => getDensities(calc_data, ilocal)
      atom_ind = getAtomIndexOfLocal(calc_data, ilocal)

      call nf_correction(ilocal)%create(mesh%irmd, atomdata%potential%lmpot)

      call local_cell(ilocal)%create(mesh%irmd, atomdata%potential%lmpot, mesh%imt)

      local_cell(ilocal)%charge_moments = densities%cmom + densities%cminst
      local_cell(ilocal)%v_intra = atomdata%potential%vons(:,:,1)
      local_cell(ilocal)%radial_points = mesh%r

      call find_near_cells(local_cell(ilocal)%near_cell_indices, &
                           local_cell(ilocal)%near_cell_dist_vec, &
                           arrays%rbasis, arrays%bravais, &
                           atom_ind, mesh%rws / alat)

      ! VERY IMPORTANT: convert distance vectors to units of alat!
      local_cell(ilocal)%near_cell_dist_vec = local_cell(ilocal)%near_cell_dist_vec * alat

    end do

    call calc_nf_correction(nf_correction, local_cell, madelung_calc%clebsch, &
                            getMySEcommunicator(my_mpi))

    do ilocal = 1, num_local_atoms
      atomdata => getAtomData(calc_data, ilocal)

      ! nf correction has to be added to each spin channel
      do ispin = 1, atomdata%potential%nspin
        atomdata%potential%vons(:,:,ispin) = atomdata%potential%vons(:,:,ispin) &
                                         + nf_correction(ilocal)%delta_potential
      end do
    end do

  end subroutine


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
  subroutine find_near_cells(near_inds, dist_vecs, rbasis, bravais, &
                             center_ind, radius_bounding)
    implicit none
    integer, allocatable, intent(out) :: near_inds(:)
    double precision, allocatable, intent(out) :: dist_vecs(:,:)

    double precision, intent(in) :: rbasis(:,:)
    double precision, intent(in) :: bravais(3,3)
    integer, intent(in) :: center_ind
    double precision, intent(in) :: radius_bounding

    integer :: ii, num
    integer :: nx, ny, nz
    integer :: count_near
    double precision :: center(3), vec(3)
    double precision :: four_rws_squared, dist_sq
    integer, parameter :: MAX_NEAR = 64 ! there should be never more than 64 near cells
    double precision :: near_inds_temp(MAX_NEAR)
    double precision :: dist_vecs_temp(3, MAX_NEAR)

    double precision, parameter :: TOL = 1.d-6

    CHECKASSERT( .not. allocated(near_inds) )
    CHECKASSERT( .not. allocated(dist_vecs) )

    center = rbasis(:, center_ind)
    four_rws_squared = 4.0d0 * radius_bounding**2

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
              dist_vecs_temp(:, count_near) = vec
            end if

          end do
        end do
      end do

    end do

    allocate(near_inds(count_near))
    near_inds = near_inds_temp(1:count_near)
    allocate(dist_vecs(3, count_near))
    dist_vecs = dist_vecs_temp(:, 1:count_near)
  end subroutine

  !----------------------------------------------------------------------------
  !> Calculate (point2 - point1) taking into account addition of a lattice vector
  !> as specified by integer indices nx, ny, nz
  function distance_vec(point1, point2, bravais, nx, ny, nz)
    implicit none
    double precision :: distance_vec(3)
    double precision, intent(in) :: point1(3)
    double precision, intent(in) :: point2(3)
    double precision, intent(in) :: bravais(3,3)
    integer, intent(in) :: nx, ny, nz
    !----------------
    double precision :: vec(3)

    vec = point2 - point1

    distance_vec = vec + nx*bravais(:, 1) + ny*bravais(:, 2) + nz*bravais(:, 3)

  end function

end module NearField_calc_mod
