!> @author Elias Rabel

#include "DebugHelpers/test_macros.h"

! Add emesh, ... ????

module CalculationData_mod

  use MadelungCalculator_mod
  use RadialMeshData_mod
  use CellData_mod
  use BasisAtom_mod
  use KKRresults_mod
  use DensityResults_mod
  use JijData_mod
  use LDAUData_mod
  use BroydenData_mod

  implicit none

  type CalculationData
    PRIVATE
    !integer :: atoms_per_proc
    integer :: num_local_atoms  ! <= atoms_per_proc
    integer, allocatable :: atom_ids(:)

    ! atom local data - different for every atom
    type (RadialMeshData), pointer     :: mesh_array(:)         => null()
    type (CellData), pointer           :: cell_array(:)         => null()
    type (BasisAtom), pointer          :: atomdata_array(:)     => null()
    type (KKRresults), pointer         :: kkr_array(:)          => null()
    type (DensityResults), pointer     :: densities_array(:)    => null()
    type (MadelungLatticeSum), pointer :: madelung_sum_array(:) => null()

    type (LDAUData), pointer           :: ldau_data_array(:)    => null()
    type (JijData), pointer            :: jij_data_array(:)     => null()
    type (BroydenData), pointer        :: broyden_array(:)      => null()

    ! global data - same for every atom
    !type (MadelungCalculator) :: madelung_calc
    !type (ShapeGauntCoefficients) :: shgaunts
    !type (GauntCoefficients) :: gaunts
    !type (EnergyMesh) :: emesh

  end type

  CONTAINS

  !----------------------------------------------------------------------------
  ! TODO: atoms_per_procs * num_procs MUST BE = naez
  ! rank = 0,1,..., num_atom_ranks-1
  subroutine createCalculationData(calc_data, dims, params, my_mpi)
    use KKRnanoParallel_mod
    use DimParams_mod
    use InputParams_mod
    implicit none

    type (CalculationData), intent(inout) :: calc_data
    type (DimParams), intent(in)  :: dims
    type (InputParams), intent(in):: params
    type (KKRnanoParallel), intent(in) :: my_mpi

    integer :: atoms_per_proc
    integer :: num_local_atoms
    integer, external :: mapblock
    integer :: ii
    integer :: atom_rank

    atoms_per_proc = dims%atoms_per_proc
    num_local_atoms = atoms_per_proc !TODO

    calc_data%num_local_atoms = num_local_atoms

    allocate(calc_data%mesh_array(num_local_atoms))
    allocate(calc_data%cell_array(num_local_atoms))
    allocate(calc_data%atomdata_array(num_local_atoms))
    allocate(calc_data%kkr_array(num_local_atoms))
    allocate(calc_data%densities_array(num_local_atoms))
    allocate(calc_data%madelung_sum_array(num_local_atoms))

    allocate(calc_data%ldau_data_array(num_local_atoms))
    allocate(calc_data%jij_data_array(num_local_atoms))
    allocate(calc_data%broyden_array(num_local_atoms))
    allocate(calc_data%atom_ids(num_local_atoms))

    atom_rank = getMyAtomRank(my_mpi)

    ! assign atom ids to processes with atom rank 'atom_rank'
    ASSERT( size(calc_data%atom_ids) == num_local_atoms)
    ASSERT ( getNumAtomRanks(my_mpi) * atoms_per_proc == dims%naez )
    ASSERT ( mod(dims%naez, atoms_per_proc) == 0 )

    do ii = 1, num_local_atoms
      calc_data%atom_ids(ii) = atom_rank * atoms_per_proc + ii
      ASSERT( calc_data%atom_ids(ii) <= dims%naez )
    end do


!    call createKKRresults(kkr, dims)
!    call createDensityResults(densities, dims)
!
!    !!!I1 = getMyAtomId(my_mpi) !assign atom number for the rest of the program
!
!    call createBasisAtom(atomdata, I1, dims%lpot, dims%nspind, dims%irmind, dims%irmd)
!    call openBasisAtomDAFile(atomdata, 37, "atoms")
!    call readBasisAtomDA(atomdata, 37, I1)
!    call closeBasisAtomDAFile(37)
!
!    if (isInMasterGroup(my_mpi)) then
!      call openBasisAtomPotentialDAFile(atomdata, 37, "vpotnew")
!      call readBasisAtomPotentialDA(atomdata, 37, I1)
!      call closeBasisAtomPotentialDAFile(37)
!    end if
!
!    call createCellData(cell, dims%irid, (2*dims%LPOT+1)**2, dims%nfund)
!    call openCellDataDAFile(cell, 37 , "cells")
!    call readCellDataDA(cell, 37, getCellIndex(atomdata))
!    call closeCellDataDAFile(37)
!
!    call associateBasisAtomCell(atomdata, cell)
!
!    call createRadialMeshData(mesh, dims%irmd, dims%ipand)
!    call openRadialMeshDataDAFile(mesh, 37 , "meshes")
!    call readRadialMeshDataDA(mesh, 37, I1)
!    call closeRadialMeshDataDAFile(37)
!
!    call associateBasisAtomMesh(atomdata, mesh)
!
!    call createLDAUData(ldau_data, params%ldau, dims%irmd, dims%lmaxd, dims%nspind)
!    call createJijData(jij_data, params%jij, params%rcutjij, dims%nxijd,dims%lmmaxd,dims%nspind)
!
!    call createBroydenData(broyden, dims%ntird, dims%itdbryd, params%imix, params%mixing)
!
!    call createMadelungCalculator(madelung_calc, dims%lmaxd, params%ALAT, params%RMAX, params%GMAX, &
!                                  arrays%BRAVAIS, dims%NMAXD, dims%ISHLD)
!
!    call createMadelungLatticeSum(madelung_sum, madelung_calc, dims%naez)
!
!    call calculateMadelungLatticeSum(madelung_sum, I1, arrays%rbasis)


  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyCalculationData(calc_data)

    implicit none

    type (CalculationData), intent(inout) :: calc_data

    deallocate(calc_data%mesh_array)
    deallocate(calc_data%cell_array)
    deallocate(calc_data%atomdata_array)
    deallocate(calc_data%kkr_array)
    deallocate(calc_data%densities_array)
    deallocate(calc_data%madelung_sum_array)

    deallocate(calc_data%ldau_data_array)
    deallocate(calc_data%jij_data_array)
    deallocate(calc_data%broyden_array)
    deallocate(calc_data%atom_ids)
  end subroutine

  !----------------------------------------------------------------------------
  integer function getNumLocalAtoms(calc_data)
    type (CalculationData), intent(in) :: calc_data

    getNumLocalAtoms = calc_data%num_local_atoms
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to atomdata for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getAtomData(calc_data, local_atom_index)
    implicit none
    type (BasisAtom), pointer :: getAtomData ! return value

    type (CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getAtomData => calc_data%atomdata_array(local_atom_index)
  end function

! ==================== Helper routines ========================================
  !subroutine

end module CalculationData_mod
