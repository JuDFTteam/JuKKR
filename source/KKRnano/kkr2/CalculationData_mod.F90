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

  use GauntCoefficients_mod
  use ShapeGauntCoefficients_mod

  implicit none

  public :: CalculationData
  public :: createCalculationData
  public :: destroyCalculationData
  private :: constructEverything

  type CalculationData
    PRIVATE
    !integer :: atoms_per_proc
    integer :: num_local_atoms  ! <= atoms_per_proc
    integer, allocatable :: atom_ids(:)

    ! atom local data - different for each atom
    type (RadialMeshData), pointer     :: mesh_array(:)         => null()
    type (CellData), pointer           :: cell_array(:)         => null()
    type (BasisAtom), pointer          :: atomdata_array(:)     => null()
    type (KKRresults), pointer         :: kkr_array(:)          => null()
    type (DensityResults), pointer     :: densities_array(:)    => null()
    type (MadelungLatticeSum), pointer :: madelung_sum_array(:) => null()

    type (LDAUData), pointer           :: ldau_data_array(:)    => null()
    type (JijData), pointer            :: jij_data_array(:)     => null()
    type (BroydenData), pointer        :: broyden_array(:)      => null()

    ! global data - same for each atom
    type (MadelungCalculator), pointer     :: madelung_calc     => null()
    type (ShapeGauntCoefficients), pointer :: shgaunts          => null()
    type (GauntCoefficients), pointer      :: gaunts            => null()
    !type (EnergyMesh), pointer :: emesh                         => null()

  end type

  CONTAINS

  !----------------------------------------------------------------------------
  ! TODO: atoms_per_procs * num_procs MUST BE = naez
  ! rank = 0,1,..., num_atom_ranks-1
  subroutine createCalculationData(calc_data, dims, params, arrays, my_mpi)
    use KKRnanoParallel_mod
    use DimParams_mod
    use InputParams_mod
    use Main2Arrays_mod
    implicit none

    type (CalculationData), intent(inout) :: calc_data
    type (DimParams), intent(in)   :: dims
    type (InputParams), intent(in) :: params
    type (Main2Arrays), intent(in) :: arrays
    type (KKRnanoParallel), intent(in) :: my_mpi

    integer :: atoms_per_proc
    integer :: num_local_atoms
    integer, external :: mapblock
    integer :: ii
    integer :: atom_rank

    atoms_per_proc = dims%atoms_per_proc
    num_local_atoms = atoms_per_proc !TODO

    calc_data%num_local_atoms = num_local_atoms

    ! one datastructure for each local atom
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

    ! These datastructures are the same for all atoms
    allocate(calc_data%madelung_calc)
    allocate(calc_data%gaunts)
    allocate(calc_data%shgaunts)

    atom_rank = getMyAtomRank(my_mpi)

    ! assign atom ids to processes with atom rank 'atom_rank'
    ! E.g. for 2 atoms per proc:
    ! process 1 treats atoms 1,2
    ! process 2 treats atoms 3,4 and so on
    ! FOR USE OF TRUNCATION THESE atoms have to be close together!!!

    ASSERT( size(calc_data%atom_ids) == num_local_atoms)
    ASSERT ( getNumAtomRanks(my_mpi) * atoms_per_proc == dims%naez )
    ASSERT ( mod(dims%naez, atoms_per_proc) == 0 )

    do ii = 1, num_local_atoms
      calc_data%atom_ids(ii) = atom_rank * atoms_per_proc + ii
      ASSERT( calc_data%atom_ids(ii) <= dims%naez )
    end do

    ! Now construct all datastructures and calculate initial data
    call constructEverything(calc_data, dims, params, arrays, my_mpi)

  end subroutine

  !----------------------------------------------------------------------------
  !> Calculate Madelung Lattice sums for all local atoms.
  subroutine prepareMadelung(calc_data, arrays)
    use Main2Arrays_mod
    implicit none

    type (CalculationData), intent(inout) :: calc_data
    type (Main2Arrays), intent(in):: arrays

    ! ----- locals ------
    integer :: I1
    integer :: ilocal
    type (MadelungLatticeSum), pointer :: madelung_sum

    do ilocal = 1, calc_data%num_local_atoms
      I1 = calc_data%atom_ids(ilocal)
      madelung_sum => calc_data%madelung_sum_array(ilocal)
      call calculateMadelungLatticeSum(madelung_sum, I1, arrays%rbasis)
    end do

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyCalculationData(calc_data)

    implicit none

    type (CalculationData), intent(inout) :: calc_data

    ! ---- locals ----
    integer :: ilocal
    integer :: I1
    type (KKRresults), pointer :: kkr
    type (DensityResults), pointer :: densities
    type (BasisAtom), pointer :: atomdata
    type (CellData), pointer :: cell
    type (RadialMeshData), pointer :: mesh
    type (LDAUData), pointer :: ldau_data
    type (JijData), pointer :: jij_data
    type (BroydenData), pointer :: broyden
    type (MadelungLatticeSum), pointer :: madelung_sum

    ! loop over all LOCAL atoms
    !--------------------------------------------------------------------------
    do ilocal = 1, calc_data%num_local_atoms
    !--------------------------------------------------------------------------

      I1 = calc_data%atom_ids(ilocal)

      kkr       => calc_data%kkr_array(ilocal)
      densities => calc_data%densities_array(ilocal)
      atomdata  => calc_data%atomdata_array(ilocal)
      cell      => calc_data%cell_array(ilocal)
      mesh      => calc_data%mesh_array(ilocal)
      ldau_data => calc_data%ldau_data_array(ilocal)
      jij_data  => calc_data%jij_data_array(ilocal)
      broyden   => calc_data%broyden_array(ilocal)
      madelung_sum   => calc_data%madelung_sum_array(ilocal)

      call destroyMadelungLatticeSum(madelung_sum)

      call destroyBasisAtom(atomdata)
      call destroyCellData(cell)
      call destroyRadialMeshData(mesh)

      call destroyBroydenData(broyden)
      call destroyLDAUData(ldau_data)
      call destroyJijData(jij_data)
      call destroyDensityResults(densities)
      call destroyKKRresults(kkr)

    end do

    call destroyMadelungCalculator(calc_data%madelung_calc)
    call destroyGauntCoefficients(calc_data%gaunts)
    call destroyShapeGauntCoefficients(calc_data%shgaunts)

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

  !----------------------------------------------------------------------------
  !> Returns reference to kkr(results) for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getKKR(calc_data, local_atom_index)
    implicit none
    type (KKRresults), pointer :: getKKR ! return value

    type (CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getKKR => calc_data%kkr_array(local_atom_index)
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to Madelung sum for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getMadelungSum(calc_data, local_atom_index)
    implicit none
    type (MadelungLatticeSum), pointer :: getMadelungSum ! return value

    type (CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getMadelungSum => calc_data%madelung_sum_array(local_atom_index)
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to mesh for atom with LOCAL atom index
  !> 'local_atom_index'.
!  function getMesh(calc_data, local_atom_index)
!    implicit none
!    type (RadialMeshData), pointer :: getMesh ! return value
!
!    type (CalculationData), intent(in) :: calc_data
!    integer, intent(in) :: local_atom_index
!
!    getMesh => calc_data%mesh_array(local_atom_index)
!  end function

  !----------------------------------------------------------------------------
  !> Returns reference to density results for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getDensities(calc_data, local_atom_index)
    implicit none
    type (DensityResults), pointer :: getDensities ! return value

    type (CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getDensities => calc_data%densities_array(local_atom_index)
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to LDA+U data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getLDAUData(calc_data, local_atom_index)
    implicit none
    type (LDAUData), pointer :: getLDAUData ! return value

    type (CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getLDAUData => calc_data%ldau_data_array(local_atom_index)
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to Jij-data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getJijData(calc_data, local_atom_index)
    implicit none
    type (JijData), pointer :: getJijData ! return value

    type (CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getJijData => calc_data%jij_data_array(local_atom_index)
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to Broyden data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getBroyden(calc_data, local_atom_index)
    implicit none
    type (BroydenData), pointer :: getBroyden ! return value

    type (CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getBroyden => calc_data%broyden_array(local_atom_index)
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to Gaunt coefficients.
  function getGaunts(calc_data)
    implicit none
    type (GauntCoefficients), pointer :: getGaunts ! return value

    type (CalculationData), intent(in) :: calc_data

    getGaunts => calc_data%gaunts
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to Shape-Gaunt coefficients.
  function getShapeGaunts(calc_data)
    implicit none
    type (ShapeGauntCoefficients), pointer :: getShapeGaunts ! return value

    type (CalculationData), intent(in) :: calc_data

    getShapeGaunts => calc_data%shgaunts
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to Madelung calculator.
  function getMadelungCalculator(calc_data)
    implicit none
    type (MadelungCalculator), pointer :: getMadelungCalculator ! return value

    type (CalculationData), intent(in) :: calc_data

    getMadelungCalculator => calc_data%madelung_calc
  end function

! ==================== Helper routines ========================================

  !----------------------------------------------------------------------------
  !> Helper routine: called by createCalculationData.
  subroutine constructEverything(calc_data, dims, params, arrays, my_mpi)
    use KKRnanoParallel_mod
    use DimParams_mod
    use InputParams_mod
    use Main2Arrays_mod
    implicit none

    type (CalculationData), intent(inout) :: calc_data
    type (DimParams), intent(in)  :: dims
    type (InputParams), intent(in):: params
    type (Main2Arrays), intent(in):: arrays
    type (KKRnanoParallel), intent(in) :: my_mpi

    ! ----- locals ------
    integer :: I1
    integer :: ilocal
    type (KKRresults), pointer :: kkr
    type (DensityResults), pointer :: densities
    type (BasisAtom), pointer :: atomdata
    type (CellData), pointer :: cell
    type (RadialMeshData), pointer :: mesh
    type (LDAUData), pointer :: ldau_data
    type (JijData), pointer :: jij_data
    type (BroydenData), pointer :: broyden
    type (MadelungLatticeSum), pointer :: madelung_sum

    if (isInMasterGroup(my_mpi)) call openBasisAtomDAFile(atomdata, 37, "atoms")
    call openCellDataDAFile(cell, 38 , "cells")
    call openRadialMeshDataDAFile(mesh, 39 , "meshes")

    ! loop over all LOCAL atoms
    !--------------------------------------------------------------------------
    do ilocal = 1, calc_data%num_local_atoms
    !--------------------------------------------------------------------------

      I1 = calc_data%atom_ids(ilocal)

      kkr       => calc_data%kkr_array(ilocal)
      densities => calc_data%densities_array(ilocal)
      atomdata  => calc_data%atomdata_array(ilocal)
      cell      => calc_data%cell_array(ilocal)
      mesh      => calc_data%mesh_array(ilocal)
      ldau_data => calc_data%ldau_data_array(ilocal)
      jij_data  => calc_data%jij_data_array(ilocal)
      broyden   => calc_data%broyden_array(ilocal)
      madelung_sum   => calc_data%madelung_sum_array(ilocal)

      call createKKRresults(kkr, dims)
      call createDensityResults(densities, dims)

      call createBasisAtom(atomdata, I1, dims%lpot, &
                           dims%nspind, dims%irmind, dims%irmd)
      call readBasisAtomDA(atomdata, 37, I1)

      ASSERT( atomdata%atom_index == I1 )

      if (isInMasterGroup(my_mpi)) then
        call readBasisAtomPotentialDA(atomdata, 37, I1)
      end if

      call createCellData(cell, dims%irid, (2*dims%LPOT+1)**2, dims%nfund)
      call readCellDataDA(cell, 38, getCellIndex(atomdata))

      call associateBasisAtomCell(atomdata, cell)

      call createRadialMeshData(mesh, dims%irmd, dims%ipand)
      call readRadialMeshDataDA(mesh, 39, I1)

      call associateBasisAtomMesh(atomdata, mesh)

      call createLDAUData(ldau_data, params%ldau, dims%irmd, dims%lmaxd, &
                          dims%nspind)
      call createJijData(jij_data, params%jij, params%rcutjij, dims%nxijd, &
                         dims%lmmaxd,dims%nspind)

      call createBroydenData(broyden, dims%ntird, dims%itdbryd, &
                             params%imix, params%mixing)

      call createMadelungCalculator(calc_data%madelung_calc, dims%lmaxd, &
                                  params%ALAT, params%RMAX, params%GMAX, &
                                  arrays%BRAVAIS, dims%NMAXD, dims%ISHLD)

      call createMadelungLatticeSum(madelung_sum, calc_data%madelung_calc, dims%naez)

    !--------------------------------------------------------------------------
    end do
    !--------------------------------------------------------------------------

    if (isInMasterGroup(my_mpi)) call closeBasisAtomPotentialDAFile(37)
    call closeCellDataDAFile(38)
    call closeRadialMeshDataDAFile(39)

    ! calculate Gaunt coefficients
    call createGauntCoefficients(calc_data%gaunts, dims%lmaxd)
    call createShapeGauntCoefficients(calc_data%shgaunts, dims%lmaxd)

  end subroutine

end module CalculationData_mod
