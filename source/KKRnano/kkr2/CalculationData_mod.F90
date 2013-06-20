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
  use EnergyResults_mod
  use RefCluster_mod
  use ClusterInfo_mod

  use GauntCoefficients_mod
  use ShapeGauntCoefficients_mod
  use TruncationZone_mod

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
    type (RefCluster), pointer         :: ref_cluster_array(:)  => null()
    type (KKRresults), pointer         :: kkr_array(:)          => null()
    type (DensityResults), pointer     :: densities_array(:)    => null()
    type (EnergyResults), pointer      :: energies_array(:)     => null()
    type (MadelungLatticeSum), pointer :: madelung_sum_array(:) => null()

    type (LDAUData), pointer           :: ldau_data_array(:)    => null()
    type (JijData), pointer            :: jij_data_array(:)     => null()
    type (BroydenData), pointer        :: broyden_array(:)      => null()

    ! global data - same for each local atom
    type (LatticeVectors), pointer         :: lattice_vectors   => null()
    type (MadelungCalculator), pointer     :: madelung_calc     => null()
    type (ShapeGauntCoefficients), pointer :: shgaunts          => null()
    type (GauntCoefficients), pointer      :: gaunts            => null()
    type (TruncationZone), pointer         :: trunc_zone        => null()
    type (ClusterInfo), pointer            :: clusters          => null()
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
    allocate(calc_data%ref_cluster_array(num_local_atoms))
    allocate(calc_data%kkr_array(num_local_atoms))
    allocate(calc_data%densities_array(num_local_atoms))
    allocate(calc_data%energies_array(num_local_atoms))
    allocate(calc_data%madelung_sum_array(num_local_atoms))

    allocate(calc_data%ldau_data_array(num_local_atoms))
    allocate(calc_data%jij_data_array(num_local_atoms))
    allocate(calc_data%broyden_array(num_local_atoms))
    allocate(calc_data%atom_ids(num_local_atoms))

    ! These datastructures are the same for all (local) atoms
    allocate(calc_data%lattice_vectors)
    allocate(calc_data%madelung_calc)
    allocate(calc_data%gaunts)
    allocate(calc_data%shgaunts)
    allocate(calc_data%trunc_zone)
    allocate(calc_data%clusters)

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
    type (EnergyResults), pointer :: energies
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
      energies  => calc_data%energies_array(ilocal)

      call destroyRefCluster(calc_data%ref_cluster_array(ilocal))

      call destroyMadelungLatticeSum(madelung_sum)

      call destroyBasisAtom(atomdata)
      call destroyCellData(cell)
      call destroyRadialMeshData(mesh)

      call destroyBroydenData(broyden)
      call destroyLDAUData(ldau_data)
      call destroyJijData(jij_data)
      call destroyDensityResults(densities)
      call destroyKKRresults(kkr)
      call destroyEnergyResults(energies)

    end do

    call destroyLatticeVectors(calc_data%lattice_vectors)

    call destroyMadelungCalculator(calc_data%madelung_calc)
    call destroyGauntCoefficients(calc_data%gaunts)
    call destroyShapeGauntCoefficients(calc_data%shgaunts)

    deallocate(calc_data%mesh_array)
    deallocate(calc_data%cell_array)
    deallocate(calc_data%atomdata_array)
    deallocate(calc_data%ref_cluster_array)
    deallocate(calc_data%kkr_array)
    deallocate(calc_data%densities_array)
    deallocate(calc_data%energies_array)
    deallocate(calc_data%madelung_sum_array)
    deallocate(calc_data%trunc_zone)
    deallocate(calc_data%clusters)

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
  integer function getAtomIndexOfLocal(calc_data, ilocal)
    type (CalculationData), intent(in) :: calc_data
    integer, intent(in) :: ilocal

    getAtomIndexOfLocal = calc_data%atom_ids(ilocal)
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
  !> Returns reference to 'reference cluster for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getRefCluster(calc_data, local_atom_index)
    implicit none
    type (RefCluster), pointer :: getRefCluster ! return value

    type (CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getRefCluster => calc_data%ref_cluster_array(local_atom_index)
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
  !> Returns reference to energy results for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getEnergies(calc_data, local_atom_index)
    implicit none
    type (EnergyResults), pointer :: getEnergies ! return value

    type (CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getEnergies => calc_data%energies_array(local_atom_index)
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

  !----------------------------------------------------------------------------
  !> Returns reference to truncation zone.
  function getTruncationZone(calc_data)
    implicit none
    type (TruncationZone), pointer :: getTruncationZone ! return value
    type (CalculationData), intent(in) :: calc_data

    getTruncationZone => calc_data%trunc_zone
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to cluster info (sparsity info).
  function getClusterInfo(calc_data)
    implicit none
    type (ClusterInfo), pointer :: getClusterInfo ! return value

    type (CalculationData), intent(in) :: calc_data

    getClusterInfo => calc_data%clusters
  end function

  !----------------------------------------------------------------------------
  !> Returns reference to lattice vector table.
  function getLatticeVectors(calc_data)
    implicit none
    type (LatticeVectors), pointer :: getLatticeVectors ! return value

    type (CalculationData), intent(in) :: calc_data

    getLatticeVectors => calc_data%lattice_vectors
  end function
! ==================== Helper routines ========================================

  !----------------------------------------------------------------------------
  !> Helper routine: called by createCalculationData.
  subroutine constructEverything(calc_data, dims, params, arrays, my_mpi)
    use KKRnanoParallel_mod
    use DimParams_mod
    use InputParams_mod
    use Main2Arrays_mod
    use TEST_lcutoff_mod
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
    type (EnergyResults), pointer :: energies
    type (BasisAtom), pointer :: atomdata
    type (CellData), pointer :: cell
    type (RadialMeshData), pointer :: mesh
    type (LDAUData), pointer :: ldau_data
    type (JijData), pointer :: jij_data
    type (BroydenData), pointer :: broyden
    type (MadelungLatticeSum), pointer :: madelung_sum

    call createLatticeVectors(calc_data%lattice_vectors, arrays%bravais)

    ! create cluster for each local atom
    !$omp parallel do private(ilocal)
    do ilocal = 1, calc_data%num_local_atoms
      call createRefCluster(calc_data%ref_cluster_array(ilocal), &
                            calc_data%lattice_vectors, arrays%rbasis, &
                            params%rclust, calc_data%atom_ids(ilocal))
      !write(*,*) "Atoms in ref. cluster: ", calc_data%ref_cluster_array(ilocal)%nacls
    end do
    !$omp end parallel do

    ! setup the truncation zone
    call initLcutoffNew(calc_data%trunc_zone, calc_data%atom_ids, arrays)


    call createClusterInfo_com(calc_data%clusters, calc_data%ref_cluster_array, &
                          calc_data%trunc_zone, getMySEcommunicator(my_mpi))

    if (isMasterRank(my_mpi)) then
      write(*,*) "Number of lattice vectors created     : ", &
                  calc_data%lattice_vectors%nrd

      write(*,*) "Max. number of reference cluster atoms: ", &
                  calc_data%clusters%naclsd

      write(*,*) "On node 0: "
      write(*,*) "Num. atoms treated with full lmax: ", num_untruncated
      write(*,*) "Num. atoms in truncation zone 1  : ", num_truncated
      write(*,*) "Num. atoms in truncation zone 2  : ", num_truncated2
    end if
    CHECKASSERT(num_truncated+num_untruncated+num_truncated2 == dims%naez)

    call createMadelungCalculator(calc_data%madelung_calc, dims%lmaxd, &
                                  params%ALAT, params%RMAX, params%GMAX, &
                                  arrays%BRAVAIS, dims%NMAXD, dims%ISHLD)

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
      energies  => calc_data%energies_array(ilocal)

      call createKKRresults(kkr, dims, calc_data%clusters%naclsd)
      call createDensityResults(densities, dims)
      call createEnergyResults(energies, dims%nspind, dims%lmaxd)

      call createBasisAtom(atomdata, I1, dims%lpot, &
                           dims%nspind, dims%irmind, dims%irmd)

      call openBasisAtomDAFile(atomdata, 37, "atoms")
      call readBasisAtomDA(atomdata, 37, I1)
      call closeBasisAtomDAFile(37)

      ASSERT( atomdata%atom_index == I1 )

      if (isInMasterGroup(my_mpi)) then
        call openBasisAtomPotentialDAFile(atomdata, 37, "vpotnew")
        call readBasisAtomPotentialDA(atomdata, 37, I1)
        call closeBasisAtomPotentialDAFile(37)
      end if

      call createCellData(cell, dims%irid, (2*dims%LPOT+1)**2, dims%nfund)
      cell%cell_index = atomdata%cell_index

      call associateBasisAtomCell(atomdata, cell)

      call createRadialMeshData(mesh, dims%irmd, dims%ipand)
      call openRadialMeshDataDAFile(mesh, 37 , "meshes")
      call readRadialMeshDataDA(mesh, 37, I1)
      call closeRadialMeshDataDAFile(37)

      call associateBasisAtomMesh(atomdata, mesh)

      call createLDAUData(ldau_data, params%ldau, dims%irmd, dims%lmaxd, &
                          dims%nspind)
      call createJijData(jij_data, params%jij, params%rcutjij, dims%nxijd, &
                         dims%lmmaxd,dims%nspind)

      call createBroydenData(broyden, dims%ntird, dims%itdbryd, &
                             params%imix, params%mixing)

      call createMadelungLatticeSum(madelung_sum, calc_data%madelung_calc, dims%naez)

      ASSERT( arrays%ZAT(I1) == atomdata%Z_nuclear )

    !--------------------------------------------------------------------------
    end do
    !--------------------------------------------------------------------------

    ! on-the-fly shapefunction generation
    call generateShapes(calc_data, dims, params, arrays)

    ! calculate Gaunt coefficients
    call createGauntCoefficients(calc_data%gaunts, dims%lmaxd)
    call createShapeGauntCoefficients(calc_data%shgaunts, dims%lmaxd)

  end subroutine

!------------------------------------------------------------------------------
  subroutine generateShapes(calc_data, dims, params, arrays)
    use KKRnanoParallel_mod
    use DimParams_mod
    use InputParams_mod
    use Main2Arrays_mod
    use ConstructShapes_mod
    use ShapefunData_mod
    implicit none

    type (CalculationData), intent(inout) :: calc_data
    type (DimParams), intent(in)  :: dims
    type (InputParams), intent(in):: params
    type (Main2Arrays), intent(in):: arrays

    !-----------------
    integer :: I1, ilocal, nfun, ii, irmd, irid
    type (InterstitialMesh) :: inter_mesh
    type (ShapefunData) :: shdata
    double precision :: new_MT_radius

    ! loop over all LOCAL atoms
    !--------------------------------------------------------------------------
    do ilocal = 1, calc_data%num_local_atoms
    !--------------------------------------------------------------------------
      I1 = calc_data%atom_ids(ilocal)

      new_MT_radius = calc_data%atomdata_array(ilocal)%RMTref / params%alat

      call construct(shdata, inter_mesh, arrays%rbasis, arrays%bravais, I1, &
                     params%rclust, 4*dims%lmaxd, dims%irid-10, 10, new_MT_radius)


      ! first test it
      nfun = calc_data%cell_array(ilocal)%shdata%nfu
      irmd = dims%irmd
      irid = dims%irid
      CHECKASSERT(irid == shdata%irid)
      CHECKASSERT(dims%nfund == shdata%nfund)

      write(*,*) "Diff xrn:   ", sum(abs(inter_mesh%xrn * params%alat - calc_data%mesh_array(ilocal)%r(irmd-irid+1:irmd)))
      write(*,*) "Diff drn:   ", sum(abs(inter_mesh%drn * params%alat - calc_data%mesh_array(ilocal)%drdi(irmd-irid+1:irmd)))

      ! then use it
      calc_data%cell_array(ilocal)%shdata = shdata ! possible in Fortran 2003

      call destroyShapefunData(shdata)
      call destroyInterstitialMesh(inter_mesh)

    end do

  end subroutine

end module CalculationData_mod
