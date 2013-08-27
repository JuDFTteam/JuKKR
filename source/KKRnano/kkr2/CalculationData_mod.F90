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
  private :: generateAtomsShapesMeshes
  private :: generateShapesTEST
  private :: recordLengths_com
  private :: writePotentialIndexFile

  type CalculationData
    PRIVATE
    !integer :: atoms_per_proc
    integer :: num_local_atoms  ! <= atoms_per_proc
    integer, allocatable :: atom_ids(:)
    integer :: max_reclen_meshes
    integer :: max_reclen_potential

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

  !----------------------------------------------------------------------------
  !> Returns record length needed for 'meshes' file.
  integer function getMaxReclenMeshes(calc_data)
    implicit none
    type (CalculationData), intent(in) :: calc_data

    getMaxReclenMeshes = calc_data%max_reclen_meshes
  end function

  !----------------------------------------------------------------------------
  !> Returns record length needed for 'meshes' file.
  integer function getMaxReclenPotential(calc_data)
    implicit none
    type (CalculationData), intent(in) :: calc_data

    getMaxReclenPotential = calc_data%max_reclen_potential
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
    type (LDAUData), pointer :: ldau_data
    type (JijData), pointer :: jij_data
    type (BroydenData), pointer :: broyden
    type (MadelungLatticeSum), pointer :: madelung_sum
    type (RadialMeshData), pointer :: mesh


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

    ! a very crucial routine
    call generateAtomsShapesMeshes(calc_data, dims, params, arrays)

    call recordLengths_com(calc_data, my_mpi)

    if (isInMasterGroup(my_mpi)) then
      call writePotentialIndexFile(calc_data)
      call writeNewMeshFiles(calc_data)
    end if

    ! loop over all LOCAL atoms
    !--------------------------------------------------------------------------
    do ilocal = 1, calc_data%num_local_atoms
    !--------------------------------------------------------------------------

      I1 = calc_data%atom_ids(ilocal)

      kkr       => calc_data%kkr_array(ilocal)
      densities => calc_data%densities_array(ilocal)
      ldau_data => calc_data%ldau_data_array(ilocal)
      jij_data  => calc_data%jij_data_array(ilocal)
      broyden   => calc_data%broyden_array(ilocal)
      madelung_sum   => calc_data%madelung_sum_array(ilocal)
      energies  => calc_data%energies_array(ilocal)
      mesh => calc_data%mesh_array(ilocal)

      call createKKRresults(kkr, dims, calc_data%clusters%naclsd)
      call createDensityResults(densities, dims, mesh%irmd)
      call createEnergyResults(energies, dims%nspind, dims%lmaxd)

      call createLDAUData(ldau_data, params%ldau, mesh%irmd, dims%lmaxd, &
                          dims%nspind)
      call createJijData(jij_data, params%jij, params%rcutjij, dims%nxijd, &
                         dims%lmmaxd,dims%nspind)

      call createBroydenData(broyden, &
      (mesh%IRMD+(mesh%IRNS+1)*(dims%LMPOTD-1))*dims%NSPIND, &  ! former NTIRD
      dims%itdbryd, params%imix, params%mixing)

      call createMadelungLatticeSum(madelung_sum, calc_data%madelung_calc, dims%naez)

      !ASSERT( arrays%ZAT(I1) == atomdata%Z_nuclear )

    !--------------------------------------------------------------------------
    end do
    !--------------------------------------------------------------------------

    ! calculate Gaunt coefficients
    call createGauntCoefficients(calc_data%gaunts, dims%lmaxd)
    call createShapeGauntCoefficients(calc_data%shgaunts, dims%lmaxd)

  end subroutine

!------------------------------------------------------------------------------
!> Generates basis atom information, radial mesh, shape-function and
!> interpolates starting potential if necessary
  subroutine generateAtomsShapesMeshes(calc_data, dims, params, arrays)
    use DimParams_mod
    use InterpolateBasisAtom_mod
    use RadialMeshData_mod
    use InputParams_mod
    use Main2Arrays_mod

    implicit none
    type (CalculationData), intent(inout) :: calc_data
    type (DimParams), intent(in)  :: dims
    type (InputParams), intent(in):: params
    type (Main2Arrays), intent(in):: arrays

    integer ilocal, I1
    type (BasisAtom), pointer :: atomdata
    type (CellData), pointer :: cell
    type (RadialMeshData), pointer :: mesh
    type (BasisAtom), pointer :: old_atom
    type (RadialMeshData), pointer :: old_mesh

    type (BasisAtom), pointer :: old_atom_array(:)
    type (RadialMeshData), pointer :: old_mesh_array(:)
    double precision, allocatable :: new_MT_radii(:)

    allocate(old_atom_array(calc_data%num_local_atoms))
    allocate(old_mesh_array(calc_data%num_local_atoms))
    allocate(new_MT_radii(calc_data%num_local_atoms))

    ! generate storage for cell information + shape-functions
    do ilocal = 1, calc_data%num_local_atoms
      cell => calc_data%cell_array(ilocal)
      call createCellData(cell, dims%irid, (2*dims%LPOT+1)**2, dims%nfund)
    end do

    ! loop over all LOCAL atoms
    !--------------------------------------------------------------------------
    do ilocal = 1, calc_data%num_local_atoms
    !--------------------------------------------------------------------------

      I1 = calc_data%atom_ids(ilocal)

      ! We want to allow the actual radial mesh to be different from the one
      ! given by the input
      ! Therefore read 'potential' and 'meshes' data into temporary data-
      ! structures
      ! Then interpolate potential to the new mesh
      old_atom  => old_atom_array(ilocal)
      old_mesh  => old_mesh_array(ilocal)

      ! load the input data
      call createBasisAtomFromFile(old_atom, "atoms", "vpotnew.0", I1)

      call createRadialMeshDataFromFile(old_mesh, "meshes.0", I1)

      call associateBasisAtomMesh(old_atom, old_mesh)

      new_MT_radii(ilocal) = old_atom%RMTref / params%alat
    !--------------------------------------------------------------------------
    end do
    !--------------------------------------------------------------------------

    ! generate shapes and meshes
    call generateShapesTEST(calc_data, dims, params, arrays, &
                            new_MT_radii, params%MT_scale)

    ! interpolate to new mesh

    !--------------------------------------------------------------------------
    do ilocal = 1, calc_data%num_local_atoms
    !--------------------------------------------------------------------------
      I1 = calc_data%atom_ids(ilocal)
      atomdata  => calc_data%atomdata_array(ilocal)
      cell      => calc_data%cell_array(ilocal)
      mesh      => calc_data%mesh_array(ilocal)
      old_atom  => old_atom_array(ilocal)
      old_mesh  => old_mesh_array(ilocal)

      ! Geometry might have changed - interpolate to new mesh
      call interpolateBasisAtom(atomdata, old_atom, mesh)

      ! set new MT radius
      atomdata%RMTref = mesh%rmt

      cell%cell_index = atomdata%cell_index
      call associateBasisAtomCell(atomdata, cell)

      CHECKASSERT(dims%IRMIND == mesh%IRMIN) !check mesh
      CHECKASSERT( atomdata%atom_index == I1 )

      call destroyBasisAtom(old_atom)
      call destroyRadialMeshData(old_mesh)
    end do

    deallocate(new_MT_radii)
    deallocate(old_atom_array)
    deallocate(old_mesh_array)

  end subroutine

!------------------------------------------------------------------------------
  subroutine generateShapesTEST(calc_data, dims, params, arrays, &
                                new_MT_radii, MT_scale)
    use KKRnanoParallel_mod
    use DimParams_mod
    use InputParams_mod
    use Main2Arrays_mod
    use ConstructShapes_mod, only: construct, InterstitialMesh, &
                                   destroyInterstitialMesh
    use ShapefunData_mod
    implicit none

    type (CalculationData), intent(inout) :: calc_data
    type (DimParams), intent(in)  :: dims
    type (InputParams), intent(in):: params
    type (Main2Arrays), intent(in):: arrays
    double precision, intent(in) :: new_MT_radii(:)
    double precision, intent(in) :: MT_scale
    !-----------------

    integer :: I1, ilocal, nfun, ii
    integer :: irmd, irid, ipand, irnsd
    type (InterstitialMesh) :: inter_mesh
    type (ShapefunData) :: shdata
    double precision :: new_MT_radius
    integer :: num_MT_points
    type (RadialMeshData), pointer :: mesh

    ! loop over all LOCAL atoms
    !--------------------------------------------------------------------------
    do ilocal = 1, calc_data%num_local_atoms
    !--------------------------------------------------------------------------
      I1 = calc_data%atom_ids(ilocal)
      mesh => calc_data%mesh_array(ilocal)

      new_MT_radius = new_MT_radii(ilocal)
      num_MT_points = params%num_MT_points

      call construct(shdata, inter_mesh, arrays%rbasis, arrays%bravais, I1, &
                     params%rclust_voronoi, 4*dims%lmaxd, &
                     dims%irid - num_MT_points, &
                     params%nmin_panel, num_MT_points, new_MT_radius, MT_scale)


      ! use it
      calc_data%cell_array(ilocal)%shdata = shdata ! possible in Fortran 2003

      irmd = dims%irmd - dims%irid + size(inter_mesh%xrn)
      irid = size(inter_mesh%xrn)
      ipand = size(inter_mesh%nm) + 1
      irnsd = irmd - (dims%irmd - dims%irnsd)

      ASSERT(inter_mesh%xrn(1) /= 0.0d0)
      !write(*,*) irmd, irid, ipand, irnsd

      call createRadialMeshData(mesh, irmd, ipand)

      call initRadialMesh(mesh, params%alat, inter_mesh%xrn, &
                          inter_mesh%drn, inter_mesh%nm, irmd - irid, irnsd)

      call destroyShapefunData(shdata)
      call destroyInterstitialMesh(inter_mesh)

    end do

  end subroutine

  !----------------------------------------------------------------------------
  !> Communicate and set record lengths.
  subroutine recordLengths_com(calc_data, my_mpi)
    use KKRnanoParallel_mod
    use RadialMeshData_mod
    use BasisAtom_mod
    implicit none
    include 'mpif.h'

    type (CalculationData), intent(inout) :: calc_data
    type (KKRnanoParallel), intent(in) :: my_mpi

    integer :: ierr
    integer :: ilocal
    integer :: sendbuf(2)
    integer :: recvbuf(2)
    type (RadialMeshData), pointer :: mesh
    type (BasisAtom), pointer :: atomdata

    sendbuf = -1
    recvbuf = -1
    do ilocal = 1, calc_data%num_local_atoms
      atomdata  => calc_data%atomdata_array(ilocal)
      mesh      => calc_data%mesh_array(ilocal)

      sendbuf(1) = max(sendbuf(1), getMinReclenBasisAtomPotential(atomdata))
      sendbuf(2) = max(sendbuf(2), getMinReclenMesh(mesh))
    end do

    call MPI_Allreduce(sendbuf, recvbuf, 2, MPI_INTEGER, &
                       MPI_MAX, getMySECommunicator(my_mpi), ierr)

    if (getMyAtomRank(my_mpi) == 0) then
      write(*,*) "Record length 'vpotnew' file: ", recvbuf(1)
      write(*,*) "Record length 'meshes'  file: ", recvbuf(2)
    end if

    calc_data%max_reclen_potential = recvbuf(1)
    calc_data%max_reclen_meshes = recvbuf(2)

    ! debug
!    atomdata  => calc_data%atomdata_array(1)
!    call openBasisAtomPotentialDAFile(atomdata, 37, "vpotnew.i", &
!                                    getMaxReclenPotential(calc_data))
!    call writeBasisAtomPotentialDA(atomdata, 37, calc_data%atom_ids(1))
!    call closeBasisAtomPotentialDAFile(37)
    ! end debug

  end subroutine

  !----------------------------------------------------------------------------
  !> Write potential index file.
  subroutine writePotentialIndexFile(calc_data)
    use BasisAtom_mod
    implicit none

    type (CalculationData), intent(in) :: calc_data

    type (BasisAtom), pointer :: atomdata
    integer :: ilocal
    integer :: I1, max_reclen

    max_reclen = getMaxReclenPotential(calc_data)

    atomdata  => calc_data%atomdata_array(1)
    call openBasisAtomPotentialIndexDAFile(atomdata, 37, 'vpotnew.idx')

    do ilocal = 1, calc_data%num_local_atoms
      atomdata  => calc_data%atomdata_array(ilocal)
      I1 = calc_data%atom_ids(ilocal)
      call writeBasisAtomPotentialIndexDA(atomdata, 37, I1, max_reclen)
    end do

    call closeBasisAtomPotentialIndexDAFile(37)
  end subroutine

  !----------------------------------------------------------------------------
  !> Write new mesh files.
  !>
  !> The mesh can deviate from the input mesh if atoms are not in an ideal
  !> position. Therefore new mesh files have to be written.
  subroutine writeNewMeshFiles(calc_data)
    use BasisAtom_mod
    implicit none
    type (CalculationData), intent(in) :: calc_data

    type (RadialMeshData), pointer :: mesh
    integer :: ilocal
    integer :: I1, max_reclen

    max_reclen = getMaxReclenMeshes(calc_data)

    mesh      => calc_data%mesh_array(1)
    call openRadialMeshDataIndexDAFile(mesh, 37, 'meshes.idx')
    call openRadialMeshDataDAFile(mesh, 38, 'meshes', max_reclen)

    do ilocal = 1, calc_data%num_local_atoms
      mesh      => calc_data%mesh_array(ilocal)
      I1 = calc_data%atom_ids(ilocal)
      call writeRadialMeshDataIndexDA(mesh, 37, I1, max_reclen)
      call writeRadialMeshDataDA(mesh, 38, I1)
    end do

    call closeRadialMeshDataDAFile(38)
    call closeRadialMeshDataIndexDAFile(37)


  end subroutine

end module CalculationData_mod
