!> @author Elias Rabel

#include "DebugHelpers/test_macros.h"

module CalculationData_mod

  use MadelungCalculator_mod, only: MadelungCalculator, create, destroy
  use MadelungCalculator_mod, only: MadelungLatticeSum, create, destroy
  use RadialMeshData_mod, only: RadialMeshData, create, destroy
  use CellData_mod, only: CellData, create, destroy
  use BasisAtom_mod, only: BasisAtom, create, destroy
  use KKRresults_mod, only: KKRresults, create, destroy
  use DensityResults_mod, only: DensityResults, create, destroy
  use JijData_mod, only: JijData, create, destroy
  use LDAUData_mod, only: LDAUData, create, destroy
  use BroydenData_mod, only: BroydenData, create, destroy
  use EnergyResults_mod, only: EnergyResults, create, destroy
  use RefCluster_mod, only: RefCluster, create, destroy
  use ClusterInfo_mod, only: ClusterInfo, create, destroy
  use GauntCoefficients_mod, only: GauntCoefficients, create, destroy
  use ShapeGauntCoefficients_mod, only: ShapeGauntCoefficients, create, destroy
  use TruncationZone_mod, only: TruncationZone, create, destroy
  use InitialGuess_mod, only: InitialGuess, create, destroy
  use RefCluster_mod, only: LatticeVectors, create, destroy
  
  implicit none
  private

  public :: CalculationData, create, destroy
  public :: createCalculationData, destroyCalculationData ! deprecated
  public :: getBroydenDim, getBroyden, getNumLocalAtoms, getAtomIndexOfLocal, getAtomData, getRefCluster, getKKR
  public :: getMadelungSum, getDensities, getEnergies, getLDAUData, getJijData, getGaunts, getShapeGaunts
  public :: getMadelungCalculator, getTruncationZone, getClusterInfo, getLatticeVectors, getInitialGuessData 
  public :: getMaxReclenMeshes, getMaxReclenPotential      
  public :: prepareMadelung, constructEverything, setup_iguess, generateAtomsShapesMeshes, generateShapesTEST         
  public :: recordLengths_com, writePotentialIndexFile, writeNewMeshFiles, print_debug_info, constructClusters          
  public :: constructTruncationZones, constructStorage           
  

  type CalculationData
    private

    integer :: num_local_atoms  ! <= atoms_per_proc
    integer, allocatable :: atom_ids(:)
    integer :: max_reclen_meshes
    integer :: max_reclen_potential

    ! atom local data - different for each atom
    type(RadialMeshData), pointer     :: mesh_array(:)         => null()
    type(CellData), pointer           :: cell_array(:)         => null()
    type(BasisAtom), pointer          :: atomdata_array(:)     => null()
    type(RefCluster), pointer         :: ref_cluster_array(:)  => null()
    type(KKRresults), pointer         :: kkr_array(:)          => null()
    type(DensityResults), pointer     :: densities_array(:)    => null()
    type(EnergyResults), pointer      :: energies_array(:)     => null()
    type(MadelungLatticeSum), pointer :: madelung_sum_array(:) => null()
    type(LDAUData), pointer           :: ldau_data_array(:)    => null()
    type(JijData), pointer            :: jij_data_array(:)     => null()

    ! global data - same for each local atom
    type(LatticeVectors), pointer         :: lattice_vectors   => null()
    type(MadelungCalculator), pointer     :: madelung_calc     => null()
    type(ShapeGauntCoefficients), pointer :: shgaunts          => null()
    type(GauntCoefficients), pointer      :: gaunts            => null()
    type(TruncationZone), pointer         :: trunc_zone        => null()
    type(ClusterInfo), pointer            :: clusters          => null()
    type(BroydenData), pointer            :: broyden           => null()

    ! storage for initial guess
    type(InitialGuess), pointer           :: iguess_data       => null()

  endtype
  
  interface create
    module procedure createCalculationData
  endinterface
  
  interface destroy
    module procedure destroyCalculationData
  endinterface

  contains

  !----------------------------------------------------------------------------
  ! TODO: atoms_per_procs * num_procs MUST BE = naez
  ! rank = 0,1,..., num_atom_ranks-1
  subroutine createCalculationData(calc_data, dims, params, arrays, my_mpi)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    use KKRnanoParallel_mod, only: getMyAtomRank, getNumAtomRanks
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays

    type(CalculationData), intent(inout) :: calc_data
    type(DimParams), intent(in)   :: dims
    type(InputParams), intent(in) :: params
    type(Main2Arrays), intent(in) :: arrays
    type(KKRnanoParallel), intent(in) :: my_mpi

    integer, external :: mapblock
    integer :: atoms_per_proc, num_local_atoms, ii, atom_rank

    atoms_per_proc = dims%naez / getNumAtomRanks(my_mpi)
    ASSERT( getNumAtomRanks(my_mpi) * atoms_per_proc == dims%naez )
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

    allocate(calc_data%atom_ids(num_local_atoms))

    ! These datastructures are the same for all (local) atoms
    allocate(calc_data%lattice_vectors)
    allocate(calc_data%madelung_calc)
    allocate(calc_data%gaunts)
    allocate(calc_data%shgaunts)
    allocate(calc_data%trunc_zone)
    allocate(calc_data%clusters)
    allocate(calc_data%broyden)
    allocate(calc_data%iguess_data)

    atom_rank = getMyAtomRank(my_mpi)

    ! assign atom ids to processes with atom rank 'atom_rank'
    ! E.g. for 2 atoms per proc:
    ! process 1 treats atoms 1,2
    ! process 2 treats atoms 3,4 and so on
    ! FOR USE OF TRUNCATION THESE atoms have to be close together!!!

    ASSERT( size(calc_data%atom_ids) == num_local_atoms )
    ASSERT( mod(dims%naez, atoms_per_proc) == 0 )

    do ii = 1, num_local_atoms
      calc_data%atom_ids(ii) = atom_rank * atoms_per_proc + ii
      ASSERT( calc_data%atom_ids(ii) <= dims%naez )
    enddo ! ii

    ! Now construct all datastructures and calculate initial data
    call constructEverything(calc_data, dims, params, arrays, my_mpi)

    !call print_debug_info(calc_data)
  endsubroutine

  !----------------------------------------------------------------------------
  !> Calculate Madelung Lattice sums for all local atoms.
  subroutine prepareMadelung(calc_data, arrays)
    use Main2Arrays_mod, only: Main2Arrays
    use MadelungCalculator_mod, only: calculateMadelungLatticeSum

    type(CalculationData), intent(inout) :: calc_data
    type(Main2Arrays), intent(in):: arrays

    ! ----- locals ------
    integer :: I1
    integer :: ilocal
    type(MadelungLatticeSum), pointer :: madelung_sum

    do ilocal = 1, calc_data%num_local_atoms
      I1 = calc_data%atom_ids(ilocal)
      madelung_sum => calc_data%madelung_sum_array(ilocal)
      call calculateMadelungLatticeSum(madelung_sum, I1, arrays%rbasis)
    enddo

  endsubroutine

  !----------------------------------------------------------------------------
  subroutine destroyCalculationData(calc_data) ! todo: text-replace calc_data by self in this routine
    ! deprecated interfaces
    use RefCluster_mod, only: destroyRefCluster
    use MadelungCalculator_mod, only: destroyMadelungLatticeSum
    use BasisAtom_mod, only: destroyBasisAtom
    use CellData_mod, only: destroyCellData
    use RadialMeshData_mod, only: destroyRadialMeshData
    use LDAUData_mod, only: destroyLDAUData
    use JijData_mod, only: destroyJijData
    use DensityResults_mod, only: destroyDensityResults
    use KKRresults_mod, only: destroyKKRresults
    use EnergyResults_mod, only: destroyEnergyResults
    use RefCluster_mod, only: destroyLatticeVectors
    use MadelungCalculator_mod, only: destroyMadelungCalculator
    use GauntCoefficients_mod, only: destroyGauntCoefficients
    use ShapeGauntCoefficients_mod, only: destroyShapeGauntCoefficients
    use BroydenData_mod, only: destroyBroydenData
    use InitialGuess_mod, only: iguess_destroy ! deprecated exceptional naming!
    
    type(CalculationData), intent(inout) :: calc_data

    ! ---- locals ----
    integer :: ilocal
    integer :: I1
    type(KKRresults), pointer :: kkr
    type(DensityResults), pointer :: densities
    type(EnergyResults), pointer :: energies
    type(BasisAtom), pointer :: atomdata
    type(CellData), pointer :: cell
    type(RadialMeshData), pointer :: mesh
    type(LDAUData), pointer :: ldau_data
    type(JijData), pointer :: jij_data
    type(MadelungLatticeSum), pointer :: madelung_sum

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
      madelung_sum   => calc_data%madelung_sum_array(ilocal)
      energies  => calc_data%energies_array(ilocal)

      call destroyRefCluster(calc_data%ref_cluster_array(ilocal))

      call destroyMadelungLatticeSum(madelung_sum)

      call destroyBasisAtom(atomdata)
      call destroyCellData(cell)
      call destroyRadialMeshData(mesh)

      call destroyLDAUData(ldau_data)
      call destroyJijData(jij_data)
      call destroyDensityResults(densities)
      call destroyKKRresults(kkr)
      call destroyEnergyResults(energies)

    enddo ! ilocal

    call destroyLatticeVectors(calc_data%lattice_vectors)
    call destroyMadelungCalculator(calc_data%madelung_calc)
    call destroyGauntCoefficients(calc_data%gaunts)
    call destroyShapeGauntCoefficients(calc_data%shgaunts)
    call destroyBroydenData(calc_data%broyden)
    call iguess_destroy(calc_data%iguess_data)

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
    deallocate(calc_data%iguess_data)

    deallocate(calc_data%ldau_data_array)
    deallocate(calc_data%jij_data_array)
    deallocate(calc_data%broyden)
    deallocate(calc_data%atom_ids)
  endsubroutine

  !----------------------------------------------------------------------------
  integer function getNumLocalAtoms(calc_data)
    type(CalculationData), intent(in) :: calc_data

    getNumLocalAtoms = calc_data%num_local_atoms
  endfunction

  !----------------------------------------------------------------------------
  integer function getAtomIndexOfLocal(calc_data, ilocal)
    type(CalculationData), intent(in) :: calc_data
    integer, intent(in) :: ilocal

    getAtomIndexOfLocal = calc_data%atom_ids(ilocal)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to atomdata for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getAtomData(calc_data, local_atom_index)
    type(BasisAtom), pointer :: getAtomData ! return value

    type(CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getAtomData => calc_data%atomdata_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to 'reference cluster for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getRefCluster(calc_data, local_atom_index)
    type(RefCluster), pointer :: getRefCluster ! return value

    type(CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getRefCluster => calc_data%ref_cluster_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to kkr(results) for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getKKR(calc_data, local_atom_index)
    type(KKRresults), pointer :: getKKR ! return value

    type(CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getKKR => calc_data%kkr_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Madelung sum for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getMadelungSum(calc_data, local_atom_index)
        
    type(MadelungLatticeSum), pointer :: getMadelungSum ! return value

    type(CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getMadelungSum => calc_data%madelung_sum_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to density results for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getDensities(calc_data, local_atom_index)
    type(DensityResults), pointer :: getDensities ! return value

    type(CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getDensities => calc_data%densities_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to energy results for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getEnergies(calc_data, local_atom_index)
    type(EnergyResults), pointer :: getEnergies ! return value

    type(CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getEnergies => calc_data%energies_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to LDA+U data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getLDAUData(calc_data, local_atom_index)
    type(LDAUData), pointer :: getLDAUData ! return value

    type(CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getLDAUData => calc_data%ldau_data_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Jij-data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getJijData(calc_data, local_atom_index)
    type(JijData), pointer :: getJijData ! return value

    type(CalculationData), intent(in) :: calc_data
    integer, intent(in) :: local_atom_index

    getJijData => calc_data%jij_data_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Broyden data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getBroyden(calc_data, local_atom_index) ! todo: remove local_atom_index from interface (although optional)
    type(BroydenData), pointer :: getBroyden ! return value

    type(CalculationData), intent(in) :: calc_data
    integer, intent(in), optional :: local_atom_index

    getBroyden => calc_data%broyden
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Gaunt coefficients.
  function getGaunts(calc_data)
    type(GauntCoefficients), pointer :: getGaunts ! return value

    type(CalculationData), intent(in) :: calc_data

    getGaunts => calc_data%gaunts
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Shape-Gaunt coefficients.
  function getShapeGaunts(calc_data)
    type(ShapeGauntCoefficients), pointer :: getShapeGaunts ! return value

    type(CalculationData), intent(in) :: calc_data

    getShapeGaunts => calc_data%shgaunts
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Madelung calculator.
  function getMadelungCalculator(calc_data)
    type(MadelungCalculator), pointer :: getMadelungCalculator ! return value

    type(CalculationData), intent(in) :: calc_data

    getMadelungCalculator => calc_data%madelung_calc
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to truncation zone.
  function getTruncationZone(calc_data)
    type(TruncationZone), pointer :: getTruncationZone ! return value
    type(CalculationData), intent(in) :: calc_data

    getTruncationZone => calc_data%trunc_zone
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to cluster info (sparsity info).
  function getClusterInfo(calc_data)
    type(ClusterInfo), pointer :: getClusterInfo ! return value

    type(CalculationData), intent(in) :: calc_data

    getClusterInfo => calc_data%clusters
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to lattice vector table.
  function getLatticeVectors(calc_data)
    type(LatticeVectors), pointer :: getLatticeVectors ! return value

    type(CalculationData), intent(in) :: calc_data

    getLatticeVectors => calc_data%lattice_vectors
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to initial guess data.
  function getInitialGuessData(calc_data)
    type(InitialGuess), pointer :: getInitialGuessData

    type(CalculationData), intent(in) :: calc_data

    getInitialGuessData => calc_data%iguess_data
  endfunction

  !----------------------------------------------------------------------------
  !> Returns record length needed for 'meshes' file.
  integer function getMaxReclenMeshes(calc_data)
    type(CalculationData), intent(in) :: calc_data

    getMaxReclenMeshes = calc_data%max_reclen_meshes
  endfunction

  !----------------------------------------------------------------------------
  !> Returns record length needed for 'meshes' file.
  integer function getMaxReclenPotential(calc_data)
    type(CalculationData), intent(in) :: calc_data

    getMaxReclenPotential = calc_data%max_reclen_potential
  endfunction

! ==================== Helper routines ========================================

  !----------------------------------------------------------------------------
  !> Helper routine: called by createCalculationData.
  subroutine constructEverything(calc_data, dims, params, arrays, my_mpi)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getMySEcommunicator, isMasterRank, isInMasterGroup
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use TEST_lcutoff_mod, only: num_untruncated, num_truncated, num_truncated2
    ! deprecated interfaces
    use RefCluster_mod, only: createLatticeVectors
    use RefCluster_mod, only: createRefCluster              
    use TEST_lcutoff_mod, only: initLcutoffNew                
    use ClusterInfo_mod, only: createClusterInfo_com         
    use MadelungCalculator_mod, only: createMadelungCalculator
    use DensityResults_mod, only: createDensityResults          
    use EnergyResults_mod, only: createEnergyResults           
    use LDAUData_mod, only: createLDAUData                
    use JijData_mod, only: createJijData                 
    use MadelungCalculator_mod, only: createMadelungLatticeSum      
    use GauntCoefficients_mod, only: createGauntCoefficients       
    use ShapeGauntCoefficients_mod, only: createShapeGauntCoefficients  
    use BroydenData_mod, only: createBroydenData             
    use KKRresults_mod, only: createKKRresults
    
    type(CalculationData), intent(inout) :: calc_data
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays
    type(KKRnanoParallel), intent(in) :: my_mpi

    ! ----- locals ------
    integer :: I1
    integer :: ilocal
    type(KKRresults), pointer :: kkr
    type(DensityResults), pointer :: densities
    type(EnergyResults), pointer :: energies
    type(LDAUData), pointer :: ldau_data
    type(JijData), pointer :: jij_data
    type(MadelungLatticeSum), pointer :: madelung_sum
    type(RadialMeshData), pointer :: mesh

    call createLatticeVectors(calc_data%lattice_vectors, arrays%bravais)

    ! create cluster for each local atom
    !$omp parallel do private(ilocal)
    do ilocal = 1, calc_data%num_local_atoms
      call createRefCluster(calc_data%ref_cluster_array(ilocal), &
                            calc_data%lattice_vectors, arrays%rbasis, &
                            params%rclust, calc_data%atom_ids(ilocal))
      !write(*,*) "Atoms in ref. cluster: ", calc_data%ref_cluster_array(ilocal)%nacls
    enddo
    !$omp endparallel do

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
    endif
    CHECKASSERT(num_truncated+num_untruncated+num_truncated2 == dims%naez)

    call createMadelungCalculator(calc_data%madelung_calc, dims%lmaxd, &
                                  params%ALAT, params%RMAX, params%GMAX, &
                                  arrays%BRAVAIS, dims%NMAXD, dims%ISHLD)

    ! a very crucial routine
    call generateAtomsShapesMeshes(calc_data, dims, params, arrays)

    call recordLengths_com(calc_data, my_mpi)

    if (isInMasterGroup(my_mpi)) then
#ifndef TASKLOCAL_FILES
      call writePotentialIndexFile(calc_data)
#endif
      call writeNewMeshFiles(calc_data)
    endif

    ! loop over all LOCAL atoms
    !--------------------------------------------------------------------------
    do ilocal = 1, calc_data%num_local_atoms
    !--------------------------------------------------------------------------

      I1 = calc_data%atom_ids(ilocal)

      kkr       => calc_data%kkr_array(ilocal)
      densities => calc_data%densities_array(ilocal)
      ldau_data => calc_data%ldau_data_array(ilocal)
      jij_data  => calc_data%jij_data_array(ilocal)
      madelung_sum   => calc_data%madelung_sum_array(ilocal)
      energies  => calc_data%energies_array(ilocal)
      mesh => calc_data%mesh_array(ilocal)

      call createKKRresults(kkr, dims, calc_data%clusters%naclsd)
      call createDensityResults(densities, dims, mesh%irmd)
      call createEnergyResults(energies, dims%nspind, dims%lmaxd)

      call createLDAUData(ldau_data, params%ldau, mesh%irmd, dims%lmaxd, dims%nspind)
      call createJijData(jij_data, .false., params%rcutjij, dims%nxijd, dims%lmmaxd,dims%nspind)

      call createMadelungLatticeSum(madelung_sum, calc_data%madelung_calc, dims%naez)

      !ASSERT( arrays%ZAT(I1) == atomdata%Z_nuclear )

    !--------------------------------------------------------------------------
    enddo ! ilocal
    !--------------------------------------------------------------------------

    ! calculate Gaunt coefficients
    call createGauntCoefficients(calc_data%gaunts, dims%lmaxd)
    call createShapeGauntCoefficients(calc_data%shgaunts, dims%lmaxd)
    call createBroydenData(calc_data%broyden, &
         getBroydenDim(calc_data), &  ! former NTIRD
         dims%itdbryd, params%imix, params%mixing)

    ! setup storage for iguess
    call setup_iguess(calc_data, dims, arrays)

  endsubroutine

  !----------------------------------------------------------------------------
  !> Initialise iguess datastructure.
  subroutine setup_iguess(calc_data, dims, arrays)
    use DimParams_mod, only: DimParams
    use Main2Arrays_mod, only: Main2Arrays
    use InitialGuess_mod, only: iguess_init ! deprecated exceptional naming

    type(CalculationData), intent(inout) :: calc_data
    type(DimParams), intent(in)  :: dims
    type(Main2Arrays), intent(in):: arrays

    integer, allocatable :: num_k_points(:)
    integer :: ii
    integer :: blocksize

    ! TODO: This is overdimensioned when l-cutoff is used!!!
    ! DO NOT USE IGUESS together with l-cutoff!!! RS-cutoff is fine
    blocksize = calc_data%trunc_zone%naez_trc * calc_data%num_local_atoms * dims%lmmaxd**2

    allocate(num_k_points(dims%iemxd))
    do ii = 1, dims%iemxd
      num_k_points(ii) = arrays%nofks(arrays%kmesh(ii))
    enddo ! ii

    ! setup storage for iguess
    if (dims%smpid == 1 .and. dims%nspind == 2) then
      ! no spin parallelisation choosen, processes must store both spin-directions
      call iguess_init(calc_data%iguess_data, num_k_points, 2, blocksize, dims%iguessd)
    else
      call iguess_init(calc_data%iguess_data, num_k_points, 1, blocksize, dims%iguessd)
    endif

  endsubroutine

!------------------------------------------------------------------------------
!> Generates basis atom information, radial mesh, shape-function and
!> interpolates starting potential if necessary
  subroutine generateAtomsShapesMeshes(calc_data, dims, params, arrays)
    use DimParams_mod, only: DimParams
    use InterpolateBasisAtom_mod, only: interpolateBasisAtom
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use BasisAtom_mod, only: destroyBasisAtom
    use RadialMeshData_mod, only: destroyRadialMeshData
    
    use CellData_mod, only: createCellData
    use BasisAtom_mod, only: createBasisAtomFromFile
    use RadialMeshData_mod, only: createRadialMeshDataFromFile
    use BasisAtom_mod, only: associateBasisAtomMesh, associateBasisAtomCell
    
    type(CalculationData), intent(inout) :: calc_data
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays

    integer ilocal, I1
    type(BasisAtom), pointer :: atomdata
    type(CellData), pointer :: cell
    type(RadialMeshData), pointer :: mesh
    type(BasisAtom), pointer :: old_atom
    type(RadialMeshData), pointer :: old_mesh

    type(BasisAtom), pointer :: old_atom_array(:)
    type(RadialMeshData), pointer :: old_mesh_array(:)
    double precision, allocatable :: new_MT_radii(:)

    allocate(old_atom_array(calc_data%num_local_atoms))
    allocate(old_mesh_array(calc_data%num_local_atoms))
    allocate(new_MT_radii(calc_data%num_local_atoms))

    ! generate storage for cell information + shape-functions
    do ilocal = 1, calc_data%num_local_atoms
      cell => calc_data%cell_array(ilocal)
      call createCellData(cell, dims%irid, (2*dims%LPOT+1)**2, (2*dims%LPOT+1)**2)
    enddo

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

      new_MT_radii(ilocal) = old_atom%radius_muffin_tin / params%alat
    !--------------------------------------------------------------------------
    enddo ! ilocal
    !--------------------------------------------------------------------------

    ! generate shapes and meshes
    call generateShapesTEST(calc_data, dims, params, arrays, new_MT_radii, params%MT_scale)

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
      call interpolateBasisAtom(atomdata, old_atom, mesh, dims%lpot)

      ! set new MT radius
      atomdata%radius_muffin_tin = mesh%rmt

      ! set radius of repulsive reference potential
      if (params%RMT_ref_scale > 0.0d0) then
        atomdata%RMTref = cell%shdata%max_muffin_tin * params%alat * params%RMT_ref_scale
      else
        atomdata%RMTref = atomdata%radius_muffin_tin ! old behaviour=Mt-radius
      endif

      cell%cell_index = atomdata%cell_index
      call associateBasisAtomCell(atomdata, cell)

      CHECKASSERT(dims%IRMIND == mesh%IRMIN) !check mesh
      CHECKASSERT( atomdata%atom_index == I1 )

      call destroyBasisAtom(old_atom)
      call destroyRadialMeshData(old_mesh)
    enddo ! ilocal

    deallocate(new_MT_radii)
    deallocate(old_atom_array)
    deallocate(old_mesh_array)

  endsubroutine

!------------------------------------------------------------------------------
  subroutine generateShapesTEST(calc_data, dims, params, arrays, new_MT_radii, MT_scale)
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use ConstructShapes_mod, only: construct, InterstitialMesh, destroyInterstitialMesh, write_shapefun_file
    use ShapefunData_mod, only: ShapefunData
    use RadialMeshData_mod, only: createRadialMeshData, initRadialMesh
    use ShapefunData_mod, only: destroyShapefunData

    type(CalculationData), intent(inout) :: calc_data
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays
    double precision, intent(in) :: new_MT_radii(:)
    double precision, intent(in) :: MT_scale

    integer :: I1, ilocal
    integer :: irmd, irid, ipand, irnsd
    type(InterstitialMesh) :: inter_mesh
    type(ShapefunData) :: shdata ! temporary shape-fun data
    double precision :: new_MT_radius
    integer :: num_MT_points
    type(RadialMeshData), pointer :: mesh

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

      ! optional output of shape functions
      if (params%write_shapes == 1) call write_shapefun_file(shdata, inter_mesh, I1)

      call destroyShapefunData(shdata)
      call destroyInterstitialMesh(inter_mesh)

    enddo ! ilocal

  endsubroutine

  !----------------------------------------------------------------------------
  !> Communicate and set record lengths.
  subroutine recordLengths_com(calc_data, my_mpi)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getMySECommunicator, getMyAtomRank
    use RadialMeshData_mod, only: getMinReclenMesh
    use BasisAtom_mod, only: getMinReclenBasisAtomPotential
    include 'mpif.h'

    type(CalculationData), intent(inout) :: calc_data
    type(KKRnanoParallel), intent(in) :: my_mpi

    integer :: ierr
    integer :: ilocal
    integer :: sendbuf(2)
    integer :: recvbuf(2)
    type(RadialMeshData), pointer :: mesh
    type(BasisAtom), pointer :: atomdata

    sendbuf = -1
    recvbuf = -1
    do ilocal = 1, calc_data%num_local_atoms
      atomdata  => calc_data%atomdata_array(ilocal)
      mesh      => calc_data%mesh_array(ilocal)

      sendbuf(1) = max(sendbuf(1), getMinReclenBasisAtomPotential(atomdata))
      sendbuf(2) = max(sendbuf(2), getMinReclenMesh(mesh))
    enddo ! ilocal

    call MPI_Allreduce(sendbuf, recvbuf, 2, MPI_INTEGER, MPI_MAX, getMySECommunicator(my_mpi), ierr)

    if (getMyAtomRank(my_mpi) == 0) then
      write(*,*) "Record length 'vpotnew' file: ", recvbuf(1)
      write(*,*) "Record length 'meshes'  file: ", recvbuf(2)
    endif

    calc_data%max_reclen_potential = recvbuf(1)
    calc_data%max_reclen_meshes = recvbuf(2)

  endsubroutine

  !----------------------------------------------------------------------------
  !> Write potential index file.
  subroutine writePotentialIndexFile(calc_data)
    use BasisAtom_mod, only: openBasisAtomPotentialIndexDAFile, writeBasisAtomPotentialIndexDA, closeBasisAtomPotentialIndexDAFile
    type(CalculationData), intent(in) :: calc_data

    type(BasisAtom), pointer :: atomdata
    integer :: ilocal
    integer :: I1, max_reclen

    max_reclen = getMaxReclenPotential(calc_data)

    atomdata  => calc_data%atomdata_array(1)
    call openBasisAtomPotentialIndexDAFile(atomdata, 37, 'vpotnew.idx')

    do ilocal = 1, calc_data%num_local_atoms
      atomdata  => calc_data%atomdata_array(ilocal)
      I1 = calc_data%atom_ids(ilocal)
      call writeBasisAtomPotentialIndexDA(atomdata, 37, I1, max_reclen)
    enddo ! ilocal

    call closeBasisAtomPotentialIndexDAFile(37)
  endsubroutine

  !----------------------------------------------------------------------------
  !> Write new mesh files.
  !>
  !> The mesh can deviate from the input mesh if atoms are not in an ideal
  !> position. Therefore new mesh files have to be written.
  subroutine writeNewMeshFiles(calc_data)
#ifndef TASKLOCAL_FILES
    use RadialMeshData_mod, only: openRadialMeshDataIndexDAFile, writeRadialMeshDataIndexDA, closeRadialMeshDataIndexDAFile
#endif
    use RadialMeshData_mod, only: openRadialMeshDataDAFile, writeRadialMeshDataDA, closeRadialMeshDataDAFile
    
    type(CalculationData), intent(in) :: calc_data

    type(RadialMeshData), pointer :: mesh
    integer :: ilocal
    integer :: I1, max_reclen

    max_reclen = getMaxReclenMeshes(calc_data)

    mesh      => calc_data%mesh_array(1)
#ifndef TASKLOCAL_FILES
    ! don't write index when using task-local files
    call openRadialMeshDataIndexDAFile(mesh, 37, 'meshes.idx')
#endif
    call openRadialMeshDataDAFile(mesh, 38, 'meshes', max_reclen)

    do ilocal = 1, calc_data%num_local_atoms
      mesh      => calc_data%mesh_array(ilocal)
      I1 = calc_data%atom_ids(ilocal)
#ifndef TASKLOCAL_FILES
      call writeRadialMeshDataIndexDA(mesh, 37, I1, max_reclen)
#endif
      call writeRadialMeshDataDA(mesh, 38, I1)
    enddo ! ilocal

    call closeRadialMeshDataDAFile(38)
#ifndef TASKLOCAL_FILES
    call closeRadialMeshDataIndexDAFile(37)
#endif

  endsubroutine

!============ Helper routines for Broyden mixing ==============================

  !----------------------------------------------------------------------------
  !> Returns the number of potential values of ALL LOCAL atoms.
  !> This is needed for dimensioning the Broyden mixing work arrays.
  integer function getBroydenDim(calc_data)
    use PotentialData_mod, only: getNumPotentialValues
    
    type(CalculationData), intent(in) :: calc_data

    integer :: ilocal
    type(BasisAtom), pointer :: atomdata

    getBroydenDim = 0
    do ilocal = 1, calc_data%num_local_atoms
      atomdata  => calc_data%atomdata_array(ilocal)
      getBroydenDim = getBroydenDim + getNumPotentialValues(atomdata%potential)
    enddo ! ilocal
    
  endfunction

  !----------------------------------------------------------------------------
  ! Print some debugging info
  subroutine print_debug_info(calc_data)
    use RadialMeshData_mod, only: repr_RadialMeshData
    use ShapefunData_mod, only: repr_ShapefunData
    use PotentialData_mod, only: repr_PotentialData

    type(CalculationData), intent(in) :: calc_data

    integer :: ilocal
    type(RadialMeshData), pointer :: mesh
    type(BasisAtom), pointer :: atomdata
    type(CellData), pointer :: cell
    character(len=:), allocatable :: str

    do ilocal = 1, calc_data%num_local_atoms
      mesh     => calc_data%mesh_array(ilocal)
      atomdata => calc_data%atomdata_array(ilocal)
      cell     => calc_data%cell_array(ilocal)
      call repr_RadialMeshData(mesh, str)
      write(*, '(A)') str
      call repr_PotentialData(atomdata%potential, str)
      write(*, '(A)') str
      call repr_ShapefunData(cell%shdata, str)
      write(*, '(A)') str
    enddo ! ilocal
    
  endsubroutine

!==============================================================================
!=             WORK in PROGRESS - not used yet                                =
!==============================================================================

  ! Factored out some routines from 'constructEverything'

  subroutine constructClusters(calc_data, params, arrays)
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use RefCluster_mod, only: createRefCluster, createLatticeVectors
    
    type(CalculationData), intent(inout) :: calc_data
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays

    integer :: ilocal

    call createLatticeVectors(calc_data%lattice_vectors, arrays%bravais)

    ! create cluster for each local atom
    !$omp parallel do private(ilocal)
    do ilocal = 1, calc_data%num_local_atoms
      call createRefCluster(calc_data%ref_cluster_array(ilocal), &
                            calc_data%lattice_vectors, arrays%rbasis, &
                            params%rclust, calc_data%atom_ids(ilocal))
    enddo ! ilocal
    !$omp endparallel do

  endsubroutine

  subroutine constructTruncationZones(calc_data, dims, arrays, my_mpi)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getMySEcommunicator, isMasterRank   
    use DimParams_mod, only: DimParams
    use Main2Arrays_mod, only: Main2Arrays
    use TEST_lcutoff_mod, only: num_untruncated, num_truncated2, num_truncated ! integer
    use TEST_lcutoff_mod, only: initLcutoffNew
    use ClusterInfo_mod, only: createClusterInfo_com

    type(CalculationData), intent(inout) :: calc_data
    type(DimParams), intent(in)  :: dims
    type(Main2Arrays), intent(in):: arrays
    type(KKRnanoParallel), intent(in) :: my_mpi

    ! setup the truncation zone
    call initLcutoffNew(calc_data%trunc_zone, calc_data%atom_ids, arrays)

    ! get information about all the reference clusters of
    ! atoms in truncation zone
    call createClusterInfo_com(calc_data%clusters, calc_data%ref_cluster_array, &
                          calc_data%trunc_zone, getMySEcommunicator(my_mpi))

    if (isMasterRank(my_mpi)) then
      write(*,*) "Number of lattice vectors created     : ", calc_data%lattice_vectors%nrd
      write(*,*) "Max. number of reference cluster atoms: ", calc_data%clusters%naclsd
      write(*,*) "On node 0: "
      write(*,*) "Num. atoms treated with full lmax: ", num_untruncated
      write(*,*) "Num. atoms in truncation zone 1  : ", num_truncated
      write(*,*) "Num. atoms in truncation zone 2  : ", num_truncated2
    endif
    CHECKASSERT(num_truncated+num_untruncated+num_truncated2 == dims%naez)
    
  endsubroutine

  subroutine constructStorage(calc_data, dims, params, arrays, my_mpi) ! todo: remove arrays, my_mpi from interface
    use KKRnanoParallel_mod, only: KKRnanoParallel
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    
    use KKRresults_mod, only: createKKRresults
    use DensityResults_mod, only: createDensityResults
    use EnergyResults_mod, only: createEnergyResults
    use LDAUData_mod, only: createLDAUData
    use JijData_mod, only: createJijData
    use MadelungCalculator_mod, only: createMadelungLatticeSum
    
    type(CalculationData), intent(inout) :: calc_data
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays
    type(KKRnanoParallel), intent(in) :: my_mpi

    integer :: I1, ilocal
    type(KKRresults), pointer :: kkr
    type(DensityResults), pointer :: densities
    type(EnergyResults), pointer :: energies
    type(LDAUData), pointer :: ldau_data
    type(JijData), pointer :: jij_data
    type(MadelungLatticeSum), pointer :: madelung_sum
    type(RadialMeshData), pointer :: mesh

    do ilocal = 1, calc_data%num_local_atoms
      I1 = calc_data%atom_ids(ilocal)

      kkr       => calc_data%kkr_array(ilocal)
      densities => calc_data%densities_array(ilocal)
      ldau_data => calc_data%ldau_data_array(ilocal)
      jij_data  => calc_data%jij_data_array(ilocal)
      madelung_sum   => calc_data%madelung_sum_array(ilocal)
      energies  => calc_data%energies_array(ilocal)
      mesh => calc_data%mesh_array(ilocal)

      call createKKRresults(kkr, dims, calc_data%clusters%naclsd)
      call createDensityResults(densities, dims, mesh%irmd)
      call createEnergyResults(energies, dims%nspind, dims%lmaxd)

      call createLDAUData(ldau_data, params%ldau, mesh%irmd, dims%lmaxd, dims%nspind)
      call createJijData(jij_data, params%jij, params%rcutjij, dims%nxijd, dims%lmmaxd,dims%nspind)

      call createMadelungLatticeSum(madelung_sum, calc_data%madelung_calc, dims%naez)
    enddo ! ilocal
    
  endsubroutine

endmodule CalculationData_mod
