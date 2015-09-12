!> @author Elias Rabel

#include "DebugHelpers/test_macros.h"

module CalculationData_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)

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
  public :: recordLengths_com, writePotentialIndexFile, writeNewMeshFiles, repr_CalculationData
  
! public :: constructTruncationZones, constructStorage, constructClusters ! not yet used


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
  subroutine createCalculationData(self, dims, params, arrays, my_mpi)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    use KKRnanoParallel_mod, only: getMyAtomRank, getNumAtomRanks
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays

    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)   :: dims
    type(InputParams), intent(in) :: params
    type(Main2Arrays), intent(in) :: arrays
    type(KKRnanoParallel), intent(in) :: my_mpi

    integer :: atoms_per_proc, num_local_atoms, ii, atom_rank

    atoms_per_proc = dims%naez / getNumAtomRanks(my_mpi)
    ASSERT( getNumAtomRanks(my_mpi) * atoms_per_proc == dims%naez )
    num_local_atoms = atoms_per_proc !TODO

    self%num_local_atoms = num_local_atoms

    ! one datastructure for each local atom
    allocate(self%mesh_array(num_local_atoms))
    allocate(self%cell_array(num_local_atoms))
    allocate(self%atomdata_array(num_local_atoms))
    allocate(self%ref_cluster_array(num_local_atoms))
    allocate(self%kkr_array(num_local_atoms))
    allocate(self%densities_array(num_local_atoms))
    allocate(self%energies_array(num_local_atoms))
    allocate(self%madelung_sum_array(num_local_atoms))
    allocate(self%ldau_data_array(num_local_atoms))
    allocate(self%jij_data_array(num_local_atoms))

    ! These datastructures are the same for all (local) atoms
    allocate(self%lattice_vectors)
    allocate(self%madelung_calc)
    allocate(self%gaunts)
    allocate(self%shgaunts)
    allocate(self%trunc_zone)
    allocate(self%clusters)
    allocate(self%broyden)
    allocate(self%iguess_data)

    atom_rank = getMyAtomRank(my_mpi)
    
    allocate(self%atom_ids(num_local_atoms))

    ! assign atom ids to processes with atom rank 'atom_rank'
    ! E.g. for 2 atoms per proc:
    ! process 1 treats atoms 1,2
    ! process 2 treats atoms 3,4 and so on
    ! FOR USE OF TRUNCATION THESE atoms have to be close together!!!

    ASSERT( size(self%atom_ids) == num_local_atoms )
    ASSERT( mod(dims%naez, atoms_per_proc) == 0 )

    do ii = 1, num_local_atoms
      self%atom_ids(ii) = atom_rank * atoms_per_proc + ii
      ASSERT( self%atom_ids(ii) <= dims%naez )
    enddo ! ii

    ! Now construct all datastructures and calculate initial data
    call constructEverything(self, dims, params, arrays, my_mpi)

    ! call repr_CalculationData(self)
  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Calculate Madelung Lattice sums for all local atoms.
  subroutine prepareMadelung(self, arrays)
    use Main2Arrays_mod, only: Main2Arrays
    use MadelungCalculator_mod, only: calculateMadelungLatticeSum

    type(CalculationData), intent(inout) :: self
    type(Main2Arrays), intent(in):: arrays

    integer :: atom_id, ila

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)
!       call calculateMadelungLatticeSum(self%madelung_sum_array(ila), atom_id, arrays%rbasis)
      call calculateMadelungLatticeSum(self%madelung_sum_array(ila), self%madelung_calc, atom_id, arrays%rbasis)
    enddo ! ila

  endsubroutine ! prepare Madelung

  !----------------------------------------------------------------------------
  subroutine destroyCalculationData(self) ! todo: text-replace self by self in this routine
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
    
    type(CalculationData), intent(inout) :: self

    integer :: ila, atom_id

    do ila = 1, self%num_local_atoms

      atom_id = self%atom_ids(ila)

      call destroyRefCluster(self%ref_cluster_array(ila))
      call destroyMadelungLatticeSum(self%madelung_sum_array(ila))
      call destroyBasisAtom(self%atomdata_array(ila))
      call destroyCellData(self%cell_array(ila))
      call destroyRadialMeshData(self%mesh_array(ila))
      call destroyLDAUData(self%ldau_data_array(ila))
      call destroyJijData(self%jij_data_array(ila))
      call destroyDensityResults(self%densities_array(ila))
      call destroyKKRresults(self%kkr_array(ila))
      call destroyEnergyResults(self%energies_array(ila))

    enddo ! ila

    call destroyLatticeVectors(self%lattice_vectors)
    call destroyMadelungCalculator(self%madelung_calc)
    call destroyGauntCoefficients(self%gaunts)
    call destroyShapeGauntCoefficients(self%shgaunts)
    call destroyBroydenData(self%broyden)
    call iguess_destroy(self%iguess_data)

    deallocate(self%mesh_array)
    deallocate(self%cell_array)
    deallocate(self%atomdata_array)
    deallocate(self%ref_cluster_array)
    deallocate(self%kkr_array)
    deallocate(self%densities_array)
    deallocate(self%energies_array)
    deallocate(self%madelung_sum_array)
    deallocate(self%trunc_zone)
    deallocate(self%clusters)
    deallocate(self%iguess_data)
    deallocate(self%ldau_data_array)
    deallocate(self%jij_data_array)
    deallocate(self%broyden)
    deallocate(self%atom_ids)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  integer function getNumLocalAtoms(self)
    type(CalculationData), intent(in) :: self

    getNumLocalAtoms = self%num_local_atoms
  endfunction

  !----------------------------------------------------------------------------
  integer function getAtomIndexOfLocal(self, ila)
    type(CalculationData), intent(in) :: self
    integer, intent(in) :: ila

    getAtomIndexOfLocal = self%atom_ids(ila)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to atomdata for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getAtomData(self, local_atom_index)
    type(BasisAtom), pointer :: getAtomData ! return value

    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getAtomData => self%atomdata_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to 'reference cluster for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getRefCluster(self, local_atom_index)
    type(RefCluster), pointer :: getRefCluster ! return value

    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getRefCluster => self%ref_cluster_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to kkr(results) for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getKKR(self, local_atom_index)
    type(KKRresults), pointer :: getKKR ! return value

    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getKKR => self%kkr_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Madelung sum for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getMadelungSum(self, local_atom_index)
        
    type(MadelungLatticeSum), pointer :: getMadelungSum ! return value

    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getMadelungSum => self%madelung_sum_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to density results for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getDensities(self, local_atom_index)
    type(DensityResults), pointer :: getDensities ! return value

    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getDensities => self%densities_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to energy results for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getEnergies(self, local_atom_index)
    type(EnergyResults), pointer :: getEnergies ! return value

    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getEnergies => self%energies_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to LDA+U data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getLDAUData(self, local_atom_index)
    type(LDAUData), pointer :: getLDAUData ! return value

    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getLDAUData => self%ldau_data_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Jij-data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getJijData(self, local_atom_index)
    type(JijData), pointer :: getJijData ! return value

    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getJijData => self%jij_data_array(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Broyden data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getBroyden(self, local_atom_index) ! todo: remove local_atom_index from interface (although optional)
    type(BroydenData), pointer :: getBroyden ! return value

    type(CalculationData), intent(in) :: self
    integer, intent(in), optional :: local_atom_index

    getBroyden => self%broyden
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Gaunt coefficients.
  function getGaunts(self)
    type(GauntCoefficients), pointer :: getGaunts ! return value

    type(CalculationData), intent(in) :: self

    getGaunts => self%gaunts
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Shape-Gaunt coefficients.
  function getShapeGaunts(self)
    type(ShapeGauntCoefficients), pointer :: getShapeGaunts ! return value

    type(CalculationData), intent(in) :: self

    getShapeGaunts => self%shgaunts
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to Madelung calculator.
  function getMadelungCalculator(self)
    type(MadelungCalculator), pointer :: getMadelungCalculator ! return value

    type(CalculationData), intent(in) :: self

    getMadelungCalculator => self%madelung_calc
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to truncation zone.
  function getTruncationZone(self)
    type(TruncationZone), pointer :: getTruncationZone ! return value
    type(CalculationData), intent(in) :: self

    getTruncationZone => self%trunc_zone
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to cluster info (sparsity info).
  function getClusterInfo(self)
    type(ClusterInfo), pointer :: getClusterInfo ! return value

    type(CalculationData), intent(in) :: self

    getClusterInfo => self%clusters
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to lattice vector table.
  function getLatticeVectors(self)
    type(LatticeVectors), pointer :: getLatticeVectors ! return value

    type(CalculationData), intent(in) :: self

    getLatticeVectors => self%lattice_vectors
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to initial guess data.
  function getInitialGuessData(self)
    type(InitialGuess), pointer :: getInitialGuessData

    type(CalculationData), intent(in) :: self

    getInitialGuessData => self%iguess_data
  endfunction

  !----------------------------------------------------------------------------
  !> Returns record length needed for 'meshes' file.
  integer function getMaxReclenMeshes(self)
    type(CalculationData), intent(in) :: self

    getMaxReclenMeshes = self%max_reclen_meshes
  endfunction

  !----------------------------------------------------------------------------
  !> Returns record length needed for 'meshes' file.
  integer function getMaxReclenPotential(self)
    type(CalculationData), intent(in) :: self

    getMaxReclenPotential = self%max_reclen_potential
  endfunction

! ==================== Helper routines ========================================

  !----------------------------------------------------------------------------
  !> Helper routine: called by createCalculationData.
  subroutine constructEverything(self, dims, params, arrays, my_mpi)
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
    
    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays
    type(KKRnanoParallel), intent(in) :: my_mpi

    integer :: atom_id, ila, irmd

    call createLatticeVectors(self%lattice_vectors, arrays%bravais)

    ! create cluster for each local atom
    !$omp parallel do private(ila)
    do ila = 1, self%num_local_atoms
      call createRefCluster(self%ref_cluster_array(ila), self%lattice_vectors, arrays%rbasis, params%rclust, self%atom_ids(ila))
      !write(*,*) "Atoms in ref. cluster: ", self%ref_cluster_array(ila)%nacls
    enddo
    !$omp endparallel do

    ! setup the truncation zone
    call initLcutoffNew(self%trunc_zone, self%atom_ids, arrays)


    call createClusterInfo_com(self%clusters, self%ref_cluster_array, self%trunc_zone, getMySEcommunicator(my_mpi))

    if (isMasterRank(my_mpi)) then
      write(*,*) "Number of lattice vectors created     : ", self%lattice_vectors%nrd
      write(*,*) "Max. number of reference cluster atoms: ", self%clusters%naclsd
      write(*,*) "On node 0: "
      write(*,*) "Num. atoms treated with full lmax: ", num_untruncated
      write(*,*) "Num. atoms in truncation zone 1  : ", num_truncated
      write(*,*) "Num. atoms in truncation zone 2  : ", num_truncated2
    endif ! master
    CHECKASSERT( num_truncated + num_untruncated + num_truncated2 == dims%naez )

    call createMadelungCalculator(self%madelung_calc, dims%lmaxd, params%alat, params%rmax, params%gmax, arrays%bravais)

    ! a very crucial routine
    call generateAtomsShapesMeshes(self, dims, params, arrays)

    call recordLengths_com(self, my_mpi)

    if (isInMasterGroup(my_mpi)) then
#ifndef TASKLOCAL_FILES
      call writePotentialIndexFile(self)
#endif
      call writeNewMeshFiles(self)
    endif

    do ila = 1, self%num_local_atoms

      atom_id = self%atom_ids(ila)
      irmd = self%mesh_array(ila)%irmd

      call createKKRresults(self%kkr_array(ila), dims, self%clusters%naclsd)
      call createDensityResults(self%densities_array(ila), dims, irmd)
      call createEnergyResults(self%energies_array(ila), dims%nspind, dims%lmaxd)

      call createLDAUData(self%ldau_data_array(ila), params%ldau, irmd, dims%lmaxd, dims%nspind)
      call createJijData(self%jij_data_array(ila), .false., params%rcutjij, dims%nxijd, dims%lmmaxd,dims%nspind)

      call createMadelungLatticeSum(self%madelung_sum_array(ila), self%madelung_calc%lmxspd, dims%naez) 

      ! ASSERT( arrays%ZAT(atom_id) == atomdata%Z_nuclear )

    enddo ! ila

    ! calculate Gaunt coefficients
    call createGauntCoefficients(self%gaunts, dims%lmaxd)
    call createShapeGauntCoefficients(self%shgaunts, dims%lmaxd)
    call createBroydenData(self%broyden, getBroydenDim(self), dims%itdbryd, params%imix, params%mixing)  ! getBroydenDim replaces former NTIRD

    call setup_iguess(self, dims, arrays) ! setup storage for iguess
    
  endsubroutine ! constructEverything

  !----------------------------------------------------------------------------
  !> Initialise iguess datastructure.
  subroutine setup_iguess(self, dims, arrays)
    use DimParams_mod, only: DimParams
    use Main2Arrays_mod, only: Main2Arrays
    use InitialGuess_mod, only: iguess_init ! deprecated exceptional naming with underscore
    use TEST_lcutoff_mod, only: num_untruncated
     
    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)  :: dims
    type(Main2Arrays), intent(in):: arrays

    integer, allocatable :: num_k_points(:)
    integer :: ii, blocksize

    ! TODO: This is overdimensioned when l-cutoff is used!!!
    if (num_untruncated /= dims%naez) &
      warn(6, "The memory proportions for iGuess are overdimensioned when l-dependent truncation is applied!")
    ! DO NOT USE IGUESS together with l-cutoff!!! RS-cutoff is fine
    blocksize = self%trunc_zone%naez_trc * self%num_local_atoms * dims%lmmaxd**2

    allocate(num_k_points(dims%iemxd))
    do ii = 1, dims%iemxd
      num_k_points(ii) = arrays%nofks(arrays%kmesh(ii))
    enddo ! ii

    ! setup storage for iguess
    if (dims%smpid == 1 .and. dims%nspind == 2) then
      ! no spin parallelisation choosen, processes must store both spin-directions
      call iguess_init(self%iguess_data, num_k_points, 2, blocksize, dims%iguessd)
    else
      call iguess_init(self%iguess_data, num_k_points, 1, blocksize, dims%iguessd)
    endif

  endsubroutine ! setup_iguess

!------------------------------------------------------------------------------
!> Generates basis atom information, radial mesh, shape-function and
!> interpolates starting potential if necessary
  subroutine generateAtomsShapesMeshes(self, dims, params, arrays)
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
    
    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays

    integer :: ila, atom_id
    type(BasisAtom), allocatable      :: old_atom_array(:)
    type(RadialMeshData), allocatable :: old_mesh_array(:)
    double precision, allocatable     :: new_MT_radii(:)

    allocate(old_atom_array(self%num_local_atoms))
    allocate(old_mesh_array(self%num_local_atoms))
    allocate(new_MT_radii(self%num_local_atoms))

    ! generate storage for cell information + shape-functions
    do ila = 1, self%num_local_atoms
      call createCellData(self%cell_array(ila), dims%irid, (2*dims%LPOT+1)**2, (2*dims%LPOT+1)**2)
    enddo ! ila

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)

      ! We want to allow the actual radial mesh to be different from the one given by the input
      ! Therefore read 'potential' and 'meshes' data into temporary data structures
      ! Then interpolate potential to the new mesh

      ! load the input data
      call createBasisAtomFromFile(old_atom_array(ila), "atoms", "vpotnew.0", atom_id)

      call createRadialMeshDataFromFile(old_mesh_array(ila), "meshes.0", atom_id)

      call associateBasisAtomMesh(old_atom_array(ila), old_mesh_array(ila))

      new_MT_radii(ila) = old_atom_array(ila)%radius_muffin_tin / params%alat
    enddo ! ila

    ! generate shapes and meshes
    call generateShapesTEST(self, dims, params, arrays, new_MT_radii, params%MT_scale)

    ! interpolate to new mesh
    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)

      ! Geometry might have changed - interpolate to new mesh
      call interpolateBasisAtom(self%atomdata_array(ila), old_atom_array(ila), self%mesh_array(ila), dims%lpot)

      ! set new MT radius
      self%atomdata_array(ila)%radius_muffin_tin = self%mesh_array(ila)%rmt

      ! set radius of repulsive reference potential
      if (params%RMT_ref_scale > 0.d0) then
        self%atomdata_array(ila)%RMTref = self%cell_array(ila)%shdata%max_muffin_tin * params%alat * params%RMT_ref_scale
      else
        self%atomdata_array(ila)%RMTref = self%atomdata_array(ila)%radius_muffin_tin ! old behaviour=Mt-radius
      endif

      self%cell_array(ila)%cell_index = self%atomdata_array(ila)%cell_index
      call associateBasisAtomCell(self%atomdata_array(ila), self%cell_array(ila))

      CHECKASSERT( dims%IRMIND == self%mesh_array(ila)%IRMIN ) !check mesh
      CHECKASSERT( self%atomdata_array(ila)%atom_index == atom_id )

      call destroyBasisAtom(old_atom_array(ila))
      call destroyRadialMeshData(old_mesh_array(ila))
    enddo ! ila

    deallocate(new_MT_radii, old_atom_array, old_mesh_array)

  endsubroutine generateAtomsShapesMeshes

!------------------------------------------------------------------------------
  subroutine generateShapesTEST(self, dims, params, arrays, new_MT_radii, MT_scale)
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use ConstructShapes_mod, only: construct, InterstitialMesh, destroyInterstitialMesh, write_shapefun_file
    use ShapefunData_mod, only: ShapefunData
    use RadialMeshData_mod, only: createRadialMeshData, initRadialMesh
    use ShapefunData_mod, only: destroyShapefunData

    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays
    double precision, intent(in) :: new_MT_radii(:)
    double precision, intent(in) :: MT_scale

    integer :: atom_id, ila, irmd, irid, ipand, irnsd
    type(InterstitialMesh) :: inter_mesh
    type(ShapefunData) :: shdata ! temporary shape-fun data
    double precision :: new_MT_radius
    integer :: num_MT_points

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)

      new_MT_radius = new_MT_radii(ila)
      num_MT_points = params%num_MT_points

      call construct(shdata, inter_mesh, arrays%rbasis, arrays%bravais, atom_id, &
                     params%rclust_voronoi, 4*dims%lmaxd, &
                     dims%irid - num_MT_points, &
                     params%nmin_panel, num_MT_points, new_MT_radius, MT_scale)

      ! use it
      self%cell_array(ila)%shdata = shdata ! possible in Fortran 2003

      irmd = dims%irmd - dims%irid + size(inter_mesh%xrn)
      irid = size(inter_mesh%xrn)
      ipand = size(inter_mesh%nm) + 1
      irnsd = irmd - (dims%irmd - dims%irnsd)

      ASSERT(inter_mesh%xrn(1) /= 0.d0)
      !write(*,*) irmd, irid, ipand, irnsd

      call createRadialMeshData(self%mesh_array(ila), irmd, ipand)

      call initRadialMesh(self%mesh_array(ila), params%alat, inter_mesh%xrn, &
                          inter_mesh%drn, inter_mesh%nm, irmd - irid, irnsd)

      ! optional output of shape functions
      if (params%write_shapes == 1) call write_shapefun_file(shdata, inter_mesh, atom_id)

      call destroyShapefunData(shdata)
      call destroyInterstitialMesh(inter_mesh)

    enddo ! ila

  endsubroutine ! generateShapesTEST

  !----------------------------------------------------------------------------
  !> Communicate and set record lengths.
  subroutine recordLengths_com(self, my_mpi)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getMySECommunicator, getMyAtomRank
    use RadialMeshData_mod, only: getMinReclenMesh
    use BasisAtom_mod, only: getMinReclenBasisAtomPotential
    include 'mpif.h'

    type(CalculationData), intent(inout) :: self
    type(KKRnanoParallel), intent(in) :: my_mpi

    integer, parameter :: ND=2
    integer :: ierr, ila, sendbuf(ND), recvbuf(ND)

    sendbuf = -1
    recvbuf = -1
    do ila = 1, self%num_local_atoms
      sendbuf(1) = max(sendbuf(1), getMinReclenBasisAtomPotential(self%atomdata_array(ila)))
      sendbuf(2) = max(sendbuf(2), getMinReclenMesh(self%mesh_array(ila)))
    enddo ! ila

    call MPI_Allreduce(sendbuf, recvbuf, ND, MPI_INTEGER, MPI_MAX, getMySECommunicator(my_mpi), ierr)

    if (getMyAtomRank(my_mpi) == 0) then
      write(*,*) "Record length 'vpotnew' file: ", recvbuf(1)
      write(*,*) "Record length 'meshes'  file: ", recvbuf(2)
    endif

    self%max_reclen_potential = recvbuf(1)
    self%max_reclen_meshes    = recvbuf(2)

  endsubroutine

  !----------------------------------------------------------------------------
  !> Write potential index file.
  subroutine writePotentialIndexFile(self)
    use BasisAtom_mod, only: openBasisAtomPotentialIndexDAFile, writeBasisAtomPotentialIndexDA, closeBasisAtomPotentialIndexDAFile
    type(CalculationData), intent(in) :: self

    integer :: ila, atom_id, max_reclen

    max_reclen = getMaxReclenPotential(self)

    ! the opening routine requires any instance of type BasisAtom
    call openBasisAtomPotentialIndexDAFile(self%atomdata_array(1), 37, 'vpotnew.idx')

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)
      call writeBasisAtomPotentialIndexDA(self%atomdata_array(ila), 37, atom_id, max_reclen)
    enddo ! ila

    call closeBasisAtomPotentialIndexDAFile(37)
  endsubroutine ! write potential file

  !----------------------------------------------------------------------------
  !> Write new mesh files.
  !>
  !> The mesh can deviate from the input mesh if atoms are not in an ideal
  !> position. Therefore new mesh files have to be written.
  subroutine writeNewMeshFiles(self)
#ifndef TASKLOCAL_FILES
    use RadialMeshData_mod, only: openRadialMeshDataIndexDAFile, writeRadialMeshDataIndexDA, closeRadialMeshDataIndexDAFile
#endif
    use RadialMeshData_mod, only: openRadialMeshDataDAFile, writeRadialMeshDataDA, closeRadialMeshDataDAFile
    
    type(CalculationData), intent(in) :: self

    integer :: ila, atom_id, max_reclen

    max_reclen = getMaxReclenMeshes(self)

    ! the opening routines require any instance of type RadialMeshData
#ifndef TASKLOCAL_FILES
    ! do not write index when using task-local files
    call openRadialMeshDataIndexDAFile(self%mesh_array(1), 37, 'meshes.idx')
#endif
    call openRadialMeshDataDAFile(self%mesh_array(1), 38, 'meshes', max_reclen)

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)
#ifndef TASKLOCAL_FILES
      call writeRadialMeshDataIndexDA(self%mesh_array(ila), 37, atom_id, max_reclen)
#endif
      call writeRadialMeshDataDA(self%mesh_array(ila), 38, atom_id)
    enddo ! ila

    call closeRadialMeshDataDAFile(38)
#ifndef TASKLOCAL_FILES
    call closeRadialMeshDataIndexDAFile(37)
#endif
  endsubroutine ! writeNewMeshFiles

!============ Helper routines for Broyden mixing ==============================

  !----------------------------------------------------------------------------
  !> Returns the number of potential values of ALL LOCAL atoms.
  !> This is needed for dimensioning the Broyden mixing work arrays.
  integer function getBroydenDim(self) result(ndof)
    use PotentialData_mod, only: getNumPotentialValues ! todo: make it an elemental function
    
    type(CalculationData), intent(in) :: self

    integer :: ila

    ndof = 0
    do ila = 1, self%num_local_atoms
      ndof = ndof + getNumPotentialValues(self%atomdata_array(ila)%potential)
    enddo ! ila
    
  endfunction ! getBroydenDim

  !----------------------------------------------------------------------------
  ! Print some debugging info
  subroutine repr_CalculationData(self)
    use RadialMeshData_mod, only: repr_RadialMeshData
    use ShapefunData_mod, only: repr_ShapefunData
    use PotentialData_mod, only: repr_PotentialData

    type(CalculationData), intent(in) :: self

    integer :: ila
    character(len=:), allocatable :: str

    do ila = 1, self%num_local_atoms
      call repr_RadialMeshData(self%mesh_array(ila), str)
      write(*, '(A)') str
      call repr_PotentialData(self%atomdata_array(ila)%potential, str)
      write(*, '(A)') str
      call repr_ShapefunData(self%cell_array(ila)%shdata, str)
      write(*, '(A)') str
    enddo ! ila
    
  endsubroutine ! represent

  
#if 0  
!==============================================================================
!=             WORK in PROGRESS - not used yet                                =
!==============================================================================

  ! Factored out some routines from 'constructEverything'

  subroutine constructClusters(self, params, arrays)
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use RefCluster_mod, only: createRefCluster, createLatticeVectors
    
    type(CalculationData), intent(inout) :: self
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays

    integer :: ila

    call createLatticeVectors(self%lattice_vectors, arrays%bravais)

    ! create cluster for each local atom
    !$omp parallel do private(ila)
    do ila = 1, self%num_local_atoms
      call createRefCluster(self%ref_cluster_array(ila), self%lattice_vectors, arrays%rbasis, params%rclust, self%atom_ids(ila))
    enddo ! ila
    !$omp endparallel do

  endsubroutine ! constructClusters

  subroutine constructTruncationZones(self, arrays, my_mpi, naez)
    use KKRnanoParallel_mod, only: KKRnanoParallel, getMySEcommunicator, isMasterRank   
    use Main2Arrays_mod, only: Main2Arrays
    use TEST_lcutoff_mod, only: num_untruncated, num_truncated2, num_truncated ! integers
    use TEST_lcutoff_mod, only: initLcutoffNew
    use ClusterInfo_mod, only: createClusterInfo_com

    type(CalculationData), intent(inout) :: self
    type(Main2Arrays), intent(in):: arrays
    type(KKRnanoParallel), intent(in) :: my_mpi
    integer, intent(in) :: naez

    ! setup the truncation zone
    call initLcutoffNew(self%trunc_zone, self%atom_ids, arrays)

    ! get information about all the reference clusters of atoms in truncation zone
    call createClusterInfo_com(self%clusters, self%ref_cluster_array, self%trunc_zone, getMySEcommunicator(my_mpi))

    if (isMasterRank(my_mpi)) then
      write(*,*) "Number of lattice vectors created     : ", self%lattice_vectors%nrd
      write(*,*) "Max. number of reference cluster atoms: ", self%clusters%naclsd
      write(*,*) "On node 0: "
      write(*,*) "Num. atoms treated with full lmax: ", num_untruncated
      write(*,*) "Num. atoms in truncation zone 1  : ", num_truncated
      write(*,*) "Num. atoms in truncation zone 2  : ", num_truncated2
    endif
    CHECKASSERT(num_truncated+num_untruncated+num_truncated2 == naez)
    
  endsubroutine ! constructTruncationZones

  subroutine constructStorage(self, dims, params)
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use KKRresults_mod, only: createKKRresults
    use DensityResults_mod, only: createDensityResults
    use EnergyResults_mod, only: createEnergyResults
    use LDAUData_mod, only: createLDAUData
    use JijData_mod, only: createJijData
    use MadelungCalculator_mod, only: createMadelungLatticeSum
    
    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params

    integer :: atom_id, ila, irmd

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)
      irmd = self%mesh_array(ila)%irmd ! abbrev.

      call createKKRresults(self%kkr_array(ila), dims, self%clusters%naclsd)
      call createDensityResults(self%densities_array(ila), dims, irmd)
      call createEnergyResults(self%energies_array(ila), dims%nspind, dims%lmaxd)
      call createLDAUData(self%ldau_data_array(ila), params%ldau, irmd, dims%lmaxd, dims%nspind)
      call createJijData(self%jij_data_array(ila), params%jij, params%rcutjij, dims%nxijd, dims%lmmaxd,dims%nspind)
      call createMadelungLatticeSum(self%madelung_sum_array(ila), self%madelung_calc, dims%naez)
    enddo ! ila
    
  endsubroutine ! constructStorage
#endif

endmodule CalculationData_mod
