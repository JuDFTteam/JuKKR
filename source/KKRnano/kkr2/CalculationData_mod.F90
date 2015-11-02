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

  public :: CalculationData, create, destroy, represent
  public :: createCalculationData, destroyCalculationData ! deprecated
  
  public :: getBroydenDim, getNumLocalAtoms, getAtomIndexOfLocal, getAtomData, getKKR
  public :: getDensities, getEnergies, getLDAUData
  public :: getMaxReclenMeshes, getMaxReclenPotential      
  
  public :: prepareMadelung         
  
! public :: constructTruncationZones, constructStorage, constructClusters ! not yet used


  type CalculationData

    integer :: num_local_atoms !< atoms in this process
    integer, allocatable :: atom_ids(:)
    integer :: max_reclen_meshes
    integer :: max_reclen_potential

    ! atom local data - different for each atom
    type(RadialMeshData), pointer     :: mesh_a(:)         => null()
    type(CellData), pointer           :: cell_a(:)         => null()
    type(BasisAtom), pointer          :: atomdata_a(:)     => null()
    type(KKRresults), pointer         :: kkr_a(:)          => null()
    type(DensityResults), pointer     :: densities_a(:)    => null()
    type(EnergyResults), pointer      :: energies_a(:)     => null()
    type(LDAUData), pointer           :: ldau_data_a(:)    => null()
    type(JijData), pointer                :: jij_data_a(:) => null()
    
    type(RefCluster),         allocatable :: ref_cluster_a(:)
    type(MadelungLatticeSum), allocatable :: madelung_sum_a(:)

    ! global data - same for each local atom
    type(LatticeVectors)                :: lattice_vectors
    type(MadelungCalculator)            :: madelung_calc
    type(ShapeGauntCoefficients)        :: shgaunts
    type(GauntCoefficients)             :: gaunts
    type(TruncationZone)                :: trunc_zone
    type(ClusterInfo)                   :: clusters
    type(BroydenData)                   :: broyden

    ! storage for initial guess
    type(InitialGuess)                  :: iguess_data

  endtype
  
  interface create
    module procedure createCalculationData
  endinterface
  
  interface destroy
    module procedure destroyCalculationData
  endinterface
  
  interface represent
    module procedure repr_CalculationData
  endinterface

  contains

  !----------------------------------------------------------------------------
  ! TODO: atoms_per_procs * num_procs MUST BE = naez, 
  !       rank = 0,1,..., num_atom_ranks-1
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

    integer :: atoms_per_proc, num_local_atoms, ila, atom_rank

    atoms_per_proc = dims%naez / getNumAtomRanks(my_mpi)
    ASSERT( getNumAtomRanks(my_mpi) * atoms_per_proc == dims%naez )
    num_local_atoms = atoms_per_proc !TODO

    self%num_local_atoms = num_local_atoms

    ! one datastructure for each local atom
    allocate(self%mesh_a(num_local_atoms))
    allocate(self%cell_a(num_local_atoms))
    allocate(self%atomdata_a(num_local_atoms))
    allocate(self%ref_cluster_a(num_local_atoms))
    allocate(self%kkr_a(num_local_atoms))
    allocate(self%densities_a(num_local_atoms))
    allocate(self%energies_a(num_local_atoms))
    allocate(self%madelung_sum_a(num_local_atoms))
    allocate(self%ldau_data_a(num_local_atoms))
    allocate(self%jij_data_a(num_local_atoms))

    atom_rank = getMyAtomRank(my_mpi)
    
    allocate(self%atom_ids(num_local_atoms))

    ! assign atom ids to processes with atom rank 'atom_rank'
    ! E.g. for 2 atoms per proc:
    ! process 1 treats atoms 1,2
    ! process 2 treats atoms 3,4 and so on
    ! FOR USE OF TRUNCATION THESE atoms have to be close together!!!

    ASSERT( size(self%atom_ids) == num_local_atoms )
    ASSERT( mod(dims%naez, atoms_per_proc) == 0 )

    do ila = 1, num_local_atoms
      self%atom_ids(ila) = atom_rank * atoms_per_proc + ila
      ASSERT( self%atom_ids(ila) <= dims%naez )
    enddo ! ila

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
      call calculateMadelungLatticeSum(self%madelung_sum_a(ila), self%madelung_calc, atom_id, arrays%rbasis)
    enddo ! ila

!   stop 'DEBUG: stop after calculateMadelungLatticeSum in CalculationData_mod.F90:166'    
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

      call destroyRefCluster(self%ref_cluster_a(ila))
      call destroyMadelungLatticeSum(self%madelung_sum_a(ila))
      call destroyBasisAtom(self%atomdata_a(ila))
      call destroyCellData(self%cell_a(ila))
      call destroyRadialMeshData(self%mesh_a(ila))
      call destroyLDAUData(self%ldau_data_a(ila))
      call destroyJijData(self%jij_data_a(ila))
      call destroyDensityResults(self%densities_a(ila))
      call destroyKKRresults(self%kkr_a(ila))
      call destroyEnergyResults(self%energies_a(ila))

    enddo ! ila

    call destroyLatticeVectors(self%lattice_vectors)
    call destroyMadelungCalculator(self%madelung_calc)
    call destroyGauntCoefficients(self%gaunts)
    call destroyShapeGauntCoefficients(self%shgaunts)
    call destroyBroydenData(self%broyden)
    call iguess_destroy(self%iguess_data)

    deallocate(self%mesh_a)
    deallocate(self%cell_a)
    deallocate(self%atomdata_a)
    deallocate(self%ref_cluster_a)
    deallocate(self%kkr_a)
    deallocate(self%densities_a)
    deallocate(self%energies_a)
    deallocate(self%madelung_sum_a)
    deallocate(self%ldau_data_a)
    deallocate(self%jij_data_a)
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

    getAtomData => self%atomdata_a(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to kkr(results) for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getKKR(self, local_atom_index)
    type(KKRresults), pointer :: getKKR ! return value
    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getKKR => self%kkr_a(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to density results for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getDensities(self, local_atom_index)
    type(DensityResults), pointer :: getDensities ! return value
    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getDensities => self%densities_a(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to energy results for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getEnergies(self, local_atom_index)
    type(EnergyResults), pointer :: getEnergies ! return value
    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getEnergies => self%energies_a(local_atom_index)
  endfunction

  !----------------------------------------------------------------------------
  !> Returns reference to LDA+U data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getLDAUData(self, local_atom_index)
    type(LDAUData), pointer :: getLDAUData ! return value
    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getLDAUData => self%ldau_data_a(local_atom_index)
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
    use TEST_lcutoff_mod, only: num_truncated
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
      call createRefCluster(self%ref_cluster_a(ila), self%lattice_vectors, arrays%rbasis, params%rclust, self%atom_ids(ila))
!     write(*,*) "Atoms in ref. cluster: ", self%ref_cluster_a(ila)%nacls
    enddo
    !$omp endparallel do

    call initLcutoffNew(self%trunc_zone, self%atom_ids, arrays) ! setup the truncation zone

    call createClusterInfo_com(self%clusters, self%ref_cluster_a, self%trunc_zone, getMySEcommunicator(my_mpi))

    if (isMasterRank(my_mpi)) then
      write(*,*) "Number of lattice vectors created     : ", self%lattice_vectors%nrd
      write(*,*) "Max. number of reference cluster atoms: ", self%clusters%naclsd
      write(*,*) "On node 0: "
      write(*,*) "Num. atoms treated with full lmax: ", num_truncated(0)
      write(*,*) "Num. atoms in truncation zone 1  : ", num_truncated(1)
      write(*,*) "Num. atoms in truncation zone 2  : ", num_truncated(2)
    endif ! master
    CHECKASSERT( sum(num_truncated) == dims%naez )

    call createMadelungCalculator(self%madelung_calc, dims%lmaxd, params%alat, params%rmax, params%gmax, arrays%bravais)

    call generateAtomsShapesMeshes(self, dims, params, arrays) ! a very crucial routine

    call recordLengths_com(self, my_mpi)

    if (isInMasterGroup(my_mpi)) then
#ifndef TASKLOCAL_FILES
      call writePotentialIndexFile(self)
#endif
      call writeNewMeshFiles(self)
    endif ! in master group

    do ila = 1, self%num_local_atoms

      atom_id = self%atom_ids(ila)
      irmd = self%mesh_a(ila)%irmd

      call createKKRresults(self%kkr_a(ila), dims, self%clusters%naclsd)
      call createDensityResults(self%densities_a(ila), dims, irmd)
      call createEnergyResults(self%energies_a(ila), dims%nspind, dims%lmaxd)

      call createLDAUData(self%ldau_data_a(ila), params%ldau, irmd, dims%lmaxd, dims%nspind)
      call createJijData(self%jij_data_a(ila), .false., params%rcutjij, dims%nxijd, dims%lmmaxd,dims%nspind)

      call createMadelungLatticeSum(self%madelung_sum_a(ila), self%madelung_calc%lmxspd, dims%naez) 

      ! ASSERT( arrays%ZAT(atom_id) == atomdata%Z_nuclear )

    enddo ! ila

    ! calculate Gaunt coefficients
    call createGauntCoefficients(self%gaunts, dims%lmaxd)
    call createShapeGauntCoefficients(self%shgaunts, dims%lmaxd)
    
!   write(*,*) __FILE__,__LINE__," createBroydenData deavtivated for DEBUG!"
    call createBroydenData(self%broyden, getBroydenDim(self), dims%itdbryd, params%imix, params%mixing)  ! getBroydenDim replaces former NTIRD

!   write(*,*) __FILE__,__LINE__," setup_iguess deavtivated for DEBUG!"
    call setup_iguess(self, dims, arrays) ! setup storage for iguess

  endsubroutine ! constructEverything

  !----------------------------------------------------------------------------
  !> Initialise iguess datastructure.
  subroutine setup_iguess(self, dims, arrays)
    use DimParams_mod, only: DimParams
    use Main2Arrays_mod, only: Main2Arrays
    use InitialGuess_mod, only: iguess_init ! deprecated exceptional naming with underscore
    use TEST_lcutoff_mod, only: num_truncated
     
    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)  :: dims
    type(Main2Arrays), intent(in):: arrays

    integer, allocatable :: num_k_points(:)
    integer :: ii, blocksize, ns

    ! TODO: This is overdimensioned when l-cutoff is used!!!
    if (num_truncated(0) /= dims%naez) & ! num_truncated(0) == former num_untruncated
      warn(6, "The memory proportions for iGuess are overdimensioned when l-dependent truncation is applied!")
    ! DO NOT USE IGUESS together with l-cutoff!!! RS-cutoff is fine
    blocksize = self%trunc_zone%naez_trc * self%num_local_atoms * dims%lmmaxd**2

    allocate(num_k_points(dims%iemxd))
    do ii = 1, dims%iemxd
      num_k_points(ii) = arrays%nofks(arrays%kmesh(ii))
    enddo ! ii

    ns = 1; if (dims%smpid == 1 .and. dims%nspind == 2) ns = 2 ! no spin parallelisation choosen, processes must store both spin-directions
    call iguess_init(self%iguess_data, num_k_points, ns, blocksize, dims%iguessd) ! setup storage for iguess

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
    type(BasisAtom), allocatable      :: old_atom_a(:)
    type(RadialMeshData), allocatable :: old_mesh_a(:)
    double precision, allocatable     :: new_MT_radii(:)

    allocate(old_atom_a(self%num_local_atoms))
    allocate(old_mesh_a(self%num_local_atoms))
    allocate(new_MT_radii(self%num_local_atoms))

    ! generate storage for cell information + shape-functions
    do ila = 1, self%num_local_atoms
      call createCellData(self%cell_a(ila), dims%irid, (2*dims%LPOT+1)**2, (2*dims%LPOT+1)**2)
    enddo ! ila

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)

      ! We want to allow the actual radial mesh to be different from the one given by the input
      ! Therefore read 'potential' and 'meshes' data into temporary data structures
      ! Then interpolate potential to the new mesh

      ! load the input data
      call createBasisAtomFromFile(old_atom_a(ila), "atoms", "vpotnew.0", atom_id)

      call createRadialMeshDataFromFile(old_mesh_a(ila), "meshes.0", atom_id)

      call associateBasisAtomMesh(old_atom_a(ila), old_mesh_a(ila))

      new_MT_radii(ila) = old_atom_a(ila)%radius_muffin_tin / params%alat
    enddo ! ila

    ! generate shapes and meshes
    call generateShapesTEST(self, dims, params, arrays, new_MT_radii, params%MT_scale)

    ! interpolate to new mesh
    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)

      ! Geometry might have changed - interpolate to new mesh
      call interpolateBasisAtom(self%atomdata_a(ila), old_atom_a(ila), self%mesh_a(ila), dims%lpot)

      ! set new MT radius
      self%atomdata_a(ila)%radius_muffin_tin = self%mesh_a(ila)%rmt

      ! set radius of repulsive reference potential
      if (params%RMT_ref_scale > 0.d0) then
        self%atomdata_a(ila)%RMTref = self%cell_a(ila)%shdata%max_muffin_tin * params%alat * params%RMT_ref_scale
      else
        self%atomdata_a(ila)%RMTref = self%atomdata_a(ila)%radius_muffin_tin ! old behaviour=Mt-radius
      endif

      self%cell_a(ila)%cell_index = self%atomdata_a(ila)%cell_index
      call associateBasisAtomCell(self%atomdata_a(ila), self%cell_a(ila))

      CHECKASSERT( dims%IRMIND == self%mesh_a(ila)%IRMIN ) ! check mesh
      CHECKASSERT( self%atomdata_a(ila)%atom_index == atom_id )

      call destroyBasisAtom(old_atom_a(ila))
      call destroyRadialMeshData(old_mesh_a(ila))
    enddo ! ila

    deallocate(new_MT_radii, old_atom_a, old_mesh_a)

  endsubroutine ! generateAtomsShapesMeshes

!------------------------------------------------------------------------------
  subroutine generateShapesTEST(self, dims, params, arrays, new_MT_radii, MT_scale)
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use ConstructShapes_mod, only: createShape, InterstitialMesh, destroyInterstitialMesh, write_shapefun_file
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

      call createShape(shdata, inter_mesh, arrays%rbasis, arrays%bravais, atom_id, &
                     params%rclust_voronoi, 4*dims%lmaxd, &
                     dims%irid - num_MT_points, &
                     params%nmin_panel, num_MT_points, new_MT_radius, MT_scale, atom_id)

      ! use it
      self%cell_a(ila)%shdata = shdata ! possible in Fortran 2003

      irmd = dims%irmd - dims%irid + size(inter_mesh%xrn)
      irid = size(inter_mesh%xrn)
      ipand = size(inter_mesh%nm) + 1
      irnsd = irmd - (dims%irmd - dims%irnsd)

      ASSERT( inter_mesh%xrn(1) /= 0.d0 ) ! write(*,*) irmd, irid, ipand, irnsd

      call createRadialMeshData(self%mesh_a(ila), irmd, ipand)

      call initRadialMesh(self%mesh_a(ila), params%alat, inter_mesh%xrn, &
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
      sendbuf(1) = max(sendbuf(1), getMinReclenBasisAtomPotential(self%atomdata_a(ila)))
      sendbuf(2) = max(sendbuf(2), getMinReclenMesh(self%mesh_a(ila)))
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
    call openBasisAtomPotentialIndexDAFile(self%atomdata_a(1), 37, 'vpotnew.idx')

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)
      call writeBasisAtomPotentialIndexDA(self%atomdata_a(ila), 37, atom_id, max_reclen)
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
    call openRadialMeshDataIndexDAFile(self%mesh_a(1), 37, 'meshes.idx')
#endif
    call openRadialMeshDataDAFile(self%mesh_a(1), 38, 'meshes', max_reclen)

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)
#ifndef TASKLOCAL_FILES
      call writeRadialMeshDataIndexDA(self%mesh_a(ila), 37, atom_id, max_reclen)
#endif
      call writeRadialMeshDataDA(self%mesh_a(ila), 38, atom_id)
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
      ndof = ndof + getNumPotentialValues(self%atomdata_a(ila)%potential)
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
      call repr_RadialMeshData(self%mesh_a(ila), str)
      write(*, '(A)') str
      call repr_PotentialData(self%atomdata_a(ila)%potential, str)
      write(*, '(A)') str
      call repr_ShapefunData(self%cell_a(ila)%shdata, str)
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
      call createRefCluster(self%ref_cluster_a(ila), self%lattice_vectors, arrays%rbasis, params%rclust, self%atom_ids(ila))
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
    call createClusterInfo_com(self%clusters, self%ref_cluster_a, self%trunc_zone, getMySEcommunicator(my_mpi))

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
      irmd = self%mesh_a(ila)%irmd ! abbrev.

      call createKKRresults(self%kkr_a(ila), dims, self%clusters%naclsd)
      call createDensityResults(self%densities_a(ila), dims, irmd)
      call createEnergyResults(self%energies_a(ila), dims%nspind, dims%lmaxd)
      call createLDAUData(self%ldau_data_a(ila), params%ldau, irmd, dims%lmaxd, dims%nspind)
      call createJijData(self%jij_data_a(ila), params%jij, params%rcutjij, dims%nxijd, dims%lmmaxd,dims%nspind)
      call createMadelungLatticeSum(self%madelung_sum_a(ila), self%madelung_calc, dims%naez)
    enddo ! ila
    
  endsubroutine ! constructStorage
#endif

endmodule CalculationData_mod
