!> @author Elias Rabel

#include "DebugHelpers/test_macros.h"

module CalculationData_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)

  use MadelungCalculator_mod, only: MadelungCalculator, create, destroy
  use MadelungCalculator_mod, only: MadelungLatticeSum, create, destroy ! todo: two types hosted in one module
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
  use LatticeVectors_mod, only: LatticeVectors, create, destroy
  
  implicit none
  private

  public :: CalculationData, create, destroy, represent
  
  public :: prepareMadelung, getBroydenDim
  public :: getDensities, getEnergies, getLDAUData, getAtomData

  type CalculationData

    integer :: num_local_atoms !< atoms in this process
    integer, allocatable :: atom_ids(:)
    integer :: max_reclen_meshes
    integer :: max_reclen_potential

    ! atom local data - different for each atom
    type(RadialMeshData), pointer     :: mesh_a(:)         => null()
    type(CellData), pointer           :: cell_a(:)         => null()
    type(BasisAtom), pointer          :: atomdata_a(:)     => null()
    type(DensityResults), pointer     :: densities_a(:)    => null()
    type(EnergyResults), pointer      :: energies_a(:)     => null()
    type(LDAUData), pointer           :: ldau_data_a(:)    => null()
    type(KKRresults),         allocatable :: kkr_a(:)
    type(JijData),            allocatable :: jij_data_a(:)
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
    type(InitialGuess)                  :: iguess_data ! storage for initial guess
  endtype ! CalculationData
  
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
  subroutine createCalculationData(self, dims, params, arrays, mp, voronano)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays

    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)   :: dims
    type(InputParams), intent(in) :: params
    type(Main2Arrays), intent(in) :: arrays
    type(KKRnanoParallel), intent(in) :: mp
    integer, intent(in) :: voronano

    integer :: atoms_per_proc, num_local_atoms, ila

    atoms_per_proc = dims%naez / mp%numAtomRanks
    ASSERT( mp%numAtomRanks * atoms_per_proc == dims%naez )
    num_local_atoms = atoms_per_proc !TODO

    self%num_local_atoms = num_local_atoms

    ! one datastructure for each local atom
    allocate(self%mesh_a(num_local_atoms)) ! only used inside this module
    allocate(self%cell_a(num_local_atoms)) ! only used inside this module
    allocate(self%ref_cluster_a(num_local_atoms)) ! only visible to this module and ScatteringCalculation_mod.F90
    allocate(self%atomdata_a(num_local_atoms))
    allocate(self%kkr_a(num_local_atoms))
    allocate(self%densities_a(num_local_atoms))
    allocate(self%energies_a(num_local_atoms))
    allocate(self%ldau_data_a(num_local_atoms))
    allocate(self%madelung_sum_a(num_local_atoms)) ! only visible to this module and MadelungPotential_mod.F90
    allocate(self%jij_data_a(num_local_atoms)) ; if(num_local_atoms > 1) warn(6, "Jij work with max. 1 atom so far!") 
    
    allocate(self%atom_ids(num_local_atoms))

    ! assign atom ids to processes with atom rank 'mp%myAtomRank'
    ! E.g. for 2 atoms per proc:
    ! process 1 treats atoms 1,2
    ! process 2 treats atoms 3,4 and so on
    ! FOR USE OF TRUNCATION THESE atoms have to be close together!!!

    ASSERT( size(self%atom_ids) == num_local_atoms )

    do ila = 1, num_local_atoms
      self%atom_ids(ila) = mp%myAtomRank * atoms_per_proc + ila
      ASSERT( self%atom_ids(ila) <= dims%naez )
    enddo ! ila

    ! Now construct all datastructures and calculate initial data
    call constructEverything(self, dims, params, arrays, mp, voronano)

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
  subroutine destroyCalculationData(self)
    type(CalculationData), intent(inout) :: self

    integer :: ila
    
    do ila = 1, self%num_local_atoms

      call destroy(self%ref_cluster_a(ila))
      call destroy(self%madelung_sum_a(ila))
      call destroy(self%atomdata_a(ila))
      call destroy(self%cell_a(ila))
      call destroy(self%mesh_a(ila))
      call destroy(self%ldau_data_a(ila))
      call destroy(self%jij_data_a(ila))
      call destroy(self%densities_a(ila))
      call destroy(self%kkr_a(ila))
      call destroy(self%energies_a(ila))

    enddo ! ila

    call destroy(self%lattice_vectors)
    call destroy(self%madelung_calc)
    call destroy(self%gaunts)
    call destroy(self%shgaunts)
    call destroy(self%broyden)
    call destroy(self%iguess_data)

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
  !> Returns reference to atomdata for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getAtomData(self, local_atom_index)
    type(BasisAtom), pointer :: getAtomData ! return value
    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getAtomData => self%atomdata_a(local_atom_index)
  endfunction ! get

  !----------------------------------------------------------------------------
  !> Returns reference to density results for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getDensities(self, local_atom_index)
    type(DensityResults), pointer :: getDensities ! return value
    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getDensities => self%densities_a(local_atom_index)
  endfunction ! get

  !----------------------------------------------------------------------------
  !> Returns reference to energy results for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getEnergies(self, local_atom_index)
    type(EnergyResults), pointer :: getEnergies ! return value
    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getEnergies => self%energies_a(local_atom_index)
  endfunction ! get

  !----------------------------------------------------------------------------
  !> Returns reference to LDA+U data for atom with LOCAL atom index
  !> 'local_atom_index'.
  function getLDAUData(self, local_atom_index)
    type(LDAUData), pointer :: getLDAUData ! return value
    type(CalculationData), intent(in) :: self
    integer, intent(in) :: local_atom_index

    getLDAUData => self%ldau_data_a(local_atom_index)
  endfunction ! get


  !----------------------------------------------------------------------------
  !> Helper routine: called by createCalculationData.
  subroutine constructEverything(self, dims, params, arrays, mp, voronano)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use TEST_lcutoff_mod, only: num_truncated
    ! deprecated interfaces
    use LatticeVectors_mod, only: createLatticeVectors
    use RefCluster_mod, only: createRefCluster              
    use TEST_lcutoff_mod, only: initLcutoffNew                
    use ClusterInfo_mod, only: createClusterInfo         
    use DensityResults_mod, only: createDensityResults          
    use EnergyResults_mod, only: createEnergyResults           
    use LDAUData_mod, only: createLDAUData                
    use JijData_mod, only: createJijData                 
    use MadelungCalculator_mod, only: createMadelungCalculator
    use MadelungCalculator_mod, only: createMadelungLatticeSum      
    use GauntCoefficients_mod, only: createGauntCoefficients       
    use ShapeGauntCoefficients_mod, only: createShapeGauntCoefficients  
    use BroydenData_mod, only: createBroydenData             
    use KKRresults_mod, only: createKKRresults
!   include 'mpif.h'
    
    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays
    type(KKRnanoParallel), intent(in) :: mp
    integer, intent(in) :: voronano

    integer :: atom_id, ila, irmd

    call createLatticeVectors(self%lattice_vectors, arrays%bravais)

    ! create cluster for each local atom
    !$omp parallel do private(ila)
    do ila = 1, self%num_local_atoms
      call createRefCluster(self%ref_cluster_a(ila), self%lattice_vectors%rr, arrays%rbasis, params%rclust, self%atom_ids(ila))
!     write(*,*) "Atoms in ref. cluster: ", self%ref_cluster_a(ila)%nacls
    enddo ! ila
    !$omp endparallel do

    call initLcutoffNew(self%trunc_zone, self%atom_ids, arrays, params%lcutoff_radii, params%cutoff_radius, params%solver) ! setup the truncation zone

    call createClusterInfo(self%clusters, self%ref_cluster_a, self%trunc_zone, mp%mySEComm)

    if (mp%isMasterRank) then
      write(*,*) "Number of lattice vectors created     : ", self%lattice_vectors%nrd
      write(*,*) "Max. number of reference cluster atoms: ", self%clusters%naclsd
      write(*,*) "On node 0: "
      write(*,*) "Num. atoms treated with full lmax:        ", num_truncated(dims%lmaxd)
      write(*,*) "Num. atoms outside of the truncation zone:", num_truncated(-1)
      write(*,*) "Num. atoms with l-dependent truncation:   ", sum(num_truncated) - num_truncated(-1) - num_truncated(dims%lmaxd)
    endif ! master
    CHECKASSERT( sum(num_truncated) == dims%naez )

    call create(self%madelung_calc, dims%lmaxd, params%alat, params%rmax, params%gmax, arrays%bravais)

    call generateAtomsShapesMeshes(self, dims, params, arrays, mp, voronano) ! a very crucial routine
    
    if (voronano == 1) call jellstart12_wrapper(self, arrays, params, dims, mp)
   
    call recordLengths_com(self, mp)

    if (mp%isInMasterGroup) then
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
    call setup_iguess(self, dims, arrays%nofks) ! setup storage for iguess

  endsubroutine ! constructEverything

  !----------------------------------------------------------------------------
  !> Initialise iguess datastructure.
  subroutine setup_iguess(self, dims, nofks)!, kmesh)
    use DimParams_mod, only: DimParams
    use InitialGuess_mod, only: create ! deprecated exceptional naming with underscore
    use TEST_lcutoff_mod, only: num_truncated

    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in) :: dims
    integer, intent(in) :: nofks(:)!, kmesh(:)

    integer, allocatable :: num_k_points(:)
    integer :: ie, blocksize, ns

    ! TODO: This is overdimensioned when l-cutoff is used!!!
    if (num_truncated(dims%lmaxd) /= dims%naez) &
      warn(6, "The memory proportions for iGuess are overdimensioned when l-dependent truncation is applied!")
      
    ! DO NOT USE IGUESS together with l-cutoff!!! RS-cutoff is fine
    blocksize = self%trunc_zone%naez_trc * self%num_local_atoms * dims%lmmaxd**2

    allocate(num_k_points(dims%iemxd))
    do ie = 1, dims%iemxd
      num_k_points(ie) = nofks(1) ! nofks(kmesh(ie)) !! cannot easily access emesh%kmesh, use the largest
    enddo ! ie

    ns = 1; if (dims%smpid == 1 .and. dims%nspind == 2) ns = 2 ! no spin parallelisation choosen, processes must store both spin-directions
    call create(self%iguess_data, num_k_points, ns, blocksize, dims%iguessd) ! setup storage for iguess

  endsubroutine ! setup_iguess

  !------------------------------------------------------------------------------
  !> Generates basis atom information, radial mesh, shape-function and
  !> interpolates starting potential if necessary
  subroutine generateAtomsShapesMeshes(self, dims, params, arrays, mp, voronano)
    use DimParams_mod, only: DimParams
    use InterpolateBasisAtom_mod, only: interpolateBasisAtom
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    
    use CellData_mod, only: createCellData
    use BasisAtom_mod, only: createBasisAtomFromFile
    use RadialMeshData_mod, only: createRadialMeshDataFromFile
    use BasisAtom_mod, only: associateBasisAtomMesh, associateBasisAtomCell
    use KKRnanoParallel_mod, only: KKRnanoParallel
     
    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)   :: dims
    type(InputParams), intent(in) :: params
    type(Main2Arrays), intent(in) :: arrays
    type(KKRnanoParallel), intent(in) :: mp
    integer, intent(in) :: voronano
    
    integer :: ila, atom_id, ist
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
    call generateShapesTEST(self, dims, params, arrays, new_MT_radii, params%MT_scale, mp, voronano)

    ! interpolate to new mesh
    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)

      ! Geometry might have changed - interpolate to new mesh
      call interpolateBasisAtom(self%atomdata_a(ila), old_atom_a(ila), self%mesh_a(ila), dims%lpot)

      ! set new MT radius
      self%atomdata_a(ila)%radius_muffin_tin = self%mesh_a(ila)%rmt

      ! set radius of repulsive reference potential
      if (params%rMT_ref_scale > 0.d0) then
        self%atomdata_a(ila)%rMTref = self%cell_a(ila)%shdata%max_muffin_tin * params%alat * params%rMT_ref_scale
      else
        self%atomdata_a(ila)%rMTref = self%atomdata_a(ila)%radius_muffin_tin ! old behaviour=Mt-radius
      endif

      self%cell_a(ila)%cell_index = self%atomdata_a(ila)%cell_index
      call associateBasisAtomCell(self%atomdata_a(ila), self%cell_a(ila))

      CHECKASSERT( dims%IRMIND == self%mesh_a(ila)%IRMIN ) ! check mesh
      CHECKASSERT( self%atomdata_a(ila)%atom_index == atom_id )

      call destroy(old_atom_a(ila))
      call destroy(old_mesh_a(ila))
    enddo ! ila

    deallocate(new_MT_radii, old_atom_a, old_mesh_a, stat=ist)
  endsubroutine ! generateAtomsShapesMeshes

  !------------------------------------------------------------------------------
  subroutine generateShapesTEST(self, dims, params, arrays, new_MT_radii, MT_scale, mp, voronano)
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use ConstructShapes_mod, only: createShape, InterstitialMesh, destroy, write_shapefun_file
    use ShapefunData_mod, only: ShapefunData, destroy
    use RadialMeshData_mod, only: createRadialMeshData, initRadialMesh
    use KKRnanoParallel_mod, only: KKRnanoParallel
    include 'mpif.h'

    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays
    double precision, intent(in) :: new_MT_radii(:)
    double precision, intent(in) :: MT_scale
    type(KKRnanoParallel), intent(in) :: mp
    integer, intent(in) :: voronano
    
    integer :: i, atom_id, ila, irmd, irid, ipand, irnsd, ierror
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
                     dims%irid-num_MT_points, &
                     params%nmin_panel, num_MT_points, new_MT_radius, MT_scale, atom_id, dims%naez)

      ! use it
      self%cell_a(ila)%shdata = shdata ! possible in Fortran 2003

      irmd = dims%irmd - dims%irid + size(inter_mesh%xrn)
      irid = size(inter_mesh%xrn)
      ipand = size(inter_mesh%nm) + 1
      irnsd = irmd - (dims%irmd - dims%irnsd)

      ASSERT( inter_mesh%xrn(1) /= 0.d0 ) ! write(*,*) irmd, irid, ipand, irnsd

      call createRadialMeshData(self%mesh_a(ila), irmd, ipand)

      call initRadialMesh(self%mesh_a(ila), params%alat, inter_mesh%xrn, &
                          inter_mesh%drn, inter_mesh%nm, irmd-irid, irnsd)

      ! optional output of shape functions
      if (params%write_shapes == 1 .or. voronano == 1) call write_shapefun_file(shdata, inter_mesh, atom_id)

      call destroy(shdata)
      call destroy(inter_mesh)

    enddo ! ila

    ! optional output of shapefunctions in single file (VORONANO)
    if (voronano == 1) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      if (mp%isMasterRank) then
        open(15, file='shapefun.header', form="formatted", status='replace', action='READWRITE')
        write(15, fmt='(i5)') dims%naez
        do i = 1, dims%naez, 4 ! write old-style scaling factor, not used any more
          write(15, fmt='(4e20.12)') 1.d0,1.d0,1.d0,1.d0
        enddo ! i
        close(15)
#ifdef HAS_EXECUTE_COMMAND_LINE       
        call EXECUTE_COMMAND_LINE("cat shapefun.header shape.* > shapefun.voronano")
        call EXECUTE_COMMAND_LINE("rm shape.*")
        call EXECUTE_COMMAND_LINE("rm shapefun.header")
#else
        warn(6, "cannot invoke shell for command ""cat shapefun.header shape.* > shapefun.voronano && rm shape.* shapefun.header""")
#endif        
      endif ! master
    endif ! VORONANO

  endsubroutine ! generateShapesTEST


  !----------------------------------------------------------------------------
  !> Wrapper for jellstart12

   subroutine jellstart12_wrapper(calc_data, arrays, params, dims, mp)
   use InputParams_mod, only: InputParams
   use Main2Arrays_mod, only: Main2Arrays
   use DimParams_mod,   only: DimParams
   use KKRnanoParallel_mod, only: KKRnanoParallel
   use JelliumPotentials_mod, only: jellstart12
   include 'mpif.h' ! imports only MPI_COMM_WORLD, maybe removed

   type(CalculationData), intent(in) :: calc_data
   type(Main2Arrays),     intent(in) :: arrays
   type(InputParams),     intent(in) :: params
   type(DimParams),       intent(in) :: dims
   type(KKRnanoParallel), intent(in) :: mp

   integer :: ila, ins, meshn(1), naez, idshape(1), irws(1), irns(1), xrn_drn_max_dimension, ierror
   double precision :: rwscl(1), rmtcl(1), qbound, z(1)
   double precision, allocatable :: xrn_2(:,:)
   double precision, allocatable :: drn_2(:,:)

    xrn_drn_max_dimension = 1
    do ila = 1, calc_data%num_local_atoms
      xrn_drn_max_dimension = max(size(calc_data%mesh_a(ila)%r), xrn_drn_max_dimension)
    enddo ! ila

    allocate(xrn_2(xrn_drn_max_dimension,1))
    allocate(drn_2(xrn_drn_max_dimension,1))

    !$omp parallel do private(ila) !! check if enough variables are declared parallel !!
    do ila = 1, calc_data%num_local_atoms

      ins         = 1                              ! use non-spherical parts
      naez        = 1                              ! 'jellstart' is called for each atom
      z(1)        = arrays%zat(calc_data%atom_ids(ila)) ! nuclear charge
      idshape(1)  = 1                              ! one shape per atom
      meshn(1)    = calc_data%cell_a(ila)%shdata%irid   ! number of interstitial mesh points 
      rmtcl(1)    = (calc_data%mesh_a(ila)%rmt)         ! MT radius
      irws(1)     = calc_data%mesh_a(ila)%irws          ! index of max. radius 
      xrn_2(1:meshn(1),1) = (calc_data%mesh_a(ila)%r((irws(1)-meshn(1)):irws(1)))/params%alat    ! radial mesh points
      drn_2(1:meshn(1),1) = (calc_data%mesh_a(ila)%drdi((irws(1)-meshn(1)):irws(1)))/params%alat ! integration weight corresponding to radial mesh point
      rwscl(1)    = xrn_2(xrn_drn_max_dimension,1) ! Wigner-Seitz radius
      irns(1)     = calc_data%mesh_a(ila)%irns
      qbound      = 1.d-7

      call jellstart12(dims%nspind,ins,naez,z,idshape,  &
            rwscl,rmtcl,meshn,xrn_2,drn_2,  &
            irws,irns,  &
            params%alat,qbound, &
            dims%lpot, dims%irmd, dims%irnsd, dims%nspind, &
            calc_data%atom_ids(ila), &
            params%elementdatabasepath)
    enddo ! ila
    !$omp endparallel do

    ! output of potential in single file
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    if (mp%isMasterRank) then
#ifdef HAS_EXECUTE_COMMAND_LINE       
        call EXECUTE_COMMAND_LINE("cat potential.0* > potential.voronano")
        call EXECUTE_COMMAND_LINE("rm potential.0*")
#else
        warn(6, "cannot invoke shell for command ""cat potential.0* > potential.voronano && rm potential.0*""")
#endif       
    endif ! master

  endsubroutine ! jellstart12_wrapper
   

  !----------------------------------------------------------------------------
  !> Communicate and set record lengths.
  subroutine recordLengths_com(self, mp)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    use RadialMeshData_mod, only: getMinReclenMesh
    use BasisAtom_mod, only: getMinReclenBasisAtomPotential
    include 'mpif.h'

    type(CalculationData), intent(inout) :: self
    type(KKRnanoParallel), intent(in) :: mp

    integer, parameter :: ND=2
    integer :: ierr, ila, sendbuf(ND), recvbuf(ND)

    sendbuf = -1
    recvbuf = -1
    do ila = 1, self%num_local_atoms
      sendbuf(1) = max(sendbuf(1), getMinReclenBasisAtomPotential(self%atomdata_a(ila)))
      sendbuf(2) = max(sendbuf(2), getMinReclenMesh(self%mesh_a(ila)))
    enddo ! ila

    call MPI_Allreduce(sendbuf, recvbuf, ND, MPI_INTEGER, MPI_MAX, mp%mySEComm, ierr)

    if (mp%myAtomRank == 0) then
      write(*,*) "Record length 'vpotnew' file: ", recvbuf(1)
      write(*,*) "Record length 'meshes'  file: ", recvbuf(2)
    endif

    self%max_reclen_potential = recvbuf(1)
    self%max_reclen_meshes    = recvbuf(2)

  endsubroutine ! recordLengths

  !----------------------------------------------------------------------------
  !> Write potential index file.
  subroutine writePotentialIndexFile(self)
    use BasisAtom_mod, only: openBasisAtomPotentialIndexDAFile, writeBasisAtomPotentialIndexDA, closeBasisAtomPotentialIndexDAFile
    type(CalculationData), intent(in) :: self

    integer :: ila, atom_id, max_reclen

    max_reclen = self%max_reclen_potential

    ! the opening routine requires any instance of type BasisAtom
    call openBasisAtomPotentialIndexDAFile(self%atomdata_a(1), 37, 'vpotnew.idx', action='write')

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

    max_reclen = self%max_reclen_meshes

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
    use PotentialData_mod, only: getNumPotentialValues ! an elemental function
    type(CalculationData), intent(in) :: self

    ndof = sum(getNumPotentialValues(self%atomdata_a(1:self%num_local_atoms)%potential))
    
  endfunction ! getBroydenDim

  !----------------------------------------------------------------------------
  !> Print some debugging info
  subroutine repr_CalculationData(self)
    use RadialMeshData_mod, only: represent
    use ShapefunData_mod, only: represent
    use PotentialData_mod, only: represent

    type(CalculationData), intent(in) :: self

    integer :: ila
    character(len=:), allocatable :: str

    do ila = 1, self%num_local_atoms
      call represent(self%mesh_a(ila), str)
      write(*, '(A)') str
      call represent(self%atomdata_a(ila)%potential, str)
      write(*, '(A)') str
      call represent(self%cell_a(ila)%shdata, str)
      write(*, '(A)') str
    enddo ! ila

  endsubroutine ! represent

endmodule ! CalculationData_mod




#if 0
! !==============================================================================
! !=             WORK in PROGRESS - not used yet                                =
! !==============================================================================
! 
!   ! Factored out some routines from 'constructEverything'
! 
!   subroutine constructClusters(self, params, arrays)
!     use InputParams_mod, only: InputParams
!     use Main2Arrays_mod, only: Main2Arrays
!     use RefCluster_mod, only: createRefCluster, createLatticeVectors
!     
!     type(CalculationData), intent(inout) :: self
!     type(InputParams), intent(in):: params
!     type(Main2Arrays), intent(in):: arrays
! 
!     integer :: ila
! 
!     call createLatticeVectors(self%lattice_vectors, arrays%bravais)
! 
!     ! create cluster for each local atom
!     !$omp parallel do private(ila)
!     do ila = 1, self%num_local_atoms
!       call createRefCluster(self%ref_cluster_a(ila), self%lattice_vectors, arrays%rbasis, params%rclust, self%atom_ids(ila))
!     enddo ! ila
!     !$omp endparallel do
! 
!   endsubroutine ! constructClusters
! 
!   subroutine constructTruncationZones(self, arrays, mp, naez)
!     use KKRnanoParallel_mod, only: KKRnanoParallel
!     use Main2Arrays_mod, only: Main2Arrays
!     use TEST_lcutoff_mod, only: num_untruncated, num_truncated2, num_truncated ! integers
!     use TEST_lcutoff_mod, only: initLcutoffNew
!     use ClusterInfo_mod, only: createClusterInfo_com
! 
!     type(CalculationData), intent(inout) :: self
!     type(Main2Arrays), intent(in):: arrays
!     type(KKRnanoParallel), intent(in) :: mp
!     integer, intent(in) :: naez
! 
!     ! setup the truncation zone
!     call initLcutoffNew(self%trunc_zone, self%atom_ids, arrays)
! 
!     ! get information about all the reference clusters of atoms in truncation zone
!     call createClusterInfo_com(self%clusters, self%ref_cluster_a, self%trunc_zone, mp%mySEComm)
! 
!     if (mp%isMasterRank) then
!       write(*,*) "Number of lattice vectors created     : ", self%lattice_vectors%nrd
!       write(*,*) "Max. number of reference cluster atoms: ", self%clusters%naclsd
!       write(*,*) "On node 0: "
!       write(*,*) "Num. atoms treated with full lmax: ", num_untruncated
!       write(*,*) "Num. atoms in truncation zone 1  : ", num_truncated
!       write(*,*) "Num. atoms in truncation zone 2  : ", num_truncated2
!     endif ! is master
!     CHECKASSERT(num_truncated+num_untruncated+num_truncated2 == naez)
!     
!   endsubroutine ! constructTruncationZones
! 
!   subroutine constructStorage(self, dims, params)
!     use DimParams_mod, only: DimParams
!     use InputParams_mod, only: InputParams
!     use KKRresults_mod, only: createKKRresults
!     use DensityResults_mod, only: createDensityResults
!     use EnergyResults_mod, only: createEnergyResults
!     use LDAUData_mod, only: createLDAUData
!     use JijData_mod, only: createJijData
!     use MadelungCalculator_mod, only: createMadelungLatticeSum
!     
!     type(CalculationData), intent(inout) :: self
!     type(DimParams), intent(in)  :: dims
!     type(InputParams), intent(in):: params
! 
!     integer :: atom_id, ila, irmd
! 
!     do ila = 1, self%num_local_atoms
!       atom_id = self%atom_ids(ila)
!       irmd = self%mesh_a(ila)%irmd ! abbrev.
! 
!       call createKKRresults(self%kkr_a(ila), dims, self%clusters%naclsd)
!       call createDensityResults(self%densities_a(ila), dims, irmd)
!       call createEnergyResults(self%energies_a(ila), dims%nspind, dims%lmaxd)
!       call createLDAUData(self%ldau_data_a(ila), params%ldau, irmd, dims%lmaxd, dims%nspind)
!       call createJijData(self%jij_data_a(ila), params%jij, params%rcutjij, dims%nxijd, dims%lmmaxd,dims%nspind)
!       call createMadelungLatticeSum(self%madelung_sum_a(ila), self%madelung_calc, dims%naez)
!     enddo ! ila
!     
!   endsubroutine ! constructStorage
#endif
