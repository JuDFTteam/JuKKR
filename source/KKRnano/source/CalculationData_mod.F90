!> @author Elias Rabel

! #include "DebugHelpers/test_macros.h"

module CalculationData_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)

  use MadelungCalculator_mod, only: MadelungCalculator, create, destroy
  use MadelungCalculator_mod, only: MadelungLatticeSum, create, destroy ! todo: two types hosted in one module
  use RadialMeshData_mod, only: RadialMeshData, create, destroy
  use ChebMeshData_mod, only: ChebMeshData ! NOCO
  use ShapefunData_mod, only: ShapefunData, create, destroy
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
  use ExchangeTable_mod, only: ExchangeTable, create, destroy
  use NonCollinearMagnetismData_mod, only: NOCOData, create, load, destroy, loadascii
  
  implicit none
  private

  public :: CalculationData, create, destroy, represent
  
  public :: prepareMadelung, getBroydenDim
  public :: getDensities, getEnergies, getLDAUData, getAtomData


  type CalculationData

    integer :: num_local_atoms !< atoms in this process
    integer :: max_local_atoms !< MPI_MAX of num_local_atoms
    integer :: max_reclen_meshes
    integer :: max_reclen_potential

!     type(LocalAtomData), allocatable :: a(:) !> dim(num_local_atoms)
    
    integer(kind=4), allocatable :: atom_ids(:) ! to be replace by a(ila)%atom_id
    
    ! atom local data - different for each atom
    type(RadialMeshData), pointer     :: mesh_a(:)         => null()
    type(ChebMeshData), pointer       :: cheb_mesh_a(:)    => null() ! NOCO
    type(ShapefunData), pointer       :: cell_a(:)         => null()
    type(BasisAtom), pointer          :: atomdata_a(:)     => null()
    type(DensityResults), pointer     :: densities_a(:)    => null()
    type(EnergyResults), pointer      :: energies_a(:)     => null()
    type(LDAUData), pointer           :: ldau_data_a(:)    => null()
    type(KKRresults),         allocatable :: kkr_a(:)
    type(JijData),            pointer     :: jij_data_a(:) => null()
    type(RefCluster),         allocatable :: ref_cluster_a(:)
    type(MadelungLatticeSum), allocatable :: madelung_sum_a(:)

    ! global data - same for each local atom
    type(LatticeVectors)                :: lattice_vectors
    type(MadelungCalculator)            :: madelung_calc
    type(ShapeGauntCoefficients)        :: shgaunts
    type(GauntCoefficients)             :: gaunts
    type(TruncationZone)                :: trunc_zone
    type(ClusterInfo)                   :: clusters
    type(ExchangeTable)                 :: xtable
    type(BroydenData)                   :: broyden
    type(InitialGuess)                  :: iguess_data ! storage for initial guess
    type(NOCOData)                      :: noco_data   ! NOCO
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
  subroutine createCalculationData(self, dims, params, arrays, mp, kmesh, voronano)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays

    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)   :: dims
    type(InputParams), intent(in) :: params
    type(Main2Arrays), intent(in) :: arrays
    type(KKRnanoParallel), intent(in) :: mp
    integer, intent(in) :: kmesh(:) !> index of k-mesh for each energy point
    integer, intent(in) :: voronano

    integer :: num_local_atoms, ila

    self%max_local_atoms = (dims%naez - 1) / mp%numAtomRanks + 1
!   num_local_atoms = min(max(0, dims%naez - self%max_local_atoms*mp%myAtomRank), self%max_local_atoms)
    num_local_atoms = min(dims%naez - self%max_local_atoms*mp%myAtomRank, self%max_local_atoms)
    if (num_local_atoms < self%max_local_atoms) &
      write(*, '(9(a,i0))') "Treat only ",num_local_atoms," instead of ",self%max_local_atoms," local atoms in rank #",mp%myAtomRank
    if (num_local_atoms < 1) die_here("The number of local atoms must be >= 1 or the rank should be inactive, found"+num_local_atoms)
    num_local_atoms = max(0, num_local_atoms) ! however, this does not work

    self%num_local_atoms = num_local_atoms ! store

!   allocate(self%a(num_local_atoms))

    ! one datastructure for each local atom
    allocate(self%mesh_a(num_local_atoms)) ! only used inside this module
    allocate(self%cheb_mesh_a(num_local_atoms)) ! only used inside this module, NOCO
    allocate(self%cell_a(num_local_atoms)) ! only used inside this module
    allocate(self%ref_cluster_a(num_local_atoms)) ! only visible to this module and ScatteringCalculation_mod.F90
    allocate(self%atomdata_a(num_local_atoms))
    allocate(self%kkr_a(num_local_atoms))
    allocate(self%densities_a(num_local_atoms))
    allocate(self%energies_a(num_local_atoms))
    allocate(self%ldau_data_a(num_local_atoms))
    allocate(self%madelung_sum_a(num_local_atoms)) ! only visible to this module and MadelungPotential_mod.F90
    allocate(self%jij_data_a(num_local_atoms)) 
    if(num_local_atoms > 1 .and. params%Jij) warn(6, "Jij work with max. 1 atom so far!") 
    
    allocate(self%atom_ids(num_local_atoms))

    ! assign atom ids to processes with atom rank 'mp%myAtomRank'
    ! E.g. for 2 atoms per proc:
    ! process 1 treats atoms 1,2
    ! process 2 treats atoms 3,4 and so on
    ! FOR USE OF TRUNCATION THESE atoms have to be close together!!!

!   write(*,*) "size(self%atom_ids) = ",size(self%atom_ids),"  num_local_atoms = ",num_local_atoms 
    assert( size(self%atom_ids) == num_local_atoms )

    do ila = 1, num_local_atoms
      self%atom_ids(ila) = mp%myAtomRank * self%max_local_atoms + ila
!     self%a(ila)%atom_id = mp%myAtomRank * self%max_local_atoms + ila
      assert( self%atom_ids(ila) <= dims%naez )
    enddo ! ila

    ! Now construct all datastructures and calculate initial data
    call constructEverything(self, dims, params, arrays, mp, kmesh, voronano)

    ! call repr_CalculationData(self)
  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Calculate Madelung Lattice sums for all local atoms.
  subroutine prepareMadelung(self, rbasis)
    use MadelungCalculator_mod, only: calculate

    type(CalculationData), intent(inout) :: self
    double precision, intent(in) :: rbasis(:,:) ! dim(3,naez_all)

    integer :: atom_id, ila

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)
      call calculate(self%madelung_sum_a(ila), self%madelung_calc, atom_id, rbasis) ! this still scales N^3
    enddo ! ila

  endsubroutine ! prepare Madelung

  !----------------------------------------------------------------------------
  elemental subroutine destroyCalculationData(self)
    type(CalculationData), intent(inout) :: self

    integer :: ist ! ignore status
    
    call destroy(self%ref_cluster_a)
    call destroy(self%madelung_sum_a)
    call destroy(self%atomdata_a)
    call destroy(self%cell_a)
    call destroy(self%mesh_a)
    call destroy(self%ldau_data_a)
    call destroy(self%jij_data_a)
    call destroy(self%densities_a)
    call destroy(self%kkr_a)
    call destroy(self%energies_a)

!   deallocate(self%a)

    call destroy(self%lattice_vectors)
    call destroy(self%madelung_calc)
    call destroy(self%shgaunts)
    call destroy(self%gaunts)
    call destroy(self%trunc_zone)
    call destroy(self%clusters)
    call destroy(self%xtable)
    call destroy(self%broyden)
    call destroy(self%iguess_data)
    call destroy(self%noco_data) ! NOCO

    deallocate(self%mesh_a, stat=ist)
    deallocate(self%cheb_mesh_a, stat=ist) ! NOCO
    deallocate(self%cell_a, stat=ist)
    deallocate(self%atomdata_a, stat=ist)
    deallocate(self%ref_cluster_a, stat=ist)
    deallocate(self%kkr_a, stat=ist)
    deallocate(self%densities_a, stat=ist)
    deallocate(self%energies_a, stat=ist)
    deallocate(self%madelung_sum_a, stat=ist)
    deallocate(self%ldau_data_a, stat=ist)
    deallocate(self%jij_data_a, stat=ist)
    deallocate(self%atom_ids, stat=ist)
    
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
  subroutine constructEverything(self, dims, params, arrays, mp, kmesh, voronano)
    use KKRnanoParallel_mod, only: KKRnanoParallel
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use Truncation_mod, only: num_truncated, initTruncation                
    
    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)  :: dims
    type(InputParams), intent(in):: params
    type(Main2Arrays), intent(in):: arrays
    type(KKRnanoParallel), intent(in) :: mp
    integer, intent(in) :: kmesh(:) !> index of k-mesh for each energy point
    integer, intent(in) :: voronano

    integer :: atom_id, ila, irmd

    call create(self%lattice_vectors, arrays%bravais) ! createLatticeVectors

    ! create cluster for each local atom
    !$omp parallel do private(ila)
    do ila = 1, self%num_local_atoms
      call create(self%ref_cluster_a(ila), self%lattice_vectors%rr, arrays%rbasis, params%rclust, self%atom_ids(ila)) ! createRefCluster
!     write(*,*) "Atoms in ref. cluster: ", self%ref_cluster_a(ila)%nacls
    enddo ! ila
    !$omp endparallel do

    call initTruncation(self%trunc_zone, self%atom_ids, dims%lmaxd, arrays%bravais, arrays%rbasis, params%lcutoff_radii, params%cutoff_radius, mp%mySEComm) ! setup the truncation zone

    call create(self%xtable, self%trunc_zone%global_atom_id, comm=mp%mySEComm, max_local_atoms=self%max_local_atoms) ! createExchangeTable

    call create(self%clusters, self%ref_cluster_a, self%trunc_zone, xTable=self%xtable) ! createClusterInfo

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
    
    if (voronano == 1) then
      call jellstart12_wrapper(self, arrays, params, dims, mp)
      return      
    endif
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
!     write(*,*) here,trim(" initialize  local atom #"-ila+"global atom_id="-atom_id)

      call create(self%kkr_a(ila), dims, self%clusters%naclsd) ! createKKRresults
      call create(self%densities_a(ila), dims%lmpotd, dims%lmaxd, dims%iemxd, dims%nspind, irmd) ! createDensityResults
      call create(self%energies_a(ila), dims%nspind, dims%lmaxd) ! createEnergyResults

      call create(self%ldau_data_a(ila), params%ldau, irmd, dims%lmaxd, dims%nspind) ! createLDAUData
      call create(self%jij_data_a(ila), .false., params%rcutjij, dims%nxijd, dims%lmmaxd,dims%nspind) ! createJijData

      call create(self%madelung_sum_a(ila), self%madelung_calc%lmxspd, dims%naez) ! createMadelungLatticeSum


      ! assert( arrays%ZAT(atom_id) == atomdata%Z_nuclear )
!     write(*,*) here,trim(" initialized local atom #"-ila+"global atom_id="-atom_id)
    enddo ! ila

    call create(self%noco_data, dims%naez) ! createNOCOData

!   in case of a NOCO calculation - read file 'nonco_angle.dat'
    if (dims%korbit == 1) then
       if(dims%nspind .NE. 2) die_here('NSPIND=2 in global.conf is mandatory for SOC calculations')
       call loadascii(self%noco_data%theta_noco, self%noco_data%phi_noco, self%noco_data%angle_fixed, dims%naez)
       !call store(noco, 'bin.noco.0')
    else
       if(dims%korbit .NE. 0) die_here('When not using NOCO: KORBIT in global.conf should be zero')    
    endif
    

!    if (dims%korbit == 1) then
!        call load(self%noco_data, 'bin.noco.0') ! every process does this!, NOCO
!    endif

    ! calculate Gaunt coefficients
    call create(self%gaunts, dims%lmaxd) ! createGauntCoefficients
    call create(self%shgaunts, dims%lmaxd) ! createShapeGauntCoefficients
    
!   write(*,*) __FILE__,__LINE__," createBroydenData deavtivated for DEBUG!"
    call create(self%broyden, getBroydenDim(self), dims%itdbryd, params%imix, params%mixing) ! createBroydenData ! getBroydenDim replaces former NTIRD

!   write(*,*) __FILE__,__LINE__," setup_iguess deavtivated for DEBUG!"
    call setup_iguess(self, dims, arrays%nofks, kmesh) ! setup storage for iguess

  endsubroutine ! constructEverything

  !----------------------------------------------------------------------------
  !> Initialise iguess datastructure.
  subroutine setup_iguess(self, dims, nofks, kmesh)
    use DimParams_mod, only: DimParams
    use InitialGuess_mod, only: create
    use Truncation_mod, only: num_elements

    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in) :: dims
    integer, intent(in) :: nofks(:) !> number of k-points for each k-mesh
    integer, intent(in) :: kmesh(:) !> index of k-mesh for each energy point

    integer, allocatable :: num_k_points(:)
    integer :: ie, blocksize, ns
    

    blocksize = sum(num_elements(:)) ! ToDo: prepare this for noco

    allocate(num_k_points(dims%iemxd))
    do ie = 1, dims%iemxd
      num_k_points(ie) = nofks(kmesh(ie))
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
    
    use ShapefunData_mod, only: create
    use BasisAtom_mod, only: load, associateBasisAtomMesh, associateBasisAtomCell, associateBasisAtomChebMesh
    use RadialMeshData_mod, only: load
    use ChebMeshData_mod, only: load, create, construct ! NOCO
    use KKRnanoParallel_mod, only: KKRnanoParallel
     
    type(CalculationData), intent(inout) :: self
    type(DimParams), intent(in)   :: dims
    type(InputParams), intent(in) :: params
    type(Main2Arrays), intent(in) :: arrays
    type(KKRnanoParallel), intent(in) :: mp
    integer, intent(in) :: voronano
    
    integer :: ila, atom_id, ist, ii
    type(BasisAtom), allocatable      :: old_atom_a(:)
    type(RadialMeshData), allocatable :: old_mesh_a(:)
    double precision, allocatable     :: new_MT_radii(:)

    integer :: npan_tot_cheb, irmd_cheb ! NOCO
    
    double precision, parameter :: tolvdist = 1.d-12

    allocate(old_atom_a(self%num_local_atoms))
    allocate(old_mesh_a(self%num_local_atoms))
    allocate(new_MT_radii(self%num_local_atoms))

#ifndef USE_OLD_MESH
    ! generate storage for cell information + shape-functions
    do ila = 1, self%num_local_atoms
      call create(self%cell_a(ila), dims%irid, (2*dims%LPOT+1)**2, (2*dims%LPOT+1)**2)
    enddo ! ila
#endif

    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)

      ! We want to allow the actual radial mesh to be different from the one given by the input
      ! Therefore read 'potential' and 'meshes' data into temporary data structures
      ! Then interpolate potential to the new mesh

      ! load the input data unless it is a voronano run
      if (voronano == 0) then
        call load(old_atom_a(ila), "bin.atoms", "bin.vpotnew.0", atom_id)

        call load(old_mesh_a(ila), "bin.meshes.0", atom_id)

        call associateBasisAtomMesh(old_atom_a(ila), old_mesh_a(ila))

        new_MT_radii(ila) = old_atom_a(ila)%radius_muffin_tin / params%alat

!        if (dims%korbit == 1) then ! NOCO
!          call load(self%cheb_mesh_a(ila),"bin.chebmeshes.0", atom_id, old_mesh_a(ila)%nfu, params)
!        endif

      else
        if(params%MT_scale < tolvdist) then
          write(*,*) 'Warning: MT_scale in input.conf is set to 0.0d0 -> MT radius is automatically set to 2.37'
        endif
        new_MT_radii(ila) = 2.37
      endif ! voronano

    enddo ! ila
#ifndef USE_OLD_MESH
    ! generate shapes and meshes
    call generateShapesTEST(self, dims, params, arrays, new_MT_radii, params%MT_scale, mp, voronano)
#endif

    ! interpolate to new mesh
    if (voronano == 0) then
      do ila = 1, self%num_local_atoms
        atom_id = self%atom_ids(ila)

#ifndef USE_OLD_MESH
        ! Geometry might have changed - interpolate to new mesh
        call interpolateBasisAtom(self%atomdata_a(ila), old_atom_a(ila), self%mesh_a(ila), dims%lpot)
#endif

#ifdef USE_OLD_MESH
        ! >>> use old mesh>>> 
        self%mesh_a(ila) = old_mesh_a(ila)
        self%atomdata_a(ila) = old_atom_a(ila)
        call associateBasisAtomMesh(self%atomdata_a(ila), self%mesh_a(ila))
        ! generate storage for cell information + shape-functions
        call create(self%cell_a(ila), self%mesh_a(ila)%meshn, (2*dims%LPOT+1)**2, (2*dims%LPOT+1)**2)
        ! >>> fill shdata with values from old mesh >>>
        do ii = 1, self%mesh_a(ila)%nfu
           self%cell_a(ila)%llmsp(ii) = self%mesh_a(ila)%llmsp(ii)
           self%cell_a(ila)%lmsp(self%mesh_a(ila)%llmsp(ii)) = 1
           self%cell_a(ila)%ifunm(self%mesh_a(ila)%llmsp(ii)) = ii
        end do
        ! copy the shape functions
        self%cell_a(ila)%theta(1:self%mesh_a(ila)%meshn, 1:self%mesh_a(ila)%nfu) = self%mesh_a(ila)%thetas(1:self%mesh_a(ila)%meshn,1:self%mesh_a(ila)%nfu)
        self%cell_a(ila)%nfu = self%mesh_a(ila)%nfu
        ! set maximum possible muffin-tin radius (ALAT units)
        self%cell_a(ila)%max_muffin_tin = self%mesh_a(ila)%rmt/params%alat
        self%cell_a(ila)%num_faces = 9999 ! not needed
        ! <<< fill shdata with values from old mesh <<<
        ! <<< use old mesh<<<
#endif
        if (dims%korbit == 1) then ! NOCO
          npan_tot_cheb = self%mesh_a(ila)%ipand-1+params%npan_eq+params%npan_log  ! number of overall intervals in Chebychev mesh  
          irmd_cheb     = npan_tot_cheb*(params%ncheb+1)                              ! number of radial mesh points in Chebychev mesh
     
          call create(self%cheb_mesh_a(ila), irmd_cheb, npan_tot_cheb, self%cell_a(ila)%nfu, params) ! create data for new radial mesh for atom iatom
        endif
        if (dims%korbit == 1) then ! NOCO
          call construct(params%r_log,params%npan_log,params%npan_eq,params%ncheb, &
                      self%cheb_mesh_a(ila)%npan_lognew,self%cheb_mesh_a(ila)%npan_eqnew, &
                      self%cheb_mesh_a(ila)%npan_tot,self%cheb_mesh_a(ila)%rnew, &
                      self%cheb_mesh_a(ila)%rpan_intervall,self%cheb_mesh_a(ila)%ipan_intervall, &
                      self%cheb_mesh_a(ila)%thetasnew,self%cell_a(ila)%theta, &
                      self%cell_a(ila)%nfu,self%mesh_a(ila))
        endif
        if (dims%korbit == 1) then ! NOCO
          call associateBasisAtomChebMesh(self%atomdata_a(ila), self%cheb_mesh_a(ila))
        endif

        ! set new MT radius
        self%atomdata_a(ila)%radius_muffin_tin = self%mesh_a(ila)%rmt

        ! set radius of repulsive reference potential
        if (params%rMT_ref_scale > 0.d0) then
          self%atomdata_a(ila)%rMTref = self%cell_a(ila)%max_muffin_tin * params%alat * params%rMT_ref_scale
        else
          self%atomdata_a(ila)%rMTref = self%atomdata_a(ila)%radius_muffin_tin ! old behaviour=Mt-radius
        endif

        self%cell_a(ila)%cell_index = self%atomdata_a(ila)%cell_index
        call associateBasisAtomCell(self%atomdata_a(ila), self%cell_a(ila))

#ifndef USE_OLD_MESH
        CHECKASSERT( dims%IRMIND == self%mesh_a(ila)%IRMIN ) ! check mesh
#endif
        CHECKASSERT( self%atomdata_a(ila)%atom_index == atom_id )

#ifndef USE_OLD_MESH
        call destroy(old_atom_a(ila))
        call destroy(old_mesh_a(ila))
#endif
      enddo ! ila
    endif

    deallocate(new_MT_radii)
#ifndef USE_OLD_MESH
    deallocate(old_atom_a, old_mesh_a, stat=ist)
#endif
  endsubroutine ! generateAtomsShapesMeshes

  !------------------------------------------------------------------------------
  subroutine generateShapesTEST(self, dims, params, arrays, new_MT_radii, MT_scale, mp, voronano)
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use ConstructShapes_mod, only: create, InterstitialMesh, destroy, store
    use ShapefunData_mod, only: ShapefunData, destroy
    use RadialMeshData_mod, only: create, initRadialMesh
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
    double precision :: new_MT_radius
    integer :: num_MT_points, nwarn

    nwarn = 0
    do ila = 1, self%num_local_atoms
      atom_id = self%atom_ids(ila)

      new_MT_radius = new_MT_radii(ila)
      num_MT_points = params%num_MT_points

      call create(self%cell_a(ila), inter_mesh, arrays%rbasis, arrays%bravais, atom_id, &
                  params%rclust_voronoi, 4*dims%lmaxd, &
                  dims%irid-num_MT_points, &
                  params%nmin_panel, num_MT_points, new_MT_radius, MT_scale, atom_id, nwarn)

      irmd = dims%irmd - dims%irid + size(inter_mesh%xrn)
      irid = size(inter_mesh%xrn)
      ipand = size(inter_mesh%nm) + 1
      irnsd = irmd - (dims%irmd - dims%irnsd)

      assert( inter_mesh%xrn(1) /= 0.d0 ) ! write(*,*) irmd, irid, ipand, irnsd

      call create(self%mesh_a(ila), irmd, ipand)

      call initRadialMesh(self%mesh_a(ila), params%alat, inter_mesh%xrn, &
                          inter_mesh%drn, inter_mesh%nm, irmd-irid, irnsd)

      ! optional output of shape functions
      if (params%write_shapes == 1 .or. voronano == 1) call store(self%cell_a(ila), inter_mesh, atom_id)

      call destroy(inter_mesh)

    enddo ! ila
    if (nwarn > 0) warn(6, "file voro_weights cannot be opened in"+nwarn+" cases, use 1.0 for all.")
    
    ! optional output of shapefunctions in single file (VORONANO)
    if (voronano == 1) then
      call MPI_Barrier(MPI_COMM_WORLD, ierror)
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
      meshn(1)    = calc_data%cell_a(ila)%irid   ! number of interstitial mesh points 
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
    call MPI_Barrier(MPI_COMM_WORLD, ierror)
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
      write(*,*) "Record length 'bin.vpotnew' file: ", recvbuf(1)
      write(*,*) "Record length 'bin.meshes'  file: ", recvbuf(2)
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
    call openBasisAtomPotentialIndexDAFile(self%atomdata_a(1), 37, 'bin.vpotnew.idx', action='write')

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
    call openRadialMeshDataIndexDAFile(self%mesh_a(1), 37, 'bin.meshes.idx')
#endif
    call openRadialMeshDataDAFile(self%mesh_a(1), 38, 'bin.meshes', max_reclen)

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
      call represent(self%cell_a(ila), str)
      write(*, '(A)') str
    enddo ! ila

  endsubroutine ! represent

endmodule ! CalculationData_mod
