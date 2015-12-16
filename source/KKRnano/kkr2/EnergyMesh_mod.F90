! TODO: Initialisation

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if ((STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if ((STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module EnergyMesh_mod
  implicit none
  private
  public :: EnergyMesh, create, destroy, load, store, update, broadcast, getEnergyMeshSize, init
  
  type EnergyMesh

    ! valence contour parameters
    double precision :: e1, e2, eFermi, tK
    integer :: npnt1, npnt2, npnt3, npol, ielast

    ! semicore contour parameters
    double precision :: ebotsemi, emusemi, fsemicore
    integer :: iesemicore, n1semi, n2semi, n3semi

    double complex, allocatable :: ez(:), wez(:), wezrn(:,:)
    integer, allocatable :: kmesh(:) ! mapping to the Brillouin zone sampling TODO, use this instead of that in Main2Arrays

  endtype
  
  interface create
    module procedure createEnergyMesh
  endinterface
  
  interface destroy
    module procedure destroyEnergyMesh
  endinterface

  interface load
    module procedure readEnergyMesh
  endinterface
  
  interface store
    module procedure writeEnergyMesh
  endinterface

  interface update
    module procedure updateEnergyMesh
  endinterface

  interface broadcast
    module procedure broadcastEnergyMesh_com
  endinterface

  interface init
    module procedure initEnergyMesh
  endinterface
  
  contains

  
  integer function getEnergyMeshSize(npol, n1val, n2val, n3val, n1sem, n2sem, n3sem) result(iemxd)
    integer, intent(in) :: npol, n1val, n2val, n3val, n1sem, n2sem, n3sem
  
    ! Calculate number of energy points
    if (n1sem > 0 .or. n2sem > 0 .or. n3sem > 0) then ! semicore contour is used

      write(*,*) "WARNING: Using semicore contour, please check your results"
      
      if (npol /= 0) then
        iemxd = npol + n1val + n2val + n3val + n1sem + n2sem + n3sem
      else ! dos-calculation
        iemxd = n2val + n2sem
      endif

    else ! semicore contour is not used

      if (npol /= 0) then
        iemxd = npol + n1val + n2val + n3val
      else ! dos-calculation
        iemxd = n2val
      endif

    endif ! semicore

  endfunction ! get

  
  subroutine initEnergyMesh(self, eFermi, emin, emax, tK, npol, npnt1, npnt2, npnt3, ebotsemi, emusemi, n1semi, n2semi, n3semi, fsemicore)
    type(EnergyMesh), intent(inout) :: self
    double precision, intent(in) :: eFermi, emin, emax, tK
    integer, intent(in) :: npol, npnt1, npnt2, npnt3
    double precision, intent(in) :: ebotsemi, emusemi
    integer, intent(in) :: n1semi, n2semi, n3semi
    double precision, intent(in) :: fsemicore

    self%eFermi = eFermi
    self%tK = tK
    self%e1 = emin
    self%e2 = emax
    self%npol = npol
    self%npnt1 = npnt1; self%npnt2 = npnt2; self%npnt3 = npnt3
    self%ebotsemi = ebotsemi
    self%emusemi = emusemi
    self%n1semi = n1semi; self%n2semi = n2semi; self%n3semi = n3semi
    self%fsemicore = fsemicore

  endsubroutine ! init
  
  !-----------------------------------------------------------------------------
  !> Constructs a EnergyMesh object.
  !> @param[inout] self    The EnergyMesh object to construct.
  !> @param[in]    ielast
  subroutine createEnergyMesh(self, ielast)
    type(EnergyMesh), intent(inout) :: self
    integer, intent(in) ::  ielast

    integer :: memory_stat

    self%ielast = ielast

    ALLOCATECHECK(self%ez(ielast))
    ALLOCATECHECK(self%wez(ielast))
    ALLOCATECHECK(self%wezrn(ielast,2)) ! for both spin channels - historical reasons
    ALLOCATECHECK(self%kmesh(ielast))
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a EnergyMesh object.
  !> @param[inout] self    The EnergyMesh object to destroy.
  elemental subroutine destroyEnergyMesh(self) ! cannot be elemental when using IO in allocatecheck
    type(EnergyMesh), intent(inout) :: self
    integer :: ist
    deallocate(self%ez, self%wez, self%wezrn, self%kmesh, stat=ist)
  endsubroutine ! destroy

  !---------------------------------------------------------------------------------
  !> Distribute EnergyMesh from master-rank to all other ranks
  subroutine broadcastEnergyMesh_com(self, comm)
    use EnergyMeshHelpers_mod, only: broadcast
    type(EnergyMesh), intent(inout) :: self
    integer, intent(in) :: comm

    call broadcast(comm, self%e1, self%e2, self%ebotsemi, self%emusemi, self%IELAST, self%ez, self%wez)

  endsubroutine ! broadcast

  !----------------------------------------------------------------------------
  !> read energy mesh data from file 'energy_mesh'
  subroutine readEnergyMesh(self)
    use EnergyMeshHelpers_mod, only: load

    type(EnergyMesh), intent(inout) :: self

    call load(self%e1, self%e2, self%eFermi, self%ez, &
              self%IELAST, self%NPNT1, self%NPNT2, self%NPNT3, &
              self%NPOL, self%tK, self%wez, self%ebotsemi, self%emusemi, &
              self%fsemicore, self%iesemicore, self%n1semi, self%n2semi, &
              self%n3semi, self%kmesh, filename='energy_mesh.0')

  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMesh(self)
    use EnergyMeshHelpers_mod, only: store

    type(EnergyMesh), intent(in) :: self

    call store(self%e1, self%e2, self%eFermi, self%ez, &
              self%IELAST, self%NPNT1, self%NPNT2, self%NPNT3, &
              self%NPOL, self%tK, self%wez, self%ebotsemi, self%emusemi, &
              self%fsemicore, self%iesemicore, self%n1semi, self%n2semi, &
              self%n3semi, self%kmesh, filename='energy_mesh')

  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Update Energy mesh. Essentially a wrapper for EPATHTB
  subroutine updateEnergyMesh(self)
    use EnergyMeshHelpers_mod, only: update

    type(EnergyMesh), intent(inout) :: self

    call update(self%ez, self%wez, self%IELAST, &
              self%e1, self%e2, self%tK, self%NPOL, &
              self%NPNT1, self%NPNT2, self%NPNT3, &
              self%ebotsemi, self%emusemi, self%iesemicore, self%fsemicore, &
              self%n1semi, self%n2semi, self%n3semi)

  endsubroutine ! update

endmodule EnergyMesh_mod
