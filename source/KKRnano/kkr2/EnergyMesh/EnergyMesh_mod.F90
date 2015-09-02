! TODO: Initialisation

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module EnergyMesh_mod
  implicit none
  private
  public :: EnergyMesh, create, destroy
  public :: createEnergyMesh, destroyEnergyMesh ! deprecated
  public :: readEnergyMesh, readEnergyMeshSemi, writeEnergyMesh, writeEnergyMeshSemi
  public :: broadcastEnergyMesh_com, updateEnergyMesh, updateEnergyMeshSemi
  
  type EnergyMesh

    ! valence contour parameters
    double precision :: E1
    double precision :: E2
    double precision :: EFERMI
    double complex, allocatable :: EZ(:)
    integer :: npnt1
    integer :: npnt2
    integer :: npnt3
    integer :: npol
    double precision :: TK
    double complex, allocatable :: WEZ(:)
    double complex, allocatable :: WEZRN(:,:)
    integer :: ielast

    ! semicore contour parameters
    double precision :: EBOTSEMI
    double precision :: EMUSEMI
    double precision :: FSEMICORE
    integer :: IESEMICORE
    integer :: N1SEMI
    integer :: N2SEMI
    integer :: N3SEMI

  endtype EnergyMesh
  
  interface create
    module procedure createEnergyMesh
  endinterface
  
  interface destroy
    module procedure destroyEnergyMesh
  endinterface

  CONTAINS

  !-----------------------------------------------------------------------------
  !> Constructs a EnergyMesh object.
  !> @param[inout] self    The EnergyMesh object to construct.
  !> @param[in]    ielast
  subroutine createEnergyMesh(self, ielast)
    type (EnergyMesh), intent(inout) :: self
    integer, intent(in) ::  ielast

    integer :: memory_stat

    self%ielast = ielast

    ALLOCATECHECK(self%EZ(ielast))
    ALLOCATECHECK(self%WEZ(ielast))
    ALLOCATECHECK(self%WEZRN(ielast, 2)) ! for both spin channels - historical reasons
  endsubroutine ! create

  !-----------------------------------------------------------------------------
  !> Destroys a EnergyMesh object.
  !> @param[inout] self    The EnergyMesh object to destroy.
  subroutine destroyEnergyMesh(self)
    type (EnergyMesh), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%EZ)
    DEALLOCATECHECK(self%WEZ)
    DEALLOCATECHECK(self%WEZRN)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> read energy mesh data from file 'energy_mesh'
  subroutine readEnergyMesh(emesh)
    use EnergyMeshHelpers_mod, only: readEnergyMeshImpl

    type (EnergyMesh), intent(inout) :: emesh

    call readEnergyMeshImpl(emesh%E1, emesh%E2, emesh%EFERMI, emesh%EZ, &
                            emesh%IELAST, emesh%NPNT1, emesh%NPNT2, emesh%NPNT3, &
                            emesh%NPOL, emesh%TK, emesh%WEZ)

  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMesh(emesh)
    use EnergyMeshHelpers_mod, only: writeEnergyMeshImpl

    type (EnergyMesh), intent(in) :: emesh

    call writeEnergyMeshImpl(emesh%E1, emesh%E2, emesh%EFERMI, emesh%EZ, &
                             emesh%IELAST, emesh%NPNT1, emesh%NPNT2, emesh%NPNT3, &
                             emesh%NPOL, emesh%TK, emesh%WEZ)

  endsubroutine ! write

  !------------------------------------------------------------------------------
  !> Update Energy mesh. Essentially a wrapper for EMESHT
  subroutine updateEnergyMesh(emesh)
    use EnergyMeshHelpers_mod, only: updateEnergyMeshImpl

    type (EnergyMesh), intent(inout) :: emesh

    call updateEnergyMeshImpl(emesh%EZ,emesh%WEZ,emesh%IELAST, &
                              emesh%E1,emesh%E2,emesh%TK,emesh%NPOL, &
                              emesh%NPNT1,emesh%NPNT2,emesh%NPNT3)

  endsubroutine ! update

  !---------------------------------------------------------------------------------
  !> Distribute EnergyMesh from master-rank to all other ranks
  subroutine broadcastEnergyMesh_com(my_mpi, emesh)
    use EnergyMeshHelpers_mod, only: broadcastEnergyMeshImpl_com
    use KKRnanoParallel_mod, only: KKRnanoParallel, getMyActiveCommunicator, getMasterRank

    type (EnergyMesh), intent(inout) :: emesh
    type (KKRnanoParallel), intent(in) :: my_mpi

    call broadcastEnergyMeshImpl_com(getMyActiveCommunicator(my_mpi), getMasterRank(my_mpi), &
         emesh%E1, emesh%E2, emesh%EZ, emesh%IELAST, emesh%WEZ)

  endsubroutine ! broadcast


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! THE FOLLOWING FUNCTIONS ARE USED ONLY IF "use_semicore=1"!!!


  ! VALENCE AND SEMICORE CONTOUR!
  !----------------------------------------------------------------------------
  !> read energy mesh data from file 'energy_mesh'
  subroutine readEnergyMeshSemi(emesh)
    use EnergyMeshHelpers_mod, only: readEnergyMeshImplSemi

    type (EnergyMesh), intent(inout) :: emesh

    call readEnergyMeshImplSemi(emesh%E1, emesh%E2, emesh%EFERMI, emesh%EZ, &
                            emesh%IELAST, emesh%NPNT1, emesh%NPNT2, emesh%NPNT3, &
                            emesh%NPOL, emesh%TK, emesh%WEZ, emesh%EBOTSEMI, emesh%EMUSEMI, &
                            emesh%FSEMICORE, emesh%IESEMICORE, emesh%N1SEMI, emesh%N2SEMI, &
                            emesh%N3SEMI)

  endsubroutine ! read

  ! VALENCE AND SEMICORE CONTOUR!
  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMeshSemi(emesh)
    use EnergyMeshHelpers_mod, only: writeEnergyMeshImplSemi

    type (EnergyMesh), intent(in) :: emesh

    call writeEnergyMeshImplSemi(emesh%E1, emesh%E2, emesh%EFERMI, emesh%EZ, &
                             emesh%IELAST, emesh%NPNT1, emesh%NPNT2, emesh%NPNT3, &
                             emesh%NPOL, emesh%TK, emesh%WEZ, emesh%EBOTSEMI, emesh%EMUSEMI, &
                             emesh%FSEMICORE, emesh%IESEMICORE, emesh%N1SEMI, emesh%N2SEMI, &
                             emesh%N3SEMI)

  endsubroutine ! write

  ! VALENCE AND SEMICORE CONTOUR!
  !----------------------------------------------------------------------------
  !> Update Energy mesh. Essentially a wrapper for EPATHTB
  subroutine updateEnergyMeshSemi(emesh)
    use EnergyMeshHelpers_mod, only: updateEnergyMeshImplSemi

    type (EnergyMesh), intent(inout) :: emesh

    call updateEnergyMeshImplSemi(emesh%EZ,emesh%WEZ,emesh%IELAST, &
                              emesh%E1,emesh%E2,emesh%TK,emesh%NPOL, &
                              emesh%NPNT1,emesh%NPNT2,emesh%NPNT3,emesh%EBOTSEMI, &
                              emesh%EMUSEMI,emesh%IESEMICORE,emesh%FSEMICORE,emesh%N1SEMI, &
                              emesh%N2SEMI,emesh%N3SEMI)

  endsubroutine ! update

endmodule EnergyMesh_mod
