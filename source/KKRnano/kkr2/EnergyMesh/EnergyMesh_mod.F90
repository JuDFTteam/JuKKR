! TODO: Initialisation

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module EnergyMesh_mod

  type EnergyMesh
    double precision  :: E1
    double precision  :: E2
    double precision  :: EFERMI
    double complex , allocatable, dimension(:)  :: EZ
    integer  :: npnt1
    integer  :: npnt2
    integer  :: npnt3
    integer  :: npol
    double precision  :: TK
    double complex , allocatable, dimension(:)  :: WEZ
    double complex , allocatable, dimension(:,:)  :: WEZRN

    integer :: ielast
  end type EnergyMesh

  CONTAINS

  !-----------------------------------------------------------------------------
  !> Constructs a EnergyMesh object.
  !> @param[inout] self    The EnergyMesh object to construct.
  !> @param[in]    ielast
  subroutine createEnergyMesh(self, ielast)
    implicit none
    type (EnergyMesh), intent(inout) :: self
    integer, intent(in) ::  ielast

    integer :: memory_stat

    self%ielast = ielast

    ALLOCATECHECK(self%EZ(ielast))
    ALLOCATECHECK(self%WEZ(ielast))
    ALLOCATECHECK(self%WEZRN(ielast, 2)) ! for both spin channels - historical reasons
  end subroutine

  !-----------------------------------------------------------------------------
  !> Destroys a EnergyMesh object.
  !> @param[inout] self    The EnergyMesh object to destroy.
  subroutine destroyEnergyMesh(self)
    implicit none
    type (EnergyMesh), intent(inout) :: self

    integer :: memory_stat

    DEALLOCATECHECK(self%EZ)
    DEALLOCATECHECK(self%WEZ)
    DEALLOCATECHECK(self%WEZRN)
  end subroutine

  !----------------------------------------------------------------------------
  !> read energy mesh data from file 'energy_mesh'
  subroutine readEnergyMesh(emesh)
    use EnergyMeshHelpers_mod
    implicit none

    type (EnergyMesh), intent(inout) :: emesh

    call readEnergyMeshImpl(emesh%E1, emesh%E2, emesh%EFERMI, emesh%EZ, &
                            emesh%IELAST, emesh%NPNT1, emesh%NPNT2, emesh%NPNT3, &
                            emesh%NPOL, emesh%TK, emesh%WEZ)

  end subroutine

  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMesh(emesh)
    use EnergyMeshHelpers_mod
    implicit none

    type (EnergyMesh), intent(in) :: emesh

    call writeEnergyMeshImpl(emesh%E1, emesh%E2, emesh%EFERMI, emesh%EZ, &
    emesh%IELAST, emesh%NPNT1, emesh%NPNT2, emesh%NPNT3, &
    emesh%NPOL, emesh%TK, emesh%WEZ)

  end subroutine

!------------------------------------------------------------------------------
!> Update Energy mesh. Essentially a wrapper for EMESHT
subroutine updateEnergyMesh(emesh)
    use EnergyMeshHelpers_mod
    implicit none

    type (EnergyMesh), intent(inout) :: emesh

    call updateEnergyMeshImpl(emesh%EZ,emesh%WEZ,emesh%IELAST, &
                              emesh%E1,emesh%E2,emesh%TK,emesh%NPOL, &
                              emesh%NPNT1,emesh%NPNT2,emesh%NPNT3)

  end subroutine

  !---------------------------------------------------------------------------------
  !> Distribute EnergyMesh from master-rank to all other ranks
  subroutine broadcastEnergyMesh_com(my_mpi, emesh)
    use EnergyMeshHelpers_mod
    use KKRnanoParallel_mod
    implicit none

    type (EnergyMesh), intent(inout) :: emesh
    type (KKRnanoParallel), intent(in) :: my_mpi

    call broadcastEnergyMeshImpl_com(getMyActiveCommunicator(my_mpi), getMasterRank(my_mpi), &
         emesh%E1, emesh%E2, emesh%EZ, emesh%IELAST, emesh%WEZ)

  end subroutine


end module
