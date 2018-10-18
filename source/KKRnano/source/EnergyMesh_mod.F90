!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module EnergyMesh_mod
!-------------------------------------------------------------------------------
!> Summary: Point mesh for the energy contour integration
!> Author: Phivos Mavropoulos, Voicu Popescu, Elias Rabel, Marcel Bornemann, Paul F Baumeister
!> Category: KKRnano, input-output, initialization, communcation
!-------------------------------------------------------------------------------
  implicit none
! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine they are used.
#define CHECKALLOC(STAT) if ((STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if ((STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)

  private
  public :: EnergyMesh, create, destroy, load, store, update, broadcast, getEnergyMeshSize, init
  
  type EnergyMesh

    ! valence contour parameters
    double precision :: e1, e2, eFermi, tK
    integer :: npnt123(3), npol, ielast

    ! semicore contour parameters
    double precision :: ebotsemi, emusemi, fsemicore
    integer :: iesemicore, npntsemi(3)

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

  
  !! Calculate number of energy points
  integer function getEnergyMeshSize(npol, n123val, n123sem) result(iemxd)
    integer, intent(in) :: npol, n123val(3), n123sem(3)
  
    if (any(n123sem > 0)) write(*,*) "WARNING: Using semicore contour, please check your results"

    if (npol /= 0) then
      iemxd = max(0, npol) + sum(max(0, n123val)) + sum(max(0, n123sem))
    else ! dos-calculation
      iemxd = max(0, n123val(2)) + max(0, n123sem(2))
    endif

  endfunction ! get

  
  subroutine initEnergyMesh(self, eFermi, emin, emax, tK, npol, npnt123, ebotsemi, emusemi, npntsemi, fsemicore)
    type(EnergyMesh), intent(inout) :: self
    double precision, intent(in) :: eFermi, emin, emax, tK
    integer, intent(in) :: npol, npnt123(3)
    double precision, intent(in) :: ebotsemi, emusemi
    integer, intent(in) :: npntsemi(3)
    double precision, intent(in) :: fsemicore

    self%eFermi = eFermi
    self%tK = tK
    self%e1 = emin
    self%e2 = emax
    self%npol = npol
    self%npnt123 = npnt123
    self%ebotsemi = ebotsemi
    self%emusemi = emusemi
    self%npntsemi = npntsemi
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
    integer :: ist ! ignore status
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
  subroutine readEnergyMesh(self, filename)
    use EnergyMeshHelpers_mod, only: load
    type(EnergyMesh), intent(inout) :: self
    character(len=*), intent(in) :: filename ! usually 'energy_mesh'

    call load(self%e1, self%e2, self%eFermi, self%ez, &
              self%IELAST, self%npnt123, &
              self%NPOL, self%tK, self%wez, self%ebotsemi, self%emusemi, &
              self%fsemicore, self%iesemicore, self%npntsemi, self%kmesh, filename)

  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMesh(self, filename)
    use EnergyMeshHelpers_mod, only: store

    type(EnergyMesh), intent(in) :: self
    character(len=*), intent(in) :: filename ! 'energy_mesh.0' or 'energy_mesh'

    call store(self%e1, self%e2, self%eFermi, self%ez, &
              self%IELAST, self%npnt123, &
              self%NPOL, self%tK, self%wez, self%ebotsemi, self%emusemi, &
              self%fsemicore, self%iesemicore, self%npntsemi, self%kmesh, filename)

  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Update Energy mesh. Essentially a wrapper for EPATHTB
  subroutine updateEnergyMesh(self)
    use EnergyMeshHelpers_mod, only: update

    type(EnergyMesh), intent(inout) :: self

    call update(self%ez, self%wez, self%IELAST, &
              self%e1, self%e2, self%tK, self%NPOL, self%npnt123, &
              self%ebotsemi, self%emusemi, self%iesemicore, self%fsemicore, self%npntsemi)

  endsubroutine ! update

endmodule ! EnergyMesh_mod
