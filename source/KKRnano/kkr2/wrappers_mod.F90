#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; endif

!------------------------------------------------------------------------------
!> This module provides wrappers for routines with enormous argument lists.
!>
!> Most of them are called in main2.
module wrappers_mod
  implicit none

  CONTAINS

!------------------------------------------------------------------------------
subroutine ESPCB_wrapper(ESPC, LCOREMAX, atomdata)
  use BasisAtom_mod
  implicit none

  double precision, intent(out) :: ESPC(:,:) !ESPC(0:3,NSPIN)
  integer, intent(out) :: LCOREMAX
  type (BasisAtom), intent(inout) :: atomdata

  !-------- locals
  integer :: nspind

  nspind = atomdata%nspin

  call ESPCB_NEW(ESPC,NSPIND,atomdata%core%ECORE,atomdata%core%LCORE(:,1:NSPIND),LCOREMAX,atomdata%core%NCORE(1:NSPIND))

end subroutine

!------------------------------------------------------------------------------
subroutine EPOTINB_wrapper(EPOTIN,RHO2NS,atomdata)
  use BasisAtom_mod
  use RadialMeshData_mod
  implicit none

  double precision, intent(out)   :: EPOTIN
  double precision, intent(inout) :: RHO2NS(:,:,:)
  type (BasisAtom), intent(inout) :: atomdata

  !-------- locals
  integer :: nspind
  integer :: irnsd
  type (RadialMeshData), pointer :: mesh

  nspind = atomdata%nspin

  mesh => atomdata%mesh_ptr

  CHECKASSERT( associated(mesh) )

  irnsd = atomdata%potential%irmd - atomdata%potential%irmind

  CHECKASSERT( atomdata%potential%irmd == mesh%irmd )

  call EPOTINB_NEW(EPOTIN,NSPIND,RHO2NS,atomdata%potential%VISP,mesh%R,mesh%DRDI, &
                   mesh%IRMIN,mesh%IRWS,atomdata%potential%LPOT,atomdata%potential%VINS, &
                   mesh%IRCUT,mesh%IPAN,atomdata%Z_nuclear, &
                   mesh%irmd, irnsd, mesh%ipand)

end subroutine

!----------------------------------------------------------------------------
subroutine ECOUB_wrapper(CMOM, ECOU, RHO2NS, shgaunts, atomdata)
  use BasisAtom_mod
  use RadialMeshData_mod
  use CellData_mod
  use ShapeGauntCoefficients_mod
  implicit none
  double precision, intent(inout) :: CMOM(:)
  double precision, intent(inout) :: ECOU(:)
  double precision, intent(inout) :: RHO2NS(:,:,:)
  type (BasisAtom), intent(inout) :: atomdata
  type (ShapeGauntCoefficients), intent(in) :: shgaunts

  !-------- locals
  integer :: nspind
  type (RadialMeshData), pointer :: mesh
  type (CellData), pointer       :: cell
  integer :: KVMAD

  nspind = atomdata%nspin

  mesh => atomdata%mesh_ptr
  cell => atomdata%cell_ptr

  CHECKASSERT( associated(mesh) )
  CHECKASSERT( associated(cell) )

  KVMAD = 0

  ! output: ECOU - l resolved Coulomb energy
  call ECOUB_NEW(CMOM,ECOU,atomdata%potential%LPOT,NSPIND,RHO2NS, &
  atomdata%potential%VONS,atomdata%Z_nuclear,mesh%R, &
  mesh%DRDI,KVMAD,mesh%IRCUT,mesh%IPAN,shgaunts%IMAXSH,cell%shdata%IFUNM, &
  shgaunts%ILM,shgaunts%GSH,cell%shdata%THETA,cell%shdata%LMSP, &
  mesh%irmd, cell%shdata%irid, cell%shdata%nfund, mesh%ipand, shgaunts%ngshd)


end subroutine

!------------------------------------------------------------------------------
!> Add exchange-correlation potential and calculate exchange correlation energy.
subroutine VXCDRV_wrapper(EXC,KXC,RHO2NS, shgaunts, atomdata)
  use BasisAtom_mod
  use RadialMeshData_mod
  use CellData_mod
  use ShapeGauntCoefficients_mod
  implicit none
  double precision, intent(inout) :: EXC(:)
  integer, intent(in)             :: KXC
  double precision, intent(inout) :: RHO2NS(:,:,:)  ! inout?
  type (BasisAtom), intent(inout) :: atomdata
  type (ShapeGauntCoefficients), intent(in) :: shgaunts

  !-------- locals
  integer :: nspind
  type (RadialMeshData), pointer :: mesh
  type (CellData), pointer       :: cell
  integer, parameter :: KTE = 1 ! always calculate exchange energy
  integer :: LPOT

  nspind = atomdata%nspin
  LPOT   = atomdata%potential%lpot

  mesh => atomdata%mesh_ptr
  cell => atomdata%cell_ptr

  CHECKASSERT( associated(mesh) )
  CHECKASSERT( associated(cell) )
  CHECKASSERT( size(EXC) == LPOT+1)

  call VXCDRV_NEW(EXC,KTE,KXC,LPOT,NSPIND,RHO2NS, &
            atomdata%potential%VONS,mesh%R,mesh%DRDI,mesh%A, &
            mesh%IRWS,mesh%IRCUT,mesh%IPAN,shgaunts%GSH,shgaunts%ILM,shgaunts%IMAXSH,cell%shdata%IFUNM, &
            cell%shdata%THETA,cell%shdata%LMSP, &
            mesh%irmd, cell%shdata%irid, cell%shdata%nfund, shgaunts%ngshd, mesh%ipand)

end subroutine

!----------------------------------------------------------------------------
!> initialise VAV0, VOL0 to 0.0d0 before calling!!!
subroutine MTZERO_wrapper(VAV0, VOL0, atomdata)
  use BasisAtom_mod
  use RadialMeshData_mod
  use CellData_mod

  implicit none
  type (BasisAtom), intent(inout) :: atomdata
  double precision, intent(in)    :: VAV0, VOL0

  !-------- locals
  integer :: nspind
  type (RadialMeshData), pointer :: mesh
  type (CellData), pointer       :: cell

  nspind = atomdata%nspin

  mesh => atomdata%mesh_ptr
  cell => atomdata%cell_ptr

  CHECKASSERT( associated(mesh) )
  CHECKASSERT( associated(cell) )

  !output: VAV0, VOL0
  call MTZERO_NEW(atomdata%potential%LMPOT,NSPIND,atomdata%potential%VONS, &
                  atomdata%Z_nuclear,mesh%R,mesh%DRDI,mesh%IMT,mesh%IRCUT, &
                  mesh%IPAN,cell%shdata%LMSP,cell%shdata%IFUNM, &
                  cell%shdata%THETA,mesh%IRWS,VAV0,VOL0, &
                  mesh%irmd, cell%shdata%irid, cell%shdata%nfund, mesh%ipand)


end subroutine


!----------------------------------------------------------------------------
subroutine CONVOL_wrapper(VBC, shgaunts, atomdata)
  use BasisAtom_mod
  use RadialMeshData_mod
  use CellData_mod
  use ShapeGauntCoefficients_mod
  implicit none
  type (BasisAtom), intent(inout) :: atomdata
  type (ShapeGauntCoefficients), intent(in) :: shgaunts
  double precision, intent(in)    :: VBC(2)

  !-------- locals
  integer :: ispin, nspind
  type (RadialMeshData), pointer :: mesh
  type (CellData), pointer       :: cell

  nspind = atomdata%nspin

  mesh => atomdata%mesh_ptr
  cell => atomdata%cell_ptr

  CHECKASSERT( associated(mesh) )
  CHECKASSERT( associated(cell) )

  do ispin = 1, nspind

    call shiftPotential(atomdata%potential%VONS(:,:,ISPIN), mesh%IRCUT(mesh%IPAN), VBC(ISPIN))

    !output: VONS (changed)
    call CONVOL_NEW(mesh%IRCUT(1),mesh%IRC, &
                    shgaunts%IMAXSH(shgaunts%LMPOTD),shgaunts%ILM,cell%shdata%IFUNM,shgaunts%LMPOTD,shgaunts%GSH, &
                    cell%shdata%THETA,atomdata%Z_nuclear, &
                    mesh%R,atomdata%potential%VONS(1,1,ISPIN),cell%shdata%LMSP, &
                    cell%shdata%irid, cell%shdata%nfund, mesh%irmd, shgaunts%ngshd)

  end do

end subroutine


!------------------------------------------------------------------------------
!> Do core relaxation for all spin directions.
!>
!> @param[in]      E1       bottom energy of valence energy contour
!> @param[in]     NSRA      flag for scalar relativistic calculation
!> @param[in,out] atomdata  basis atom - changed on output
subroutine RHOCORE_wrapper(E1, NSRA, atomdata)
  use BasisAtom_mod
  use RadialMeshData_mod
  implicit none

  double precision, intent(in)    :: E1
  integer, intent(in)             :: NSRA
  type (BasisAtom), intent(inout) :: atomdata

  !-------- locals
  integer :: ispin, nspind
  type (RadialMeshData), pointer :: mesh_ptr

  nspind = atomdata%nspin

  mesh_ptr => atomdata%mesh_ptr

  CHECKASSERT( associated(mesh_ptr) )

  do ispin = 1, nspind

    ! output: ECORE, NCORE, LCORE, RHOCAT?, QC
    call RHOCORE(E1,NSRA,ISPIN,NSPIND,atomdata%atom_index, &  ! atom_index is used only for debugging output
                 mesh_ptr%DRDI,mesh_ptr%R,atomdata%potential%VISP(:,ISPIN), &
                 mesh_ptr%A,mesh_ptr%B,atomdata%Z_nuclear, &
                 mesh_ptr%IRCUT,atomdata%core%RHOCAT,atomdata%core%QC_corecharge, &
                 atomdata%core%ECORE(:,ISPIN),atomdata%core%NCORE(ispin),atomdata%core%LCORE(:,ispin), &
                 mesh_ptr%irmd, mesh_ptr%ipand)

  end do

end subroutine

!-------------------------------------------------------------------------------
!> A wrapper for the subroutine RHOMOM_NEW.
subroutine RHOMOM_NEW_wrapper(CMOM,CMINST,RHO2NS, cell, mesh, shgaunts)

  use CellData_mod
  use RadialMeshData_mod
  use ShapeGauntCoefficients_mod
  implicit none

  type (CellData), intent(in) :: cell
  type (RadialMeshData), intent(in) :: mesh
  type (ShapeGauntCoefficients), intent(in) :: shgaunts

  DOUBLE PRECISION, dimension(:) :: CMINST
  DOUBLE PRECISION, dimension(:) :: CMOM
  DOUBLE PRECISION, dimension(:,:) :: RHO2NS

  ! locals
  integer :: lpot

  lpot = shgaunts%lmax * 2

  CHECKASSERT( size(CMOM) == (LPOT+1)**2)
  CHECKASSERT( size(CMINST) == (LPOT+1)**2)
  CHECKASSERT( size(RHO2NS, 1) == mesh%irmd)
  CHECKASSERT( size(RHO2NS, 2) == (LPOT+1)**2)
  CHECKASSERT( cell%shdata%lmmax_shape == (2*LPOT + 1) ** 2)

  !call RHOMOM_NEW(CMOM,CMINST,LPOT,RHO2NS,R,DRDI,IRCUT,IPAN,ILM,IFUNM,IMAXSH,GSH,THETAS,LMSP,irmd,irid,nfund,ipand,ngshd)

  !output: CMOM, CMINST
  call RHOMOM_NEW(CMOM,CMINST,LPOT,RHO2NS, &
            mesh%R,mesh%DRDI,mesh%IRCUT,mesh%IPAN,shgaunts%ILM, &
            cell%shdata%IFUNM,shgaunts%IMAXSH,shgaunts%GSH, &
            cell%shdata%THETA,cell%shdata%LMSP, &
            mesh%irmd, cell%shdata%irid, cell%shdata%nfund, mesh%ipand, shgaunts%ngshd)

end subroutine

!==============================================================================
! Some helper routines for the wrappers
!==============================================================================

  !----------------------------------------------------------------------------
  !> Shift lm-decomposed potential be a constant sqrt(4 * pi) * VBC.
  subroutine shiftPotential(VONS_ISPIN, index_rmax, VBC)
    implicit none
    double precision, dimension(:,:), intent(inout) :: VONS_ISPIN
    integer, intent(in) :: index_rmax
    double precision, intent(in) :: VBC
    !----------------------------------

    integer :: IR
    double precision :: RFPI

    !do IR = 1,IRCUT(IPAN(I1),I1)
    !  VONS_ISPIN(IR,1,ISPIN) = VONS_ISPIN(IR,1,ISPIN) + RFPI*VBC(ISPIN)
    !end do

    RFPI = SQRT(16.0D0*ATAN(1.0D0))

    do IR = 1, index_rmax
      ! a constant potential shift affects only lm=1 component
      VONS_ISPIN(IR,1) = VONS_ISPIN(IR,1) + RFPI*VBC
    end do

  end subroutine

end module wrappers_mod
