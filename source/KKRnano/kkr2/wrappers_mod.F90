#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; endif

!> This module provides wrappers for routines with enormous argument lists.
module wrappers_mod
  implicit none

  CONTAINS

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

end module wrappers_mod
