#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; endif

!> This module provides wrappers for routines with enormous argument lists.
module wrappers_mod
  implicit none

  CONTAINS

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
