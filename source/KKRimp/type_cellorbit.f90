!------------------------------------------------------------------------------------
!> Summary: Cell spin-orbit coupling type
!> Author: David Bauer
!> Category: KKRimp, geometry, spin-orbit-coupling
!> Deprecated: False 
!> Cell information about the usage of spin-orbit coupling
!------------------------------------------------------------------------------------
module type_cellorbit
use nrtype
 TYPE CELL_TYPEORBIT
 integer,allocatable              :: use_spinorbit(:)       !! spin-orbit coupling used for the 0->no, 1->yes
 END TYPE CELL_TYPEORBIT
end module type_cellorbit
