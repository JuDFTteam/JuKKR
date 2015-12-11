!
! cell information using intervals containing expansions in Chebyshev polynomials
!
module type_cellorbit
use nrtype
 TYPE CELL_TYPEORBIT
 integer,allocatable              :: use_spinorbit(:)       ! spin-orbit coupling used for the 0->no, 1->yes
 END TYPE CELL_TYPEORBIT
end module type_cellorbit
