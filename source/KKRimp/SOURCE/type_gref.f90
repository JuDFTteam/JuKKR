module type_gref
use nrtype
 TYPE                              ::  GREF_TYPE
   COMPLEX(KIND=DP),allocatable    ::  MAT(:,:)  ! max radius, core radius, muffin tin radius
!    INTEGER                         ::   

 END TYPE GREF_TYPE

end module type_gref