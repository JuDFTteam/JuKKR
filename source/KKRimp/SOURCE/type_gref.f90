!-------------------------------------------------------------------------------
!> Summary: 
!> Author: 
!> Category: KKRimp, 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module type_gref
use nrtype, only: dp
 TYPE                              ::  GREF_TYPE
   COMPLEX(KIND=DP),allocatable    ::  MAT(:,:)  ! max radius, core radius, muffin tin radius
 END TYPE GREF_TYPE

end module type_gref