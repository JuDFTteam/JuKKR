!-------------------------------------------------------------------------------
!> Summary: 
!> Author: 
!> Category: KKRimp, 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module type_tmat

 TYPE                              ::  TMAT_TYPE
   DOUBLE COMPLEX,ALLOCATABLE              ::  TMAT(:,:)
   DOUBLE COMPLEX,ALLOCATABLE              ::  deltaT_Jij(:,:,:)

 END TYPE TMAT_TYPE

end module type_tmat