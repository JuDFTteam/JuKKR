!-------------------------------------------------------------------------------
!> Summary: Type holding onsite gmat
!> Author: 
!> Category: KKRimp, structural-greensfunction
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module type_gmatonsite

 type ::  gmatonsite_type
   double complex, allocatable ::  gmat(:,:)
 end type gmatonsite_type

end module type_gmatonsite