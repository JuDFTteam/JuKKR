!-------------------------------------------------------------------------------
!> Summary: 
!> Author: 
!> Category: KKRimp, 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module type_shapefun
use nrtype, only: dp
 TYPE :: SHAPEFUN_TYPE
   integer,allocatable :: index2lm(:) !! shape function : array index -> lm value
   integer,allocatable :: lmused(:) !! shape function : 1 -> lm-shape/=0, 0 -> lm-shape=0 
   integer,allocatable :: lm2index(:) !! shape function : lm value -> array index

   real(kind=dp),allocatable :: thetas(:,:) !! shape function
   integer :: nrshaped, nlmshaped
   INTEGER :: NRSHAPE, NLMSHAPE

 END TYPE SHAPEFUN_TYPE

end module type_shapefun