!-------------------------------------------------------------------------------
!> Summary: Type holing shape functions and helper index arrays
!> Author: 
!> Category: KKRimp, shape-functions
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module type_shapefun
use nrtype, only: dp
 type :: shapefun_type
   integer,allocatable :: index2lm(:) !! shape function : array index -> lm value
   integer,allocatable :: lmused(:)   !! shape function : 1 -> lm-shape/=0, 0 -> lm-shape=0 
   integer,allocatable :: lm2index(:) !! shape function : lm value -> array index

   real(kind=dp),allocatable :: thetas(:,:) !! shape function
   integer :: nrshaped, nlmshaped
   integer :: nrshape, nlmshape

 end type shapefun_type

end module type_shapefun