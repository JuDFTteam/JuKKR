!------------------------------------------------------------------------------------
!> Summary: Green function matrix type
!> Author: 
!> Category: KKRimp 
!> Deprecated: False 
!> Contains informations about the Green function matrix of the host and the Green function of the impurity cluster
!------------------------------------------------------------------------------------
module type_gmat
use nrtype
 type                              ::  gmat_type
   double complex,allocatable               ::  gmathost(:,:)  				!! host Green function matrix of size (ntotatom*lmsizehost, ntotatom*lmsizehost)
   integer                                  :: gmathostdim               		!! dimension of the host Green function 
   integer                                  :: gmathost_lmmax            		!! maximal angular momentum lmmax of the host
   integer                                  ::  kgrefsoc		 		!! spin-orbit coupling in the host
   double complex,allocatable               ::  gmat(:,:)		 		!! Green function matrix
   integer                                  ::  gmatdim			 		!! dimension of the Green funtion matrix
   integer,allocatable                      ::  iatom2nlmindex(:,:)	 		!! some pointer array		
   integer,allocatable                      ::  iatom2nlmindexhost(:,:)	 		!! some pointer array
   integer,allocatable                      ::  nlmindex2iatom(:)	 		!! some pointer array

 end type gmat_type

end module type_gmat
