module type_gmat
use nrtype
 type                              ::  gmat_type
   double complex,allocatable               ::  gmathost(:,:)
   integer                                  :: gmathostdim
   integer                                  :: gmathost_lmmax
   integer                                  ::  kgrefsoc
   double complex,allocatable               ::  gmat(:,:)
   integer                                  ::  gmatdim
   integer,allocatable                      ::  iatom2nlmindex(:,:)
   integer,allocatable                      ::  iatom2nlmindexhost(:,:)
   integer,allocatable                      ::  nlmindex2iatom(:)

 end type gmat_type

end module type_gmat