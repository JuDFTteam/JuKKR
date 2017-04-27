module type_gmatbulk
use nrtype
 type                              ::  gmatbulk_type
   integer                                  :: hostdim
   integer                                  :: lmax
   integer                                  :: lmmax
   integer                                  :: nspin
   integer                                  :: natom
 end type gmatbulk_type

end module type_gmatbulk