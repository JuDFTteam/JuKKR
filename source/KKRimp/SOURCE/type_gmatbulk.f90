!------------------------------------------------------------------------------------
!> Summary: Type for the Green function matrix of the host 
!> Author: 
!> Category: KKRimp
!> Deprecated: False 
!> Contains information about the host Green function
!------------------------------------------------------------------------------------
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
