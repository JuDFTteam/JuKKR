module type_corestate
use nrtype
type                                ::  corestate_type
   integer                          ::  ncorestated                ! maximum number of core states (for allocating)
   integer                          ::  ncore                 ! number of core states for atom iatom
   integer                          ::  lcoremax               ! number of core states for atom iatom
   integer,allocatable              ::  lcore(:,:)          ! angular momentum of core state
   real(kind=dp),allocatable        ::  ecore(:,:)          ! energy of core state
end type corestate_type

end module type_corestate