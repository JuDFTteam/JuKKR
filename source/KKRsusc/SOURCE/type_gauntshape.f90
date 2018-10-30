module type_gauntshape
use nrtype
  type                               :: gauntshape_type
     integer,allocatable             :: ilm(:,:)!(NGSHD,3)
     real*8,allocatable              :: gsh(:) !(NGSHD) 
     integer,allocatable             :: imaxsh(:) !(0:LMPOTD)
      integer                        :: NGSHD=13079
  end type gauntshape_type

end module type_gauntshape