module type_gauntcoeff
use nrtype
  type                               :: gauntcoeff_type
     integer,allocatable             :: icleb(:,:)           !: pointer array
     integer,allocatable             :: loflm(:)              !: l of lm=(l,m) (gaunt)
     real(kind=dp),allocatable       :: cleb(:,:)            !: gaunt coefficients (gaunt)
     integer                         :: iend                     !: number of nonzero gaunt coeffizients
     integer,allocatable             :: jend(:,:,:)              !: pointer array for icleb()
     integer                         :: ncleb                     
     real(kind=dp),allocatable       :: wg(:)
     real(kind=dp),allocatable       :: yrg(:,:,:)
  end type gauntcoeff_type


end module type_gauntcoeff
