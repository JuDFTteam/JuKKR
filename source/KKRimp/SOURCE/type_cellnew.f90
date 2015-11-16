!
! cell information using intervals containing expansions in Chebyshev polynomials
!
module type_cellnew
use nrtype
 TYPE                             ::  CELL_TYPENEW
 real(kind=dp),allocatable        ::  vpotnew(:,:,:)       ! potential for the Cheb. mesh
 real(kind=dp),allocatable        ::  shapefun(:,:)        ! shape function for Cheb. mesh
 integer,allocatable              ::  shapefun_lm2index(:) ! shape function : lm value -> array index
 !integer                          :: use_spinorbit=1       ! spin-orbit coupling used for the 0->no, 1->yes
 integer                          :: use_spinorbit         ! spin-orbit coupling used for the 0->no, 1->yes
 integer                          :: ncheb                 ! maximum number of Chebyshev expansion fn used
 real(kind=dp),allocatable        :: rpan_intervall(:)     ! larger boundary value for panel m 
                                                           ! rpan_intervall(0) smaller boundary value for 1st panel
 integer,allocatable              :: ipan_intervall(:)     ! larger index for the potential and shape fn array
                                                           ! vpotnew(ipan_intervall(1)) last potential value for
                                                           ! first panel
 real(kind=dp),allocatable        :: rmeshnew(:)           ! radial mesh containing mesh points for each panel
 integer                          :: npan_log,    &        ! number of panel in the logarithmic region
                                     npan_eq,     &        !                        equidistant region
                                     npan_inst,   &        !                        interstitial region 
                                     npan_tot              ! total number of panel (npan_* summed up)
 integer                          :: nrmaxnew              ! number of mesh points
 END TYPE CELL_TYPENEW

end module type_cellnew
