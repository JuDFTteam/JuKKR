module type_density
use nrtype
 TYPE                              ::  DENSITY_TYPE
   real(kind=dp),allocatable       ::  rho2ns(:,:,:)              ! density
   double complex,allocatable       ::  rho2ns_complex(:,:,:)     ! density including imag. part
   double complex,allocatable       ::  rho2ns_complexnew(:,:,:)  ! same for new mesh
   double complex,allocatable       ::  gflle(:,:,:) ! Energy- and lms,lms'-resolved GF lms1,lms2,ie       ! lda+u
   double complex,allocatable       ::  gfint(:,:)   ! Integrated GF resolved in lms1,lms2 (integr. grlle) ! lda+u
   double complex        :: rho2ns_integrated(4)                  ! density integrated over r for all
                                                                  ! for all spin components
   double complex        :: rho2ns_integrated_scattering(4)       ! integrated density for the multiple
                                                                  ! scattering part of the potential (R G R)
                                                                  ! used to speed up angular convergence
                                                                  ! as explained in Bauer, PhD page 113
   double complex        :: orbitalmom(3)                         ! orbital magnetic moment
   double complex        :: orbitalmom_lm(0:9,3)
   double complex        :: orbitalmom_sp(2,3)
   double complex        :: orbitalmom_ns(3)
   real(kind=dp),allocatable   ::   ncharge(:,:) 
   complex(kind=dpc),allocatable   ::   den(:,:,:)                 ! density of states 
   complex(kind=dpc),allocatable   ::   denlm(:,:,:) 
   real(kind=DP),allocatable        ::  RHOC(:,:) 

   real(kind=DP)        ::  theta, phi !(CELL%NRMAX,NSPIN)         ! angle of the mag. moment towards z-axis
   real(kind=DP)        ::  theta2, phi2 !(CELL%NRMAX,NSPIN)       !
   real(kind=DP)        ::  thetaold, phiold !(CELL%NRMAX,NSPIN)   ! theta of prev. iteration
   integer              :: magmomentfixed=1                       ! keep the angle fixed within sc-cycle
   real(kind=DP)        ::  magmoment(3) !(CELL%NRMAX,NSPIN)        ! magnetic moment
   real(kind=DP)        ::  magmoment2(3) !(CELL%NRMAX,NSPIN)
   real(kind=DP)        ::  magmomentold(3) !(CELL%NRMAX,NSPIN)
   real(kind=DP)        ::  magmomentold2(3) !(CELL%NRMAX,NSPIN)
   integer             :: nangleconfigur    ! in mode testflag calctmatfirstIter
                                         ! angle configurations are read in and set in each interation
                                         ! nangleconfigur are the number of different energy configuarations


 END TYPE DENSITY_TYPE

end module type_density