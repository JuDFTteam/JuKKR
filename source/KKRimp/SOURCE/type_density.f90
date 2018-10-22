!------------------------------------------------------------------------------------
!> Summary: Density type
!> Author: 
!> Category: KKRimp, physical-observables
!> Deprecated: False 
!> Contains information of the density of each atom like e.g.:
!> the density matrix for both meshes
!> spin and orbital moments 
!> directions of the spin moment, which defines the local frame
!------------------------------------------------------------------------------------
module type_density
use nrtype
 TYPE                              ::  DENSITY_TYPE
   real(kind=dp),allocatable       ::  rho2ns(:,:,:)              !! diagonal of the density matrix (1:nrmax,1:lmmax,2), in energyloop this is always weighted by the energy weight
   double complex,allocatable       ::  rho2ns_complex(:,:,:)     !! complex density matrix (1:nrmax,1:lmmax,nspinden), in energyloop this is always weighted by the energy weight
   double complex,allocatable       ::  rho2ns_complexnew(:,:,:)  !! complex density matrix in the new radial mesh (1:nrmaxnew,1:lmmax,nspinden), in energyloop this is always weighted by the energy weight
   double complex,allocatable       ::  gflle(:,:,:) 		  !! Energy- and lms,lms'-resolved GF lms1,lms2,ie       ! lda+u
   double complex,allocatable       ::  gfint(:,:)   		  !! Integrated GF resolved in lms1,lms2 (integr. grlle) ! lda+u
   double complex        :: rho2ns_integrated(4)                  !! density integrated over r for all for all spin components
   double complex        :: rho2ns_integrated_scattering(4)       !! integrated density for the multiple scattering part of the potential (R G R) used to speed up angular convergence as explained in Bauer, PhD page 113 orbital magnetic moment
   double complex        :: orbitalmom_lm(0:9,3)		  !! orbital moment lm expansion
   double complex        :: orbitalmom_sp(2,3)			  !! spherical part of the orbital moment	
   double complex        :: orbitalmom_ns(3)			  !! non-spherical part of the orbital moment
   real(kind=dp),allocatable   ::   ncharge(:,:) 		  !! integrated denisty of states l-expansion (1:lmax,1:2)
   complex(kind=dpc),allocatable   ::   den(:,:,:)                !! complex density of states l=expansion (1:lmax,1:2,1:ne)
   complex(kind=dpc),allocatable   ::   denlm(:,:,:) 		  !! complex density of states lm-expansion for the diagonals of the density matrix (1:lmmax,1:2,1:ne)
   real(kind=DP),allocatable        ::  RHOC(:,:) 		  !! probably core density
   real(kind=DP)        ::  theta                                 !! polar angle of the magnetic moment towards z-axis, defines the local frame
   real(kind=DP)        ::  phi                                   !! azimuthal angle of the magnetic moment towards z-axis, defines the local frame
   real(kind=DP)        ::  theta2                                !! polar angle from the scattering part 
   real(kind=DP)        ::  phi2                                  !! azimuthal angle from the scattering part
   real(kind=DP)        ::  thetaold                              !! polar angle of prev. iteration
   real(kind=DP)        ::  phiold                                !! azimuthal angle of prev. iteration
   integer              ::  magmomentfixed=1                      !! keep the angle fixed within sc-cycle
   real(kind=DP)        ::  magmoment(3)                          !! magnetic moment in the global frame
   real(kind=DP)        ::  magmoment2(3)                   	  !! magnetic moment calculated from the scattering part
   real(kind=DP)        ::  magmomentold(3)                    	  !! magnetic moment of the previous iteration in the global frame
   real(kind=DP)        ::  magmomentold2(3)                      !! magnetic moment calculated from the scattering part of the previous iteration in the global frame
   integer             :: nangleconfigur    			  !! in mode testflag calctmatfirstIte angle configurations are read in and set in each interationr nangleconfigur are the number of different energy configuarations
 END TYPE DENSITY_TYPE

end module type_density
