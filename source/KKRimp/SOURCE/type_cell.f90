!------------------------------------------------------------------------------------
!> Summary: Cell type
!> Author: 
!> Category: KKRimp, geometry, old-mesh
!> Deprecated: False 
!> Contains information of each individual cell in the old radial mesh.
!------------------------------------------------------------------------------------
module type_cell
use nrtype
 TYPE                               ::  CELL_TYPE
   REAL(KIND=DP)                    ::  RMAX                      !! maximal radius
   REAL(KIND=DP)                    ::  RCORE                     !! core radius
   REAL(KIND=DP)                    ::  RMT                       !! muffin tin radius
   INTEGER                          ::  NRMAX        		  !! total radial mesh points
   INTEGER                          ::  NRCORE        		  !! core gp
   INTEGER                          ::  NRNS        		  !! non-spherical gp
   INTEGER                          ::  NRMAXD                    !! maximum number of radial point for all cells
   INTEGER,ALLOCATABLE              ::  NRCUT(:)                  !! array index of the last position of each panel
   REAL(KIND=DP)                    ::  LOGPARAMS(2)              !! a,b for logarithmic mesh : r(i) = b*(exp(a*(i-1))-1)
   REAL(KIND=DP),ALLOCATABLE        ::  RMESH(:)                  !! logarithmic mesh
   REAL(KIND=DP),ALLOCATABLE        ::  DRMESHDI(:)               !! dR/dI spacing between the radial mesh points
   REAL(KIND=DP),ALLOCATABLE        ::  DRMESHOR(:)               !! dR/R 
   INTEGER                          ::  NRMIN_NS                  !! index of the most inner data for a non-spherical treatment
   INTEGER                          ::  KXC                       !! type of exchange correlation functional 

   INTEGER                          ::  NPAN 	                  !! Number of panels
   INTEGER                          ::  NPAND 			  !! Maximum number of panels for all cells
   INTEGER,ALLOCATABLE              ::  NMESHPAN(:)	          !! probably number of radial points in each panel read from the potential file (array size npand)
   CHARACTER(LEN=28)                ::  VPOT_NAME(2) 		  !! header in the potential file for each spin channel
 END TYPE CELL_TYPE

end module type_cell
