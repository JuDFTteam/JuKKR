!
! Contains information of each individual cell
!
module type_cell
use nrtype
 TYPE                              ::  CELL_TYPE
   REAL(KIND=DP)                    ::   RMAX, RCORE, RMT         ! max radius, core radius, muffin tin radius
   INTEGER                          ::  NRMAX,NRCORE, NRNS        ! total grid points, core gp, non-sph gp 
   INTEGER                          ::  NRMAXD                    ! maximum NRMAX (for all cells)
   INTEGER,ALLOCATABLE              ::  NRCUT(:)                  ! array index of last position of a panel
   REAL(KIND=DP)                    ::  LOGPARAMS(2)              ! a,b for log. mesh : r(i) = b*(exp(a*(i-1))-1)
   REAL(KIND=DP),ALLOCATABLE        ::  RMESH(:)                  ! log. mesh
   REAL(KIND=DP),ALLOCATABLE        ::  DRMESHDI(:)               ! dR/dI
   REAL(KIND=DP),ALLOCATABLE        ::  DRMESHOR(:)               ! dR/R
   INTEGER                          ::  NRMIN_NS                  ! index of the most inner data for a non-
    															  ! spherical treatment
   INTEGER                          ::  KXC                       ! type of exchange corr. functional 

   INTEGER                          ::  NPAN 
   INTEGER                          ::  NPAND 
   INTEGER,ALLOCATABLE              ::  NMESHPAN(:)
   CHARACTER(LEN=28)                ::  VPOT_NAME(2)
 END TYPE CELL_TYPE

end module type_cell