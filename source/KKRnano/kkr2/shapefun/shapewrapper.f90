!------------------------------------------------------------------------------
SUBROUTINE shapewrapper(NPOI,AFACE,BFACE,CFACE,DFACE, &
NMIN, &
NVERTICES,XVERT,YVERT,ZVERT,NFACE,LMAX, DLT, &
NPAN, NM, XRN, DRN, MESHN, &  ! radial mesh ! output parameters
THETAS_S, LMIFUN_S, NFUN, & ! shape function
IBMAXD,MESHND, NPAND,NFACED, NVERTD)

  use ShapeFunctions_mod, only: SHAPEF
  implicit none

  integer, parameter :: DP = 8
  integer :: NVERTD !< maximal number of cell vertices
  integer :: NFACED !< maximal number of cell faces


  integer,intent(in) :: NPOI
  real(kind=DP), parameter :: TOLVDIST = 1.d-12
  real(kind=DP), parameter :: TOLEULER = 1.d-10
  integer, intent(in) :: NMIN
  integer, intent(in) :: NFACE
  integer, intent(in) :: LMAX
  integer, parameter :: KEYPAN = 0
  real(kind=DP), intent(in) :: DLT

  integer, intent(in) :: IBMAXD
  integer, intent(in) :: MESHND
  integer, intent(in) :: NPAND

  integer ::   NVERTICES(NFACED)
  real(kind=DP) ::    AFACE(NFACED),BFACE(NFACED),CFACE(NFACED),DFACE(NFACED)
  real(kind=DP) ::    XVERT(NVERTD,NFACED),YVERT(NVERTD,NFACED), &
  ZVERT(NVERTD,NFACED)

  ! output
  integer, intent(out) ::   NM(NPAND)
  real(kind=DP), intent(out)  ::   XRN(MESHND)
  real(kind=DP), intent(out)  ::   DRN(MESHND)
  integer, intent(out) ::   NPAN
  integer, intent(out) ::   MESHN

  real(kind=DP), intent(out) ::  THETAS_S(MESHND,IBMAXD)
  integer, intent(out) :: LMIFUN_S(IBMAXD)
  integer, intent(out) :: NFUN

  call SHAPEF(NPOI,AFACE,BFACE,CFACE,DFACE, &
    TOLVDIST, &
    TOLEULER, &
    NMIN, &
    NVERTICES,XVERT,YVERT,ZVERT,NFACE,LMAX, &
    KEYPAN, DLT, &
    NPAN, NM, XRN, DRN, MESHN, &  ! radial mesh ! output parameters
    THETAS_S, LMIFUN_S, NFUN, & ! shape function
    IBMAXD,MESHND, NPAND)


end subroutine
