!------------------------------------------------------------------------------
!> wrapper for shape-function generation from voronoi data using sensible
!> defaults for some parameters.

subroutine shapewrapper(npoi,aface,bface,cface,dface, &
  nmin, &
  nvertices,xvert,yvert,zvert,nface,lmax, dlt, &
  npan, nm, xrn, drn, meshn, &  ! radial mesh ! output parameters
  thetas_s, lmifun_s, nfun, & ! shape function
  ibmaxd,meshnd, npand,nfaced, nvertd)

  use shapefunctions_mod, only: shapef
  implicit none

  integer :: nvertd !< maximal number of cell vertices
  integer :: nfaced !< maximal number of cell faces


  integer, intent(in) :: npoi
  integer, intent(in) :: nmin
  integer, intent(in) :: nface
  integer, intent(in) :: lmax
  double precision, intent(in) :: dlt
  integer, intent(in) :: ibmaxd
  integer, intent(in) :: meshnd
  integer, intent(in) :: npand

  integer :: nvertices(nfaced)
  double precision :: aface(nfaced), bface(nfaced), cface(nfaced), dface(nfaced)
  double precision :: xvert(nvertd,nfaced), yvert(nvertd,nfaced), zvert(nvertd,nfaced)

  ! output
  integer, intent(out) ::   nm(npand)
  double precision, intent(out)  ::   xrn(meshnd)
  double precision, intent(out)  ::   drn(meshnd)
  integer, intent(out) ::   npan
  integer, intent(out) ::   meshn

  double precision, intent(out) ::  thetas_s(meshnd,ibmaxd)
  integer, intent(out) :: lmifun_s(ibmaxd)
  integer, intent(out) :: nfun

  double precision, parameter :: tolvdist = 1.d-12, toleuler = 1.d-10
  integer, parameter :: keypan = 0
  
  call shapef(npoi,aface,bface,cface,dface, &
    tolvdist, &
    toleuler, &
    nmin, &
    nvertices,xvert,yvert,zvert,nface,lmax, &
    keypan, dlt, &
    npan, nm, xrn, drn, meshn, &  ! radial mesh ! output parameters
    thetas_s, lmifun_s, nfun, & ! shape function
    ibmaxd,meshnd, npand)


end subroutine
