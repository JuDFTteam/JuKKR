!------------------------------------------------------------------------------
!> Module that provides a simple interface for doing shape-function construction
!> from positions only.
!>
!> This effectively hides the Voronoi construction.
!> Reuses code from RefCluster_mod.f90
!>
!> @author Elias Rabel
!>
module ConstructShapes_mod

!> A datastructure containing the corresponding interstitial mesh.
type InterstitialMesh
  double precision, allocatable, dimension(:) :: xrn !< radial mesh points r(i)
  double precision, allocatable, dimension(:) :: drn !< integration weights dr/di (i)
  integer :: npan !< number of panels
  integer, allocatable, dimension(:) :: nm !< positions of panels
end type

CONTAINS

!------------------------------------------------------------------------------
!> Construct shape functions and interstitial mesh.
!>
!> Creates ShapefunData datastructure -> user has to deal with deallocation!
!> Creates InterstitialMesh datastructure -> user has to deal with deallocation!
!> (use destroyInterstitialMesh)
!>
!> MT_scale > 0.0 overrides new_MT_radius!!!
subroutine construct(shdata, inter_mesh, rbasis, bravais, center_ind, &
                     rcluster, lmax_shape, npoints_min, nmin_panel, &
                     num_MT_points, new_MT_radius, MT_scale)
  use RefCluster_mod
  use ShapefunData_mod
  implicit none

  ! Output (shape-functions and interstitial mesh):
  type (ShapefunData), intent(inout) :: shdata
  type (InterstitialMesh), intent(inout) :: inter_mesh

  ! Input:
  double precision, intent(in) :: rbasis(:,:)
  double precision, intent(in) :: bravais(3,3)
  integer, intent(in)          :: center_ind
  double precision, intent(in) :: rcluster
  integer, intent(in) :: lmax_shape
  integer, intent(in) :: npoints_min
  integer, intent(in) :: nmin_panel
  integer, intent(in) :: num_MT_points
  double precision, intent(in) :: new_MT_radius
  double precision, intent(in) :: MT_scale

  !----------
  type (LatticeVectors) :: lattice_vectors
  type (RefCluster) :: ref_cluster

  call createLatticeVectors(lattice_vectors, bravais)
  call createRefCluster(ref_cluster, lattice_vectors, rbasis, rcluster, center_ind)

  ! the cluster positions are in ref_cluster%rcls
  ! they are sorted by distance from center
  ! the voronoi routine expects array without position (0,0,0) -> pass rcls(:,2:)

  if (sum(abs(ref_cluster%rcls(:,1))) > 1.d-8) then
    write(*,*) "Expected origin in ref_cluster%rcls(:,1)"
    STOP
  end if

  call constructFromCluster(shdata, inter_mesh, ref_cluster%rcls(:,2:), lmax_shape, &
                            npoints_min, nmin_panel, &
                            num_MT_points, new_MT_radius, MT_scale)

  call destroyLatticeVectors(lattice_vectors)
  call destroyRefCluster(ref_cluster)

end subroutine


!------------------------------------------------------------------------------
!> @param num_MT_mesh add 'num_MT_mesh' radial points of MT-region to
!>        shape-function mesh -> non-touching MT-spheres
!>        = 0 to not use this feature
! TODO: return 'NM' -> panel info!!!
subroutine constructFromCluster(shdata, inter_mesh, rvec, lmax_shape, &
                                npoints_min, nmin, &
                                num_MT_points, new_MT_radius, MT_scale)
  use ShapefunData_mod
  implicit none

  type (ShapefunData), intent(inout) :: shdata
  type (InterstitialMesh), intent(inout) :: inter_mesh

  double precision, intent(in) :: rvec(:, :)
  integer, intent(in) :: lmax_shape
  integer, intent(in) :: npoints_min
  integer, intent(in) :: nmin !< minimum number of points in panel
  integer, intent(in) :: num_MT_points
  double precision, intent(in) :: new_MT_radius
  double precision, intent(in) :: MT_scale

  integer, parameter :: NVERTMAX = 30  ! hoping for at most 30 vertices for each face
  logical, parameter :: OUTPUT = .false.
  
  double precision, parameter :: TOLVDIST = 1.d-12
  double precision, parameter :: TOLVAREA = 1.d-10
  double precision, parameter :: DLT = 0.05d0 ! step-size angular integration

  double precision :: rmt, rout, volume
  integer :: ibmaxd
  integer :: npand
  integer :: nface
  integer :: meshnd
  integer :: meshn
  integer :: nfaced
  integer :: nfun
  double precision, allocatable, dimension(:) :: aface, bface, cface, dface
  double precision, allocatable, dimension(:) :: weight
  integer, allocatable, dimension(:) :: nm
  integer, allocatable, dimension(:) :: nvertices
  double precision, allocatable, dimension(:) :: xrn, drn
  double precision, allocatable, dimension(:, :, :) :: vertices 
  double precision, allocatable, dimension(:, :) :: thetas_s
  integer :: npan
  integer, allocatable, dimension(:) :: lmifun_s
  integer :: ii

  double precision, parameter :: weight0 = 1.0d0
  double precision :: radius

  nfaced = size(rvec,2)

  ibmaxd = (lmax_shape + 1)**2
  allocate(lmifun_s(ibmaxd))
  lmifun_s = 0
  
  allocate( aface(nfaced) )
  allocate( bface(nfaced) )
  allocate( cface(nfaced) )
  allocate( dface(nfaced) )
  allocate( weight(nfaced) )
  ! support only unweighted voronoi diagrams for now
  weight = weight0
  allocate( nvertices(nfaced) )
  allocate( vertices(NVERTMAX, nfaced, 3) )
  nvertices = 0


  call voronoi08( &
       nfaced,rvec,NVERTMAX,nfaced,weight0,weight,TOLVDIST,TOLVAREA, &
       rmt,rout,volume,nface,aface,bface,cface,dface,nvertices, &
       vertices(:,:,1),vertices(:,:,2),vertices(:,:,3), &
       OUTPUT)

  npand = sum(nvertices) + nface + 1  ! +1 for possible muffin-tinisation
  meshnd = max(npoints_min, npand * NMIN) + num_MT_points

  allocate(nm(npand))
  allocate(thetas_s(meshnd, ibmaxd))
  allocate(xrn(meshnd))
  allocate(drn(meshnd))

  call shapewrapper(npoints_min,aface,bface,cface,dface, &
                    NMIN, &
                    nvertices,vertices(:,:,1),vertices(:,:,2),vertices(:,:,3), &
                    nface,lmax_shape, DLT, &
                    npan, nm, xrn, drn, meshn, & 
                    thetas_s, lmifun_s, nfun, & 
                    ibmaxd,meshnd, npand,nfaced, NVERTMAX)

  ! muffin-tinization
  radius = new_MT_radius
  if (MT_scale > TOLVDIST) then    ! MT_scale > 0.0 overrides new_MT_radius
    radius = min(rmt * MT_scale, rmt)
  end if

  if (num_MT_points > 0) then
    call mtmesh(num_MT_points,npan,meshn,nm,xrn,drn,nfun,thetas_s,lmifun_s, radius)
  end if

  call createShapefunData(shdata, meshn, ibmaxd, nfun)

  shdata%theta = thetas_s(1:meshn, 1:nfun)
  shdata%nfu = nfun
  shdata%ifunm = 0
  shdata%lmsp = 0

  do ii = 1, nfun
    shdata%llmsp(ii) = lmifun_s(ii)
    shdata%lmsp(lmifun_s(ii)) = 1
    shdata%ifunm(lmifun_s(ii)) = ii
  end do

  allocate(inter_mesh%xrn(meshn))
  allocate(inter_mesh%drn(meshn))
  allocate(inter_mesh%nm(npan))

  inter_mesh%xrn = xrn(1:meshn)
  inter_mesh%drn = drn(1:meshn)
  inter_mesh%nm  = nm(1:npan)
  inter_mesh%npan = npan

  !TODO: check that nm is always >= 5

end subroutine

!------------------------------------------------------------------------------
subroutine destroyInterstitialMesh(inter_mesh)
  implicit none
  type(InterstitialMesh), intent(inout) :: inter_mesh

  deallocate(inter_mesh%xrn)
  deallocate(inter_mesh%drn)
  deallocate(inter_mesh%nm)
  inter_mesh%npan = 0

end subroutine

!******************************************************************************
SUBROUTINE MTMESH(NRAD,NPAN,MESHN,NM,XRN,DRN, &
                  NFU,THETAS,LMIFUN,MTRADIUS)
!
IMPLICIT NONE
!
!
!
! Program  mtmesh.f adds one extra pannel inside the
! muffin-tin sphere to allow lattice relaxations.
! stores the mt-nized shapes in unit 15 as shapefun
!     .. Parameters ..
!     nrad : number of points added inside the MT radius
Integer ibmaxd, irid, npand
INTEGER NRAD
!      Parameter (nrad=20)
!     ..
!     .. Local Scalars ..
REAL*8           dist,dn1,mtradius,pi,scale
Integer ifun,ipan1,ir,lm,meshn,nfu,npan,number
integer ip,i
!     ..
!     .. Arrays ..
!REAL*8           drn(IRID),thetas(IRID,ibmaxd),xrn(IRID)
!Integer nm(npand),lmIFUN(ibmaxd)

REAL*8           drn(:),thetas(:,:),xrn(:)
Integer nm(:),lmIFUN(:)

!     ..
!     .. Local Arrays ..      
!REAL*8           drn1(IRID),thetas1(IRID,ibmaxd),xrn1(IRID)
!Integer nm1(npand)

REAL*8  drn1(meshn+nrad),thetas1(meshn+nrad,size(thetas,2)),xrn1(meshn+nrad)
Integer nm1(npan+1)
integer meshn1,npan1
!     ..
!     .. Intrinsic Functions ..
INTRINSIC ABS,DATAN,DSQRT,SQRT
!     ..

ibmaxd = size(thetas, 2)
irid   = size(thetas, 1)
npand  = size(nm)

PI = 4.D0*DATAN(1.D0)

NPAN1 = NPAN + 1
MESHN1 = MESHN + NRAD
NM1(1) = NRAD

DO IP=2,NPAN1
   NM1(IP) = NM(IP-1)
END DO

IF (NPAN1.GT.NPAND) THEN
  WRITE (6,FMT=*) ' NPAN , NPAND ',NPAN1,NPAND
  STOP
END IF

IF (MESHN1.GT.IRID) THEN
  WRITE (6,FMT=*) ' MESHN , IRID ',MESHN1,IRID
  STOP
END IF

DIST = XRN(1) - MTRADIUS

If (dist.lt.1.0d-5) Then
   Write (6,FMT=*) 'Error from MTMESH '
   write (6,*) &
 'Your MT-radious is biger that the minimum shape radious ' 
   write(6,*) 'Your MT-Radious .....',mtradius
   write(6,*) 'Shape Radious .......',xrn(1)    
  Stop
End if
     
dn1 = dist/ (nrad-1)
xrn1(nrad) = xrn(1)
drn1(nrad) = dn1
Do ir = 1,nrad-1
  xrn1(ir) = mtradius + dn1* (ir-1)
  drn1(ir) = dn1
End do
do i=1,meshn1-nrad
   xrn1(i+nrad) = xrn(i)
   drn1(i+nrad) = drn(i) 
end do

Do ir = 1,nrad
  thetas1(ir,1) = dsqrt(4.d0*pi)
  Do ifun = 2,ibmaxd
    thetas1(ir,ifun) = 0.0d0
  End do
End do

Do 10 ifun = 1,nfu
  do ir=1,meshn1-nrad
   thetas1(nrad+ir,ifun) = thetas(ir,ifun)
  end do
10 Continue
!
!  Now map back and return. 
!
npan = npan1
meshn = meshn1

do ip=1,npan
   nm(ip) = nm1(ip)
end do  
do ir=1,meshn
   xrn(ir) = xrn1(ir)
   drn(ir) = drn1(ir)
end do

Do ifun = 1,nfu
  do ir=1,meshn
   thetas(ir,ifun) = thetas1(ir,ifun)
  end do
end do
!
! Now store on disk
!
!WRITE(15,FMT=9000) NPAN,MESHN
!WRITE(15,FMT=9000) (nm(ipan1),ipan1=1,npan)
!WRITE(15,FMT=9010) (xrn(ir),drn(ir),ir=1,meshn)
!WRITE(15,FMT=9000) nfu
!DO IFUN = 1,NFU
!  WRITE(15,FMT=9000) lmifun(ifun)
!!         Write (6,FMT=*)     lmifun(ifun)
!  WRITE(15,FMT=9010) (thetas(ir,ifun),ir=1,meshn)
!ENDDO
!RETURN
9000 FORMAT (16i5)
9010 FORMAT (4d20.12)
End subroutine


end module
