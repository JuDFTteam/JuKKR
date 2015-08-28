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
  implicit none
  private
  public :: InterstitialMesh
  public :: destroyInterstitialMesh
  public :: construct
  public :: write_shapefun_file

  !> A datastructure containing the corresponding interstitial mesh.
  type InterstitialMesh
    integer :: npan !< number of panels
    integer, allocatable :: nm(:) !< positions of panels
    double precision, allocatable :: xrn(:) !< radial mesh points r(i)
    double precision, allocatable :: drn(:) !< integration weights dr/di (i)
  end type

  CONTAINS

!------------------------------------------------------------------------------
!> Reads Voronoi weights from file 'voro_weights'
!> this is kind of a hack and does not scale well.
subroutine read_voro_weights(weights, ATOM, num_atoms)
  double precision, intent(inout) :: weights(:)
  integer, intent(in) :: ATOM(:)
  integer, intent(in) :: num_atoms

  integer :: ii
  double precision, allocatable :: weight_table(:)

  allocate(weight_table(num_atoms))

  open(32, file='voro_weights', form='formatted')
    do ii = 1, num_atoms
      read(32, *) weight_table(ii)
    end do
  close(32)

  do ii = 1, size(ATOM)
    weights(ii) = weight_table(ATOM(ii))
  enddo ! ii
endsubroutine ! read_voro_weights

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
  use RefCluster_mod, only: LatticeVectors, RefCluster
  use RefCluster_mod, only: createLatticeVectors, createRefCluster, destroyLatticeVectors, destroyRefCluster
  use ShapefunData_mod, only: ShapefunData

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
  type(LatticeVectors) :: lattice_vectors
  type(RefCluster) :: ref_cluster
  double precision, allocatable :: weights(:)
  integer :: ist

  call createLatticeVectors(lattice_vectors, bravais)
  call createRefCluster(ref_cluster, lattice_vectors, rbasis, rcluster, center_ind)

  allocate(weights(size(ref_cluster%rcls, 2)))
  weights = 1.0d0

#ifdef USE_VOROWEIGHTS
  ! HACK: Read voronoi weights from file - this does not scale well!
  call read_voro_weights(weights, ref_cluster%atom, size(rbasis, 2))
#endif

  ! the cluster positions are in ref_cluster%rcls
  ! they are sorted by distance from center
  ! the voronoi routine expects array without position (0,0,0) -> pass rcls(:,2:)
  ! therfore check if the first index of %rcls is the origin

  if (any(abs(ref_cluster%rcls(:,1)) > 1.d-8)) then
    write(*,*) "Expected origin in ref_cluster%rcls(:,1)"
    STOP
  endif

  call constructFromCluster(shdata, inter_mesh, ref_cluster%rcls(:,2:), weights, &
                            lmax_shape, npoints_min, nmin_panel, &
                            num_MT_points, new_MT_radius, MT_scale)

  call destroyLatticeVectors(lattice_vectors)
  call destroyRefCluster(ref_cluster)
  
  deallocate(weights, stat=ist)

endsubroutine ! construct


!------------------------------------------------------------------------------
!> @param num_MT_mesh add 'num_MT_mesh' radial points of MT-region to
!>        shape-function mesh -> non-touching MT-spheres
!>        = 0 to not use this feature
! TODO: return 'NM' -> panel info!!!
!
!> @param weights weights for weighted Voronoi diagram (power diagram)
subroutine constructFromCluster(shdata, inter_mesh, rvec, weights, &
                                lmax_shape, npoints_min, nmin, &
                                num_MT_points, new_MT_radius, MT_scale)
  use ShapefunData_mod, only: ShapefunData, createShapefunData
  use voronoi08_mod, only: voronoi08
  use shapefunctions_mod, only: shapef


  type (ShapefunData), intent(inout) :: shdata
  type (InterstitialMesh), intent(inout) :: inter_mesh

  double precision, intent(in) :: rvec(:,:)
  double precision, intent(in) :: weights(:)
  integer, intent(in) :: lmax_shape
  integer, intent(in) :: npoints_min
  integer, intent(in) :: nmin !< minimum number of points in panel
  integer, intent(in) :: num_MT_points
  double precision, intent(in) :: new_MT_radius
  double precision, intent(in) :: MT_scale
  
  integer, parameter :: nvertmax = 30  ! hoping for at most 30 vertices for each face
  logical, parameter :: output = .false.
  
  double precision, parameter :: tolvdist = 1.d-12
  double precision, parameter :: tolvarea = 1.d-10, toleuler = 1.d-10
  double precision, parameter :: dlt = 0.05d0 ! step-size angular integration

  double precision :: rmt, rout, volume
  integer :: ibmaxd
  integer :: npand
  integer :: nface
  integer :: meshnd
  integer :: meshn
  integer :: nfaced
  integer :: nfun
  double precision, allocatable :: planes(:,:)
  integer, allocatable :: nm(:), nvertices(:), lmifun_s(:)
  double precision, allocatable :: xrn(:), drn(:)
  double precision, allocatable :: vert(:,:,:)
  double precision, allocatable :: thetas_s(:,:)
  integer :: npan, ii

  double precision :: radius

  integer, parameter :: keypan = 0
  
  nfaced = size(rvec,2)

  ibmaxd = (lmax_shape+1)**2
  allocate(lmifun_s(ibmaxd))
  lmifun_s = 0
  
  allocate(planes(0:3,nfaced)) ! 0: distance, 1..3:normal vector

  allocate(nvertices(nfaced), vert(3,nvertmax,nfaced))
  nvertices = 0

  call voronoi08(nfaced,rvec,nvertmax,nfaced, weights(1), weights(2:),tolvdist,tolvarea, rmt,rout,volume,nface, planes, nvertices, vert, output)

  npand = sum(nvertices) + nface + 1  ! +1 for possible muffin-tinisation
  meshnd = max(npoints_min, npand*nmin) + num_mt_points

  ! increase meshnd, since in some cases the above empirical formula was
  ! not enough
  meshnd = meshnd + npoints_min

  ! these arrays are overdimensioned - but used only in this routine
  allocate(nm(npand))
  allocate(thetas_s(meshnd, ibmaxd))
  allocate(xrn(meshnd))
  allocate(drn(meshnd))

!   call shapewrapper(npoints_min,planes, &
!                     nmin, &
!                     nvertices, vert, &
!                     nface, lmax_shape, dlt, &
!                     npan, nm, xrn, drn, meshn, & 
!                     thetas_s, lmifun_s, nfun, & 
!                     ibmaxd, meshnd, npand, nfaced, nvertmax)
  
  call shapef(npoints_min, planes, &
    tolvdist, toleuler, &
    nmin, &
    nvertices, vert, nface, lmax_shape, &
    keypan, dlt, &
    npan, nm, xrn, drn, meshn, &  ! radial mesh ! output parameters
    thetas_s, lmifun_s, nfun, & ! shape function
    ibmaxd, meshnd, npand)
                   

  ! muffin-tinization
  radius = new_mt_radius
  if (mt_scale > tolvdist) radius = min(rmt*MT_scale, rmt) ! MT_scale > 0.0 overrides new_MT_radius

  if (num_MT_points > 0) call mtmesh(num_MT_points,npan,meshn,nm,xrn,drn,nfun,thetas_s,lmifun_s, radius)

  ! Construct shape-fun datastructure
  call createShapefunData(shdata, meshn, ibmaxd, nfun)

  shdata%theta = thetas_s(1:meshn,1:nfun) ! store the shape functions
  
  shdata%max_muffin_tin = rmt ! set maximum possible muffin-tin radius (ALAT units)
  shdata%num_faces = nface ! store number of faces of Voronoi cell

  shdata%nfu = nfun
  shdata%ifunm = 0
  shdata%lmsp = 0
  do ii = 1, nfun
    shdata%llmsp(ii) = lmifun_s(ii)
    shdata%lmsp(lmifun_s(ii)) = 1
    shdata%ifunm(lmifun_s(ii)) = ii
  enddo ! ii

  allocate(inter_mesh%xrn(meshn), inter_mesh%drn(meshn), inter_mesh%nm(npan))

  inter_mesh%xrn  = xrn(1:meshn)
  inter_mesh%drn  = drn(1:meshn)
  inter_mesh%nm   = nm(1:npan)
  inter_mesh%npan = npan

  !TODO: check that nm is always >= 5

  !DEBUG
  ! uncommenting the next line would replace shape function with atomic sphere shape function
  !call replace_with_PseudoASA(shdata, inter_mesh, volume)

endsubroutine

!------------------------------------------------------------------------------
subroutine destroyInterstitialMesh(inter_mesh)
  type(InterstitialMesh), intent(inout) :: inter_mesh

  deallocate(inter_mesh%xrn, inter_mesh%drn, inter_mesh%nm)
  inter_mesh%npan = 0

endsubroutine ! destroy

!******************************************************************************
  subroutine mtmesh(nrad, npan, meshn, nm, xrn, drn, nfu, thetas, lmifun, mtradius)
    use shape_constants_mod, only: pi
    ! program  mtmesh.f adds one extra pannel inside the
    ! muffin-tin sphere to allow lattice relaxations.
    ! stores the mt-nized shapes in unit 15 as shapefun
    !     nrad : number of points added inside the mt radius
    integer, intent(in) :: nrad, nfu
    
    integer, intent(inout) :: npan, meshn
    double precision, intent(inout) :: drn(:), thetas(:,:), xrn(:)
    integer, intent(inout) :: nm(:)

    integer, intent(in) :: lmifun(:) ! not used, please remove from interface


    double precision :: drn1(meshn+nrad), thetas1(meshn+nrad,size(thetas,2)), xrn1(meshn+nrad)
    integer :: ibmaxd, irid, npand, meshn1, npan1, nm1(npan+1), ifun, ir, ip
    double precision :: dist, dn1, mtradius

    ibmaxd = size(thetas, 2)
    irid   = size(thetas, 1)
    npand  = size(nm)

    npan1 = npan + 1
    meshn1 = meshn + nrad
    nm1(1) = nrad

    nm1(2:npan1) = nm(1:npan1-1)

    if (npan1 > npand) then
      write (6,fmt=*) ' npan , npand ',npan1,npand
      stop
    endif

    if (meshn1 > irid) then
      write (6,fmt=*) ' meshn , irid ',meshn1,irid
      stop
    endif

    dist = xrn(1) - mtradius

    if (dist < 1.d-5) then
      write(6,fmt=*) 'error from mtmesh '
      write(6,*) 'your mt-radius is bigger that the minimum shape radius ' 
      write(6,*) 'your mt-radius .....',mtradius
      write(6,*) 'shape radius .......',xrn(1)    
      stop
    endif
        
    dn1 = dist/(nrad-1)
    xrn1(nrad) = xrn(1)
    drn1(nrad) = dn1
    do ir = 1, nrad-1
      xrn1(ir) = mtradius + dn1*(ir-1)
      drn1(ir) = dn1
    enddo ! ir

    xrn1(1+nrad:meshn1) = xrn(1:meshn1-nrad)
    drn1(1+nrad:meshn1) = drn(1:meshn1-nrad) 

    thetas1(1:nrad,1) = sqrt(4.d0*pi)
    do ifun = 2, ibmaxd
      thetas1(1:nrad,ifun) = 0.d0
    enddo ! ifun

    do ifun = 1, nfu
      thetas1(1+nrad:meshn1,ifun) = thetas(1:meshn1-nrad,ifun)
    enddo ! ifun
    !
    !  now map back and return. 
    !
    npan = npan1
    meshn = meshn1

    do ip = 1, npan
      nm(ip) = nm1(ip)
    enddo ! ip

    xrn(1:meshn) = xrn1(1:meshn)
    drn(1:meshn) = drn1(1:meshn)
      
    thetas(1:meshn,1:nfu) = thetas1(1:meshn,1:nfu)

  endsubroutine

!------------------------------------------------------------------------------
!> Write shape function, interstitial mesh + panels in a format compatible
!> to what the Juelich KKR programs expect.
!>
!> The name of the file written is shape.<shape_index>
subroutine write_shapefun_file(shdata, inter_mesh, shape_index)
  use ShapefunData_mod, only: ShapefunData

  type(ShapefunData), intent(in) :: shdata
  type(InterstitialMesh), intent(in) :: inter_mesh
  integer, intent(in) :: shape_index

  integer :: ir, ifun, npan, meshn
  character(len=16) :: filename
  character(len=*), parameter :: f9000="(16i5)", f9010="(4d20.12)" 

  write(filename, fmt="(a,i7.7)") "shape.",shape_index

  npan = inter_mesh%npan
  meshn = size(inter_mesh%xrn)

  open(15, file=filename, form="formatted")
  write(15,fmt=f9000) npan,meshn
  write(15,fmt=f9000) inter_mesh%nm(1:npan)
  write(15,fmt=f9010) (inter_mesh%xrn(ir), inter_mesh%drn(ir), ir=1,meshn) ! interleaved pairs (r, dr)
  write(15,fmt=f9000) shdata%nfu
  do ifun = 1, shdata%nfu
    write(15,fmt=f9000) shdata%llmsp(ifun)
    write(15,fmt=f9010) shdata%theta(1:meshn,ifun)
  enddo ! ifun
  close(15)

endsubroutine write_shapefun_file

!------------------------------------------------------------------------------
!> Replace actual shape function with shape function of the atom sphere.
!>
!> Only the LM=1 component is nonzero
!> Atomic sphere has same volume as Voronoi cell
subroutine replace_with_PseudoASA(shdata, inter_mesh, volume)
  use ShapefunData_mod, only: ShapefunData
  use shape_constants_mod, only: pi

  type(ShapefunData), intent(inout) :: shdata
  type(InterstitialMesh), intent(in) :: inter_mesh
  double precision, intent(in) :: volume

  double precision :: radius
! double precision, parameter :: PI = 3.1415926535897932d0
  integer :: ir

  shdata%lmsp = 0
  shdata%lmsp(1) = 1
  shdata%nfu = 1

  radius = (3*volume/(4*pi))**(1.d0/3.d0)

  shdata%theta = 0.d0
  do ir = 1, size(inter_mesh%xrn)
    if (inter_mesh%xrn(ir) <= radius) shdata%theta(ir,1) = sqrt(4.d0*pi)
  enddo ! ir

endsubroutine ! replace_with_PseudoASA

endmodule ! ConstructShapes_mod
