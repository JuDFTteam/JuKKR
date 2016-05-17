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
  public :: destroy
  public :: createShape
  public :: write_shapefun_file

  !> A datastructure containing the corresponding interstitial mesh.
  type InterstitialMesh
    integer :: npan !< number of panels
    integer, allocatable :: nm(:) !< positions of panels
    double precision, allocatable :: xrn(:) !< radial mesh points r(i)
    double precision, allocatable :: drn(:) !< integration weights dr/di (i)
  endtype
  
  interface destroy
    module procedure destroyInterstitialMesh
  endinterface

  contains

  !------------------------------------------------------------------------------
  !> Reads Voronoi weights from file 'voro_weights'
  !> this is kind of a hack and does not scale well.
  subroutine read_voro_weights(weights, atom_index, num_atoms)
    double precision, intent(inout) :: weights(:)
    integer, intent(in) :: atom_index(:) ! index table
    integer, intent(in) :: num_atoms

    integer :: ii, ios
    double precision :: weight_table(num_atoms)
 
    open(32, file='voro_weights', form='formatted', action='read', status='old', iostat=ios)
    if (ios /= 0) then
      write(0,*) "Warning! file voro_weights cannot be opened, use 1.0 for all."
      write(*,*) "Warning! file voro_weights cannot be opened, use 1.0 for all."
      weights(:) = 1.d0
      return
    endif

    do ii = 1, num_atoms
      read(32, *) weight_table(ii)
    enddo ! ii
    close(32, iostat=ios)

    do ii = 1, size(atom_index)
      weights(ii) = weight_table(atom_index(ii))
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
  subroutine createShape(self, inter_mesh, rbasis, bravais, center_ind, &
                      rcluster, lmax_shape, npoints_min, nmin_panel, &
                      num_MT_points, new_MT_radius, MT_scale, atom_id, num_atoms)
    use LatticeVectors_mod, only: LatticeVectors, create, destroy
    use RefCluster_mod, only: RefCluster, create, destroy
    use ShapefunData_mod, only: ShapefunData

    ! Output (shape-functions and interstitial mesh):
    type(ShapefunData), intent(inout) :: self
    type(InterstitialMesh), intent(inout) :: inter_mesh

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
    integer, intent(in) :: atom_id
    integer, intent(in) :: num_atoms

    type(LatticeVectors) :: lattice_vectors
    type(RefCluster) :: ref_cluster
    double precision, allocatable :: weights(:)
    integer :: ist

    call create(lattice_vectors, bravais)
    call create(ref_cluster, lattice_vectors%rr, rbasis, rcluster, center_ind)

    allocate(weights(size(ref_cluster%rcls, 2)))
    weights = 1.d0 ! all weights are the same, so their differences are zero 
    ! --> planes cut space in the exact center between two atomic cores

#ifdef USE_VOROWEIGHTS
    ! HACK: Read voronoi weights from file - this does not scale well!
    call read_voro_weights(weights, ref_cluster%atom, size(rbasis, 2))
#endif

    ! the cluster positions are in ref_cluster%rcls
    ! they are sorted by distance from center
    ! the voronoi routine expects array without position (0,0,0) -> pass rcls(:,2:)
    ! therfore check if the first index of %rcls is the origin

    if (any(abs(ref_cluster%rcls(:,1)) > 1.d-8)) stop "Expected origin in ref_cluster%rcls(:,1)"

    call constructFromCluster(self, inter_mesh, ref_cluster%rcls(:,2:), weights, &
                              lmax_shape, npoints_min, nmin_panel, &
                              num_MT_points, new_MT_radius, MT_scale, atom_id, num_atoms)

    call destroy(lattice_vectors)
    call destroy(ref_cluster)
    
    deallocate(weights, stat=ist)

  endsubroutine ! construct


  !------------------------------------------------------------------------------
  !> @param num_MT_mesh add 'num_MT_mesh' radial points of MT-region to
  !>        shape-function mesh -> non-touching MT-spheres
  !>        = 0 to not use this feature
  ! TODO: return 'NM' -> panel info!!!
  !
  !> @param weights weights for weighted Voronoi diagram (power diagram)
  subroutine constructFromCluster(self, inter_mesh, rvec, weights, &
                                  lmax_shape, npoints_min, nmin, &
                                  num_MT_points, new_MT_radius, MT_scale, atom_id, num_atoms)
    use ShapefunData_mod, only: ShapefunData, createShapefunData
    use Voronoi_mod, only: Voronoi_construction
    use ShapeFunctions_mod, only: shapef

    type(ShapefunData), intent(out) :: self
    type(InterstitialMesh), intent(out) :: inter_mesh

    double precision, intent(in) :: rvec(:,:)
    double precision, intent(in) :: weights(:)
    integer, intent(in) :: lmax_shape
    integer, intent(in) :: npoints_min
    integer, intent(in) :: nmin !< minimum number of points in panel
    integer, intent(in) :: num_MT_points
    double precision, intent(in) :: new_MT_radius
    double precision, intent(in) :: MT_scale
    integer, intent(in) :: atom_id
    integer, intent(in) :: num_atoms
    
    integer, parameter :: nvertmax = 32 ! hoping for at most 32 vertices for each face
    
    double precision, parameter :: tolvdist = 1.d-12
    double precision, parameter :: tolvarea = 1.d-10
    double precision, parameter :: toleuler = 1.d-10
    double precision, parameter :: dlt = 0.05d0 ! step-size angular integration

    double precision :: rmt, rout, volume
    integer :: ibmaxd, npand, nface, meshnd, meshn, nfaced, nfun
    double precision, allocatable :: planes(:,:)
    double precision, allocatable :: vert(:,:,:)
    double precision, allocatable :: xrn(:), drn(:)
    double precision, allocatable :: thetas_s(:,:)
    integer, allocatable :: nm(:), nvertices(:), lmifun_s(:)
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

    call Voronoi_construction(nfaced, rvec, nvertmax, nfaced, weights(1), weights(2:), tolvdist, tolvarea, rmt, rout, volume, nface, planes, nvertices, vert, &
#ifdef DEBUGSHAPEFUNCTIONS
         atom_id, output=.true.) ! for debug
#else
         atom_id, output=.false.)
#endif
         
    npand = sum(nvertices) + nface + 1  ! +1 for possible muffin-tinisation
    meshnd = max(npoints_min, npand*nmin) + num_mt_points

    ! increase meshnd, since in some cases the above empirical formula was not enough
    meshnd = meshnd + npoints_min

    ! these arrays are overdimensioned - but used only in this routine
    allocate(nm(npand))
    allocate(thetas_s(meshnd, ibmaxd))
    allocate(xrn(meshnd))
    allocate(drn(meshnd))
    
    call shapef(npoints_min, planes, tolvdist, toleuler, nmin, nvertices, vert, nface, lmax_shape, keypan, dlt, & 
      npan, nm, xrn, drn, meshn, &  ! radial mesh (output)
      thetas_s, lmifun_s, nfun, & ! shape function (output)
      ibmaxd, meshnd, npand, atom_id)

    ! muffin-tinization
    radius = new_MT_radius; if (MT_scale > tolvdist) radius = min(rmt*MT_scale, rmt) ! MT_scale > 0.0 overrides new_MT_radius

    if (num_MT_points > 0) call mtmesh(num_MT_points, npan, meshn, nm, xrn, drn, nfun, thetas_s, lmifun_s, radius, atom_id)

    ! Construct shape-fun datastructure
    call createShapefunData(self, meshn, ibmaxd, nfun, num_atoms)

    self%theta = thetas_s(1:meshn,1:nfun) ! store the shape functions
    
    self%max_muffin_tin = rmt ! set maximum possible muffin-tin radius (ALAT units)
    self%num_faces = nface ! store number of faces of Voronoi cell

    self%nfu = nfun
    self%ifunm = 0
    self%lmsp = 0
    do ii = 1, nfun
      self%llmsp(ii) = lmifun_s(ii)
      self%lmsp(lmifun_s(ii)) = 1
      self%ifunm(lmifun_s(ii)) = ii
    enddo ! ii

    ! createInterstitialMesh
    allocate(inter_mesh%xrn(meshn), inter_mesh%drn(meshn), inter_mesh%nm(npan))

    inter_mesh%xrn  = xrn(1:meshn)
    inter_mesh%drn  = drn(1:meshn)
    inter_mesh%nm   = nm(1:npan)
    inter_mesh%npan = npan

    !TODO: check that nm is always >= 5

    !call replace_with_PseudoASA(self, inter_mesh, volume) ! uncommenting this line would replace shape function with atomic sphere shape function

  endsubroutine ! construct

  !------------------------------------------------------------------------------
  elemental subroutine destroyInterstitialMesh(inter_mesh)
    type(InterstitialMesh), intent(inout) :: inter_mesh
    integer :: ist
    
    deallocate(inter_mesh%xrn, inter_mesh%drn, inter_mesh%nm, stat=ist)
    inter_mesh%npan = 0
  endsubroutine ! destroy

  !------------------------------------------------------------------------------
  subroutine mtmesh(nrad, npan, meshn, nm, xrn, drn, nfu, thetas, lmifun, mtradius, atom_id)
    use Constants_mod, only: pi
    ! program  mtmesh.f adds one extra pannel inside the muffin-tin sphere to allow lattice relaxations.
    ! stores the mt-nized shapes in unit 15 as shapefun
    !     nrad : number of points added inside the mt radius
    integer, intent(in) :: nrad, nfu
    integer, intent(inout) :: npan, meshn
    integer, intent(inout) :: nm(:)
    double precision, intent(inout) :: xrn(:), drn(:)
    double precision, intent(inout) :: thetas(:,:)
    integer, intent(inout) :: lmifun(:) ! unchanged
    double precision, intent(in) :: mtradius
    integer, intent(in) :: atom_id

    double precision :: drn1(meshn+nrad), thetas1(meshn+nrad,size(thetas,2)), xrn1(meshn+nrad)
    integer :: ibmaxd, irid, npand, meshn1, npan1, nm1(npan+1), ir
    double precision :: dist, dn1

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
      write(6,*) 'error from mtmesh for atom #',atom_id
      write(6,*) 'your mt-radius is bigger than the minimum shape radius ' 
      write(6,*) 'your mt-radius .....',mtradius
      write(6,*) 'shape radius .......',xrn(1)    
!       stop
      write(6,*) 'WARNING! stop statement deactivated:',__FILE__,":",__LINE__
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

    thetas1(1:nrad,:) = 0.d0
    thetas1(1:nrad,1) = sqrt(4.d0*pi)

    thetas1(nrad+1:meshn1,1:nfu) = thetas(1:meshn1-nrad,1:nfu)
    !
    !  now map back and return. 
    !
    npan = npan1
    meshn = meshn1

    nm(1:npan) = nm1(1:npan)

    xrn(1:meshn) = xrn1(1:meshn)
    drn(1:meshn) = drn1(1:meshn)
      
    thetas(1:meshn,1:nfu) = thetas1(1:meshn,1:nfu)

  endsubroutine ! mtmesh

  !------------------------------------------------------------------------------
  !> Write shape function, interstitial mesh + panels in a format compatible
  !> to what the Juelich KKR programs expect.
  !>
  !> The name of the file written is shape.<shape_index>
  subroutine write_shapefun_file(self, inter_mesh, shape_index)
    use ShapefunData_mod, only: ShapefunData

    type(ShapefunData), intent(in) :: self
    type(InterstitialMesh), intent(in) :: inter_mesh
    integer, intent(in) :: shape_index

    integer :: ir, ifun, npan, meshn
    character(len=16) :: filename
    character(len=*), parameter :: F9000="(16i5)", F9010="(4d20.12)" ! todo: define these formats with meaningful names in ShapefunData_mod

    write(filename, fmt="(a,i7.7)") "shape.",shape_index

    npan = inter_mesh%npan
    meshn = size(inter_mesh%xrn)

    open(15, file=filename, form="formatted", action='write')
    write(15, fmt=F9000) npan, meshn
    write(15, fmt=F9000) inter_mesh%nm(1:npan)
    write(15, fmt=F9010) (inter_mesh%xrn(ir), inter_mesh%drn(ir), ir=1,meshn) ! interleaved pairs (r, dr)
    write(15, fmt=F9000) self%nfu
    do ifun = 1, self%nfu
      write(15, fmt=F9000) self%llmsp(ifun)
      write(15, fmt=F9010) self%theta(1:meshn,ifun)
    enddo ! ifun
    close(15)

  endsubroutine ! write_shapefun_file

  !------------------------------------------------------------------------------
  !> Replace actual shape function with shape function of the atom sphere.
  !>
  !> Only the LM=1 component is nonzero
  !> Atomic sphere has same volume as Voronoi cell
  subroutine replace_with_PseudoASA(self, inter_mesh, volume)
    use ShapefunData_mod, only: ShapefunData
    use Constants_mod, only: pi

    type(ShapefunData), intent(inout) :: self
    type(InterstitialMesh), intent(in) :: inter_mesh
    double precision, intent(in) :: volume

    double precision :: radius
    integer :: ir

    self%lmsp = 0
    self%lmsp(1) = 1
    self%nfu = 1

    radius = (3*volume/(4*pi))**(1.d0/3.d0)

    self%theta = 0.d0
    do ir = 1, size(inter_mesh%xrn)
      if (inter_mesh%xrn(ir) <= radius) self%theta(ir,1) = sqrt(4.d0*pi)
    enddo ! ir

  endsubroutine ! replace_with_PseudoASA

endmodule ! ConstructShapes_mod
