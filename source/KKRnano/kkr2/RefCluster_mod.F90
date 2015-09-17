!> Definitions
!> RefCluster: a cluster of atoms with a distance smaller than rcut to a
!>             central atom
!>             periodic boundary conditions are taken into account:
!>             mirror images of atoms can be included in the cluster
!>
!> LatticeVectors: linear combinations of Bravais vectors with integer coefficients
!>
!> Note: if N (=number of atoms) reference clusters are created, then this algorithm scales as O(N**2)
! #define DEBUG

module RefCluster_mod
  implicit none
  private
  
  public :: RefCluster, LatticeVectors, create, destroy
  public :: createLatticeVectors, destroyLatticeVectors, createRefCluster, destroyRefCluster

  type LatticeVectors
    integer :: nrd
    double precision, allocatable :: rr(:,:) ! dim(3,0:nrd)
  endtype

  type RefCluster
    !> reference to LatticeVectors datastructure
    double precision, allocatable :: rcls(:,:) !< dim(3,nacls) positions relative to center
    integer, allocatable :: atom(:)  !< dim(nacls) basis atom indices of cluster atoms
    integer, allocatable :: ezoa(:)  !< dim(nacls) points into lattice_vectors%rr
    integer, allocatable :: indn0(:) !< dim(numn0) indices of inequivalent cluster atoms
    integer :: numn0 !< number of inequivalent cluster atoms
    integer :: nacls !< number of cluster atoms
    integer :: atom_index !< basis atom index of central cluster atom
  endtype

  interface create
    module procedure createLatticeVectors, createRefCluster
  endinterface
  
  interface destroy
    module procedure destroyLatticeVectors, destroyRefCluster
  endinterface
  
  contains

  !------------------------------------------------------------------------------
  !> Creates a table of lattice vectors.
  !>
  !> Before any reference clusters can be created, a table of
  !> lattice vectors has to be created.
  !> Several different reference clusters can (and should) share
  !> the same lattice vectors.
  subroutine createLatticeVectors(lattice_vectors, bravais)
    type(LatticeVectors), intent(inout) :: lattice_vectors
    double precision, intent(in) :: bravais(3,3)

    call rrgen(bravais, lattice_vectors%rr, lattice_vectors%nrd)
#ifdef DEBUG  
    write(*,*) "Num. real space vectors: ", lattice_vectors%nrd
#endif
  endsubroutine

  !------------------------------------------------------------------------------
  !> Destroys table of lattice vectors.
  !>
  subroutine destroyLatticeVectors(lattice_vectors)
    type(LatticeVectors), intent(inout) :: lattice_vectors

    integer :: ist
    deallocate(lattice_vectors%rr, stat=ist)
    lattice_vectors%nrd = 0
  endsubroutine

  !------------------------------------------------------------------------------
  !> Creates reference cluster.
  !>
  !> Usage example:
  !> type(LatticeVectors) :: lattice_vectors
  !> type(RefCluster) :: ref_cluster
  !> ! ...
  !> call createLatticeVectors(lattice_vectors, bravais)
  !> call createRefCluster(ref_cluster, lattice_vectors, rbasis, rcut, center_ind)
  !> ! ... some code
  !> call destroyRefCluster(ref_cluster)
  !> call destroyLatticeVectors(lattice_vectors)
  !
  subroutine createRefCluster(self, lattice_vectors, rbasis, rcut, center_ind)
    type(RefCluster), intent(inout) :: self

    type(LatticeVectors), intent(in) :: lattice_vectors
    double precision, intent(in) :: rbasis(:,:) !< dim(3,naez) all atomic positions
    double precision, intent(in) :: rcut
    integer, intent(in) :: center_ind

    integer :: naez
    naez = size(rbasis, 2)

    self%atom_index = center_ind

    call clsgen99(center_ind, rcut, lattice_vectors%nrd, lattice_vectors%rr, naez, rbasis, self%nacls, self%rcls, self%atom, self%ezoa, self%indn0, self%numn0)
#ifdef DEBUG
    write(*,*)  "Number of atoms in cluster:            ", self%nacls
    write(*,*)  "Number of inequivalent cluster atoms: :", self%numn0
#endif
  endsubroutine ! create

  !------------------------------------------------------------------------------
  !> Destroys reference cluster.
  subroutine destroyRefCluster(self)
    type(RefCluster), intent(inout) :: self

    integer :: ist
    deallocate(self%atom, self%ezoa, self%indn0, self%rcls, stat=ist)
    self%nacls = 0
    self%numn0 = 0
    self%atom_index = 0
  endsubroutine ! destroy


  ! ************************************************************************
  subroutine clsgen99(jatom, rcut, nr, rr, naez, rbasis, nacls, rcls, atom, ezoa, indn0, numn0)
    use Sorting_mod, only: dsort
    
    ! this subroutine is used to create the clusters around each atom 
    ! where repulsive potentials will be positioned.
    !
    ! strategy : 
    ! calculate the cluster of each atom by the lattice
    ! parameters avaliable. sort the atoms in a unique way :big r, big z, big y
    ! compare the positions with the previous clusters to see if there is 
    ! a difference. if not keep only previous clusters and make indexing if
    ! a new cluster is found then check dimensions and continue for the new
    ! atom.  
    !
    integer, intent(in) :: jatom  !< central atom index atom
    double precision, intent(in) :: rcut
    integer, intent(in) :: naez                 !< number of atoms in ez
    double precision, intent(in) :: rbasis(:,:) !< dim(3,naez) pos. of basis atoms in ez
    integer, intent(in) :: nr                 !< number of lattice vectors rr
    double precision, intent(in) :: rr(:,0:)  !< dim(3,0:nr) set of lattice vectors

    integer, intent(out) :: nacls
    double precision, allocatable, intent(out) :: rcls(:,:) !< dim(3,nacls) real space position of atom in cluster
    integer, intent(out) :: numn0
    integer, allocatable, intent(out) :: indn0(:) !< dim(numn0)
    integer, allocatable, intent(out) :: atom(:)  !< dim(nacls) index to atom in elem/cell at site in cluster
    integer, allocatable, intent(out) :: ezoa(:)  !< dim(nacls) index to bravais lattice  at site in cluster

    !     .. locals
    double precision, parameter :: epsshl=1.d-4
    integer :: ia, ib, iat, n, ist, macls ! macls==naez*(nr+1)
    logical :: couplmat
    double precision :: r2, tmp(3), rcut2
    
    integer(kind=4), allocatable :: indni(:) ! dim(naez)
    double precision, allocatable :: rsort(:), rg(:,:) ! dim(macls) and dim(3,macls)
    integer, allocatable :: iatom(:), iezoa(:), isort(:) ! dim(macls)
#ifdef ZYLINDRICAL_CLUSTERS
    double precision :: rcutxy2, rcutxy, rxy2
#endif

    rcut2 = (rcut + epsshl)**2
#ifdef ZYLINDRICAL_CLUSTERS
    rcutxy = rcut ! with this configuration, rcutxy == rcut, the cluster becomes spherical anyway
    rcutxy2 = (rcutxy + epsshl)**2
#endif
    
    macls = naez*(nr+1)
    allocate(rsort(macls), rg(3,macls), iatom(macls), iezoa(macls), isort(macls), indni(naez), stat=ist)
    if (ist /= 0) stop 'RefCluster: gen: allocation of [rsort, rg, iatom, iezoa, isort, indni] failed!'
    
    numn0 = 0 ! counter for unique atoms in cluster
    nacls = 0 ! counter for atoms in cluster
    do iat = 1, naez  ! loop in all atoms
      couplmat = .false.
      
      do n = 0, nr ! loop in all selected periodic images
        tmp(1:3) = rr(1:3,n) + rbasis(1:3,iat) - rbasis(1:3,jatom)
        
#ifdef ZYLINDRICAL_CLUSTERS
        rxy2 = tmp(1)*tmp(1) + tmp(2)*tmp(2) ! prepare for a cylindrical cluster construction with the cylinder in z-direction
        r2   = tmp(3)*tmp(3) + rxy2
        if (rxy2 <= rcutxy2 .and. r2 <= rcut2) then
#else
        r2 = tmp(1)*tmp(1) + tmp(2)*tmp(2) + tmp(3)*tmp(3)
        if (r2 <= rcut2) then
#endif
          nacls = nacls + 1
          iatom(nacls) = iat ! store the atom in elem cell
          iezoa(nacls) = n ! store the bravais vector
          rg(1:3,nacls) = tmp(1:3)
          rsort(nacls) = 1.d9*sqrt(r2) + dot_product([1.d6, 1.d3, 1.d0], tmp)
          couplmat = .true.
        endif ! inside
      enddo ! n loop in bravais
      
      if (couplmat) then ! atom iat is in the cluster with at least one periodic image 
        numn0 = numn0 + 1
        indni(numn0) = iat ! indn0 is an array of ordered global atom indices
        ! write(6,*) 'jatom,numn0,indn0',jatom,numn0(jatom),iat
      endif
      
    enddo ! iat loop in naez
    
    call dsort(rsort(1:nacls), isort(1:nacls), nacls)
    
    deallocate(rsort, rcls, atom, indn0, ezoa, stat=ist) ! ignore status
    
    allocate(indn0(numn0), stat=ist) ! todo: check status
    if (ist /= 0) stop 'RefCluster: gen: allocation of [indn0] failed!'

    indn0(:) = indni(1:numn0)
    deallocate(indni, stat=ist) ! ignore status
    
    allocate(rcls(3,nacls), atom(nacls), ezoa(nacls), stat=ist) ! todo: check status
    if (ist /= 0) stop 'RefCluster: gen: allocation of [rcls, atom, ezoa] failed!'
    
    ! now apply the correct order
    do ia = 1, nacls
      ib = isort(ia)
      rcls(1:3,ia) = rg(1:3,ib)
      atom(ia) = iatom(ib)
      ezoa(ia) = iezoa(ib)
    enddo ! ia
    
    deallocate(rg, iatom, iezoa, isort, stat=ist) ! ignore status

#ifdef DEBUG  
    write(*,'(a,2(i0,a),9999(" ",i0))') 'for atom #',jatom,'  indn0(1:',numn0,') =',indn0(1:numn0)
    write(*,'(a,3F16.12)') 'rcls(:,1) = ',rcls(:,1)
#endif
    ! todo: display some statistics
  endsubroutine clsgen99



  !*==rrgen.f    processed by SPAG 6.05Rc at 20:37 on 17 May 2004
  ! 02.08.95 *************************************************************
  subroutine rrgen(bv, rr, nr)
    !> generates a number of real space vectors to construct the clusters representing the local surrounding of the atoms in routine CLSGEN99
    use Sorting_mod, only: dsort
    
    double precision, intent(in) :: bv(3,3) ! bravais vectors or matrix
    double precision, allocatable, intent(out) :: rr(:,:) ! (3,0:nr)
    integer, intent(out) :: nr
    
! #define ORIGINAL_VERSION  
  
    !    .. Locals
    double precision :: r, rmax, rr2, rs
    integer :: i, i1, i2, i3, nn(1:3), iprint, mr
    double precision :: v(3)
#ifdef  ORIGINAL_VERSION  
    double precision :: vx(3), vy(3), vz(3), vx0(3), vy0(3), vz0(3)
#endif
    double precision, allocatable :: rabs2(:), rr1(:,:) 
    integer, allocatable :: ind(:)!(nr) ! permutation for sorting
    double precision, parameter :: epsshl = 1.d-5
#ifdef DEBUG
    character(len=*), parameter :: F8='(10x,i6,3f24.12,f15.4)' ! DEBUG
    iprint = 1
#else
    character(len=*), parameter :: F8='(10x,i6,3f12.3,f15.4)'
    iprint = 0
#endif
    ! write (6,'(5X,A,/)') '< RRGEN > : generation of real space mesh RR(NR)'

    

    do i = 1, 3 
#ifdef  ORIGINAL_VERSION  
      call scalpr(bv(1:3,i), bv(1:3,i), v(i)) ! v(i) = bv(1,i)*bv(1,i) + bv(2,i)*bv(2,i) + bv(3,i)*bv(3,i)
#else
      v(i) = sum(bv(1:3,i)**2)
#endif
    enddo ! i
    v(1:3) = sqrt(v(1:3))
    
    rmax = 5.d0
    r = 1.5d0*rmax + sqrt(sum(v(1:3)**2)) + epsshl
    rs = r*r
    
    nn = min(max(2, nint(r/v)), 12)
    
    mr = product(nn+1+nn) ! preliminary and about a factor 2 too large

    if (iprint > 0) then
      write(6, fmt="(10x,'Radius R        : ',f15.6,' (ALAT    units)')") r
      write(6, fmt="(10x,'       R**2     : ',f15.6,' (ALAT**2 units)')") rs
      write(6, fmt="(10x,'mesh divisions  : ',3i5)") nn
#ifdef DEBUG
      write(6, fmt="(10x,'max N vectors   : ',3i10)") mr
#endif
    endif ! output

    allocate(rabs2(mr), rr1(3,mr))

#ifdef  ORIGINAL_VERSION  
    call vmul(bv(1,1),dble(-nn(1)-1),vx0(1)) ! vx0(1:3) = (-nn(1)-1)*bv(1:3,1)
    call vmul(bv(1,2),dble(-nn(2)-1),vy0(1)) ! vy0(1:3) = (-nn(2)-1)*bv(1:3,2)
    call vmul(bv(1,3),dble(-nn(3)-1),vz0(1)) ! vz0(1:3) = (-nn(3)-1)*bv(1:3,3)
    call veq(vx0,vx) ! vx(1:3) = vx0(1:3)
#endif

    nr = 0
    do i1 = -nn(1), nn(1)
#ifdef  ORIGINAL_VERSION  
      call vadd(vx,bv(1,1),vx) ! vx(1:3) = vx(1:3) + bv(1:3,1)
      call veq(vy0,vy) ! vy(1:3) = vy0(1:3)
#endif
      do i2 = -nn(2), nn(2)
#ifdef  ORIGINAL_VERSION  
        call vadd(vy,bv(1,2),vy) ! vy(1:3) = vy(1:3) + bv(1:3,2)
        call veq(vz0,vz) ! vz(1:3) = vz0(1:3)
#endif    
        do i3 = -nn(3), nn(3)
#ifdef  ORIGINAL_VERSION  
          call vadd(vz,bv(1,3),vz) ! vz(1:3) = vz(1:3) + bv(1:3,3)
          call vadd(vx,vy,v) ! v(1:3) = vx(1:3) + vy(1:3)
          call vadd(v,vz,v) ! v(1:3) = v(1:3) + vz(1:3)
          call scalpr(v,v,rr2) ! rr2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
#else
          v(1:3) = i1*bv(1:3,1) + i2*bv(1:3,2) + i3*bv(1:3,3)
          rr2 = sum(v(1:3)**2)
#endif
                  
          ! todo: verify that v(1:3) == i1*bv(1:3,1) + i2*bv(1:3,2) + i3*bv(1:3,3) 
          !       ... and then get rid of the vector-subroutines vmul, vadd, veq
        
          if ((rr2 <= rs .or. abs(i1)+abs(i2)+abs(i3) <= 6) .and. rr2 > epsshl) then
            nr = nr + 1
            rr1(1:3,nr) = v(1:3)
            rabs2(nr) = rr2
          endif
        enddo ! i3
      enddo ! i2
    enddo ! i1

#ifdef DEBUG
    write(6, fmt="(10x,'found N vectors : ',3i10)") nr
#endif
    
    deallocate(rr, stat=i) ! ignore status
    allocate(rr(1:3,0:nr), ind(nr))
    
    rr = 0.d0 ! init
    rr(1:3,0) = 0.d0 ! init

    !!! for an OpenMP implementation, one could not treat the origin separatly, i.e. we can skip the (rr2 > epsshl) condition
    !!! and simply create a list of all nrd vectors (we also do need to increase nr automically since we can compute nr from nn and [i1,i2,i3]
    !!! so the vector setup region could be thread-parallel, furthermore, the sorting and reordering also can be threaded
    
    if (iprint > 0) then
      write(6, fmt="(/,10x,60('+'),/,18x,'generated real-space mesh-points (ALAT units)',/,10x,60('+'))")
      write(6, fmt="(13x,'index      x           y           z          distance  ',/,10x,60('-'))")
      write(6, fmt=F8) 0, rr(1:3,0), 0.d0
    endif

    call dsort(rabs2, ind, nr)
! #ifdef DEBUG
!     write(6, fmt='(1(a,i0),a,9999(" ",i0))'  ) __FILE__,__LINE__,' ind(:) =', ind
!     write(6, fmt='(1(a,i0),a,9999(" ",f0.3))') __FILE__,__LINE__,' rr2(:) =',rabs2(ind)
! #endif
    
    do i = 1, nr
      rr(1:3,i) = rr1(1:3,ind(i)) ! store
      if (iprint > 0 .and. (iprint > 1 .or. 9 > i .or. i > nr-9)) &
      write(6, fmt=F8) i, rr(1:3,i), sqrt(rabs2(ind(i)))
    enddo ! i

    deallocate(ind, rr1, rabs2, stat=i) ! ignore status
    
    if (iprint > 0)  write(6, fmt="(10x,60('+'))")
    
  endsubroutine rrgen



#ifdef ORIGINAL_VERSION  
  !+ORIGINAL_VERSION

  !todo: inline these very simple functions
  subroutine scalpr(x,y,z)
    ! scalsp computes the scalar product of x and y returning it into z.
    double precision :: x(*), y(*), z
    z = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)
  endsubroutine

  subroutine vadd(a,b,c)
    double precision :: a(*), b(*), c(*)
    c(1:3) = a(1:3) + b(1:3)
  endsubroutine

  subroutine veq(a,b)
    double precision :: a(*), b(*)
    b(1:3) = a(1:3)
  endsubroutine

  subroutine vmul(a,b,c)
    double precision :: a(*), b, c(*)
    c(1:3) = b*a(1:3)
  endsubroutine

  !-ORIGINAL_VERSION
#endif

! subroutine rrgencount(bv1, nrd)
!   ! **********************************************************************
!   ! * COUNTS      number of real space vectors to construct the          *
!   ! * clusters representing the local surrounding of the atoms in        *
!   ! * routine CLSGEN99                                                   *
!   ! **********************************************************************
!   integer, intent(out) :: nrd
!   double precision, intent(in) :: bv1(3,3)
!   
!   double precision :: r,r1,r2,r3,rmax,rr2,rs
!   integer :: i,j,k,n1,n2,n3
!   double precision :: v(3),vx(3),vy(3),vz(3), vx0(3),vy0(3),vz0(3)
!   double precision, parameter :: epsshl = 1.d-5
! 
!   
!   call scalpr(bv1(1,1),bv1(1,1),r1)
!   call scalpr(bv1(1,2),bv1(1,2),r2)
!   call scalpr(bv1(1,3),bv1(1,3),r3)
!   rmax = 5.d0
! 
!   r1 = sqrt(r1)
!   r2 = sqrt(r2)
!   r3 = sqrt(r3)
!   r = 1.5d0*rmax + sqrt(r1*r1 + r2*r2 + r3*r3) + epsshl
!   
!   rs = r*r
!   n1 = nint(r/r1)
!   n2 = nint(r/r2)
!   n3 = nint(r/r3)
!   
!   n1 = min(12,n1)
!   n2 = min(12,n2)
!   n3 = min(12,n3)
! 
!   n1 = max(2,n1)
!   n2 = max(2,n2)
!   n3 = max(2,n3)
!   
!   nrd = 0
! 
!   call vmul(bv1(1,1),dble(-n1-1),vx0(1))
!   call vmul(bv1(1,2),dble(-n2-1),vy0(1))
!   call vmul(bv1(1,3),dble(-n3-1),vz0(1))
!   call veq(vx0,vx)
! 
!   do i = -n1, n1
!      call vadd(vx,bv1(1,1),vx)
!      call veq(vy0,vy)
! 
!      do j = -n2, n2
!         call vadd(vy,bv1(1,2),vy)
!         call veq(vz0,vz)
! 
!         do k = -n3, n3
!            call vadd(vz,bv1(1,3),vz)
!            call vadd(vx,vy,v)
!            call vadd(v,vz,v)
!            call scalpr(v,v,rr2)
! 
!            if ((rr2 <= rs .or. abs(i)+abs(j)+abs(k) <= 6) .and. rr2 > epsshl) nrd = nrd + 1
!         enddo ! k
!      enddo ! j
!   enddo ! i
! 
! endsubroutine rrgencount


! subroutine clsgen99count(center_ind, naez, rr, rbasis, rcut, nrd, nacls)
!   ! This subroutine is used to create the clusters around each atom 
!   ! where repulsive potentials will be positioned.
!   !
!   ! STRATEGY : 
!   ! Calculate the cluster of each atom by the lattice
!   ! parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
!   ! compare the positions with the previous clusters to see if there is 
!   ! a difference. If not keep only previous clusters and make indexing if
!   ! a new cluster is found then check dimensions and continue for the new
!   ! atom.  
!   !
!   !     .. arguments
!   !
!   integer, intent(in) :: center_ind 
!   integer, intent(in) :: nrd
!   integer, intent(out) :: nacls
! 
!   double precision, intent(in) :: rcut
!   double precision, intent(in) :: &
!        rbasis(3,*), &      ! pos. of basis atoms in EZ 
!   rr(3,0:nrd)              ! set of lattice vectors
! 
!   integer, intent(in) :: naez !< number of atoms in EZ
!   
!   !     .. locals
!   integer :: nr                       ! number of lattice vectors RR
!   integer i, &
!           na,number,n, &
!           jatom
! 
!   double precision r2,epsshl, tmp(3)
!   double precision rcut2,rcutxy2,rxy2,rcutxy
!   !
!   data     epsshl   / 1.0d-4 /
! 
!   nr = nrd
! 
!   rcutxy = rcut ! only spherical clusters allowed
! 
!   rcutxy2 = (rcutxy + epsshl)*(rcutxy+epsshl)
!   rcut2   = (rcut+epsshl)*(rcut+epsshl)
! 
!   !do jatom = 1,naez       ! loop in all atoms
!   jatom = center_ind       ! modification to original: treat only one atom
!      number = 0           ! counter for atoms in cluster
!      do na = 1,naez  ! loop in all atoms
!         do n=0,nr    ! loop in all bravais vectors    
!            do i=1,3
!               tmp(i) = rr(i,n)+rbasis(i,na)-rbasis(i,jatom)
!            enddo
!            rxy2 =  tmp(1)**2+tmp(2)**2
!            r2   =  tmp(3)**2 + tmp(1)**2+tmp(2)**2
! 
!            if ((rxy2 <= rcutxy2).and.(r2 <= rcut2))  then
!               number = number + 1
!               !atom(number) = na ! store the atom in elem cell
!               !ezoa(number) = n ! store the bravais vector
!            endif
!         enddo              ! N loop in bravais
! 
!      enddo                 ! NA loop in NAEZ
! 
!   nacls = number ! return number of atoms in cluster
! 
! endsubroutine clsgen99count



endmodule

#ifdef TEST_CLUSTERS__
program test_clusters
  use RefCluster_mod, only: RefCluster, LatticeVectors!, create, destroy
  use RefCluster_mod, only: createLatticeVectors, destroyLatticeVectors ! deprecated
  use RefCluster_mod, only: createRefCluster, destroyRefCluster ! deprecated

  double precision bravais(3,3)
  double precision rbasis(3,4)

  type(RefCluster) :: ref_cluster
  type(LatticeVectors), target :: lattice_vectors

  integer ii
  double precision :: rcut = 1.5
  integer :: center_ind = 4
  double precision :: dist

  bravais = reshape([ 1., 0., 0., &
                      0., 1., 0., &
                      0., 0., 1. ], [3, 3])

  rbasis = reshape([ .0, .0, .0, &
                     .5, .5, .0, &
                     .5, .0, .5, &
                     .0, .5, .5  ], [3, 4])

  call createLatticeVectors(lattice_vectors, bravais)
  call createRefCluster(ref_cluster, lattice_vectors, rbasis, rcut, center_ind)

  write(*,*)  "Number of lattice vectors:             ", ref_cluster%lattice_vectors%nrd
  write(*,*)  "Number of atoms in cluster:            ", ref_cluster%nacls
  write(*,*)  "Number of inequivalent cluster atoms: :", ref_cluster%numn0

  do ii = 1, ref_cluster%nacls
    dist = sqrt(ref_cluster%rcls(1,ii)**2 + ref_cluster%rcls(2,ii)**2 + ref_cluster%rcls(3,ii)**2)
    write(*,"(79('-'))")
    write(*, 92000) "Cluster atom", ii
    write(*, 91000) "Basis atom / dist. ", ref_cluster%atom(ii), dist
    write(*, 90000) "Lattice vector", ref_cluster%lattice_vectors%rr(:, ref_cluster%ezoa(ii))
    write(*, 90000) "Position in basis", rbasis(:, ref_cluster%atom(ii))
    write(*, 90000) "Position from center", ref_cluster%rcls(:, ii)
    if (dist > rcut) then
      write(*,*) "ERROR in cluster generation."
      STOP
    endif
  enddo

  call destroyLatticeVectors(lattice_vectors)
  call destroyRefCluster(ref_cluster)

90000 format (A20, X, 3(F13.6, X))
91000 format (A20, X, I7, X, F13.6)
92000 format (A20, X, I7)
endprogram
#endif

