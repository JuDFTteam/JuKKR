!> Definitions
!> RefCluster: a cluster of atoms with a distance smaller than rcut to a
!>             central atom
!>             periodic boundary conditions are taken into account:
!>             mirror images of atoms can be included in the cluster
!>
!> LatticeVectors: linear combinations of Bravais vectors with integer
!>                 coefficients
!>
!> Note: if N (=number of atoms) reference clusters are created,
!>       then this algorithm scales as O(N**2)


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
    type(LatticeVectors), pointer :: lattice_vectors
    double precision, allocatable :: rcls(:,:) !< positions relative to center
    integer, allocatable :: atom(:)  !< basis atom indices of cluster atoms
    integer, allocatable :: ezoa(:)  !< points into lattice_vectors%rr
    integer, allocatable :: indn0(:) !< indices of inequivalent cluster atoms
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
  !write(*,*) "Num. real space vectors: ", lattice_vectors%nrd

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
!> type(LatticeVectors), target :: lattice_vectors
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

  type(LatticeVectors), target, intent(in) :: lattice_vectors
  double precision, intent(in) :: rbasis(:,:)
  double precision, intent(in) :: rcut
  integer, intent(in) :: center_ind


  integer nacls
  integer num_atoms

  num_atoms = size(rbasis,2)

  self%lattice_vectors => lattice_vectors
  self%atom_index = center_ind

  call clsgen99count(center_ind, num_atoms ,lattice_vectors%rr,rbasis,rcut, &
                     lattice_vectors%nrd, nacls)


  allocate (self%atom(nacls), self%ezoa(nacls), self%indn0(nacls))
  allocate (self%rcls(3, nacls))

  call clsgen99(center_ind, num_atoms , lattice_vectors%rr , rbasis, &
                self%atom,self%ezoa, &
                self%rcls,rcut, &
                self%numn0,self%indn0, &
                lattice_vectors%nrd, nacls)

  !write(*,*)  "Number of atoms in cluster:            ", nacls
  !write(*,*)  "Number of inequivalent cluster atoms: :", numn0

  self%nacls = nacls
endsubroutine

!------------------------------------------------------------------------------
!> Destroys reference cluster.
subroutine destroyRefCluster(self)
  type(RefCluster), intent(inout) :: self

  integer :: ist
  deallocate(self%atom, self%ezoa, self%indn0, self%rcls, stat=ist)
  self%nacls = 0
endsubroutine

!todo: inline these very simple functions
subroutine scalpr(x,y,z)
  ! scalsp computes the scalar product of x and y returning it into z.
  double precision x(*), y(*), z
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

subroutine clsgen99count(center_ind, naez,rr,rbasis,rcut, nrd, nacls)
  ! This subroutine is used to create the clusters around each atom 
  ! where repulsive potentials will be positioned.
  !
  ! STRATEGY : 
  ! Calculate the cluster of each atom by the lattice
  ! parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
  ! compare the positions with the previous clusters to see if there is 
  ! a difference. If not keep only previous clusters and make indexing if
  ! a new cluster is found then check dimensions and continue for the new
  ! atom.  
  !
  !     .. arguments
  !
  integer, intent(in) :: center_ind 
  integer, intent(in) :: nrd
  integer, intent(out) :: nacls

  double precision rcut,rcutxy
  double precision &
       rbasis(3,*), &      ! pos. of basis atoms in EZ 
  rr(3,0:nrd)              ! set of lattice vectors
  !
  integer naez, &          ! number of atoms in EZ
  nr                       ! number of lattice vectors RR
  !     .. locals
  !
  integer i, &
          na,number,n, &
          jatom

  double precision r2,epsshl, tmp(3)
  double precision rcut2,rcutxy2,rxy2 
  !
  data     epsshl   / 1.0d-4 /

  nr = nrd

  rcutxy = rcut ! only spherical clusters allowed

  rcutxy2 = (rcutxy+epsshl)*(rcutxy+epsshl)
  rcut2   = (rcut+epsshl)*(rcut+epsshl)

  !do jatom = 1,naez       ! loop in all atoms
  jatom = center_ind       ! modification to original: treat only one atom
     number = 0           ! counter for atoms in cluster
     do na = 1,naez  ! loop in all atoms
        do n=0,nr    ! loop in all bravais vectors    
           do i=1,3
              tmp(i) = rr(i,n)+rbasis(i,na)-rbasis(i,jatom)
           enddo
           rxy2 =  tmp(1)**2+tmp(2)**2
           r2   =  tmp(3)**2 + tmp(1)**2+tmp(2)**2

           if ((rxy2 <= rcutxy2).and.(r2 <= rcut2))  then
              number = number + 1
              !atom(number) = na ! store the atom in elem cell
              !ezoa(number) = n ! store the bravais vector
           endif
        enddo              ! N loop in bravais

     enddo                 ! NA loop in NAEZ

  nacls = number ! return number of atoms in cluster

endsubroutine clsgen99count


! ************************************************************************
subroutine clsgen99(jatom, naez,rr,rbasis, atom,ezoa, rcls, rcut, numn0, indn0, nrd, nacls)
  use Sorting_mod, only: dsort
  
  ! This subroutine is used to create the clusters around each atom 
  ! where repulsive potentials will be positioned.
  !
  ! STRATEGY : 
  ! Calculate the cluster of each atom by the lattice
  ! parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
  ! compare the positions with the previous clusters to see if there is 
  ! a difference. If not keep only previous clusters and make indexing if
  ! a new cluster is found then check dimensions and continue for the new
  ! atom.  
  !
  !     .. arguments
  !
  integer, intent(in) :: jatom  !< modification to original: treat only one atom
  integer, intent(in) :: nrd
  integer, intent(in) :: nacls
  double precision rcut
  double precision rbasis(3,*), &      ! pos. of basis atoms in EZ 
  rcls(3,nacls), &        ! real space position of atom in cluster
  rr(3,0:nrd)              ! set of lattice vectors
  integer naez, &          ! number of atoms in EZ
  nr                       ! number of lattice vectors RR
  integer numn0, indn0(nacls), &
  atom(nacls), &         ! index to atom in elem/cell at site in cluster
  ezoa(nacls)           ! index to bravais lattice  at site in cluster

  !     .. locals
  integer :: i, na, number, n, iat!, icluster
! integer, parameter  :: nclsd = 1
! integer iatcls(nclsd)
  logical(kind=1) :: icouplmat(naez) ! this array could be large
  double precision :: r2, tmp(3), rcut2, rcutxy2, rxy2, rcutxy
  double precision rsort(nacls), rg(3,nacls) ! for sorting
  integer ia, ib, pos, iatom(nacls), iezoa(nacls), isort(nacls)
  double precision, parameter :: epsshl=1.d-4

  nr = nrd
  indn0(:) = -1

  rcutxy = rcut ! only spherical clusters allowed

! icluster = 1
! iatcls(:) = 0

  rcutxy2 = (rcutxy + epsshl)**2
  rcut2   = (rcut + epsshl)**2
  
  number = 0           ! counter for atoms in cluster
  do na = 1, naez  ! loop in all atoms
    do n = 0, nr    ! loop in all bravais vectors    
      tmp(1:3) = rr(1:3,n) + rbasis(1:3,na) - rbasis(1:3,jatom)
      rxy2 = tmp(1)*tmp(1) + tmp(2)*tmp(2)
      r2   = tmp(3)*tmp(3) + rxy2

      if (rxy2 <= rcutxy2 .and. r2 <= rcut2)  then
        number = number + 1
        if (number > nacls) stop '   < CLSGEN99 > Dimension NACLSD too small'
        atom(number) = na ! store the atom in elem cell
        ezoa(number) = n ! store the bravais vector
        rcls(1:3,number) = tmp(1:3)
      endif
    enddo ! N loop in bravais
  enddo ! NA loop in NAEZ

  !     sort the atoms of the cluster in increasing order. First by distance
  !     Then by z then by y
  !
  do ia = 1, number
    rsort(ia) = sqrt(sum(rcls(1:3,ia)**2))
!       rsort(ia) = 100 000 000.d0*rsort(ia) + 100 000.d0*rcls(3,ia) + 100.d0*rcls(2,ia) + 0.1d0*rcls(1,ia)
    rsort(ia) = 1.d9*rsort(ia) + 1.d6*rcls(3,ia) + 1.d3*rcls(2,ia) + rcls(1,ia)
  enddo ! ia

  call dsort(rsort, isort, number, pos)
  !     Rearange exchange ia with ib
  ! MAP temporarily to another array
  do ia = 1, number
    rg(1:3,ia) = rcls(1:3,ia)
    iatom(ia) = atom(ia)
    iezoa(ia) = ezoa(ia)
  enddo ! ia
  
  ! Now use correct order
  do ia = 1, number
    ib = isort(ia)
    rcls(1:3,ia) = rg(1:3,ib)
    atom(ia) = iatom(ib)
    ezoa(ia) = iezoa(ib)
  enddo ! ia

  do iat = 1, naez
    icouplmat(iat) = .false.
    do i = 1, number 
      if (atom(i) == iat) icouplmat(iat) = .true.
    enddo ! i
  enddo ! iat

  numn0 = 0
  do iat = 1, naez
    if (icouplmat(iat)) then
      numn0 = numn0 + 1
      indn0(numn0) = iat
      ! write(6,*) 'jatom,numn0,indn0',jatom,numn0(jatom),iat
    endif
  enddo ! iat

endsubroutine clsgen99


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


!*==rrgen.f    processed by SPAG 6.05Rc at 20:37 on 17 May 2004
! 02.08.95 *************************************************************
subroutine rrgen(bv1, rr, nrd)
  use Sorting_mod, only: dsort
  ! **********************************************************************
  ! * generates a number of real space vectors to construct the          *
  ! * clusters representing the local surrounding of the atoms in        *
  ! * routine CLSGEN99                                                   *
  ! **********************************************************************
  double precision, intent(in) :: bv1(3,3)
  double precision, allocatable, intent(out) :: rr(:,:) ! (3,0:nrd)
  integer, intent(out) :: nrd
  
  !    .. Locals
  integer :: nr
  double precision :: r, rmax, rr2, rs
  integer :: i, i1, i2, i3, nn(1:3), iprint
  double precision :: v(3), vx(3), vy(3), vz(3), vx0(3), vy0(3), vz0(3)
  double precision, allocatable :: rabs2(:), rr1(:,:) 
  integer, allocatable :: ind(:)!(nrd) ! permutation for sorting
  double precision, parameter :: epsshl = 1.d-5
  character(len=*), parameter :: F8='(10x,i6,3f12.3,f15.4)'
  
  ! write (6,'(5X,A,/)') '< RRGEN > : generation of real space mesh RR(NR)'
  iprint = 0

  do i = 1, 3 
!   call scalpr(bv1(1:3,i), bv1(1:3,i), v(i))
    v(i) = dot_product(bv1(1:3,i), bv1(1:3,i))
  enddo ! i
  rmax = 5.d0

  v = sqrt(v)
  r = 1.5d0*rmax + sqrt(sum(v(1:3)**2)) + epsshl
  rs = r*r
  
  nn = min(max(2, nint(r/v)), 12)

  if (iprint == 1) then
    write(6, fmt="(10x,'Radius R        : ',f15.6,' (ALAT    units)')") r
    write(6, fmt="(10x,'       R**2     : ',f15.6,' (ALAT**2 units)')") rs
    write(6, fmt="(10x,'mesh divisions  : ',3i5)") nn
  endif ! output

  nrd = product(nn+1+nn) ! preliminary and about a factor 2 too large
  allocate(rabs2(nrd), rr1(3,nrd))
  
  nr = 0

  call vmul(bv1(1,1),dble(-nn(1)-1),vx0(1))
  call vmul(bv1(1,2),dble(-nn(2)-1),vy0(1))
  call vmul(bv1(1,3),dble(-nn(3)-1),vz0(1))
  call veq(vx0,vx)
  
  do i1 = -nn(1), nn(1)
    call vadd(vx,bv1(1,1),vx)
    call veq(vy0,vy)

    do i2 = -nn(2), nn(2)
      call vadd(vy,bv1(1,2),vy)
      call veq(vz0,vz)

      do i3 = -nn(3), nn(3)
        call vadd(vz,bv1(1,3),vz)
        call vadd(vx,vy,v)
        call vadd(v,vz,v)
        call scalpr(v,v,rr2)
        
        ! todo: verify that v == i1*bv1(1:3,1) + i2*bv1(1:3,2) + i3*bv1(1:3,3) 
        !       ... and then get rid of the vector-subroutines vmul, vadd, veq
      
        if ((rr2 <= rs .or. abs(i1)+abs(i2)+abs(i3) <= 6) .and. rr2 > epsshl) then
          nr = nr + 1
          rr1(1:3,nr) = v(1:3)
          rabs2(nr) = rr2
        endif
      enddo ! i3
    enddo ! i2
  enddo ! i1

  deallocate(rr, stat=i) ! ignore status
  allocate(rr(1:3,0:nr), ind(nr))
  nrd = nr ! store
  rr = 0.d0 ! init
  rr(1:3,0) = 0.d0 ! init

  if (iprint > 0) then
    write(6, fmt="(/,10x,60('+'),/,18x,'generated real-space mesh-points (ALAT units)',/,10x,60('+'))")
    write(6, fmt="(13x,'index      x           y           z          distance  ',/,10x,60('-'))")
    write(6, fmt=F8) 0, rr(1:3,0), 0.d0
  endif

  call dsort(rabs2, ind, nr)
  do i = 1, nr
    rr(1:3,i) = rr1(1:3,ind(i)) ! store
    if (iprint > 0) write(6, fmt=F8) i, rr(1:3,i), sqrt(rabs2(ind(i)))
  enddo ! i

  deallocate(ind, rr1, rabs2, stat=i) ! ignore status
  
  if (iprint > 0)  write(6, fmt="(10x,60('+'))")
  
endsubroutine rrgen

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

