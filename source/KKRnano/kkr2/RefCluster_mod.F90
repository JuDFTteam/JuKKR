!> Definitions
!> RefCluster: a cluster of atoms with a distance smaller than rcut to a
!>             central atom
!>             periodic boundary conditions are taken into account:
!>             mirror images of atoms can be included in the cluster
!>
!> Note: if N (=number of atoms) reference clusters are created, then this algorithm scales as O(N**2)
! #define DEBUG

module RefCluster_mod
  implicit none
  private
  
  public :: RefCluster, create, destroy

  type RefCluster
    integer :: numn0 !< number of inequivalent cluster atoms
    integer :: nacls !< number of cluster atoms
    integer :: atom_index !< basis atom index of central cluster atom
    !> reference to LatticeVectors datastructure
    integer, allocatable :: atom(:)  !< dim(nacls) basis atom indices of cluster atoms
    integer, allocatable :: indn0(:) !< dim(numn0) indices of inequivalent cluster atoms
    integer, allocatable :: ezoa(:)  !< dim(nacls) points into lattice_vectors%rr, which periodic image
    double precision, allocatable :: rcls(:,:) !< dim(3,nacls) positions relative to center
  endtype

  interface create
    module procedure createRefCluster
  endinterface
  
  interface destroy
    module procedure destroyRefCluster
  endinterface
  
  contains

  !------------------------------------------------------------------------------
  !> Creates reference cluster.
  !>
  !> Usage example:
  !> type(LatticeVectors) :: lattice_vectors
  !> type(RefCluster) :: ref_cluster
  !> ! ...
  !> call createLatticeVectors(lattice_vectors, bravais)
  !> call createRefCluster(ref_cluster, lattice_vectors, rbasis, rcut, global_atom_id)
  !> ! ... some code
  !> call destroyRefCluster(ref_cluster)
  !> call destroyLatticeVectors(lattice_vectors)
  !
  subroutine createRefCluster(self, lattice_vectors, rbasis, rcut, global_atom_id)
    type(RefCluster), intent(inout) :: self

    double precision, intent(in) :: lattice_vectors(:,0:) !< dim(3,0:nrd)
    double precision, intent(in) :: rbasis(:,:) !< dim(3,naez) all atomic positions
    double precision, intent(in) :: rcut
    integer, intent(in) :: global_atom_id !< global atom index

    self%atom_index = global_atom_id

    call clsgen99(global_atom_id, rcut, lattice_vectors, rbasis, self%nacls, self%rcls, self%atom, self%ezoa, self%indn0, self%numn0)
#ifdef DEBUG
    write(*,*)  "Number of atoms in cluster:            ", self%nacls
    write(*,*)  "Number of inequivalent cluster atoms: :", self%numn0
#endif
  endsubroutine ! create

  !------------------------------------------------------------------------------
  !> Destroys reference cluster.
  elemental subroutine destroyRefCluster(self)
    type(RefCluster), intent(inout) :: self

    integer :: ist ! ignore status
    deallocate(self%atom, self%ezoa, self%indn0, self%rcls, stat=ist)
    self%nacls = 0
    self%numn0 = 0
    self%atom_index = 0
  endsubroutine ! destroy


  !------------------------------------------------------------------------------
  subroutine clsgen99(jatom, rcut, rr, rbasis, nacls, rcls, atom, ezoa, indn0, numn0)
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
    integer, intent(in) :: jatom !< central atom index atom into rbasis(1:3,:)
    double precision, intent(in) :: rcut
    double precision, intent(in) :: rbasis(:,:) !< dim(3,naez) pos. of basis atoms in ez
    double precision, intent(in) :: rr(:,0:)  !< dim(3,0:nr-1) set of lattice vectors

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
    integer :: naez  !< number of atoms in ez
    integer :: nr    !< number of lattice vectors rr
    
    integer(kind=4), allocatable :: indni(:) ! dim(naez)
    double precision, allocatable :: rsort(:), rg(:,:) ! dim(macls) and dim(3,macls)
    integer, allocatable :: iatom(:), iezoa(:), isort(:) ! dim(macls)
#ifdef CYLINDRICAL_CLUSTERS
    double precision :: rcutxy2, rcutxy, rxy2
#endif

    naez = size(rbasis, 2)
    nr = size(rr, 2)

    rcut2 = (rcut + epsshl)**2
#ifdef CYLINDRICAL_CLUSTERS
    rcutxy = rcut ! with this configuration, rcutxy == rcut, the cluster becomes spherical anyway
    rcutxy2 = (rcutxy + epsshl)**2
#endif

    macls = naez*nr
    allocate(rsort(macls), rg(3,macls), iatom(macls), iezoa(macls), isort(macls), indni(naez), stat=ist)
    if (ist /= 0) stop 'RefCluster: gen: allocation of [rsort, rg, iatom, iezoa, isort, indni] failed!'
    
    numn0 = 0 ! counter for unique atoms in cluster
    nacls = 0 ! counter for atoms in cluster
    do iat = 1, naez  ! loop in all atoms
      couplmat = .false.
      
      do n = 0, nr-1 ! loop in all selected periodic images
        tmp(1:3) = rr(1:3,n) + rbasis(1:3,iat) - rbasis(1:3,jatom)
        
#ifdef CYLINDRICAL_CLUSTERS
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
      endif ! couplmat
      
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
  endsubroutine ! clsgen99

endmodule ! RefCluster_mod

#ifdef TEST_CLUSTERS__
program test_clusters
  use LatticeVectors_mod, only: LatticeVectors, create, destroy
  use RefCluster_mod, only: RefCluster, create, destroy

  double precision :: bravais(3,3), rbasis(3,4), dist, rcut = 1.5
  integer :: ii, global_atom_id = 4

  type(RefCluster) :: ref_cluster
  type(LatticeVectors) :: lattice_vectors

  bravais = reshape([ 1., 0., 0., &
                      0., 1., 0., &
                      0., 0., 1. ], [3, 3])

  rbasis = reshape([ .0, .0, .0, &
                     .5, .5, .0, &
                     .5, .0, .5, &
                     .0, .5, .5  ], [3, 4])

  call create(lattice_vectors, bravais)
  call create(ref_cluster, lattice_vectors%rr, rbasis, rcut, global_atom_id)

  write(*,*) "Number of lattice vectors:             ", ref_cluster%lattice_vectors%nrd
  write(*,*) "Number of atoms in cluster:            ", ref_cluster%nacls
  write(*,*) "Number of inequivalent cluster atoms: :", ref_cluster%numn0

  do ii = 1, ref_cluster%nacls
    dist = sqrt(ref_cluster%rcls(1,ii)**2 + ref_cluster%rcls(2,ii)**2 + ref_cluster%rcls(3,ii)**2)
    write(*,"(79('-'))")
    write(*, '(a20,x,i7)')         "Cluster atom", ii
    write(*, '(a20,x,i7,x,f13.6)') "Basis atom / dist. ", ref_cluster%atom(ii), dist
    write(*, '(a20,x,3(f13.6,x))') "Lattice vector", ref_cluster%lattice_vectors%rr(:, ref_cluster%ezoa(ii))
    write(*, '(a20,x,3(f13.6,x))') "Position in basis", rbasis(:, ref_cluster%atom(ii))
    write(*, '(a20,x,3(f13.6,x))') "Position from center", ref_cluster%rcls(:, ii)
    if (dist > rcut) stop "ERROR in cluster generation."
  enddo ! ii

  call destroy(lattice_vectors)
  call destroy(ref_cluster)

endprogram ! test
#endif

