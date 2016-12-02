!> Definitions
!> RefCluster: a cluster of atoms with a distance smaller than radius to a
!>             central atom
!>             periodic boundary conditions are taken into account:
!>             mirror images of atoms can be included in the cluster
!>
!> Note: if periodic_image_index (=number of atoms) reference clusters are created, then this algorithm scales as O(periodic_image_index**2)
! #define DEBUG

module RefCluster_mod
  implicit none
  private
  
  public :: RefCluster, create, destroy

  type RefCluster
    integer(kind=4)               :: source_atom_index !< global atom index of central cluster atom
    integer                       :: nacls !< number of target atoms in the cluster (periodic images included)
    integer(kind=4), allocatable  :: atom(:)  !< ddim(nacls) atom index of the target atom in cluster
    !> reference to LatticeVectors datastructure
    integer(kind=2), allocatable  :: ezoa(:)  !< dim(nacls) index to bravais lattice vector at site in cluster
    integer                       :: numn0 !< number of inequivalent target atoms in the cluster
    integer(kind=4), allocatable  :: indn0(:) !< dim(numn0) atom indices of the set of inequivalent target atoms
    double precision, allocatable :: rcls(:,:) !< dim(3,nacls) real space position of atoms in the cluster relative to center
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
  subroutine createRefCluster(self, lattice_vectors, atom_positions, radius, source_atom_index)
    use Sorting_mod, only: dsort
    type(RefCluster), intent(inout) :: self

    double precision, intent(in) :: lattice_vectors(:,0:) !< dim(1:3,0:nr-1)
    double precision, intent(in) :: atom_positions(:,:) !< dim(1:3,natoms) all atomic positions
    double precision, intent(in) :: radius
    integer, intent(in)          :: source_atom_index !< global atom index
    
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

    !  .. locals
    double precision, parameter :: epsshl=1.d-4
    integer :: ia, ib, target_atom_index, periodic_image_index, ist, macls ! macls==natoms*nr
    logical :: couplmat
    double precision :: r2, tmp(3), rcut2
#ifdef CYLINDRICAL_CLUSTERS
    double precision :: rcutxy2, rcutxy, rxy2
#endif
    integer :: natoms  !< number of atoms in ez
    integer :: nr    !< number of lattice vectors
    
    integer(kind=4),   allocatable :: t_indn0(:) ! dim(natoms)
    integer(kind=4),   allocatable :: t_atom(:)  ! dim(macls)
    integer(kind=4),   allocatable :: iperm(:)   ! dim(macls)
    integer(kind=2),   allocatable :: t_ezoa(:)  ! dim(macls)
    real(kind=8),      allocatable :: rsort(:)   ! dim(macls) ! this could be float32 if sorting routine available, ToDo: check influence
    integer, parameter :: ConstructVectorsOnce = 0 ! 1 or 0, 0:saves memory
    double precision, allocatable :: rg(:,:) ! dim(3,macls)
    double precision, parameter   :: LiftDegeneracy(3) = [1.d5, 1.d2, 1.d0] ! equivalent to JM code

    call destroy(self)
  
    self%source_atom_index = source_atom_index

    if (ConstructVectorsOnce**2 /= ConstructVectorsOnce) stop 'RefCluster: internal error!' ! must be 0 or 1
    
    natoms = size(atom_positions, 2)
    nr = size(lattice_vectors, 2)
    
    rcut2 = (radius + epsshl)**2
#ifdef CYLINDRICAL_CLUSTERS
    rcutxy = radius ! with this configuration, rcutxy == radius, the cluster becomes spherical anyway
    rcutxy2 = (rcutxy + epsshl)**2
#endif

    macls = natoms*nr ! this can be very large, can we use less? ToDo: user may propose a reduction factor
    allocate(rsort(macls), rg(3,1+(macls-1)*ConstructVectorsOnce), t_atom(macls), t_ezoa(macls), t_indn0(natoms), stat=ist)
    if (ist /= 0) stop 'RefCluster: gen: allocation of [rsort, rg, t_atom, t_ezoa, t_indn0] failed!'

    self%numn0 = 0 ! counter for inequivalent atoms in cluster
    self%nacls = 0 ! counter for atoms in cluster
    do target_atom_index = 1, natoms ! loop in all global atom ids of target atoms
      couplmat = .false.
      
      do periodic_image_index = 0, nr-1 ! loop in all selected periodic images
        tmp(1:3) = position_vector(periodic_image_index, target_atom_index, source_atom_index)

#ifndef CYLINDRICAL_CLUSTERS
        r2 = tmp(1)*tmp(1) + tmp(2)*tmp(2) + tmp(3)*tmp(3)
        if (r2 <= rcut2) then
#else
        rxy2 = tmp(1)*tmp(1) + tmp(2)*tmp(2) ! prepare for a cylindrical cluster construction with the cylinder in z-direction
        r2   = tmp(3)*tmp(3) + rxy2
        if (rxy2 <= rcutxy2 .and. r2 <= rcut2) then
#endif
          self%nacls = self%nacls + 1
          t_atom(self%nacls) = target_atom_index ! store the index of the target atom in the unit cell
          t_ezoa(self%nacls) = periodic_image_index ! store the index of the periodic image
          if (ConstructVectorsOnce == 1) rg(1:3,self%nacls) = tmp(1:3) ! store for later usage
          rsort(self%nacls) = 1.d9*sqrt(r2) + dot_product(LiftDegeneracy, tmp) ! the sqrt is not needed for sorting, ToDo: remove and check influence
          couplmat = .true. ! yes, target_atom_index is close enough to source_atom_index
        endif ! inside
      enddo ! periodic_image_index loop in bravais

      if (couplmat) then ! atom target_atom_index is in the cluster with at least one periodic image 
        self%numn0 = self%numn0 + 1 ! add to the set of inequivalent target atoms
        t_indn0(self%numn0) = target_atom_index ! indn0 becomes an array of ordered global atom indices of inequivalent atoms
      endif ! couplmat

    enddo ! target_atom_index loop in natoms

    allocate(iperm(self%nacls), stat=ist) ; if (ist /= 0) stop 'RefCluster: gen: allocation of [iperm] failed!'

    call dsort(rsort(1:self%nacls), iperm, self%nacls)

    deallocate(rsort, stat=ist) ! ignore status

    allocate(self%indn0(self%numn0), stat=ist) ; if (ist /= 0) stop 'RefCluster: gen: allocation of [indn0] failed!'
    self%indn0(:) = t_indn0(:self%numn0)
    deallocate(t_indn0, stat=ist) ! ignore status

    allocate(self%rcls(3,self%nacls), self%atom(self%nacls), self%ezoa(self%nacls), stat=ist)
    if (ist /= 0) stop 'RefCluster: gen: allocation of [rcls, atom, ezoa] failed!'

    ! now apply the correct order
    do ia = 1, self%nacls
      ib = iperm(ia)
      self%atom(ia) = t_atom(ib)
      self%ezoa(ia) = t_ezoa(ib)
      if (ConstructVectorsOnce == 1) then 
        self%rcls(1:3,ia) = rg(1:3,ib) ! take the preconstructed and stored vector
      else
        periodic_image_index = self%ezoa(ia)
        self%rcls(1:3,ia) = position_vector(periodic_image_index, self%atom(ia), source_atom_index)
      endif
    enddo ! ia

    deallocate(rg, t_atom, t_ezoa, iperm, stat=ist) ! ignore status

#ifdef DEBUG  
    write(*,'(a,2(i0,a),9999(" ",i0))') 'for atom #',source_atom_index,'  indn0(1:',self%numn0,') =',self%indn0(:)
    write(*,'(a,3F16.12)') 'rcls(:,1) = ',rcls(1:3,1)
#endif
    ! todo: display some statistics
#ifdef DEBUG
    write(*,*)  "Number of atoms in cluster:            ", self%nacls
    write(*,*)  "Number of inequivalent cluster atoms: :", self%numn0
#endif
    
  contains 
  
    function position_vector(peridic_image, target_atom, source_atom) result(v)
      double precision :: v(1:3)
      integer, intent(in) :: peridic_image, target_atom, source_atom
      
      v(1:3) = lattice_vectors(1:3,peridic_image) + atom_positions(1:3,target_atom) - atom_positions(1:3,source_atom) ! construct the vector
    endfunction ! position_vector

  endsubroutine ! create

  !------------------------------------------------------------------------------
  !> Destroys reference cluster.
  elemental subroutine destroyRefCluster(self)
    type(RefCluster), intent(inout) :: self

    integer :: ist ! ignore status
    deallocate(self%atom, self%ezoa, self%indn0, self%rcls, stat=ist)
    self%nacls = 0
    self%numn0 = 0
    self%source_atom_index = 0
  endsubroutine ! destroy
  
endmodule ! RefCluster_mod

#ifdef TEST_CLUSTERS__
program test_clusters
  use LatticeVectors_mod, only: LatticeVectors, create, destroy
  use RefCluster_mod, only: RefCluster, create, destroy

  double precision :: bravais(3,3), atom_positions(3,4), dist2, radius = 1.5
  integer :: ii, source_atom_index = 4

  type(RefCluster) :: ref_cluster
  type(LatticeVectors) :: lattice_vectors

  bravais = reshape([ 1., 0., 0., &
                      0., 1., 0., &
                      0., 0., 1. ], [3, 3])

  atom_positions = reshape([ .0, .0, .0, &
                             .5, .5, .0, &
                             .5, .0, .5, &
                             .0, .5, .5  ], [3, 4]) ! fcc positions

  call create(lattice_vectors, bravais)
  call create(ref_cluster, lattice_vectors%rr, atom_positions, radius, source_atom_index)

  write(*,*) "Number of lattice vectors:             ", ref_cluster%lattice_vectors%nrd
  write(*,*) "Number of atoms in cluster:            ", ref_cluster%nacls
  write(*,*) "Number of inequivalent cluster atoms: :", ref_cluster%numn0

  do ii = 1, ref_cluster%nacls
    dist2 = sum(ref_cluster%rcls(1:3,ii)**2)
    write(*,"(79('-'))")
    write(*, '(a20,x,i7)')         "Cluster atom", ii
    write(*, '(a20,x,i7,x,f13.6)') "Basis atom / dist. ", ref_cluster%atom(ii), sqrt(dist2)
    write(*, '(a20,x,3(f13.6,x))') "Lattice vector", ref_cluster%lattice_vectors%rr(:,ref_cluster%ezoa(ii))
    write(*, '(a20,x,3(f13.6,x))') "Position in basis", atom_positions(:,ref_cluster%atom(ii))
    write(*, '(a20,x,3(f13.6,x))') "Position from center", ref_cluster%rcls(1:3,ii)
    if (dist2 > radius*radius) stop "ERROR in cluster generation."
  enddo ! ii

  call destroy(lattice_vectors)
  call destroy(ref_cluster)

endprogram ! test
#endif

