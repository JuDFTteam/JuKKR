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

public

type LatticeVectors
  integer nrd
  double precision, allocatable :: rr(:,:)   !dim 0:nrd
end type

type RefCluster
  !> reference to LatticeVectors datastructure
  type (LatticeVectors), pointer :: lattice_vectors
  double precision, allocatable :: rcls(:,:) !< positions relative to center
  integer, allocatable :: atom(:)  !< basis atom indices of cluster atoms
  integer, allocatable :: ezoa(:)  !< points into lattice_vectors%rr
  integer, allocatable :: indn0(:) !< indices of inequivalent cluster atoms
  integer :: numn0 !< number of inequivalent cluster atoms
  integer :: nacls !< number of cluster atoms

  integer :: atom_index !< basis atom index of central cluster atom
end type

!public createRefCluster

private vmul, veq, vadd, scalpr, xsort
private rrgen, rrgencount, clsgen99, clsgen99count

CONTAINS

!------------------------------------------------------------------------------
!> Creates a table of lattice vectors.
!>
!> Before any reference clusters can be created, a table of
!> lattice vectors has to be created.
!> Several different reference clusters can (and should) share
!> the same lattice vectors.
subroutine createLatticeVectors(lattice_vectors, bravais)
  implicit none
  type (LatticeVectors), intent(inout) :: lattice_vectors
  double precision, intent(in) :: bravais(3,3)

  integer nrd

  call rrgencount(bravais, nrd)
  !write(*,*) "Num. real space vectors: ", nrd
  allocate( lattice_vectors%rr(3,0:nrd) )
  lattice_vectors%nrd = nrd
  call rrgen (bravais, lattice_vectors%rr, lattice_vectors%nrd)

end subroutine

!------------------------------------------------------------------------------
!> Destroys table of lattice vectors.
!>
subroutine destroyLatticeVectors(lattice_vectors)
  implicit none
  type (LatticeVectors), intent(inout) :: lattice_vectors

  deallocate( lattice_vectors%rr )
  lattice_vectors%nrd = 0
end subroutine

!------------------------------------------------------------------------------
!> Creates reference cluster.
!>
!> Usage example:
!> type (LatticeVectors), target :: lattice_vectors
!> type (RefCluster) :: ref_cluster
!> ! ...
!> call createLatticeVectors(lattice_vectors, bravais)
!> call createRefCluster(ref_cluster, lattice_vectors, rbasis, rcut, center_ind)
!> ! ... some code
!> call destroyRefCluster(ref_cluster)
!> call destroyLatticeVectors(lattice_vectors)
!
subroutine createRefCluster(self, lattice_vectors, rbasis, rcut, center_ind)
  implicit none
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


  allocate ( self%atom(nacls), self%ezoa(nacls), self%indn0(nacls) )
  allocate ( self%rcls(3, nacls) )

  call clsgen99(center_ind, num_atoms , lattice_vectors%rr , rbasis, &
                self%atom,self%ezoa, &
                self%rcls,rcut, &
                self%numn0,self%indn0, &
                lattice_vectors%nrd, nacls)

  !write(*,*)  "Number of atoms in cluster:            ", nacls
  !write(*,*)  "Number of inequivalent cluster atoms: :", numn0

  self%nacls = nacls
end subroutine

!------------------------------------------------------------------------------
!> Destroys reference cluster.
subroutine destroyRefCluster(self)
  implicit none
  type(RefCluster), intent(inout) :: self

  deallocate ( self%atom, self%ezoa, self%indn0)
  deallocate ( self%rcls )
  self%nacls = 0
end subroutine

! ************************************************************************
! HELPER ROUTINES
! ************************************************************************

! ************************************************************************
Subroutine SCALPR(x,y,z)
! ************************************************************************
!     SCALSP COMPUTES THE scalar PRODUCT OF X AND Y RETURNING
!     IT INTO Z.
! ------------------------------------------------------------------------
DOUBLE PRECISION X(*), Y(*), Z
z= x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
return
END SUBROUTINE

! ************************************************************************
SUBROUTINE VADD(a,b,c)
! ************************************************************************
double precision a(*),b(*),c(*)
integer i
do 1 i=1,3
  c(i)=a(i)+b(i)
1 continue
return
END SUBROUTINE

! ************************************************************************
SUBROUTINE VEQ(a,b)
! ************************************************************************
double precision a(*),b(*)
integer i
do 1 i=1,3
  b(i)=a(i)
1 continue
return
END SUBROUTINE

! ************************************************************************
SUBROUTINE VMUL(a,b,c)
! ************************************************************************
double precision a(*),b,c(*)
integer i
do 1 i=1,3
  c(i)=b*a(i)
1 continue
return
END SUBROUTINE

! ************************************************************************
subroutine xsort (w,ind,max,pos)
  implicit none
  ! ************************************************************************
  !     p.zahn, april 96
  !     W   is the original array returned unchanged
  !     IND is an array that holds the new positions
  !     max number of ellements to be sorted
  !     pos the position where the first element is found
  ! ------------------------------------------------------------------------
  integer max,pos
  double precision  w(*), bound, diff
  integer ind(*)

  integer i,ii,j,jj,k
  data bound /1.0d-12/
  ! ------------------------------------------------------------------------
  do 10 i = 1,max
     ind(i) = i
10 end do

  j = max
  j = 1
  do 60 while (j.lt.max/3)
     j = 3*j+1
60 end do

  do 20 while (j.gt.1)
     j = j/3
     jj = 1
     do 30 while (jj.eq.1)
        jj = 0
        do 40 k=1,max-j
           diff = abs( w(ind(k)) - w(ind(k+j)) )
           if ( w(ind(k)) .gt. w(ind(k+j)) .and. &
                diff.gt.bound ) then
              ii       = ind(k)
              ind(k)   = ind(k+j)
              ind(k+j) = ii
              jj = 1
           end if
40      end do                    ! K=1,MAX-J
30   end do                      ! WHILE (JJ.EQ.1)
20 end do

  do 50 i=1,max
     if (ind(i) .eq. 1) pos=i
50 end do

  return
end subroutine xsort

! ************************************************************************
subroutine clsgen99count(center_ind, naez,rr,rbasis,rcut, &
                         nrd, nacls)
  implicit none
  ! ************************************************************************
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
           end do
           rxy2 =  tmp(1)**2+tmp(2)**2
           r2   =  tmp(3)**2 + tmp(1)**2+tmp(2)**2

           if ( (rxy2.le.rcutxy2).and.(r2.le.rcut2) )  then
              number = number + 1
              !atom(number) = na ! store the atom in elem cell
              !ezoa(number) = n ! store the bravais vector
           end if
        end do              ! N loop in bravais

     end do                 ! NA loop in NAEZ

  nacls = number ! return number of atoms in cluster

end subroutine clsgen99count


! ************************************************************************
subroutine clsgen99(center_ind, naez,rr,rbasis, &
     atom,ezoa, &
     rcls,rcut, &
     numn0,indn0, &
     nrd, nacls)
  implicit none
  ! ************************************************************************
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
  integer, parameter :: nclsd = 1
  integer, intent(in) :: nacls

  double precision rcut,rcutxy
  double precision &
       rbasis(3,*), &      ! pos. of basis atoms in EZ 
  rcls(3,nacls), &        ! real space position of atom in cluster
  rr(3,0:nrd)              ! set of lattice vectors
  !
  integer naez, &          ! number of atoms in EZ
  nr                       ! number of lattice vectors RR
  !
  integer numn0,indn0(nacls), &
  atom(nacls), &         ! index to atom in elem/cell at site in cluster
  ezoa(nacls)           ! index to bravais lattice  at site in cluster
  !
  !     .. locals
  !
  integer &
       i, &
       na,number,n, &
       jatom,iat,icluster
  integer iatcls(nclsd)
  logical icouplmat(naez)
  !
  double precision r2,epsshl, tmp(3)
  double precision rcut2,rcutxy2,rxy2 
  !
  external xsort,clustcomp
  intrinsic min,sqrt

  ! for sorting
  double precision rsort(nacls), rg(3,nacls)
  integer ia, ib, pos, iatom(nacls),iezoa(nacls),isort(nacls)
  !
  data     epsshl   / 1.0d-4 /

  nr = nrd
  indn0 = -1

  rcutxy = rcut ! only spherical clusters allowed

  icluster = 1

  do n = 1,nclsd
     iatcls(n) = 0
  end do

  rcutxy2 = (rcutxy+epsshl)*(rcutxy+epsshl)
  rcut2   = (rcut+epsshl)*(rcut+epsshl)

  !do jatom = 1,naez       ! loop in all atoms
  jatom = center_ind       ! modification to original: treat only one atom
     number = 0           ! counter for atoms in cluster
     do na = 1,naez  ! loop in all atoms
        do n=0,nr    ! loop in all bravais vectors    
           do i=1,3
              tmp(i) = rr(i,n)+rbasis(i,na)-rbasis(i,jatom)
           end do
           rxy2 =  tmp(1)**2+tmp(2)**2
           r2   =  tmp(3)**2 + tmp(1)**2+tmp(2)**2

           if ( (rxy2.le.rcutxy2).and.(r2.le.rcut2) )  then
              number = number + 1
              if (number.gt.nacls) then 
                 write (6,*) &
                      ' ERROR: Dimension NACLSD in inc.cls too small', &
                      number, nacls
                 stop '   < CLSGEN99 >'
              end if
              !
              atom(number) = na ! store the atom in elem cell
              ezoa(number) = n ! store the bravais vector
              do i=1,3
                 rcls(i,number) = tmp(i)
              end do
           end if
        end do              ! N loop in bravais

     end do                 ! NA loop in NAEZ

     !     sort the atoms of the cluster in increasing order. First by distance
     !     Then by z then by y
     !
     do ia=1,number
        rsort(ia) = sqrt(rcls(1,ia)**2+ &
             rcls(2,ia)**2+ &
             rcls(3,ia)**2)
        rsort(ia) = 100000000.d0*rsort(ia)+ &
             100000.d0*rcls(3,ia)+ &
             100.d0*rcls(2,ia)+ &
             0.1d0*rcls(1,ia)
     end do
     !
     call xsort(rsort,isort,number,pos)
     !     Rearange exchange ia with ib
     ! MAP temporarily to another array
     do ia=1,number
        do i=1,3
           rg(i,ia)    = rcls(i,ia)
        end do
        iatom(ia) = atom(ia)
        iezoa(ia) = ezoa(ia)
     end do
     ! Now use correct order
     do ia =1,number
        ib = isort(ia)
        do i=1,3
           rcls(i,ia) = rg(i,ib)
        end do
        atom(ia) = iatom(ib)
        ezoa(ia) = iezoa(ib)
     end do

 ! do jatom = 1,naez

     do iat = 1,naez
        icouplmat(iat) = .false.
        do i=1,number 
           if (atom(i).eq.iat) then
              icouplmat(iat) = .true.
           end if
        end do
     end do

     numn0 = 0
     do iat = 1,naez
        if(icouplmat(iat)) then
           numn0 = numn0 + 1
           indn0(numn0) = iat
           !     WRITE(6,*) 'JATOM,NUMN0,INDN0',JATOM,NUMN0(JATOM),IAT
        endif
     enddo

  return
end subroutine clsgen99


subroutine rrgencount (bv1,nrd)
  ! **********************************************************************
  ! *                                                                    *
  ! * COUNTS      number of real space vectors to construct the          *
  ! * clusters representing the local surrounding of the atoms in        *
  ! * routine CLSGEN99                                                   *
  ! *                                                                    *
  ! **********************************************************************
  implicit none
  !     ..
  integer nr,nrd
  !    ..
  double precision bv1(3,3)
  !    ..
  double precision epsshl,r,r1,r2,r3,rmax,rr2,rs
  integer i,j,k,n1,n2,n3,iprint
  integer nint
  double precision dble
  !     ..
  !     .. Local arrays
  double precision &
       v(3),vx(3),vy(3),vz(3), &
       vx0(3),vy0(3),vz0(3)

  !     .. Data Statements ..
  data  epsshl /1.0d-5/
  !     ..................................................................
  iprint = 0
  !
  call scalpr(bv1(1,1),bv1(1,1),r1)
  call scalpr(bv1(1,2),bv1(1,2),r2)
  call scalpr(bv1(1,3),bv1(1,3),r3)
  rmax = 5.d0
  !
  r1 = sqrt(r1)
  r2 = sqrt(r2)
  r3 = sqrt(r3)
  r = 1.5d0*rmax + sqrt(r1*r1+r2*r2+r3*r3) + epsshl
  rs = r*r
  n1 = nint(r/r1)
  n2 = nint(r/r2)
  n3 = nint(r/r3)
  !
  n1 = min(12,n1)
  n2 = min(12,n2)
  n3 = min(12,n3)
  !
  n1 = max(2,n1)
  n2 = max(2,n2)
  n3 = max(2,n3)
  !
  !write (6,99001) r
  !write (6,99002) rs
  !write (6,99004) n1,n2,n3
  !
  nr = 0
  !
  call vmul(bv1(1,1),dble(-n1-1),vx0(1))
  call vmul(bv1(1,2),dble(-n2-1),vy0(1))
  call vmul(bv1(1,3),dble(-n3-1),vz0(1))
  call veq(vx0,vx)
  ! **********************************************************************
  do i = -n1,n1
     call vadd(vx,bv1(1,1),vx)
     call veq(vy0,vy)
     ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     do j = -n2,n2
        call vadd(vy,bv1(1,2),vy)
        call veq(vz0,vz)
        ! ----------------------------------------------------------------------
        do k = -n3,n3
           call vadd(vz,bv1(1,3),vz)
           call vadd(vx,vy,v)
           call vadd(v,vz,v)
           call scalpr(v,v,rr2)
           !
           if ( ((rr2.le.rs) .or. (abs(i)+abs(j)+abs(k).le.6)) &
                .and. (rr2.gt.epsshl) ) then

              nr = nr + 1

           end if
        end do
        ! ----------------------------------------------------------------------
     end do
     ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do

  nrd = nr

99001 format (10x,'Radius R        : ',f15.6,' (ALAT    units)')
99002 format (10x,'       R**2     : ',f15.6,' (ALAT**2 units)')
99003 format (10x,'mesh divisions  : ',5x,2i5)
99004 format (10x,'mesh divisions  : ',3i5)
99005 format (10x,'vectors created : ',i15)
99006 format (/,10x,60('+'),/,18x, &
       'generated real-space mesh-points (ALAT units)',/, &
       10x,60('+'),/,13x, &
       'index      x           y           z          distance  ' &
       ,/,10x,60('-'))
99007 format (10x,60('+'))
99008 format (10x,i6,3f12.3,f15.4)
end subroutine rrgencount


!*==rrgen.f    processed by SPAG 6.05Rc at 20:37 on 17 May 2004
! 02.08.95 *************************************************************
subroutine rrgen (bv1,rr,nrd)
  ! **********************************************************************
  ! *                                                                    *
  ! * generates a number of real space vectors to construct the          *
  ! * clusters representing the local surrounding of the atoms in        *
  ! * routine CLSGEN99                                                   *
  ! *                                                                    *
  ! **********************************************************************
  implicit none
  !     ..
  !     .. Scalar arguments ..
  integer nr,nrd
  !    ..
  !    .. Array arguments ..
  double precision bv1(3,3),rr(3,0:nrd)
  !    ..
  !    .. Local scalars ..
  double precision epsshl,r,r1,r2,r3,rmax,rr2,rs
  integer i,j,k,n1,n2,n3,pos,iprint
  integer nint
  double precision dble
  !     ..
  !     .. Local arrays
  double precision rabs(nrd),rr1(3,nrd), &
       v(3),vx(3),vy(3),vz(3), &
       vx0(3),vy0(3),vz0(3)
  integer ind(nrd)
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic abs,min,sqrt,nint
  !     ..

  !     .. Data Statements ..
  data  epsshl /1.0d-5/
  !     ..................................................................
  !write (6,'(5X,A,/)') &
  !     '< RRGEN > : generation of real space mesh RR(NR)'
  !
  iprint = 0
  !
  call scalpr(bv1(1,1),bv1(1,1),r1)
  call scalpr(bv1(1,2),bv1(1,2),r2)
  call scalpr(bv1(1,3),bv1(1,3),r3)
  rmax = 5.d0
  !
  r1 = sqrt(r1)
  r2 = sqrt(r2)
  r3 = sqrt(r3)
  r = 1.5d0*rmax + sqrt(r1*r1+r2*r2+r3*r3) + epsshl
  rs = r*r
  n1 = nint(r/r1)
  n2 = nint(r/r2)
  n3 = nint(r/r3)
  !
  n1 = min(12,n1)
  n2 = min(12,n2)
  n3 = min(12,n3)
  !
  n1 = max(2,n1)
  n2 = max(2,n2)
  n3 = max(2,n3)
  !
  if (iprint == 1) then
      write (6,99001) r
      write (6,99002) rs
      write (6,99004) n1,n2,n3
  end if
  !
  nr = 0
  rr(1,0) = 0.0d0
  rr(2,0) = 0.0d0
  rr(3,0) = 0.0d0
  !
  call vmul(bv1(1,1),dble(-n1-1),vx0(1))
  call vmul(bv1(1,2),dble(-n2-1),vy0(1))
  call vmul(bv1(1,3),dble(-n3-1),vz0(1))
  call veq(vx0,vx)
  ! **********************************************************************
  do i = -n1,n1
     call vadd(vx,bv1(1,1),vx)
     call veq(vy0,vy)
     ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     do j = -n2,n2
        call vadd(vy,bv1(1,2),vy)
        call veq(vz0,vz)
        ! ----------------------------------------------------------------------
        do k = -n3,n3
           call vadd(vz,bv1(1,3),vz)
           call vadd(vx,vy,v)
           call vadd(v,vz,v)
           call scalpr(v,v,rr2)
           !
           if ( ((rr2.le.rs) .or. (abs(i)+abs(j)+abs(k).le.6)) &
                .and. (rr2.gt.epsshl) ) then
              nr = nr + 1
              !
              if ( nr.gt.nrd ) then
                 write (6,*) 'Dimension ERROR. Please, change the ', &
                      'parameter NRD in inc.p to ',nr
                 stop
              end if
              !
              rr1(1,nr) = v(1)
              rr1(2,nr) = v(2)
              rr1(3,nr) = v(3)
              rabs(nr) = sqrt(rr2)
           end if
        end do
        ! ----------------------------------------------------------------------
     end do
     ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do
  ! **********************************************************************

  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  if ( iprint.gt.0 ) then
     write (6,99006)
     write (6,99008) 0,0.0,0.0,0.0,0.0
  end if
  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

  call xsort(rabs,ind,nr,pos)

  do i = 1,nr
     pos = ind(i)
     rr(1,i) = rr1(1,pos)
     rr(2,i) = rr1(2,pos)
     rr(3,i) = rr1(3,pos)
     ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
     if ( iprint.gt.0 ) write (6,99008) i,rr(1,i),rr(2,i), &
          rr(3,i),rabs(pos)
     ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  end do

  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  if ( iprint.gt.0 )  write (6,99007)
  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

99001 format (10x,'Radius R        : ',f15.6,' (ALAT    units)')
99002 format (10x,'       R**2     : ',f15.6,' (ALAT**2 units)')
99003 format (10x,'mesh divisions  : ',5x,2i5)
99004 format (10x,'mesh divisions  : ',3i5)
99005 format (10x,'vectors created : ',i15)
99006 format (/,10x,60('+'),/,18x, &
       'generated real-space mesh-points (ALAT units)',/, &
       10x,60('+'),/,13x, &
       'index      x           y           z          distance  ' &
       ,/,10x,60('-'))
99007 format (10x,60('+'))
99008 format (10x,i6,3f12.3,f15.4)
end subroutine rrgen

end module

#ifdef TEST_CLUSTERS__
program test_clusters
  use ClusterHelpers_mod
  implicit none

  double precision bravais(3,3)
  double precision rbasis(3,4)

  type (RefCluster) :: ref_cluster
  type (LatticeVectors), target :: lattice_vectors

  integer ii
  double precision :: rcut = 1.5
  integer :: center_ind = 4
  double precision :: dist

  bravais = reshape( (/ 1.0, 0.0, 0.0, &
                        0.0, 1.0, 0.0, &
                        0.0, 0.0, 1.0 /), (/ 3, 3 /) )

  rbasis = reshape( (/ 0.0, 0.0, 0.0, &
                       0.5, 0.5, 0.0, &
                       0.5, 0.0, 0.5, &
                       0.0, 0.5, 0.5  /), (/ 3, 4 /) )

  call createLatticeVectors(lattice_vectors, bravais)
  call createRefCluster(ref_cluster, lattice_vectors, rbasis, rcut, center_ind)

  write(*,*)  "Number of lattice vectors:             ", ref_cluster%lattice_vectors%nrd
  write(*,*)  "Number of atoms in cluster:            ", ref_cluster%nacls
  write(*,*)  "Number of inequivalent cluster atoms: :", ref_cluster%numn0

  do ii = 1, ref_cluster%nacls
    dist = sqrt(ref_cluster%rcls(1,ii)**2 + ref_cluster%rcls(2,ii)**2 + ref_cluster%rcls(3,ii)**2 )
    write(*,"(79('-'))")
    write(*, 92000) "Cluster atom", ii
    write(*, 91000) "Basis atom / dist. ", ref_cluster%atom(ii), dist
    write(*, 90000) "Lattice vector", ref_cluster%lattice_vectors%rr(:, ref_cluster%ezoa(ii))
    write(*, 90000) "Position in basis", rbasis(:, ref_cluster%atom(ii))
    write(*, 90000) "Position from center", ref_cluster%rcls(:, ii)
    if (dist > rcut) then
      write(*,*) "ERROR in cluster generation."
      STOP
    end if
  end do

  call destroyLatticeVectors(lattice_vectors)
  call destroyRefCluster(ref_cluster)

90000 format (A20, X, 3(F13.6, X))
91000 format (A20, X, I7, X, F13.6)
92000 format (A20, X, I7)
end program
#endif TEST_CLUSTERS__

