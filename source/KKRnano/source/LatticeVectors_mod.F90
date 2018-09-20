module LatticeVectors_mod
!-------------------------------------------------------------------------------
!> Summary: Linear combinations of Bravais vectors with integer coefficients
!> Author: Paul F Baumeister, Elias Rabel
!> Category: KKRnano, geometry, memory-management, initialization
!>
!> Note: if N (=number of atoms) reference clusters are created, then this algorithm scales as O(N**2)
!-------------------------------------------------------------------------------
  implicit none
  private
  
  public :: LatticeVectors, create, destroy

  type LatticeVectors
    integer :: nrd
    double precision, allocatable :: rr(:,:) !< dim(3,0:nrd)
  endtype

  interface create
    module procedure createLatticeVectors
  endinterface
  
  interface destroy
    module procedure destroyLatticeVectors
  endinterface
  
  contains

  !------------------------------------------------------------------------------
  !> Creates a table of lattice vectors.
  !>
  !> Before any reference clusters can be created, a table of
  !> lattice vectors has to be created.
  !> Several different reference clusters can (and should) share
  !> the same lattice vectors.
  subroutine createLatticeVectors(self, bravais, rmax)
    type(LatticeVectors), intent(inout) :: self
    double precision, intent(in) :: bravais(3,3)
    double precision, intent(in), optional :: rmax
    
    double precision :: rmx
    rmx = 5.d0; if (present(rmax)) rmx = rmax
     
    self%nrd = rrgen(bravais, rmx, self%rr)
#ifdef DEBUG  
    write(*,*) "Num. real space vectors: ", self%nrd
#endif
  endsubroutine ! create

  !------------------------------------------------------------------------------
  !> Destroys table of lattice vectors.
  elemental subroutine destroyLatticeVectors(self)
    type(LatticeVectors), intent(inout) :: self

    integer :: ist
    deallocate(self%rr, stat=ist)
    self%nrd = 0
  endsubroutine ! destroy

  !------------------------------------------------------------------------------
  integer function rrgen(bv, rmax, rr) result(nr)
    !> generates a number of real space vectors to construct the clusters representing the local surrounding of the atoms in routine clsgen99
    use Sorting_mod, only: dsort
    
    double precision, intent(in) :: bv(3,3) ! bravais vectors or matrix
    double precision, intent(in) :: rmax ! radius
    double precision, allocatable, intent(out) :: rr(:,:) ! (3,0:nr)
    
    !    .. Locals
    double precision :: r, rr2, rs
    integer :: i, i1, i2, i3, nn(1:3), mr
    double precision :: v(3)
    double precision, allocatable :: rabs2(:), rr1(:,:) 
    integer, allocatable :: ind(:) ! permutation for sorting
    double precision, parameter :: epsshl = 1.d-5
#ifdef DEBUG
    character(len=*), parameter :: F8='(10x,i6,3f24.12,f15.4)' ! DEBUG
    integer, parameter :: iprint = 1
#else
    character(len=*), parameter :: F8='(10x,i6,3f12.3,f15.4)'
    integer, parameter :: iprint = 0
#endif
    ! write (6,'(5X,A,/)') '< RRGEN > : generation of real space mesh RR(NR)'

    do i = 1, 3 
      v(i) = sum(bv(1:3,i)**2)
    enddo ! i
    v(1:3) = sqrt(v(1:3))
    
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

    nr = 0
    do i1 = -nn(1), nn(1)
      do i2 = -nn(2), nn(2)
        do i3 = -nn(3), nn(3)
        
          v(1:3) = i1*bv(1:3,1) + i2*bv(1:3,2) + i3*bv(1:3,3)
          rr2 = sum(v(1:3)*v(1:3))
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
    
    rr = 0.d0 ! init all
    rr(1:3,0) = 0.d0 ! init 1st entry

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
    
  endfunction ! rrgen

endmodule ! LatticeVectors_mod
