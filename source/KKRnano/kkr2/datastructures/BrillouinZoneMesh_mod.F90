module BrillouinZoneMesh_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  
  public :: BrillouinZoneMesh, create, load, store, destroy

  type BrillouinZoneMesh
    integer :: nofks = 0 !< number of kpoints in this mesh
    integer :: nks(3) = 0 !< numbers of kpoint in the mesh before symmetry reduction
    double precision :: volBz = 0.d0 !< volume of the Brillouin zone
    double precision, allocatable :: kwxyz(:,:) ! (0:3,nofks), 0: weight (prev. called VOLCUB)
  endtype

  interface create
    module procedure createBrillouinZoneMesh_n, createBrillouinZoneMesh_k
  endinterface
  
  interface load
    module procedure loadBrillouinZoneMesh, loadBrillouinZoneMeshes
  endinterface
  
  interface store
    module procedure storeBrillouinZoneMesh, storeBrillouinZoneMeshes
  endinterface
  
  interface destroy
    module procedure destroyBrillouinZoneMesh
  endinterface

  contains

  subroutine createBrillouinZoneMesh_n(self, nofks, volBz, nks)
    type(BrillouinZoneMesh), intent(inout) :: self
    integer, intent(in) :: nofks
    double precision, intent(in), optional :: volBz !<
    integer, intent(in), optional :: nks(3)

    integer :: ist
    
    if (nofks < 1) warn(6, "try to initialize BrillouinZoneMesh with less than 1 k-point, found"+nofks)
    self%nofks = max(1, nofks)
    deallocate(self%kwxyz, stat=ist) ! ignore status
    allocate(self%kwxyz(0:3,self%nofks), stat=ist)
    assert(ist == 0)
    self%kwxyz(:,:) = 0.d0; self%kwxyz(0,1) = 1.d0 ! set Gamma point
    self%volBz = 1.d0
    
    if (present(volBz)) self%volBz = volBz
    if (self%volBz < 1e-6) warn(6, "volume of the Brillouin zone seems small, found"+self%volBz)
    if (present(nks)) self%nks = nks
    
  endsubroutine ! create

  
  subroutine createBrillouinZoneMesh_k(self, kwxyz, volBz, nks)
    type(BrillouinZoneMesh), intent(inout) :: self
    double precision, intent(in) :: kwxyz(0:,:) !< dim(0:3,nofks)
    double precision, intent(in), optional :: volBz !<
    integer, intent(in), optional :: nks(3)

    call create(self, size(kwxyz, 2), volBz)
    
    assert(ubound(kwxyz, 1) == 3)
    assert(size(kwxyz, 2) == size(self%kwxyz, 2))
    self%kwxyz = kwxyz ! copy
    if (present(nks)) self%nks = nks
    
  endsubroutine ! create
  
  subroutine loadBrillouinZoneMesh(self, fu, comm, rank)
    type(BrillouinZoneMesh), intent(inout) :: self
    integer, intent(in) :: fu !< file unit
    integer, intent(in) :: comm, rank

    include 'mpif.h'
    integer :: ist, nk, ik, jk
    double precision, allocatable :: wxyz(:,:)
    character(len=96) :: line
    double precision :: vol, xyzw(1:4)

    if (rank == 0) then

      line = ""; read(unit=fu, fmt='(a)', iostat=ist) line
      nk = 0; vol = 0.d0
      read(unit=line, fmt=*, iostat=ist) nk, vol
      if (ist /= 0) warn(6, "unable to read a Brillouin zone mesh from file unit"+fu-", expected int double but found <"-line-">")

    endif ! root     

    call MPI_Bcast(nk, 1, MPI_INTEGER, 0, comm, ist)

    if (nk < 1) then
      call create(self, nofks=0, volBz=0.d0)
      warn(6, "create an empty k-point mesh")
      return
    endif

    allocate(wxyz(0:3,1:nk), stat=ist)

    if (rank == 0) then

      line = ""; read(unit=fu, fmt='(a)', iostat=ist) line ! read comment line

      jk = 0
      do ik = 1, nk
        line = ""; read(unit=fu, fmt='(a)', iostat=ist) line
        if (ist == 0) then
          read(unit=line, fmt=*, iostat=ist) xyzw
          if (ist == 0) then
            jk = jk + 1
            wxyz(0,jk) = xyzw(4); wxyz(1:3,jk) = xyzw(1:3)
          else
            warn(6, "unable to parse the"+ik-"th k-point but found <"-line-">")
          endif
        else
          warn(6, "unable to read line expecting the"+ik-"th k-point")
        endif
      enddo ! ik
      if (jk /= nk) die_here("unable to read find"+nk+"valid k-points in file unit"+fu-", found only"+jk)

    endif ! root
    
    call MPI_Bcast(wxyz, 4*nk, MPI_DOUBLE_PRECISION, 0, comm, ist)
    call MPI_Bcast(vol, 1, MPI_DOUBLE_PRECISION, 0, comm, ist)
    
    call create(self, kwxyz=wxyz, volBz=vol)

    deallocate(wxyz, stat=ist)
  endsubroutine ! load

  
  subroutine loadBrillouinZoneMeshes(self, filename, comm)
    type(BrillouinZoneMesh), intent(inout) :: self(1:)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: comm

    include 'mpif.h'
    integer, parameter :: fu = 654 !< file unit
    integer :: i, rank, n, ios
    character(len=16) :: comment
    call MPI_Comm_rank(comm, rank, ios)
    if (rank == 0) open(unit=fu, file=filename, action='read', status='old', iostat=ios)
    call MPI_Bcast(ios, 1, MPI_INTEGER, 0, comm, ios)
    if (ios /= 0) return
    n = 0
    if (rank == 0) read(unit=fu, iostat=ios) n ! 1st line of the file must contain the number of meshes
    call MPI_Bcast(n, 1, MPI_INTEGER, 0, comm, ios)
    if (rank == 0) read(unit=fu, iostat=ios) comment
    do i = 1, min(size(self, 1), n)
      call load(self(i), fu, comm, rank)
    enddo ! i
    if (rank == 0) close(unit=fu, iostat=ios)
  endsubroutine ! load
  
  
  subroutine storeBrillouinZoneMesh(self, fu)
    type(BrillouinZoneMesh), intent(in) :: self
    integer, intent(in) :: fu !< file unit
    integer :: ist, ik
    if (self%nofks < 1 .or. size(self%kwxyz, 2) < 1) return
    write(unit=fu, fmt='(i8,f15.10)', iostat=ist) self%nofks, self%volBz
    write(unit=fu, fmt='(9(a,i0))', iostat=ist) '# ',self%nks(1),'x',self%nks(2),'x',self%nks(3),' k-mesh, ', &
                                                     self%nofks,' irreducible points' ! comment line
    write(unit=fu, fmt='(3f12.8,d20.10)', iostat=ist) (self%kwxyz(1:3,ik), self%kwxyz(0,ik), ik=1,self%nofks)
  endsubroutine ! store

  
  integer function storeBrillouinZoneMeshes(self, filename) result(ios)
    type(BrillouinZoneMesh), intent(in) :: self(1:)
    character(len=*), intent(in) :: filename
    integer, parameter :: fu = 654 !< file unit
    integer :: i
    open(unit=fu, file=filename, action='write', iostat=ios)
    if (ios /= 0) return
    write(unit=fu, fmt='(9(a,i0))', iostat=ios) ' ',size(self, 1),'   =NumberOfMeshes'
    do i = 1, size(self, 1)
      call store(self(i), fu)
    enddo ! i
    close(unit=fu, iostat=ios)
  endfunction ! store
  
  elemental subroutine destroyBrillouinZoneMesh(self)
    type(BrillouinZoneMesh), intent(inout) :: self
    integer :: ist
    deallocate(self%kwxyz, stat=ist) ! ignore status
  endsubroutine ! destroy
  
endmodule ! BrillouinZoneMesh_mod
