module BrillouinZoneMesh_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  
  public :: BrillouinZoneMesh, create, load, store, destroy

  type BrillouinZoneMesh
    integer :: nofks = 0 !< number of kpoints in this mesh
    double precision :: volBz = 0.d0 !< volume of the Brillouin zone
    double precision, allocatable :: kwxyz(:,:) ! (0:3,nofks), 0: weight (prev. VOLCUB)
  endtype

  interface create
    module procedure createBrillouinZoneMesh
  endinterface
  
  interface load
    module procedure loadBrillouinZoneMesh
  endinterface
  
  interface store
    module procedure storeBrillouinZoneMesh
  endinterface
  
  interface destroy
    module procedure destroyBrillouinZoneMesh
  endinterface

  contains
  
  subroutine createBrillouinZoneMesh(self, nofks, kwxyz, volBz)
    type(BrillouinZoneMesh), intent(inout) :: self
    integer, intent(in) :: nofks
    double precision, intent(in), optional :: kwxyz(0:,:) !< dim(0:3,nofks)
    double precision, intent(in), optional :: volBz !<

    integer :: ist
    
    if (nofks < 1) warn(6, "try to initialize BrillouinZoneMesh with less than 1 k-point, found"+nofks)
    self%nofks = max(1, nofks)
    deallocate(self%kwxyz, stat=ist) ! ignore status
    allocate(self%kwxyz(0:3,self%nofks), stat=ist)
    assert(ist == 0)
    self%kwxyz(:,:) = 0.d0; self%kwxyz(0,1) = 1.d0 ! set Gamma point
    self%volBz = 1.d0
    
    if (present(kwxyz)) then
      assert(ubound(kwxyz, 1) == 3)
      assert(size(kwxyz, 2) == size(self%kwxyz, 2))
      self%kwxyz = kwxyz ! copy
    endif
    
    if (present(volBz)) self%volBz = volBz
    if (self%volBz < 1e-6) warn(6, "volume of the Brillouin zone seems small, found"+self%volBz)
    
  endsubroutine ! create

  subroutine loadBrillouinZoneMesh(self, fu)
    type(BrillouinZoneMesh), intent(inout) :: self
    integer, intent(in) :: fu !< file unit

    integer :: ist, nk, ik, jk
    double precision, allocatable :: wxyz(:,:)
    character(len=96) :: line
    double precision :: vol, xyzw(1:4)
    
    line = ""; read(unit=fu, fmt='(A)', iostat=ist) line
    nk = 0; vol = 0.d0
    read(unit=line, fmt=*, iostat=ist) nk, vol
    if (ist /= 0) warn(6, "unable to read a Brillouin zone mesh from file unit"+fu-", expected int double but found <"-line-">")
    
    if (nk < 1) then
      call createBrillouinZoneMesh(self, nofks=0, volBz=0.d0)
      return
    endif
    
    allocate(wxyz(0:3,1:nk), stat=ist)

    jk = 0
    do ik = 1, nk
      line = ""; read(unit=fu, fmt='(A)', iostat=ist) line
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
    
    call createBrillouinZoneMesh(self, nofks=nk, kwxyz=wxyz, volBz=vol)
    
    deallocate(wxyz, stat=ist)
  endsubroutine ! load
  
  subroutine storeBrillouinZoneMesh(self, fu)
    type(BrillouinZoneMesh), intent(inout) :: self
    integer, intent(in) :: fu !< file unit
    integer :: ist, ik
    write(unit=fu, fmt='(i8,f15.10,/,(3f12.8,d20.10))', iostat=ist) &
      self%nofks, self%volBz, (self%kwxyz(1:3,ik), self%kwxyz(0,ik), ik=1,self%nofks)
  endsubroutine ! store
  
  elemental subroutine destroyBrillouinZoneMesh(self)
    type(BrillouinZoneMesh), intent(inout) :: self
    integer :: ist
    deallocate(self%kwxyz, stat=ist) ! ignore status
  endsubroutine ! destroy
  
endmodule ! BrillouinZoneMesh_mod