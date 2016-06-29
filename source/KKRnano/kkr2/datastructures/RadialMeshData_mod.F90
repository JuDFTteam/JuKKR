#ifdef TASKLOCAL_FILES
#define FILEWRITE(X,Y) write(X)
#define FILEREAD(X,Y) read(X)
#else
#define FILEWRITE(X,Y) write(X,Y)
#define FILEREAD(X,Y) read(X,Y)
#endif

!> This module defines a datatype that contains data related to the radial mesh
!> @author Elias Rabel
module RadialMeshData_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: RadialMeshData, create, destroy, represent, load
  public :: getMinReclenMesh
  public :: initRadialMesh, initMuffinTinMesh, initInterstitialMesh, writeRadialMeshDataDA           
  public :: readRadialMeshDataDA, openRadialMeshDataDAFile, closeRadialMeshDataDAFile, writeRadialMeshDataIndexDA      
  public :: readRadialMeshDataIndexDA, openRadialMeshDataIndexDAFile, closeRadialMeshDataIndexDAFile, repr_RadialMeshData             
  public :: readRadialMeshDataHeader, test_mesh_consistency           

  type RadialMeshData

    ! dimensioning parameters
    integer :: irmd   !< number of mesh points
    integer :: ipand  !< dimension variable panels for historical reasons

    ! scalars
    double precision :: A    !< logarithmic mesh parameter A
    double precision :: B    !< logarithmic mesh parameter A
    double precision :: RWS  !< maximal radius
    double precision :: RMT  !< muffin-tin radius
    integer :: IPAN   !< number of mesh panels
    integer :: IRC    !< ?
    integer :: IMT    !< endof muffin-tin region
    integer :: IRNS
    integer :: IRWS   !< index of max. radius
    integer :: IRMIN
#ifdef USE_OLD_MESH
    integer :: meshn    !< number of interstitial mesh points 
    integer :: nfu    !< number of non-zero shape function components 
#endif
    ! arrays
    double precision, dimension(:), allocatable :: R
    double precision, dimension(:), allocatable :: DRDI
    integer, dimension(:), allocatable :: IRCUT  !< panel locations
#ifdef USE_OLD_MESH
    integer, dimension(:), allocatable :: llmsp  !<  LM index of non-zero shape function component
    double precision, dimension(:,:), allocatable :: THETAS !< radial values of shape function component \Theta_L MESHN values
#endif
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  ! meshn and nfu are optional!
  subroutine createRadialMeshData(meshdata, irmd, ipand, meshn, nfu)
    implicit none
    type (RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: irmd
    integer, intent(in) :: ipand
    integer, intent(in), optional :: meshn
    integer, intent(in), optional :: nfu

    meshdata%irmd = irmd
    meshdata%ipand = ipand

    allocate(meshdata%R(irmd))
    allocate(meshdata%DRDI(irmd))
    allocate(meshdata%IRCUT(0:ipand))
#ifdef USE_OLD_MESH
    if(present(meshn) .AND. present(nfu)) then
      allocate(meshdata%LLMSP(nfu))
      allocate(meshdata%THETAS(meshn,nfu))
    endif
#endif
    meshdata%R = 0.0d0
    meshdata%DRDI = 0.0d0
    meshdata%IRCUT = 0
#ifdef USE_OLD_MESH
    if(present(meshn) .AND. present(nfu)) then
      meshdata%LLMSP = 0.0d0
      meshdata%THETAS = 0.0d0
    endif
#endif
    meshdata%A = 0.0d0
    meshdata%B = 0.0d0

    meshdata%IPAN = 0
    meshdata%IRC = 0
    meshdata%IMT = 0
    meshdata%IRNS = 0
    meshdata%IRWS = 0
    meshdata%IRMIN = 0
#ifdef USE_OLD_MESH
    if(present(meshn) .AND. present(nfu)) then
      meshdata%MESHN = 0
      meshdata%NFU = 0
    endif
#endif
    meshdata%RWS = 0.0d0
    meshdata%RMT = 0.0d0

    self%RWS = 0.d0
    self%RMT = 0.d0

  endsubroutine ! create

  !----------------------------------------------------------------------------
  subroutine destroyRadialMeshData(self)
    type(RadialMeshData), intent(inout) :: self

    deallocate(meshdata%R)
    deallocate(meshdata%DRDI)
    deallocate(meshdata%IRCUT)
#ifdef USE_OLD_MESH
    deallocate(meshdata%LLMSP)
    deallocate(meshdata%THETAS)
#endif
  end subroutine

  !----------------------------------------------------------------------------
  ! meshn, nfu, llmsp and thetas are optional!
  subroutine initRadialMesh(meshdata, alat, xrn, drn, nm, imt, irns, meshn, nfu, llmsp, thetas)
    implicit none
    type (RadialMeshData), intent(inout) :: meshdata
    double precision, intent(in) :: alat
    double precision, intent(in) :: xrn(:), drn(:)
    integer, intent(in) :: nm(:)
    integer, intent(in) :: imt
    integer, intent(in) :: irns
    integer, intent(in), optional :: meshn
    integer, intent(in), optional :: nfu
    integer, intent(in), optional :: llmsp (:)
    double precision, intent(in), optional :: thetas(:,:)

    call initInterstitialMesh(meshdata, alat, xrn, drn, nm, imt, irns, meshn, nfu, llmsp, thetas)
    ! note radius_mt = xrn(1) * alat
    call initMuffinTinMesh(self, imt, xrn(1)*alat)

    call test_mesh_consistency(self)

  endsubroutine ! init

  !---------------------------------------------------------------------------
  ! note radius_mt = xrn(1) * alat - in units of Bohr
  subroutine initMuffinTinMesh(self, imt, radius_mt)
    type(RadialMeshData), intent(inout) :: self
    integer, intent(in) :: imt
    double precision, intent(in) :: radius_mt

    double precision, parameter :: A = 0.025d0
    integer :: ii

    self%A = A
    self%B = radius_mt / (exp(A * (imt - 1)) - 1.d0)

    do ii = 1, imt
      self%r(ii)    = self%B * (exp(A * (ii - 1)) - 1.d0)
      self%drdi(ii) = A * self%B * exp(A * (ii - 1))
    enddo ! ii

    self%rmt = radius_MT

  endsubroutine ! init

  !---------------------------------------------------------------------------
  ! note radius_mt = xrn(1) * alat
  ! meshn, nfu, llmsp and thetas are optional!
  subroutine initInterstitialMesh(meshdata, alat, xrn, drn, nm, imt, irns, meshn, nfu, llmsp, thetas)
    implicit none
    type (RadialMeshData), intent(inout) :: meshdata
    double precision, intent(in) :: alat
    double precision, intent(in) :: xrn(:), drn(:)
    integer, intent(in) :: nm(:)
    integer, intent(in) :: imt
    integer, intent(in) :: irns
    integer, intent(in), optional :: meshn
    integer, intent(in), optional :: nfu
    integer, intent(in), optional :: llmsp (:)
    double precision, intent(in), optional :: thetas(:,:)

    integer :: isum, ii, ipan

    ipan = size(nm) + 1 ! +1 for muffin-tin region 1..imt
    meshdata%ipan = ipan
    meshdata%imt = imt
#ifdef USE_OLD_MESH
    if(present(meshn) .AND. present(nfu)) then
      meshdata%meshn = meshn
      meshdata%nfu = nfu
    endif
#endif
    ! ircut(0) has to be 0, integrations start at ircut(i)+1
    self%ircut(0) = 0
    self%ircut(1) = imt

    isum = imt
    do ii = 2, ipan
      isum = isum + nm(ii - 1)
      self%ircut(ii) = isum
    enddo ! ii
    self%irc = self%ircut(ipan)

    self%irws = isum
    if (isum > self%irmd) &
      die_here("Error creating mesh, irmd too small, isum ="+isum-", self%irmd ="+self%irmd)

    meshdata%irmin = meshdata%irws - irns
    meshdata%irns = irns
    meshdata%rws  = meshdata%r(meshdata%irws)
#ifdef USE_OLD_MESH
    if(present(llmsp) .AND. present(thetas)) then
      meshdata%llmsp = llmsp
      meshdata%thetas = thetas
    endif
#endif
    end subroutine

!==============================================================================
!=                                    I/O                                     =
!==============================================================================

  !----------------------------------------------------------------------------
  !> Read and create radial mesh data from files <filename> and <filename>.idx
  !>
  !> The index file with extension .idx is necessary because different
  !> atoms can have a different number of radial mesh points
  subroutine createRadialMeshDataFromFile(self, filename, recnr)
    type(RadialMeshData), intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: recnr

    integer :: FILEUNIT = 37
    integer :: irmd, ipand, max_reclen, meshn, nfu

    ! index file has extension .idx

    call openRadialMeshDataIndexDAFile(meshdata, FILEUNIT, &
                                       filename // ".idx") !ignored for task-local files
    call readRadialMeshDataIndexDA(meshdata, FILEUNIT, recnr, &
                                       irmd, ipand, max_reclen &
#ifdef USE_OLD_MESH    
                                       ,meshn, nfu &
#endif
                                       )
    call closeRadialMeshDataIndexDAFile(FILEUNIT)

    call createRadialMeshData(meshdata, irmd, ipand &
#ifdef USE_OLD_MESH    
                              , meshn, nfu &
#endif
                              )

#ifdef USE_OLD_MESH
    meshdata%meshn = meshn
    meshdata%nfu   = nfu
#endif

    call openRadialMeshDataDAFile(self, fu, filename, max_reclen)
    call readRadialMeshDataDA(self, fu, recnr)
    call closeRadialMeshDataDAFile(fu)

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Write mesh data to direct access file 'fileunit' at record 'recnr'
  subroutine writeRadialMeshDataDA(self, fileunit, recnr)
    type(RadialMeshData), intent(in) :: self
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

#ifdef TASKLOCAL_FILES
    character(len=16) :: num
    write(unit=num, fmt='(a,i7.7)') "mesh.",recnr
    open(fileunit, file=num, form='unformatted', action='write')

    call writeRadialMeshDataIndexDA(self, fileunit, recnr, 0)
#endif

    FILEWRITE (fileunit, rec=recnr) MAGIC_NUMBER + recnr, &
                                meshdata%A, &
                                meshdata%B, &
                                meshdata%RWS, &
                                meshdata%RMT, &
                                meshdata%IPAN, &
                                meshdata%IRC, &
                                meshdata%IMT, &
                                meshdata%IRNS, &
                                meshdata%IRWS, &
                                meshdata%IRMIN, &
                                meshdata%R, &
                                meshdata%DRDI, &
                                meshdata%IRCUT, &
#ifdef USE_OLD_MESH
                                meshdata%LLMSP, &
                                meshdata%THETAS, &
#endif
                                MAGIC_NUMBER + recnr

#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Read mesh data from direct access file 'fileunit' at record 'recnr'
  subroutine readRadialMeshDataDA(self, fileunit, recnr)
    type(RadialMeshData), intent(inout) :: self
    integer, intent(in) :: fileunit, recnr

    integer :: magic, magic2, checkmagic

#ifdef TASKLOCAL_FILES
    character(len=7) :: num
    integer :: ipand, irmd, max_reclen, meshn, nfu    

    write(unit=num, fmt='(a,i7.7)') "mesh.",recnr
    open(fileunit, file=num, form='unformatted', action='read', status='old')

    ! read header at beginning of file
    call readRadialMeshDataHeader(meshdata, FILEUNIT, recnr, irmd, ipand, max_reclen &
#ifdef USE_OLD_MESH
                                  , meshn, nfu &
#endif
                                  )
#endif

    checkmagic = MAGIC_NUMBER + recnr

    FILEREAD (fileunit, rec=recnr) magic, &
                                meshdata%A, &
                                meshdata%B, &
                                meshdata%RWS, &
                                meshdata%RMT, &
                                meshdata%IPAN, &
                                meshdata%IRC, &
                                meshdata%IMT, &
                                meshdata%IRNS, &
                                meshdata%IRWS, &
                                meshdata%IRMIN, &
                                meshdata%R, &
                                meshdata%DRDI, &
                                meshdata%IRCUT, &
#ifdef USE_OLD_MESH
                                meshdata%LLMSP, &
                                meshdata%THETAS, &
#endif
                                magic2
#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif

    if (magic /= checkmagic .or. magic2 /= checkmagic) die_here("Invalid mesh data read.")

    call test_mesh_consistency(self)

  endsubroutine ! read

  !---------------------------------------------------------------------------
  !> Return MINIMUM record length needed to store mesh of this atom
  !>
  !> Note: for a file containing meshes of several atoms, the maximum
  !> of all their record lengths has to be determined
  endfunction ! get
  integer function getMinReclenMesh(meshdata) result(reclen)
    implicit none

    type (RadialMeshData), intent(in) :: meshdata
    !------

    integer, parameter :: MAGIC_NUMBER = -889271554

    inquire (iolength = reclen) MAGIC_NUMBER, &
                                meshdata%A, &
                                meshdata%B, &
                                meshdata%RWS, &
                                meshdata%RMT, &
                                meshdata%IPAN, &
                                meshdata%IRC, &
                                meshdata%IMT, &
                                meshdata%IRNS, &
                                meshdata%IRWS, &
                                meshdata%IRMIN, &
                                meshdata%R, &
                                meshdata%DRDI, &
                                meshdata%IRCUT, &
#ifdef USE_OLD_MESH
                                meshdata%LLMSP, &
                                meshdata%THETAS, &
#endif
                                MAGIC_NUMBER

  end function

  !----------------------------------------------------------------------------
  !> Opens RadialMeshData direct access file.
  subroutine openRadialMeshDataDAFile(self, fileunit, filename, max_reclen)
    type(RadialMeshData), intent(in) :: self
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: max_reclen

#ifndef TASKLOCAL_FILES
    integer :: reclen
    
    if (present(max_reclen)) then
      reclen = max_reclen
    else
      reclen = getMinReclenMesh(self)
    endif

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')
#endif
  endsubroutine ! open

  !----------------------------------------------------------------------------
  !> Closes RadialMeshData direct access file.
  subroutine closeRadialMeshDataDAFile(fileunit)
    integer, intent(in) :: fileunit
#ifndef TASKLOCAL_FILES
    close(fileunit)
#endif
  endsubroutine ! close

!===========  Index file ======================================================

  !----------------------------------------------------------------------------
  !> Write mesh dimension data to direct access file 'fileunit' at record 'recnr'
  subroutine writeRadialMeshDataIndexDA(self, fileunit, recnr, max_reclen)
    type(RadialMeshData), intent(in) :: self
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr
    integer, intent(in) :: max_reclen

    integer, parameter :: MAGIC_NUMBER = -889271554

    FILEWRITE (fileunit, rec=recnr) meshdata%irmd, &
                                meshdata%ipand, &
#ifdef USE_OLD_MESH
                                meshdata%meshn, &
                                meshdata%NFU, &
#endif
                                max_reclen, &
                                MAGIC_NUMBER + recnr

  end subroutine

  !----------------------------------------------------------------------------
  !> Read mesh dimension data from direct access file 'fileunit' at record 'recnr'
  !>
  !> Returns dimensions irmd and ipand
  subroutine readRadialMeshDataIndexDA(meshdata, fileunit, recnr, &
                                       irmd, ipand, max_reclen, meshn, nfu)
    implicit none

    type (RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr
    integer, intent(out) :: irmd
    integer, intent(out) :: ipand
    integer, intent(out) :: max_reclen
    integer, intent(out), optional :: meshn
    integer, intent(out), optional :: nfu

#ifdef TASKLOCAL_FILES
    character(len=16) :: num

    write(unit=num, fmt='(a,i7.7)') "mesh.",recnr
    open(fileunit, file=num, form='unformatted', action='read', status='old')
#endif

    call readRadialMeshDataHeader(meshdata, fileunit, recnr, &
                                       irmd, ipand, max_reclen &
#ifdef USE_OLD_MESH
                                       , meshn, nfu &
#endif
                                       )

#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Opens RadialMeshData index file.
  subroutine openRadialMeshDataIndexDAFile(self, fileunit, filename)
    type(RadialMeshData), intent(in) :: self
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename

#ifndef TASKLOCAL_FILES
    integer :: reclen, max_reclen

    max_reclen = 0
    inquire (iolength = reclen) meshdata%irmd, &
                                meshdata%ipand, &
#ifdef USE_OLD_MESH
                                meshdata%meshn, &
                                meshdata%NFU, &
#endif
                                max_reclen, &
                                MAGIC_NUMBER

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')
#endif
  endsubroutine ! open

  !----------------------------------------------------------------------------
  !> Closes RadialMeshData index file.
  subroutine closeRadialMeshDataIndexDAFile(fileunit)
    integer, intent(in) :: fileunit

#ifndef TASKLOCAL_FILES
    close(fileunit)
#endif

  endsubroutine ! close

  !----------------------------------------------------------------------------
  !> Returns a string representation of RadialMeshData.
  !>
  !> Example usage:
  !>
  !>    character(len=:), allocatable :: str ! needs variable length string
  !>
  !>    call repr_RadialMeshData(old_mesh, str)
  !>    write(*,'(A)') str  ! format (A) necessary to avoid messy output
  subroutine repr_RadialMeshData(self, str)
    class(RadialMeshData), intent(in) :: self
    character(len=:), allocatable, intent(inout) :: str

    character :: nl
    character(len=80) :: buffer
    integer :: ind

    nl = new_line(' ')

    str = ''
    write(buffer, *) "irmd  = ", self%irmd   !< number of mesh points
    str = str // trim(buffer) // nl
    write(buffer, *) "ipand = ", self%ipand  !< dimension variable panels for historical reasons
    str = str // trim(buffer) // nl
    write(buffer, *) "A     = ", self%A    !< logarithmic mesh parameter A
    str = str // trim(buffer) // nl
    write(buffer, *) "B     = ", self%B    !< logarithmic mesh parameter A
    str = str // trim(buffer) // nl
    write(buffer, *) "RWS   = ", self%RWS  !< maximal radius
    str = str // trim(buffer) // nl
    write(buffer, *) "RMT   = ", self%RMT !< muffin-tin radius
    str = str // trim(buffer) // nl
    write(buffer, *) "IPAN  = ", self%IPAN !< number of mesh panels
    str = str // trim(buffer) // nl
    write(buffer, *) "IRC   = ", self%IRC
    str = str // trim(buffer) // nl
    write(buffer, *) "IMT   = ", self%IMT    !< endof muffin-tin region
    str = str // trim(buffer) // nl
    write(buffer, *) "IRNS  = ", self%IRNS
    str = str // trim(buffer) // nl
    write(buffer, *) "IRWS  = ", self%IRWS  !< index of max. radius
    str = str // trim(buffer) // nl
    write(buffer, *) "IRMIN = ", self%IRMIN
    str = str // trim(buffer) // nl // nl
    write(buffer, *) "nr.    R                          DRDI"
    str = str // trim(buffer) // nl
    write(buffer, '(79("="))')
    str = str // trim(buffer) // nl

    do ind = 1, size(self%R)
      write(buffer, '(I5, 2X, E23.16,2X,E23.16)') ind, self%R(ind), self%DRDI(ind)
      str = str // trim(buffer) // nl
    enddo ! ind
    str = str // nl

    write(buffer, *) "IRCUT = "
    str = str // trim(buffer) // nl

    do ind = 0, size(self%IRCUT) - 1
      write(buffer, *) self%IRCUT(ind)
      str = str // trim(buffer) // nl
    enddo ! ind
    
  endsubroutine ! represent

!================= private helper functions ===================================
  subroutine readRadialMeshDataHeader(meshdata, fileunit, recnr, &
                                       irmd, ipand, max_reclen, meshn, nfu)
    implicit none

    type (RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr
    integer, intent(out) :: irmd
    integer, intent(out) :: ipand
    integer, intent(out) :: max_reclen
    integer, intent(out), optional :: meshn
    integer, intent(out), optional :: nfu

    integer :: magic, checkmagic

    checkmagic = MAGIC_NUMBER + recnr

    FILEREAD  (fileunit, rec=recnr) irmd, &
                                ipand, &
#ifdef USE_OLD_MESH                                
                                meshn, &
                                nfu, &
#endif
                                max_reclen, &
                                magic

    if (magic /= checkmagic) then
      write (*,*) "ERROR: Invalid mesh index data read. ", __FILE__, __LINE__
      STOP
    end if

    if (magic /= checkmagic) die_here("Invalid mesh index data read.") 
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Test consistency of mesh parameters.
  subroutine test_mesh_consistency(self)
    class(RadialMeshData), intent(in) :: self

    double precision :: rmt_test
    double precision, parameter :: TOLERANCE = 1.d-8

    rmt_test = self%B * (exp(self%A * (self%IMT - 1)) - 1.0d0)

    if (abs(rmt_test - self%rmt) > TOLERANCE) &
      die_here("Mesh parameters A, B, IMT inconsistent with Rmt, rmt ="+rmt_test+"but self%rmt ="+self%rmt)

    if (abs(self%r(self%IMT) - self%rmt) > TOLERANCE) &
      die_here("radial value at self%imt not consistent with Rmt, r(imt) ="+self%r(self%IMT)+"but rmt ="+self%rmt)

    if (self%irc /= self%irmd) &
      die_here("self%irc ="+self%irc+"not equal to self%irmd ="+self%irmd)

    if (self%irws /= self%irmd) &
      die_here("self%irws ="+self%irws+"not equal to self%irmd ="+self%irmd)

    if (abs(self%r(self%IRWS) - self%rws) > TOLERANCE) &
      die_here("radial value at self%irws ="+self%r(self%IRWS)+"not consistent with rWS ="+self%rws)

    if (self%ircut(0) /= 0) &
      die_here("self%ircut(0) ="+self%ircut(0)+"is not 0.")

    if (self%ircut(1) /= self%imt) &
      die_here("self%ircut(1) ="+self%ircut(1)+" is not equal to self%imt ="+self%imt)

    if (self%ircut(self%ipan) /= self%irmd) &
      die_here("self%ircut(self%ipan) ="+self%ircut(self%ipan)+"is not equal to self%irmd ="+self%irmd)

  endsubroutine ! test

endmodule ! RadialMeshData_mod
