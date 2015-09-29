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
  public :: RadialMeshData, create, destroy, represent
  public :: createRadialMeshData, destroyRadialMeshData ! deprecated
  public :: getMinReclenMesh
  public :: initRadialMesh, initMuffinTinMesh, initInterstitialMesh, createRadialMeshDataFromFile, writeRadialMeshDataDA           
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

    ! arrays
    double precision, allocatable :: R(:)
    double precision, allocatable :: DRDI(:)
    integer, allocatable :: IRCUT(:)  !< panel locations
  endtype

  interface create
    module procedure createRadialMeshData
  endinterface
  
  interface destroy
    module procedure destroyRadialMeshData
  endinterface
  
  interface represent
    module procedure repr_RadialMeshData
  endinterface

  integer, parameter :: MAGIC_NUMBER = -889271554
  
  contains

  !----------------------------------------------------------------------------
  subroutine createRadialMeshData(meshdata, irmd, ipand)
    type(RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: irmd, ipand

    meshdata%irmd = irmd
    meshdata%ipand = ipand

    allocate(meshdata%R(irmd))
    allocate(meshdata%DRDI(irmd))
    allocate(meshdata%IRCUT(0:ipand))

    meshdata%R = 0.d0
    meshdata%DRDI = 0.d0
    meshdata%IRCUT = 0

    meshdata%A = 0.d0
    meshdata%B = 0.d0

    meshdata%IPAN = 0
    meshdata%IRC = 0
    meshdata%IMT = 0
    meshdata%IRNS = 0
    meshdata%IRWS = 0
    meshdata%IRMIN = 0

    meshdata%RWS = 0.d0
    meshdata%RMT = 0.d0

  endsubroutine ! create

  !----------------------------------------------------------------------------
  subroutine destroyRadialMeshData(meshdata)
    type(RadialMeshData), intent(inout) :: meshdata

    deallocate(meshdata%R)
    deallocate(meshdata%DRDI)
    deallocate(meshdata%IRCUT)

  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  subroutine initRadialMesh(meshdata, alat, xrn, drn, nm, imt, irns)
    type(RadialMeshData), intent(inout) :: meshdata
    double precision, intent(in) :: alat
    double precision, intent(in) :: xrn(:), drn(:)
    integer, intent(in) :: nm(:)
    integer, intent(in) :: imt, irns

    call initInterstitialMesh(meshdata, alat, xrn, drn, nm, imt, irns)
    ! note radius_mt = xrn(1) * alat
    call initMuffinTinMesh(meshdata, imt, xrn(1)*alat)

    call test_mesh_consistency(meshdata)

  endsubroutine ! init

  !---------------------------------------------------------------------------
  ! note radius_mt = xrn(1) * alat - in units of Bohr
  subroutine initMuffinTinMesh(meshdata, imt, radius_mt)
    type(RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: imt
    double precision, intent(in) :: radius_mt

    double precision, parameter :: A = 0.025d0
    integer :: ii

    meshdata%A = A
    meshdata%B = radius_mt / (exp(A * (imt - 1)) - 1.d0)

    do ii = 1, imt
      meshdata%r(ii)    = meshdata%B * (exp(A * (ii - 1)) - 1.d0)
      meshdata%drdi(ii) = A * meshdata%B * exp(A * (ii - 1))
    enddo ! ii

    meshdata%rmt = radius_MT

  endsubroutine ! init

  !---------------------------------------------------------------------------
  ! note radius_mt = xrn(1) * alat
  subroutine initInterstitialMesh(meshdata, alat, xrn, drn, nm, imt, irns)
    type(RadialMeshData), intent(inout) :: meshdata
    double precision, intent(in) :: alat
    double precision, intent(in) :: xrn(:), drn(:)
    integer, intent(in) :: nm(:)
    integer, intent(in) :: imt, irns

    integer :: isum, ii, ipan

    ipan = size(nm) + 1 ! +1 for muffin-tin region 1..imt
    meshdata%ipan = ipan
    meshdata%imt = imt

    ! ircut(0) has to be 0, integrations start at ircut(i)+1
    meshdata%ircut(0) = 0
    meshdata%ircut(1) = imt

    isum = imt
    do ii = 2, ipan
      isum = isum + nm(ii - 1)
      meshdata%ircut(ii) = isum
    enddo ! ii
    meshdata%irc = meshdata%ircut(ipan)

    meshdata%irws = isum
    if (isum > meshdata%irmd) &
      die_here("Error creating mesh, irmd too small, isum ="+isum-", meshdata%irmd ="+meshdata%irmd)

    do ii = 1, meshdata%irws - imt
      meshdata%r(ii + imt) = xrn(ii) * alat
      meshdata%drdi(ii + imt) = drn(ii) * alat
    enddo ! ii

    meshdata%irmin = meshdata%irws - irns
    meshdata%irns = irns
    meshdata%rws  = meshdata%r(meshdata%irws)

  endsubroutine ! init

!==============================================================================
!=                                    I/O                                     =
!==============================================================================

  !----------------------------------------------------------------------------
  !> Read and create radial mesh data from files <filename> and <filename>.idx
  !>
  !> The index file with extension .idx is necessary because different
  !> atoms can have a different number of radial mesh points
  subroutine createRadialMeshDataFromFile(meshdata, filename, recnr)
    type(RadialMeshData), intent(inout) :: meshdata
    character(len=*), intent(in) :: filename
    integer, intent(in) :: recnr

    integer, parameter :: fu = 37
    integer :: irmd, ipand, max_reclen

    ! index file has extension .idx but the filename is ignored for task-local files
    call openRadialMeshDataIndexDAFile(meshdata, fu, filename-".idx")
    call readRadialMeshDataIndexDA(meshdata, fu, recnr, irmd, ipand, max_reclen)
    call closeRadialMeshDataIndexDAFile(fu)

    call createRadialMeshData(meshdata, irmd, ipand)

    call openRadialMeshDataDAFile(meshdata, fu, filename, max_reclen)
    call readRadialMeshDataDA(meshdata, fu, recnr)
    call closeRadialMeshDataDAFile(fu)

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Write mesh data to direct access file 'fileunit' at record 'recnr'
  subroutine writeRadialMeshDataDA(meshdata, fileunit, recnr)
    type(RadialMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

#ifdef TASKLOCAL_FILES
    character(len=16) :: num
    write(unit=num, fmt='(a,i7.7)') "mesh.",recnr
    open(fileunit, file=num, form='unformatted', action='write')

    call writeRadialMeshDataIndexDA(meshdata, fileunit, recnr, 0)
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
                                    MAGIC_NUMBER + recnr
#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Read mesh data from direct access file 'fileunit' at record 'recnr'
  subroutine readRadialMeshDataDA(meshdata, fileunit, recnr)
    type(RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: fileunit, recnr

    integer :: magic, magic2, checkmagic

#ifdef TASKLOCAL_FILES
    character(len=16) :: num
    integer :: ipand, irmd, max_reclen    

    write(unit=num, fmt='(a,i7.7)') "mesh.",recnr
    open(fileunit, file=num, form='unformatted', action='read', status='old')

    ! read header at beginning of file
    call readRadialMeshDataHeader(meshdata, fileunit, recnr, irmd, ipand, max_reclen)
#endif

    checkmagic = MAGIC_NUMBER + recnr

    FILEREAD (fileunit, rec=recnr)  magic, &
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
                                    magic2
#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif

    if (magic /= checkmagic .or. magic2 /= checkmagic) die_here("Invalid mesh data read.")

    call test_mesh_consistency(meshdata)

  endsubroutine ! read

  !---------------------------------------------------------------------------
  !> Return MINIMUM record length needed to store mesh of this atom
  !>
  !> Note: for a file containing meshes of several atoms, the maximum
  !> of all their record lengths has to be determined
  integer function getMinReclenMesh(meshdata) result(reclen)
    type(RadialMeshData), intent(in) :: meshdata

    inquire (iolength=reclen) MAGIC_NUMBER, &
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
                              MAGIC_NUMBER

  endfunction ! get

  !----------------------------------------------------------------------------
  !> Opens RadialMeshData direct access file.
  subroutine openRadialMeshDataDAFile(meshdata, fileunit, filename, max_reclen)
    type(RadialMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: max_reclen

#ifndef TASKLOCAL_FILES
    integer :: reclen
    
    if (present(max_reclen)) then
      reclen = max_reclen
    else
      reclen = getMinReclenMesh(meshdata)
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
  subroutine writeRadialMeshDataIndexDA(meshdata, fileunit, recnr, max_reclen)
    type(RadialMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr
    integer, intent(in) :: max_reclen

    FILEWRITE (fileunit, rec=recnr) meshdata%irmd, &
                                    meshdata%ipand, &
                                    max_reclen, &
                                    MAGIC_NUMBER + recnr
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Read mesh dimension data from direct access file 'fileunit' at record 'recnr'
  !>
  !> Returns dimensions irmd and ipand
  subroutine readRadialMeshDataIndexDA(meshdata, fileunit, recnr, irmd, ipand, max_reclen)
    type(RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr
    integer, intent(out) :: irmd
    integer, intent(out) :: ipand
    integer, intent(out) :: max_reclen

#ifdef TASKLOCAL_FILES
    character(len=16) :: num

    write(unit=num, fmt='(a,i7.7)') "mesh.",recnr
    open(fileunit, file=num, form='unformatted', action='read', status='old')
#endif

    call readRadialMeshDataHeader(meshdata, fileunit, recnr, irmd, ipand, max_reclen)

#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Opens RadialMeshData index file.
  subroutine openRadialMeshDataIndexDAFile(meshdata, fileunit, filename)
    type(RadialMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename

#ifndef TASKLOCAL_FILES
    integer :: reclen, max_reclen

    max_reclen = 0
    inquire (iolength=reclen) meshdata%irmd, &
                              meshdata%ipand, &
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
  subroutine repr_RadialMeshData(meshdata, str)
    class(RadialMeshData), intent(in) :: meshdata
    character(len=:), allocatable, intent(inout) :: str

    character :: nl
    character(len=80) :: buffer
    integer :: ind

    nl = new_line(' ')

    str = ''
    write(buffer, *) "irmd  = ", meshdata%irmd   !< number of mesh points
    str = str // trim(buffer) // nl
    write(buffer, *) "ipand = ", meshdata%ipand  !< dimension variable panels for historical reasons
    str = str // trim(buffer) // nl
    write(buffer, *) "A     = ", meshdata%A    !< logarithmic mesh parameter A
    str = str // trim(buffer) // nl
    write(buffer, *) "B     = ", meshdata%B    !< logarithmic mesh parameter A
    str = str // trim(buffer) // nl
    write(buffer, *) "RWS   = ", meshdata%RWS  !< maximal radius
    str = str // trim(buffer) // nl
    write(buffer, *) "RMT   = ", meshdata%RMT !< muffin-tin radius
    str = str // trim(buffer) // nl
    write(buffer, *) "IPAN  = ", meshdata%IPAN !< number of mesh panels
    str = str // trim(buffer) // nl
    write(buffer, *) "IRC   = ", meshdata%IRC
    str = str // trim(buffer) // nl
    write(buffer, *) "IMT   = ", meshdata%IMT    !< endof muffin-tin region
    str = str // trim(buffer) // nl
    write(buffer, *) "IRNS  = ", meshdata%IRNS
    str = str // trim(buffer) // nl
    write(buffer, *) "IRWS  = ", meshdata%IRWS  !< index of max. radius
    str = str // trim(buffer) // nl
    write(buffer, *) "IRMIN = ", meshdata%IRMIN
    str = str // trim(buffer) // nl // nl
    write(buffer, *) "nr.    R                          DRDI"
    str = str // trim(buffer) // nl
    write(buffer, '(79("="))')
    str = str // trim(buffer) // nl

    do ind = 1, size(meshdata%R)
      write(buffer, '(I5, 2X, E23.16,2X,E23.16)') ind, meshdata%R(ind), meshdata%DRDI(ind)
      str = str // trim(buffer) // nl
    enddo ! ind
    str = str // nl

    write(buffer, *) "IRCUT = "
    str = str // trim(buffer) // nl

    do ind = 0, size(meshdata%IRCUT) - 1
      write(buffer, *) meshdata%IRCUT(ind)
      str = str // trim(buffer) // nl
    enddo ! ind
    
  endsubroutine ! represent

!================= private helper functions ===================================
  subroutine readRadialMeshDataHeader(meshdata, fileunit, recnr, irmd, ipand, max_reclen)
    type(RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: fileunit, recnr
    integer, intent(out) :: irmd, ipand, max_reclen

    integer :: magic, checkmagic

    checkmagic = MAGIC_NUMBER + recnr

    FILEREAD  (fileunit, rec=recnr) irmd, &
                                    ipand, &
                                    max_reclen, &
                                    magic

    if (magic /= checkmagic) die_here("Invalid mesh index data read.") 
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Test consistency of mesh parameters.
  subroutine test_mesh_consistency(meshdata)
    class(RadialMeshData), intent(in) :: meshdata

    double precision :: rmt_test
    double precision, parameter :: TOLERANCE = 1.d-8

    rmt_test = meshdata%B * (exp(meshdata%A * (meshdata%IMT - 1)) - 1.0d0)

    if (abs(rmt_test - meshdata%rmt) > TOLERANCE) &
      die_here("Mesh parameters A, B, IMT inconsistent with Rmt, rmt ="+rmt_test+"but meshdata%rmt ="+meshdata%rmt)

    if (abs(meshdata%r(meshdata%IMT) - meshdata%rmt) > TOLERANCE) &
      die_here("radial value at meshdata%imt not consistent with Rmt, r(imt) ="+meshdata%r(meshdata%IMT)+"but rmt ="+meshdata%rmt)

    if (meshdata%irc /= meshdata%irmd) &
      die_here("meshdata%irc ="+meshdata%irc+"not equal to meshdata%irmd ="+meshdata%irmd)

    if (meshdata%irws /= meshdata%irmd) &
      die_here("meshdata%irws ="+meshdata%irws+"not equal to meshdata%irmd ="+meshdata%irmd)

    if (abs(meshdata%r(meshdata%IRWS) - meshdata%rws) > TOLERANCE) &
      die_here("radial value at meshdata%irws ="+meshdata%r(meshdata%IRWS)+"not consistent with rWS ="+meshdata%rws)

    if (meshdata%ircut(0) /= 0) &
      die_here("meshdata%ircut(0) ="+meshdata%ircut(0)+"is not 0.")

    if (meshdata%ircut(1) /= meshdata%imt) &
      die_here("meshdata%ircut(1) ="+meshdata%ircut(1)+" is not equal to meshdata%imt ="+meshdata%imt)

    if (meshdata%ircut(meshdata%ipan) /= meshdata%irmd) &
      die_here("meshdata%ircut(meshdata%ipan) ="+meshdata%ircut(meshdata%ipan)+"is not equal to meshdata%irmd ="+meshdata%irmd)

  endsubroutine ! test

endmodule ! RadialMeshData_mod
