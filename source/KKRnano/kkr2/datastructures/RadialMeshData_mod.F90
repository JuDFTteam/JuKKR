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
  implicit none

  !private :: initMuffinTinMesh
  private :: initInterstitialMesh

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
    integer :: IMT    !< end of muffin-tin region
    integer :: IRNS
    integer :: IRWS   !< index of max. radius
    integer :: IRMIN

    ! arrays
    double precision, dimension(:), allocatable :: R
    double precision, dimension(:), allocatable :: DRDI
    integer, dimension(:), allocatable :: IRCUT  !< panel locations
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createRadialMeshData(meshdata, irmd, ipand)
    implicit none
    type (RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: irmd
    integer, intent(in) :: ipand

    meshdata%irmd = irmd
    meshdata%ipand = ipand

    allocate(meshdata%R(irmd))
    allocate(meshdata%DRDI(irmd))
    allocate(meshdata%IRCUT(0:ipand))

    meshdata%R = 0.0d0
    meshdata%DRDI = 0.0d0
    meshdata%IRCUT = 0

    meshdata%A = 0.0d0
    meshdata%B = 0.0d0

    meshdata%IPAN = 0
    meshdata%IRC = 0
    meshdata%IMT = 0
    meshdata%IRNS = 0
    meshdata%IRWS = 0
    meshdata%IRMIN = 0

    meshdata%RWS = 0.0d0
    meshdata%RMT = 0.0d0

  end subroutine


  !----------------------------------------------------------------------------
  subroutine destroyRadialMeshData(meshdata)
    implicit none
    type (RadialMeshData), intent(inout) :: meshdata

    deallocate(meshdata%R)
    deallocate(meshdata%DRDI)
    deallocate(meshdata%IRCUT)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine initRadialMesh(meshdata, alat, xrn, drn, nm, imt, irns)
    implicit none
    type (RadialMeshData), intent(inout) :: meshdata
    double precision, intent(in) :: alat
    double precision, intent(in) :: xrn(:)
    double precision, intent(in) :: drn(:)
    integer, intent(in) :: nm(:)
    integer, intent(in) :: imt
    integer, intent(in) :: irns

    call initInterstitialMesh(meshdata, alat, xrn, drn, nm, imt, irns)
    ! note radius_mt = xrn(1) * alat
    call initMuffinTinMesh(meshdata, imt, xrn(1)*alat)

  end subroutine

  !---------------------------------------------------------------------------
  ! note radius_mt = xrn(1) * alat - in units of Bohr
  subroutine initMuffinTinMesh(meshdata, imt, radius_mt)
    implicit none
    type (RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: imt
    double precision, intent(in) :: radius_mt

    double precision, parameter :: A = 0.025d0
    integer :: ii

    meshdata%A = A
    meshdata%B = radius_mt / (exp(A * (imt - 1)) - 1.d0)

    do ii = 1, imt
      meshdata%r(ii)    = meshdata%B * (exp(A * (ii - 1)) - 1.d0)
      meshdata%drdi(ii) = A * meshdata%B * exp(A * (ii - 1))
    end do

    meshdata%rmt = radius_MT

  end subroutine

  !---------------------------------------------------------------------------
  ! note radius_mt = xrn(1) * alat
  subroutine initInterstitialMesh(meshdata, alat, xrn, drn, nm, imt, irns)
    implicit none
    type (RadialMeshData), intent(inout) :: meshdata
    double precision, intent(in) :: alat
    double precision, intent(in) :: xrn(:)
    double precision, intent(in) :: drn(:)
    integer, intent(in) :: nm(:)
    integer, intent(in) :: imt
    integer, intent(in) :: irns

    integer :: isum, ii
    integer :: ipan

    ipan = size(nm) + 1 ! +1 for muffin-tin region 1..imt
    meshdata%ipan = ipan
    meshdata%imt = imt

    ! ircut(0) has to be 0, integrations start at ircut(i)+1
    meshdata%ircut(0) = 0
    meshdata%ircut(1) = imt

    isum = imt
    do ii = 2,ipan
      isum = isum + nm(ii - 1)
      meshdata%ircut(ii) = isum
    end do
    meshdata%irc = meshdata%ircut(ipan)

    meshdata%irws = isum
    if (isum > meshdata%irmd) then
      write(*,*) "Error creating mesh, irmd too small", isum, meshdata%irmd
      STOP
    endif

    do ii = 1, meshdata%irws - imt
      meshdata%r(ii + imt) = xrn(ii) * alat
      meshdata%drdi(ii + imt) = drn(ii) * alat
    end do

    meshdata%irmin = meshdata%irws - irns
    meshdata%irns = irns
    meshdata%rws  = meshdata%r(meshdata%irws)

  end subroutine

!==============================================================================
!=                                    I/O                                     =
!==============================================================================

  !----------------------------------------------------------------------------
  !> Read and create radial mesh data from files <filename> and <filename>.idx
  !>
  !> The index file with extension .idx is necessary because different
  !> atoms can have a different number of radial mesh points
  subroutine createRadialMeshDataFromFile(meshdata, filename, recnr)
    implicit none
    type (RadialMeshData), intent(inout) :: meshdata
    character(len=*), intent(in) :: filename
    integer, intent(in) :: recnr

    integer :: FILEUNIT = 37
    integer :: irmd, ipand, max_reclen

    ! index file has extension .idx

    call openRadialMeshDataIndexDAFile(meshdata, FILEUNIT, &
                                       filename // ".idx") !ignored for task-local files
    call readRadialMeshDataIndexDA(meshdata, FILEUNIT, recnr, &
                                       irmd, ipand, max_reclen)
    call closeRadialMeshDataIndexDAFile(FILEUNIT)

    call createRadialMeshData(meshdata, irmd, ipand)

    call openRadialMeshDataDAFile(meshdata, FILEUNIT, filename, max_reclen)
    call readRadialMeshDataDA(meshdata, FILEUNIT, recnr)
    call closeRadialMeshDataDAFile(FILEUNIT)

  end subroutine

  !----------------------------------------------------------------------------
  !> Write mesh data to direct access file 'fileunit' at record 'recnr'
  subroutine writeRadialMeshDataDA(meshdata, fileunit, recnr)

    implicit none
    type (RadialMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = -889271554

#ifdef TASKLOCAL_FILES
    character(len=7) :: num
    write(num, '(I7.7)') recnr
    open(fileunit, file="mesh." // num, form='unformatted')

    call writeRadialMeshDataIndexDA(meshdata, FILEUNIT, recnr, 0)
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

  end subroutine

  !----------------------------------------------------------------------------
  !> Read mesh data from direct access file 'fileunit' at record 'recnr'
  subroutine readRadialMeshDataDA(meshdata, fileunit, recnr)
    implicit none

    type (RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = -889271554
    integer :: magic, magic2
    integer :: checkmagic

#ifdef TASKLOCAL_FILES
    character(len=7) :: num
    integer :: ipand, irmd, max_reclen    

    write(num, '(I7.7)') recnr
    open(fileunit, file="mesh." // num, form='unformatted')

    ! read header at beginning of file
    call readRadialMeshDataHeader(meshdata, FILEUNIT, recnr, irmd, ipand, max_reclen)
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
                                magic2
#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif

    if (magic /= checkmagic .or. magic2 /= checkmagic) then
      write (*,*) "ERROR: Invalid mesh data read. ", __FILE__, __LINE__
      STOP
    end if

  end subroutine

  !---------------------------------------------------------------------------
  !> Return MINIMUM record length needed to store mesh of this atom
  !>
  !> Note: for a file containing meshes of several atoms, the maximum
  !> of all their record lengths has to be determined
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
                                MAGIC_NUMBER

  end function

  !----------------------------------------------------------------------------
  !> Opens RadialMeshData direct access file.
  subroutine openRadialMeshDataDAFile(meshdata, fileunit, filename, max_reclen)
    implicit none

    type (RadialMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: max_reclen
    !------
    integer :: reclen

#ifndef TASKLOCAL_FILES
    if (present(max_reclen)) then
      reclen = max_reclen
    else
      reclen = getMinReclenMesh(meshdata)
    end if

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')
#endif

  end subroutine

  !----------------------------------------------------------------------------
  !> Closes RadialMeshData direct access file.
  subroutine closeRadialMeshDataDAFile(fileunit)
    implicit none
    integer, intent(in) :: fileunit
#ifndef TASKLOCAL_FILES
    close(fileunit)
#endif
  end subroutine

!===========  Index file ======================================================

  !----------------------------------------------------------------------------
  !> Write mesh dimension data to direct access file 'fileunit' at record 'recnr'
  subroutine writeRadialMeshDataIndexDA(meshdata, fileunit, recnr, max_reclen)

    implicit none
    type (RadialMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr
    integer, intent(in) :: max_reclen

    integer, parameter :: MAGIC_NUMBER = -889271554

    FILEWRITE (fileunit, rec=recnr) meshdata%irmd, &
                                meshdata%ipand, &
                                max_reclen, &
                                MAGIC_NUMBER + recnr

  end subroutine

  !----------------------------------------------------------------------------
  !> Read mesh dimension data from direct access file 'fileunit' at record 'recnr'
  !>
  !> Returns dimensions irmd and ipand
  subroutine readRadialMeshDataIndexDA(meshdata, fileunit, recnr, &
                                       irmd, ipand, max_reclen)
    implicit none

    type (RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr
    integer, intent(out) :: irmd
    integer, intent(out) :: ipand
    integer, intent(out) :: max_reclen

#ifdef TASKLOCAL_FILES
    character(len=7) :: num

    write(num, '(I7.7)') recnr
    open(fileunit, file="mesh." // num, form='unformatted')
#endif

    call readRadialMeshDataHeader(meshdata, fileunit, recnr, &
                                       irmd, ipand, max_reclen)

#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif

  end subroutine

  !----------------------------------------------------------------------------
  !> Opens RadialMeshData index file.
  subroutine openRadialMeshDataIndexDAFile(meshdata, fileunit, filename)
    implicit none

    type (RadialMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename
    !------
    integer :: reclen
    integer :: max_reclen
    integer, parameter :: MAGIC_NUMBER = -889271554

#ifndef TASKLOCAL_FILES
    max_reclen = 0
    inquire (iolength = reclen) meshdata%irmd, &
                                meshdata%ipand, &
                                max_reclen, &
                                MAGIC_NUMBER

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')
#endif

  end subroutine

  !----------------------------------------------------------------------------
  !> Closes RadialMeshData index file.
  subroutine closeRadialMeshDataIndexDAFile(fileunit)
    implicit none
    integer, intent(in) :: fileunit

#ifndef TASKLOCAL_FILES
    close(fileunit)
#endif

  end subroutine

  !----------------------------------------------------------------------------
  !> Returns a string representation of RadialMeshData.
  subroutine repr_RadialMeshData(meshdata, str)
    implicit none
    class (RadialMeshData), intent(in) :: meshdata
    character(len=:), allocatable, intent(inout) :: str

    character :: nl
    character(80) :: buffer
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
    write(buffer, *) "IPAN  = ", meshdata%IPAN   !< number of mesh panels
    str = str // trim(buffer) // nl
    write(buffer, *) "IRC   = ", meshdata%IRC
    str = str // trim(buffer) // nl
    write(buffer, *) "IMT   = ", meshdata%IMT    !< end of muffin-tin region
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
    end do
    str = str // nl

    write(buffer, *) "IRCUT = "
    str = str // trim(buffer) // nl

    do ind = 0, size(meshdata%IRCUT) - 1
      write(buffer, *) meshdata%IRCUT(ind)
      str = str // trim(buffer) // nl
    end do
  end subroutine

!================= private helper functions ===================================
  subroutine readRadialMeshDataHeader(meshdata, fileunit, recnr, &
                                       irmd, ipand, max_reclen)
    implicit none

    type (RadialMeshData), intent(inout) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr
    integer, intent(out) :: irmd
    integer, intent(out) :: ipand
    integer, intent(out) :: max_reclen

    integer, parameter :: MAGIC_NUMBER = -889271554
    integer :: magic
    integer :: checkmagic

    checkmagic = MAGIC_NUMBER + recnr

    FILEREAD  (fileunit, rec=recnr) irmd, &
                                ipand, &
                                max_reclen, &
                                magic

    if (magic /= checkmagic) then
      write (*,*) "ERROR: Invalid mesh index data read. ", __FILE__, __LINE__
      STOP
    end if

  end subroutine

end module
