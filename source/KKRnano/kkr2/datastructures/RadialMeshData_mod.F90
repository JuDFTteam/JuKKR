!> This module defines a datatype that contains data related to the radial mesh
!> @author Elias Rabel
module RadialMeshData_mod
  implicit none

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
  !> Write mesh data to direct access file 'fileunit' at record 'recnr'
  subroutine writeRadialMeshDataDA(meshdata, fileunit, recnr)

    implicit none
    type (RadialMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = -889271554

    write (fileunit, rec=recnr) MAGIC_NUMBER, &
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

    read  (fileunit, rec=recnr) magic, &
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

    if (magic /= MAGIC_NUMBER .or. magic2 /= MAGIC_NUMBER) then
      write (*,*) "ERROR: Invalid cell data read. ", __FILE__, __LINE__
      STOP
    end if

  end subroutine

  !----------------------------------------------------------------------------
  !> Opens RadialMeshData direct access file.
  subroutine openRadialMeshDataDAFile(meshdata, fileunit, filename)
    implicit none

    type (RadialMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename
    !------
    integer :: reclen
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

    write (*,*) reclen

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')

  end subroutine

  !----------------------------------------------------------------------------
  !> Closes RadialMeshData direct access file.
  subroutine closeRadialMeshDataDAFile(fileunit)
    implicit none
    integer, intent(in) :: fileunit

    close(fileunit)

  end subroutine


end module
