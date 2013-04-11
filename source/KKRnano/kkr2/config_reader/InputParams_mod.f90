
!------------------------------------------------------------------------------
! Automatically generated source file. Do not edit by hand.
! To add/remove/modify input parameters:
! Edit InputParamsNew.txt and run 
! 'inputgenerator.py InputParams InputParamsNew.txt > source.f90'
! to generate source code.
!------------------------------------------------------------------------------


module InputParams_mod
type InputParams
  integer :: icst
  integer :: kpre
  double precision :: bravais_b (3)
  double precision :: bravais_c (3)
  double precision :: gmax
  double precision :: qmrbound
  integer :: nsra
  integer :: kxc
  double precision :: rmax
  double precision :: tempr
  double precision :: alat
  double precision :: rclust
  double precision :: mixing
  double precision :: fcm
  integer :: scfsteps
  integer :: kforce
  double precision :: emin
  double precision :: bravais_a (3)
  logical :: jij
  logical :: ldau
  double precision :: rcutjij
  double precision :: basisscale (3)
  double precision :: emax
  integer :: kte
  integer :: imix
  integer :: bzdivide (3)
  logical :: cartesian
  integer :: npol
  integer :: npnt1
  integer :: npnt2
  integer :: npnt3
end type InputParams

CONTAINS
!-------------------------------------------------------------------------------
integer function getInputParamsValues(filename, confvalues) result(ierror)
  use Config_Reader
  implicit none

  character(len=*) :: filename
  type (ConfigReader) :: conf
  type (InputParams), intent(inout) :: confvalues

  ierror = 0
  call createConfigReader(conf)
  call parseFile(conf, filename, ierror)
  if (ierror /= 0) then
    write(*,*) "Error reading configfile", filename
    call destroyConfigReader(conf)
    return
  end if


  call getValueInteger(conf, "icst", confvalues%icst, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for icst."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "kpre", confvalues%kpre, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for kpre."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDoubleVector(conf, "bravais_b", confvalues%bravais_b, 3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_b."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDoubleVector(conf, "bravais_c", confvalues%bravais_c, 3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_c."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "gmax", confvalues%gmax, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for gmax."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "qmrbound", confvalues%qmrbound, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for qmrbound."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "nsra", confvalues%nsra, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for nsra."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "kxc", confvalues%kxc, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for kxc."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "rmax", confvalues%rmax, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rmax."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "tempr", confvalues%tempr, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for tempr."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "alat", confvalues%alat, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for alat."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "rclust", confvalues%rclust, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rclust."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "mixing", confvalues%mixing, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for mixing."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "fcm", confvalues%fcm, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for fcm."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "scfsteps", confvalues%scfsteps, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for scfsteps."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "kforce", confvalues%kforce, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for kforce."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "emin", confvalues%emin, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for emin."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDoubleVector(conf, "bravais_a", confvalues%bravais_a, 3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_a."
    call destroyConfigReader(conf)
    return
  end if
  call getValueLogical(conf, "jij", confvalues%jij, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for jij."
    call destroyConfigReader(conf)
    return
  end if
  call getValueLogical(conf, "ldau", confvalues%ldau, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for ldau."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "rcutjij", confvalues%rcutjij, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rcutjij."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDoubleVector(conf, "basisscale", confvalues%basisscale, 3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for basisscale."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "emax", confvalues%emax, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for emax."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "kte", confvalues%kte, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for kte."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "imix", confvalues%imix, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for imix."
    call destroyConfigReader(conf)
    return
  end if
  call getValueIntVector(conf, "bzdivide", confvalues%bzdivide, 3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bzdivide."
    call destroyConfigReader(conf)
    return
  end if
  call getValueLogical(conf, "cartesian", confvalues%cartesian, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for cartesian."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "npol", confvalues%npol, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npol."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "npnt1", confvalues%npnt1, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt1."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "npnt2", confvalues%npnt2, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt2."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "npnt3", confvalues%npnt3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt3."
    call destroyConfigReader(conf)
    return
  end if
  call destroyConfigReader(conf)
end function

!-------------------------------------------------------------------------------
integer function readInputParamsFromFile(filename, confvalues) result(ierror)
  implicit none
  character(len=*), intent(in) :: filename
  type (InputParams), intent(inout) :: confvalues

  integer, parameter :: FILEHANDLE = 67

  ierror = 0
  open(FILEHANDLE, file=filename, form="unformatted")
  read(FILEHANDLE) confvalues%icst
  read(FILEHANDLE) confvalues%kpre
  read(FILEHANDLE) confvalues%bravais_b
  read(FILEHANDLE) confvalues%bravais_c
  read(FILEHANDLE) confvalues%gmax
  read(FILEHANDLE) confvalues%qmrbound
  read(FILEHANDLE) confvalues%nsra
  read(FILEHANDLE) confvalues%kxc
  read(FILEHANDLE) confvalues%rmax
  read(FILEHANDLE) confvalues%tempr
  read(FILEHANDLE) confvalues%alat
  read(FILEHANDLE) confvalues%rclust
  read(FILEHANDLE) confvalues%mixing
  read(FILEHANDLE) confvalues%fcm
  read(FILEHANDLE) confvalues%scfsteps
  read(FILEHANDLE) confvalues%kforce
  read(FILEHANDLE) confvalues%emin
  read(FILEHANDLE) confvalues%bravais_a
  read(FILEHANDLE) confvalues%jij
  read(FILEHANDLE) confvalues%ldau
  read(FILEHANDLE) confvalues%rcutjij
  read(FILEHANDLE) confvalues%basisscale
  read(FILEHANDLE) confvalues%emax
  read(FILEHANDLE) confvalues%kte
  read(FILEHANDLE) confvalues%imix
  read(FILEHANDLE) confvalues%bzdivide
  read(FILEHANDLE) confvalues%cartesian
  read(FILEHANDLE) confvalues%npol
  read(FILEHANDLE) confvalues%npnt1
  read(FILEHANDLE) confvalues%npnt2
  read(FILEHANDLE) confvalues%npnt3
  close(FILEHANDLE)
end function

!-------------------------------------------------------------------------------
integer function writeInputParamsToFile(filename, confvalues) result(ierror)
  implicit none
  character(len=*), intent(in) :: filename
  type (InputParams), intent(inout) :: confvalues

  integer, parameter :: FILEHANDLE = 67

  ierror = 0
  open(FILEHANDLE, file=filename, form="unformatted")
  write(FILEHANDLE) confvalues%icst
  write(FILEHANDLE) confvalues%kpre
  write(FILEHANDLE) confvalues%bravais_b
  write(FILEHANDLE) confvalues%bravais_c
  write(FILEHANDLE) confvalues%gmax
  write(FILEHANDLE) confvalues%qmrbound
  write(FILEHANDLE) confvalues%nsra
  write(FILEHANDLE) confvalues%kxc
  write(FILEHANDLE) confvalues%rmax
  write(FILEHANDLE) confvalues%tempr
  write(FILEHANDLE) confvalues%alat
  write(FILEHANDLE) confvalues%rclust
  write(FILEHANDLE) confvalues%mixing
  write(FILEHANDLE) confvalues%fcm
  write(FILEHANDLE) confvalues%scfsteps
  write(FILEHANDLE) confvalues%kforce
  write(FILEHANDLE) confvalues%emin
  write(FILEHANDLE) confvalues%bravais_a
  write(FILEHANDLE) confvalues%jij
  write(FILEHANDLE) confvalues%ldau
  write(FILEHANDLE) confvalues%rcutjij
  write(FILEHANDLE) confvalues%basisscale
  write(FILEHANDLE) confvalues%emax
  write(FILEHANDLE) confvalues%kte
  write(FILEHANDLE) confvalues%imix
  write(FILEHANDLE) confvalues%bzdivide
  write(FILEHANDLE) confvalues%cartesian
  write(FILEHANDLE) confvalues%npol
  write(FILEHANDLE) confvalues%npnt1
  write(FILEHANDLE) confvalues%npnt2
  write(FILEHANDLE) confvalues%npnt3
  close(FILEHANDLE)
end function

end module
