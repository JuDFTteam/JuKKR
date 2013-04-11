
!------------------------------------------------------------------------------
! Automatically generated source file. Do not edit by hand.
! To add/remove/modify input parameters:
! Edit dimensions.txt and run 
! 'inputgenerator.py Dimensions dimensions.txt > source.f90'
! to generate source code.
!------------------------------------------------------------------------------


module Dimensions_mod
type Dimensions
  integer :: LMAXD
  integer :: NSPIND
  integer :: NAEZD
  integer :: IRNSD
  integer :: IRMD
  integer :: NREFD
  integer :: NRD
  integer :: IRID
  integer :: NFUND
  integer :: NCELLD
  integer :: NACLSD
  integer :: NCLSD
  integer :: IPAND
  integer :: NXIJD
  integer :: KPOIBZ
  integer :: IGUESSD
  integer :: BCPD
  integer :: NMAXD
  integer :: ISHLD
  integer :: LLY
  integer :: SMPID
  integer :: EMPID
  integer :: NTHRDS
  integer :: XDIM
  integer :: YDIM
  integer :: ZDIM
  integer :: NATBLD
  integer :: ITDBRYD
  integer :: num_atom_procs
end type Dimensions

CONTAINS
!-------------------------------------------------------------------------------
integer function getDimensionsValues(filename, confvalues) result(ierror)
  use Config_Reader
  implicit none

  character(len=*) :: filename
  type (ConfigReader) :: conf
  type (Dimensions), intent(inout) :: confvalues

  ierror = 0
  call createConfigReader(conf)
  call parseFile(conf, filename, ierror)
  if (ierror /= 0) then
    write(*,*) "Error reading configfile", filename
    call destroyConfigReader(conf)
    return
  end if


  call getValueInteger(conf, "LMAXD", confvalues%LMAXD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for LMAXD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NSPIND", confvalues%NSPIND, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NSPIND."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NAEZD", confvalues%NAEZD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NAEZD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "IRNSD", confvalues%IRNSD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for IRNSD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "IRMD", confvalues%IRMD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for IRMD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NREFD", confvalues%NREFD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NREFD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NRD", confvalues%NRD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NRD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "IRID", confvalues%IRID, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for IRID."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NFUND", confvalues%NFUND, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NFUND."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NCELLD", confvalues%NCELLD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NCELLD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NACLSD", confvalues%NACLSD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NACLSD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NCLSD", confvalues%NCLSD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NCLSD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "IPAND", confvalues%IPAND, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for IPAND."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NXIJD", confvalues%NXIJD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NXIJD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "KPOIBZ", confvalues%KPOIBZ, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for KPOIBZ."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "IGUESSD", confvalues%IGUESSD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for IGUESSD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "BCPD", confvalues%BCPD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for BCPD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NMAXD", confvalues%NMAXD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NMAXD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "ISHLD", confvalues%ISHLD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for ISHLD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "LLY", confvalues%LLY, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for LLY."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "SMPID", confvalues%SMPID, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for SMPID."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "EMPID", confvalues%EMPID, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for EMPID."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NTHRDS", confvalues%NTHRDS, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NTHRDS."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "XDIM", confvalues%XDIM, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for XDIM."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "YDIM", confvalues%YDIM, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for YDIM."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "ZDIM", confvalues%ZDIM, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for ZDIM."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "NATBLD", confvalues%NATBLD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for NATBLD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "ITDBRYD", confvalues%ITDBRYD, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for ITDBRYD."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "num_atom_procs", confvalues%num_atom_procs, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for num_atom_procs."
    call destroyConfigReader(conf)
    return
  end if
  call destroyConfigReader(conf)
end function

!-------------------------------------------------------------------------------
integer function readDimensionsFromFile(filename, confvalues) result(ierror)
  implicit none
  character(len=*), intent(in) :: filename
  type (Dimensions), intent(inout) :: confvalues

  integer, parameter :: FILEHANDLE = 67

  ierror = 0
  open(FILEHANDLE, file=filename, form="unformatted")
  read(FILEHANDLE) confvalues%LMAXD
  read(FILEHANDLE) confvalues%NSPIND
  read(FILEHANDLE) confvalues%NAEZD
  read(FILEHANDLE) confvalues%IRNSD
  read(FILEHANDLE) confvalues%IRMD
  read(FILEHANDLE) confvalues%NREFD
  read(FILEHANDLE) confvalues%NRD
  read(FILEHANDLE) confvalues%IRID
  read(FILEHANDLE) confvalues%NFUND
  read(FILEHANDLE) confvalues%NCELLD
  read(FILEHANDLE) confvalues%NACLSD
  read(FILEHANDLE) confvalues%NCLSD
  read(FILEHANDLE) confvalues%IPAND
  read(FILEHANDLE) confvalues%NXIJD
  read(FILEHANDLE) confvalues%KPOIBZ
  read(FILEHANDLE) confvalues%IGUESSD
  read(FILEHANDLE) confvalues%BCPD
  read(FILEHANDLE) confvalues%NMAXD
  read(FILEHANDLE) confvalues%ISHLD
  read(FILEHANDLE) confvalues%LLY
  read(FILEHANDLE) confvalues%SMPID
  read(FILEHANDLE) confvalues%EMPID
  read(FILEHANDLE) confvalues%NTHRDS
  read(FILEHANDLE) confvalues%XDIM
  read(FILEHANDLE) confvalues%YDIM
  read(FILEHANDLE) confvalues%ZDIM
  read(FILEHANDLE) confvalues%NATBLD
  read(FILEHANDLE) confvalues%ITDBRYD
  read(FILEHANDLE) confvalues%num_atom_procs
  close(FILEHANDLE)
end function

!-------------------------------------------------------------------------------
integer function writeDimensionsToFile(filename, confvalues) result(ierror)
  implicit none
  character(len=*), intent(in) :: filename
  type (Dimensions), intent(inout) :: confvalues

  integer, parameter :: FILEHANDLE = 67

  ierror = 0
  open(FILEHANDLE, file=filename, form="unformatted")
  write(FILEHANDLE) confvalues%LMAXD
  write(FILEHANDLE) confvalues%NSPIND
  write(FILEHANDLE) confvalues%NAEZD
  write(FILEHANDLE) confvalues%IRNSD
  write(FILEHANDLE) confvalues%IRMD
  write(FILEHANDLE) confvalues%NREFD
  write(FILEHANDLE) confvalues%NRD
  write(FILEHANDLE) confvalues%IRID
  write(FILEHANDLE) confvalues%NFUND
  write(FILEHANDLE) confvalues%NCELLD
  write(FILEHANDLE) confvalues%NACLSD
  write(FILEHANDLE) confvalues%NCLSD
  write(FILEHANDLE) confvalues%IPAND
  write(FILEHANDLE) confvalues%NXIJD
  write(FILEHANDLE) confvalues%KPOIBZ
  write(FILEHANDLE) confvalues%IGUESSD
  write(FILEHANDLE) confvalues%BCPD
  write(FILEHANDLE) confvalues%NMAXD
  write(FILEHANDLE) confvalues%ISHLD
  write(FILEHANDLE) confvalues%LLY
  write(FILEHANDLE) confvalues%SMPID
  write(FILEHANDLE) confvalues%EMPID
  write(FILEHANDLE) confvalues%NTHRDS
  write(FILEHANDLE) confvalues%XDIM
  write(FILEHANDLE) confvalues%YDIM
  write(FILEHANDLE) confvalues%ZDIM
  write(FILEHANDLE) confvalues%NATBLD
  write(FILEHANDLE) confvalues%ITDBRYD
  write(FILEHANDLE) confvalues%num_atom_procs
  close(FILEHANDLE)
end function

end module
