!---------------------------------------------------------------------
! Writes array dimension parameters (previously in inc.p and inc.cls)
! to file inp0.unf
!---------------------------------------------------------------------
! missing: TRC (not used in kkr2)
! could be removed: KREL ?
subroutine write_dimension_parameters( &
LMAXD, &
NSPIND, &
NAEZD, &
IRNSD, &
IRMD, &
NREFD, &
NRD, &
IRID, &
NFUND, &
NCELLD, &
NACLSD, &
NCLSD, &
IPAND, &
NXIJD, &
KPOIBZ, &
IGUESSD, &
BCPD, &
NMAXD, &
ISHLD, &
LLY, &
SMPID, &
EMPID, &
NTHRDS, &
XDIM, &
YDIM, &
ZDIM, &
NATBLD, &
ITDBRYD, &
IEMXD, &
EKMD)

  implicit none

  integer, intent(in) :: LMAXD
  integer, intent(in) :: NSPIND
  integer, intent(in) :: NAEZD
  integer, intent(in) :: IRNSD
  integer, intent(in) :: IRMD
  integer, intent(in) :: NREFD
  integer, intent(in) :: NRD
  integer, intent(in) :: IRID
  integer, intent(in) :: NFUND
  integer, intent(in) :: NCELLD
  integer, intent(in) :: NACLSD
  integer, intent(in) :: NCLSD
  integer, intent(in) :: IPAND
  integer, intent(in) :: NXIJD
  integer, intent(in) :: KPOIBZ
  integer, intent(in) :: IGUESSD
  integer, intent(in) :: BCPD
  integer, intent(in) :: NMAXD
  integer, intent(in) :: ISHLD
  integer, intent(in) :: LLY
  integer, intent(in) :: SMPID
  integer, intent(in) :: EMPID
  integer, intent(in) :: NTHRDS
  integer, intent(in) :: XDIM
  integer, intent(in) :: YDIM
  integer, intent(in) :: ZDIM
  integer, intent(in) :: NATBLD
  integer, intent(in) :: ITDBRYD
  integer, intent(in) :: IEMXD
  integer, intent(in) :: EKMD

  integer, parameter :: FILEHANDLE = 67

  open (FILEHANDLE, FILE='inp0.unf', FORM='unformatted')

  write(FILEHANDLE) LMAXD
  write(FILEHANDLE) NSPIND
  write(FILEHANDLE) NAEZD
  write(FILEHANDLE) IRNSD
  write(FILEHANDLE) IRMD
  write(FILEHANDLE) NREFD
  write(FILEHANDLE) NRD
  write(FILEHANDLE) IRID
  write(FILEHANDLE) NFUND
  write(FILEHANDLE) NCELLD
  write(FILEHANDLE) NACLSD
  write(FILEHANDLE) NCLSD
  write(FILEHANDLE) IPAND
  write(FILEHANDLE) NXIJD
  write(FILEHANDLE) KPOIBZ
  write(FILEHANDLE) IGUESSD
  write(FILEHANDLE) BCPD
  write(FILEHANDLE) NMAXD
  write(FILEHANDLE) ISHLD
  write(FILEHANDLE) LLY
  write(FILEHANDLE) SMPID
  write(FILEHANDLE) EMPID
  write(FILEHANDLE) NTHRDS
  write(FILEHANDLE) XDIM
  write(FILEHANDLE) YDIM
  write(FILEHANDLE) ZDIM
  write(FILEHANDLE) NATBLD
  write(FILEHANDLE) ITDBRYD
  write(FILEHANDLE) IEMXD
  write(FILEHANDLE) EKMD

  close(FILEHANDLE)

end subroutine write_dimension_parameters
