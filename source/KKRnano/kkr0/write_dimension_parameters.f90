!---------------------------------------------------------------------
! Writes array dimension parameters (previously in inc.p and inc.cls)
! to file inp0.unf
!---------------------------------------------------------------------
! missing: TRC (not used in kkr2)
! could be removed: KREL ?
subroutine write_dimension_parameters( &
     NAEZD, &
     LMAXD, &
     NREFD, &
     IRID, &
     BCPD, &
     NACLSD, &
     NCLEB, &
     IRMD, &
     IEXMD, &
     NGSHD, &
     IGUESSD, &
     IPAND, &
     ISHLD, &
     IRNSD, &
     KPOIBZ, &
     KREL, &
     NFUND, &
     NATRCD, &
     NCLSD, &
     NMAXD, &
     NRD, &
     NSPIND, &
     NUTRCD, &
     NXIJD, &
     LLY, &
     EKMD)

  implicit none

  integer, intent(in) :: NAEZD
  integer, intent(in) :: LMAXD
  integer, intent(in) :: NREFD
  integer, intent(in) :: IRID
  integer, intent(in) :: BCPD
  integer, intent(in) :: NACLSD
  integer, intent(in) :: NCLEB
  integer, intent(in) :: IRMD
  integer, intent(in) :: IEXMD
  integer, intent(in) :: NGSHD
  integer, intent(in) :: IGUESSD
  integer, intent(in) :: IPAND
  integer, intent(in) :: ISHLD
  integer, intent(in) :: IRNSD
  integer, intent(in) :: KPOIBZ
  integer, intent(in) :: KREL
  integer, intent(in) :: NFUND
  integer, intent(in) :: NATRCD
  integer, intent(in) :: NCLSD
  integer, intent(in) :: NMAXD
  integer, intent(in) :: NRD
  integer, intent(in) :: NSPIND
  integer, intent(in) :: NUTRCD
  integer, intent(in) :: NXIJD
  integer, intent(in) :: LLY
  integer, intent(in) :: EKMD

  integer, parameter :: FILEHANDLE = 67

  open (FILEHANDLE, FILE='inp0.unf', FORM='unformatted')

  write(FILEHANDLE) NAEZD
  write(FILEHANDLE) LMAXD
  write(FILEHANDLE) NREFD
  write(FILEHANDLE) IRID
  write(FILEHANDLE) BCPD
  write(FILEHANDLE) NACLSD
  write(FILEHANDLE) NCLEB
  write(FILEHANDLE) IRMD
  write(FILEHANDLE) IEXMD
  write(FILEHANDLE) NGSHD
  write(FILEHANDLE) IGUESSD
  write(FILEHANDLE) IPAND
  write(FILEHANDLE) ISHLD
  write(FILEHANDLE) IRNSD
  write(FILEHANDLE) KPOIBZ
  write(FILEHANDLE) KREL
  write(FILEHANDLE) NFUND
  write(FILEHANDLE) NATRCD
  write(FILEHANDLE) NCLSD
  write(FILEHANDLE) NMAXD
  write(FILEHANDLE) NRD
  write(FILEHANDLE) NSPIND
  write(FILEHANDLE) NUTRCD
  write(FILEHANDLE) NXIJD
  write(FILEHANDLE) LLY
  write(FILEHANDLE) EKMD

  close(FILEHANDLE)

end subroutine write_dimension_parameters
