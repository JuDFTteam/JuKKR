!---------------------------------------------------------------------
! Reads array dimension parameters (previously in inc.p and inc.cls)
! from file inp0.unf
!---------------------------------------------------------------------
! missing: TRC (not used in kkr2)
! could be removed: KREL ?
subroutine read_dimension_parameters( &
LMAXD, &
NSPIND, &
NAEZD, &
IRNSD, &
TRC, &
IRMD, &
NREFD, &
NRD, &
IRID, &
NFUND, &
NCELLD, &
NGSHD, &
NACLSD, &
NCLSD, &
IPAND, &
NXIJD, &
NATRCD, &
NUTRCD, &
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

  integer, intent(out) :: LMAXD
  integer, intent(out) :: NSPIND
  integer, intent(out) :: NAEZD
  integer, intent(out) :: IRNSD
  integer, intent(out) :: TRC
  integer, intent(out) :: IRMD
  integer, intent(out) :: NREFD
  integer, intent(out) :: NRD
  integer, intent(out) :: IRID
  integer, intent(out) :: NFUND
  integer, intent(out) :: NCELLD
  integer, intent(out) :: NGSHD
  integer, intent(out) :: NACLSD
  integer, intent(out) :: NCLSD
  integer, intent(out) :: IPAND
  integer, intent(out) :: NXIJD
  integer, intent(out) :: NATRCD
  integer, intent(out) :: NUTRCD
  integer, intent(out) :: KPOIBZ
  integer, intent(out) :: IGUESSD
  integer, intent(out) :: BCPD
  integer, intent(out) :: NMAXD
  integer, intent(out) :: ISHLD
  integer, intent(out) :: LLY
  integer, intent(out) :: SMPID
  integer, intent(out) :: EMPID
  integer, intent(out) :: NTHRDS
  integer, intent(out) :: XDIM
  integer, intent(out) :: YDIM
  integer, intent(out) :: ZDIM
  integer, intent(out) :: NATBLD
  integer, intent(out) :: ITDBRYD
  integer, intent(out) :: IEMXD
  integer, intent(out) :: EKMD

  integer :: FILEHANDLE = 67

  open (FILEHANDLE, FILE='inp0.unf', FORM='unformatted')

  read(FILEHANDLE) LMAXD
  read(FILEHANDLE) NSPIND
  read(FILEHANDLE) NAEZD
  read(FILEHANDLE) IRNSD
  read(FILEHANDLE) TRC
  read(FILEHANDLE) IRMD
  read(FILEHANDLE) NREFD
  read(FILEHANDLE) NRD
  read(FILEHANDLE) IRID
  read(FILEHANDLE) NFUND
  read(FILEHANDLE) NCELLD
  read(FILEHANDLE) NGSHD
  read(FILEHANDLE) NACLSD
  read(FILEHANDLE) NCLSD
  read(FILEHANDLE) IPAND
  read(FILEHANDLE) NXIJD
  read(FILEHANDLE) NATRCD
  read(FILEHANDLE) NUTRCD
  read(FILEHANDLE) KPOIBZ
  read(FILEHANDLE) IGUESSD
  read(FILEHANDLE) BCPD
  read(FILEHANDLE) NMAXD
  read(FILEHANDLE) ISHLD
  read(FILEHANDLE) LLY
  read(FILEHANDLE) SMPID
  read(FILEHANDLE) EMPID
  read(FILEHANDLE) NTHRDS
  read(FILEHANDLE) XDIM
  read(FILEHANDLE) YDIM
  read(FILEHANDLE) ZDIM
  read(FILEHANDLE) NATBLD
  read(FILEHANDLE) ITDBRYD
  read(FILEHANDLE) IEMXD
  read(FILEHANDLE) EKMD

  close(FILEHANDLE)

!  FILEHANDLE = 6
!  write(FILEHANDLE, *) "# -- KKRnano: global options --"
!  write(FILEHANDLE, *) "# Replacement for old inc.p and inc.cls files"
!  write(FILEHANDLE, *) "# Rerun kkr0.exe if a parameter of this file is changed."
!  write(FILEHANDLE, *)
!  write(FILEHANDLE, *) "#     general settings"
!  write(FILEHANDLE, *) "NAEZD  = ", NAEZD, " ! number of atoms in unit cell"
!  write(FILEHANDLE, *) "NSPIND = ", NSPIND
!  write(FILEHANDLE, *) "LMAXD  = ", LMAXD
!  write(FILEHANDLE, *) "LLY    = ", LLY, " ! LLY = 1(0), (in)active"
!  write(FILEHANDLE, *) "# ---------------------------------------------------------------------"
!  write(FILEHANDLE, *) "IRMD   = ", IRMD
!  write(FILEHANDLE, *) "IRNSD  = ", IRNSD
!  write(FILEHANDLE, *) "NRD    = ", NRD
!  write(FILEHANDLE, *) "KPOIBZ = ", KPOIBZ
!  write(FILEHANDLE, *) "NMAXD  = ", NMAXD
!  write(FILEHANDLE, *) "ISHLD  = ", ISHLD
!  write(FILEHANDLE, *)
!  write(FILEHANDLE, *) "# ---------------------------------------------------------------------"
!  write(FILEHANDLE, *) "#     structure-dependent"
!  write(FILEHANDLE, *) "NREFD  = ", NREFD
!  write(FILEHANDLE, *) "NXIJD  = ", NXIJD
!  write(FILEHANDLE, *) "TRC    = ", TRC, " ! truncation for order-N method"
!  write(FILEHANDLE, *) "NATRCD = ", NATRCD
!  write(FILEHANDLE, *) "NUTRCD = ", NUTRCD
!  write(FILEHANDLE, *)
!  write(FILEHANDLE, *) "# ---------------------------------------------------------------------"
!  write(FILEHANDLE, *) "#     Parallelization"
!  write(FILEHANDLE, *) "SMPID  = ", SMPID
!  write(FILEHANDLE, *) "EMPID  = ", EMPID
!  write(FILEHANDLE, *) "NTHRDS = ", NTHRDS
!  write(FILEHANDLE, *)
!  write(FILEHANDLE, *) "# ---------------------------------------------------------------------"
!  write(FILEHANDLE, *) "#     non-spherical potential, shape-functions"
!  write(FILEHANDLE, *) "NCELLD = ", NCELLD
!  write(FILEHANDLE, *) "IPAND  = ", IPAND
!  write(FILEHANDLE, *) "NFUND  = ", NFUND
!  write(FILEHANDLE, *) "IRID   = ", IRID
!  write(FILEHANDLE, *) "NGSHD  = ", NGSHD
!  write(FILEHANDLE, *)
!  write(FILEHANDLE, *) "# ---------------------------------------------------------------------"
!  write(FILEHANDLE, *) "#     mixing"
!  write(FILEHANDLE, *) "ITDBRYD = ", ITDBRYD, " ! Broyden history"
!  write(FILEHANDLE, *)
!  write(FILEHANDLE, *) "# ---------------------------------------------------------------------"
!  write(FILEHANDLE, *) "#     PRECONDITIONING"
!  write(FILEHANDLE, *) "BCPD    = ", BCPD
!  write(FILEHANDLE, *) "IGUESSD = ", IGUESSD
!  write(FILEHANDLE, *) "NATBLD  = ", NATBLD, " ! number of atoms per preconditioning block"
!  write(FILEHANDLE, *) "XDIM    = ", XDIM
!  write(FILEHANDLE, *) "YDIM    = ", YDIM
!  write(FILEHANDLE, *) "ZDIM    = ", ZDIM, " ! number of blocks in x,y,z direction"
!  write(FILEHANDLE, *)
!  write(FILEHANDLE, *) "# --------------------------------------------------------------------"
!  write(FILEHANDLE, *) "# Reference Cluster parameters (from inc.cls)"
!  write(FILEHANDLE, *) "NCLSD  = ", NCLSD
!  write(FILEHANDLE, *) "NACLSD = ", NACLSD
!
!  write(FILEHANDLE, *) "EKMD = ", EKMD
!  write(FILEHANDLE, *) "IEMXD = ", IEMXD

end subroutine read_dimension_parameters
