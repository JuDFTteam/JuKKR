!> Module that defines a datatype, which contains
!> the parameters used to dimension important arrays
!> E.R.

module DimParams_mod

  type DimParams
    integer  :: NSYMAXD
    integer  :: NAEZ
    integer  :: LMAXD
    integer  :: NREFD
    integer  :: IRID
    integer  :: BCPD
    integer  :: NACLSD
    integer  :: IRMD
    integer  :: IEMXD
    integer  :: IGUESSD
    integer  :: IPAND
    integer  :: ISHLD
    integer  :: IRNSD
    integer  :: KPOIBZ
    integer  :: NFUND
    integer  :: NCLSD
    integer  :: NMAXD
    integer  :: NRD
    integer  :: NSPIND
    integer  :: NXIJD
    integer  :: LLY
    integer  :: EKMD
    integer  :: NCELLD
    integer  :: XDIM
    integer  :: YDIM
    integer  :: ZDIM
    integer  :: NATBLD
    integer  :: ITDBRYD
    integer  :: LMMAXD
    integer  :: LMAXD1
    integer  :: MMAXD
    integer  :: LMXSPD
    integer  :: LMPOTD
    integer  :: IRMIND
    integer  :: LRECRES2
    integer  :: NTIRD
    integer  :: LPOT
    integer  :: NGUESSD
    integer  :: SMPID
    integer  :: EMPID
    integer  :: NTHRDS
    integer  :: MAXMSHD

    integer :: atoms_per_proc

  end type DimParams

  CONTAINS

  !-----------------------------------------------------------------------------
  !> Constructs a DimParams object from inp0.unf file
  !> @param[inout] self    The DimParams object to construct.
  subroutine createDimParams(self)
    implicit none
    type (DimParams), intent(inout) :: self

    integer, parameter :: FILEHANDLE = 67

    open (FILEHANDLE, FILE='inp0.unf', FORM='unformatted')

    read(FILEHANDLE) self%LMAXD
    read(FILEHANDLE) self%NSPIND
    read(FILEHANDLE) self%NAEZ
    read(FILEHANDLE) self%IRNSD
    read(FILEHANDLE) self%IRMD
    read(FILEHANDLE) self%NREFD
    read(FILEHANDLE) self%NRD
    read(FILEHANDLE) self%IRID
    read(FILEHANDLE) self%NFUND
    read(FILEHANDLE) self%NCELLD
    read(FILEHANDLE) self%NACLSD
    read(FILEHANDLE) self%NCLSD
    read(FILEHANDLE) self%IPAND
    read(FILEHANDLE) self%NXIJD
    read(FILEHANDLE) self%KPOIBZ
    read(FILEHANDLE) self%IGUESSD
    read(FILEHANDLE) self%BCPD
    read(FILEHANDLE) self%NMAXD
    read(FILEHANDLE) self%ISHLD
    read(FILEHANDLE) self%LLY
    read(FILEHANDLE) self%SMPID
    read(FILEHANDLE) self%EMPID
    read(FILEHANDLE) self%NTHRDS
    read(FILEHANDLE) self%XDIM
    read(FILEHANDLE) self%YDIM
    read(FILEHANDLE) self%ZDIM
    read(FILEHANDLE) self%NATBLD
    read(FILEHANDLE) self%ITDBRYD
    read(FILEHANDLE) self%IEMXD
    read(FILEHANDLE) self%EKMD

    close(FILEHANDLE)

    ! derived parameters
    self%MAXMSHD = 8
    self%NSYMAXD  = 48
    self%LPOT = 2*self%LMAXD

    ! derived dimension parameters
    self%LMMAXD= (self%LMAXD+1)**2

    self%LMAXD1=self%LMAXD+1
    self%MMAXD  = 2*self%LMAXD + 1
    self%LMXSPD= (2*self%LPOT+1)**2

    self%LMPOTD= (self%LPOT+1)**2
    self%NTIRD=(self%IRMD+(self%IRNSD+1)*(self%LMPOTD-1))*self%NSPIND
    self%IRMIND=self%IRMD-self%IRNSD

    self%NGUESSD = 1 + self%IGUESSD * ( self%NAEZ * (self%LMAXD+1)**2 - 1 )

    ! Record lengths
    self%LRECRES2=4+8*(self%NSPIND*(self%LMAXD+7)+2*self%LPOT+4+2)

    ! Only 1 atom per MPI process supported (for now)
    self%atoms_per_proc = 4

    call consistencyCheck01(self%IEMXD, self%LMAXD, self%NSPIND, self%SMPID)

  end subroutine

  !-----------------------------------------------------------------------------
  !> Destroys a DimParams object.
  !> @param[inout] self    The DimParams object to destroy.
  subroutine destroyDimParams(self)
    implicit none
    type (DimParams), intent(inout) :: self

    integer :: memory_stat

    continue ! Nothing to do.

  end subroutine
!-------------------------------------------------------------------------------
!> Returns NAEZ.
!> @param[in] self DimParams object
integer  function getNAEZ(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNaez = self%NAEZ
end function

!-------------------------------------------------------------------------------
!> Returns LMAXD.
!> @param[in] self DimParams object
integer  function getLMAXD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getLmaxd = self%LMAXD
end function

!-------------------------------------------------------------------------------
!> Returns NREFD.
!> @param[in] self DimParams object
integer  function getNREFD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNrefd = self%NREFD
end function

!-------------------------------------------------------------------------------
!> Returns IRID.
!> @param[in] self DimParams object
integer  function getIRID(self)
  implicit none
  type(DimParams), intent(in) :: self

  getIrid = self%IRID
end function

!-------------------------------------------------------------------------------
!> Returns BCPD.
!> @param[in] self DimParams object
integer  function getBCPD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getBcpd = self%BCPD
end function

!-------------------------------------------------------------------------------
!> Returns NACLSD.
!> @param[in] self DimParams object
integer  function getNACLSD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNaclsd = self%NACLSD
end function

!-------------------------------------------------------------------------------
!> Returns IRMD.
!> @param[in] self DimParams object
integer  function getIRMD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getIrmd = self%IRMD
end function

!-------------------------------------------------------------------------------
!> Returns IEMXD.
!> @param[in] self DimParams object
integer  function getIEMXD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getIemxd = self%IEMXD
end function

!-------------------------------------------------------------------------------
!> Returns IGUESSD.
!> @param[in] self DimParams object
integer  function getIGUESSD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getIguessd = self%IGUESSD
end function

!-------------------------------------------------------------------------------
!> Returns IPAND.
!> @param[in] self DimParams object
integer  function getIPAND(self)
  implicit none
  type(DimParams), intent(in) :: self

  getIpand = self%IPAND
end function

!-------------------------------------------------------------------------------
!> Returns ISHLD.
!> @param[in] self DimParams object
integer  function getISHLD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getIshld = self%ISHLD
end function

!-------------------------------------------------------------------------------
!> Returns IRNSD.
!> @param[in] self DimParams object
integer  function getIRNSD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getIrnsd = self%IRNSD
end function

!-------------------------------------------------------------------------------
!> Returns KPOIBZ.
!> @param[in] self DimParams object
integer  function getKPOIBZ(self)
  implicit none
  type(DimParams), intent(in) :: self

  getKpoibz = self%KPOIBZ
end function

!-------------------------------------------------------------------------------
!> Returns NFUND.
!> @param[in] self DimParams object
integer  function getNFUND(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNfund = self%NFUND
end function

!-------------------------------------------------------------------------------
!> Returns NCLSD.
!> @param[in] self DimParams object
integer  function getNCLSD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNclsd = self%NCLSD
end function

!-------------------------------------------------------------------------------
!> Returns NMAXD.
!> @param[in] self DimParams object
integer  function getNMAXD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNmaxd = self%NMAXD
end function

!-------------------------------------------------------------------------------
!> Returns NRD.
!> @param[in] self DimParams object
integer  function getNRD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNrd = self%NRD
end function

!-------------------------------------------------------------------------------
!> Returns NSPIND.
!> @param[in] self DimParams object
integer  function getNSPIND(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNspind = self%NSPIND
end function

!-------------------------------------------------------------------------------
!> Returns NXIJD.
!> @param[in] self DimParams object
integer  function getNXIJD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNxijd = self%NXIJD
end function

!-------------------------------------------------------------------------------
!> Returns LLY.
!> @param[in] self DimParams object
integer  function getLLY(self)
  implicit none
  type(DimParams), intent(in) :: self

  getLly = self%LLY
end function

!-------------------------------------------------------------------------------
!> Returns EKMD.
!> @param[in] self DimParams object
integer  function getEKMD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getEkmd = self%EKMD
end function

!-------------------------------------------------------------------------------
!> Returns NCELLD.
!> @param[in] self DimParams object
integer  function getNCELLD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNcelld = self%NCELLD
end function

!-------------------------------------------------------------------------------
!> Returns XDIM.
!> @param[in] self DimParams object
integer  function getXDIM(self)
  implicit none
  type(DimParams), intent(in) :: self

  getXdim = self%XDIM
end function

!-------------------------------------------------------------------------------
!> Returns YDIM.
!> @param[in] self DimParams object
integer  function getYDIM(self)
  implicit none
  type(DimParams), intent(in) :: self

  getYdim = self%YDIM
end function

!-------------------------------------------------------------------------------
!> Returns ZDIM.
!> @param[in] self DimParams object
integer  function getZDIM(self)
  implicit none
  type(DimParams), intent(in) :: self

  getZdim = self%ZDIM
end function

!-------------------------------------------------------------------------------
!> Returns NATBLD.
!> @param[in] self DimParams object
integer  function getNATBLD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNatbld = self%NATBLD
end function

!-------------------------------------------------------------------------------
!> Returns ITDBRYD.
!> @param[in] self DimParams object
integer  function getITDBRYD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getItdbryd = self%ITDBRYD
end function

!-------------------------------------------------------------------------------
!> Returns LMMAXD.
!> @param[in] self DimParams object
integer  function getLMMAXD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getLmmaxd = self%LMMAXD
end function

!-------------------------------------------------------------------------------
!> Returns LMAXD1.
!> @param[in] self DimParams object
integer  function getLMAXD1(self)
  implicit none
  type(DimParams), intent(in) :: self

  getLmaxd1 = self%LMAXD1
end function

!-------------------------------------------------------------------------------
!> Returns MMAXD.
!> @param[in] self DimParams object
integer  function getMMAXD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getMmaxd = self%MMAXD
end function

!-------------------------------------------------------------------------------
!> Returns LMXSPD.
!> @param[in] self DimParams object
integer  function getLMXSPD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getLmxspd = self%LMXSPD
end function

!-------------------------------------------------------------------------------
!> Returns LMPOTD.
!> @param[in] self DimParams object
integer  function getLMPOTD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getLmpotd = self%LMPOTD
end function

!-------------------------------------------------------------------------------
!> Returns IRMIND.
!> @param[in] self DimParams object
integer  function getIRMIND(self)
  implicit none
  type(DimParams), intent(in) :: self

  getIrmind = self%IRMIND
end function

!-------------------------------------------------------------------------------
!> Returns LRECRES2.
!> @param[in] self DimParams object
integer  function getLRECRES2(self)
  implicit none
  type(DimParams), intent(in) :: self

  getLrecres2 = self%LRECRES2
end function

!-------------------------------------------------------------------------------
!> Returns NTIRD.
!> @param[in] self DimParams object
integer  function getNTIRD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNtird = self%NTIRD
end function

!-------------------------------------------------------------------------------
!> Returns LPOT.
!> @param[in] self DimParams object
integer  function getLPOT(self)
  implicit none
  type(DimParams), intent(in) :: self

  getLpot = self%LPOT
end function

!-------------------------------------------------------------------------------
!> Returns NGUESSD.
!> @param[in] self DimParams object
integer  function getNGUESSD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNguessd = self%NGUESSD
end function

!-------------------------------------------------------------------------------
!> Returns SMPID.
!> @param[in] self DimParams object
integer  function getSMPID(self)
  implicit none
  type(DimParams), intent(in) :: self

  getSmpid = self%SMPID
end function

!-------------------------------------------------------------------------------
!> Returns EMPID.
!> @param[in] self DimParams object
integer  function getEMPID(self)
  implicit none
  type(DimParams), intent(in) :: self

  getEmpid = self%EMPID
end function

!-------------------------------------------------------------------------------
!> Returns NTHRDS.
!> @param[in] self DimParams object
integer  function getNTHRDS(self)
  implicit none
  type(DimParams), intent(in) :: self

  getNthrds = self%NTHRDS
end function

!-------------------------------------------------------------------------------
!> Returns MAXMESHD.
!> @param[in] self DimParams object
integer  function getMAXMSHD(self)
  implicit none
  type(DimParams), intent(in) :: self

  getMaxmshd = self%MAXMSHD
end function

  ! Consistency checks
  !----------------------------------------------------------------------------
  subroutine consistencyCheck01(IEMXD, LMAXD, NSPIND, SMPID)
    implicit none
    integer :: IEMXD
    integer :: LMAXD
    integer :: NSPIND
    integer :: SMPID

    ! -------------------------------------------------------------------------
    ! consistency checks
    if (IEMXD < 1) then
      write (*,*) "main2: IEMXD must be >= 1"
      stop
    end if

    if (LMAXD < 0) then
      write (*,*) "main2: LMAXD must be >= 0"
      stop
    end if

    if (SMPID /= 1 .and. SMPID /=2) then
      write (*,*) "main2: SMPID must be 1 or 2"
      stop
    end if

    if (NSPIND /= 1 .and. NSPIND /=2) then
      write (*,*) "main2: NSPIND must be 1 or 2"
      stop
    end if

    if ((SMPID == 2) .and. (NSPIND /= 2)) then
      write (*,*) "main2: Spin parallelism is only possible if NSPIND=2."
      write (*,*) "Set SMPID=1"
      stop
    end if
  end subroutine


end module
