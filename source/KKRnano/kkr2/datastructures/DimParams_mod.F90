!> Module that defines a datatype, which contains the parameters used to dimension important arrays
!> E.R.

! Dependencies: ConfigReader_mod

module DimParams_mod
implicit none
  private
  public :: DimParams, create, destroy
  public :: createDimParams, destroyDimParams ! deprecated
  public :: writeDimParams, createDimParamsFromFile

  type DimParams
    integer  :: NSYMAXD
    integer  :: NAEZ
    integer  :: LMAXD
    integer  :: IRID
    integer  :: BCPD
    integer  :: IRMD
    integer  :: IEMXD
    integer  :: IGUESSD
    integer  :: ISHLD
    integer  :: IRNSD
    integer  :: KPOIBZ
    integer  :: NMAXD
    integer  :: NSPIND
    integer  :: NXIJD
    integer  :: LLY
    integer  :: EKMD
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
    integer  :: LPOT
    integer  :: NGUESSD
    integer  :: SMPID
    integer  :: EMPID
    integer  :: NTHRDS
    integer  :: MAXMSHD

    integer :: num_atom_procs !< atom-parallelisation number of MPI-procs
    integer :: atoms_per_proc

  end type DimParams

  
  interface create
    module procedure createDimParams
  endinterface
  
  interface destroy
    module procedure destroyDimParams
  endinterface
  
  CONTAINS


  !-----------------------------------------------------------------------------
  !> Constructs a DimParams object from FORMATTED global.conf file
  !> @param[in,out] self    The DimParams object to construct.
  subroutine createDimParamsFromFile(self, filename)
    use ConfigReader_mod, only: ConfigReader, getValueInteger, getUnreadVariable, parseFile
    use ConfigReader_mod, only: createConfigReader, destroyConfigReader ! deprecated
    
    type (DimParams), intent(inout) :: self
    character(len=*), intent(in) :: filename ! usually 'global.conf'

    type (ConfigReader) :: conf
    integer :: ierror
    character(len=40) :: variable
    integer :: next_ptr

    ! these two parameters have to be determined later
    self%IEMXD = 0
    self%EKMD = 0

    call createConfigReader(conf)
    call parseFile(conf, filename, ierror)
    if (ierror /= 0) stop

    call getValueInteger(conf, "LMAXD", self%LMAXD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NSPIND", self%NSPIND, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NAEZD", self%NAEZ, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "IRNSD", self%IRNSD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "IRMD", self%IRMD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "IRID", self%IRID, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NXIJD", self%NXIJD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "KPOIBZ", self%KPOIBZ, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "IGUESSD", self%IGUESSD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "BCPD", self%BCPD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NMAXD", self%NMAXD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "ISHLD", self%ISHLD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "LLY", self%LLY, ierror)
    if (ierror /= 0) stop

    call getValueInteger(conf, "SMPID", self%SMPID, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "EMPID", self%EMPID, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NTHRDS", self%NTHRDS, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "XDIM", self%XDIM, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "YDIM", self%YDIM, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "ZDIM", self%ZDIM, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "NATBLD", self%NATBLD, ierror)
    if (ierror /= 0) stop
    call getValueInteger(conf, "ITDBRYD", self%ITDBRYD, ierror)
    if (ierror /= 0) stop

    call getValueInteger(conf, "num_atom_procs", self%num_atom_procs, ierror)
    if (ierror /= 0) then
      write(*,*) "WARNING: num_atom_procs not specified, using default = NAEZD."
      self%num_atom_procs = self%naez
    end if

    write(*,*) "The following variables have not been read from global.conf:"
    next_ptr = 1
    do
      call getUnreadVariable(conf, variable, next_ptr, ierror)
      if (ierror /= 0) exit
      write (*,*) variable
    enddo

    call destroyConfigReader(conf)

    ! deal with derived parameters
    call calculateDerivedParameters(self)

  endsubroutine

  
  !-----------------------------------------------------------------------------
  !> Constructs a DimParams object from UNFORMATTED inp0.unf file
  !> @param[in,out] self    The DimParams object to construct.
  subroutine createDimParams(self, filename)
    type (DimParams), intent(inout) :: self
    character(len=*), intent(in) :: filename ! usually 'inp0.unf' binary file

    integer, parameter :: FILEHANDLE = 67

    open (FILEHANDLE, FILE=filename, FORM='unformatted')

    read(FILEHANDLE) self%LMAXD
    read(FILEHANDLE) self%NSPIND
    read(FILEHANDLE) self%NAEZ
    read(FILEHANDLE) self%IRNSD
    read(FILEHANDLE) self%IRMD
    read(FILEHANDLE) self%IRID
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
    read(FILEHANDLE) self%num_atom_procs

    close(FILEHANDLE)

    ! deal with derived parameters
    call calculateDerivedParameters(self)

  endsubroutine
  
  !-----------------------------------------------------------------------------
  !> Write DimParams object to inp0.unf file
  !> @param[in,out] self    The DimParams object to construct.
  subroutine writeDimParams(self, filename)
    type (DimParams), intent(in) :: self
    character(len=*), intent(in) :: filename ! usually 'inp0.unf'

    integer, parameter :: FILEHANDLE = 67

    open (FILEHANDLE, FILE=filename, FORM='unformatted')

    write(FILEHANDLE) self%LMAXD
    write(FILEHANDLE) self%NSPIND
    write(FILEHANDLE) self%NAEZ
    write(FILEHANDLE) self%IRNSD
    write(FILEHANDLE) self%IRMD
    write(FILEHANDLE) self%IRID
    write(FILEHANDLE) self%NXIJD
    write(FILEHANDLE) self%KPOIBZ
    write(FILEHANDLE) self%IGUESSD
    write(FILEHANDLE) self%BCPD
    write(FILEHANDLE) self%NMAXD
    write(FILEHANDLE) self%ISHLD
    write(FILEHANDLE) self%LLY
    write(FILEHANDLE) self%SMPID
    write(FILEHANDLE) self%EMPID
    write(FILEHANDLE) self%NTHRDS
    write(FILEHANDLE) self%XDIM
    write(FILEHANDLE) self%YDIM
    write(FILEHANDLE) self%ZDIM
    write(FILEHANDLE) self%NATBLD
    write(FILEHANDLE) self%ITDBRYD
    write(FILEHANDLE) self%IEMXD
    write(FILEHANDLE) self%EKMD
    write(FILEHANDLE) self%num_atom_procs

    close(FILEHANDLE)

  end subroutine

  !-----------------------------------------------------------------------------
  !> Destroys a DimParams object.
  !> @param[in,out] self    The DimParams object to destroy.
  subroutine destroyDimParams(self)
    type (DimParams), intent(inout) :: self

!   integer :: memory_stat
    ! Nothing to do.
  endsubroutine

!------------------------------------------------------------------------------
!  HELPER ROUTINES
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !> Helper routine to initialise some parameters derived from others
  !> correctly.
  subroutine calculateDerivedParameters(self)
    type (DimParams), intent(inout) :: self

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
    !self%NTIRD=(self%IRMD+(self%IRNSD+1)*(self%LMPOTD-1))*self%NSPIND
    self%IRMIND=self%IRMD-self%IRNSD
    self%NGUESSD = 1 + self%IGUESSD * ( self%NAEZ * (self%LMAXD+1)**2 - 1 )

    ! Record lengths
    self%LRECRES2=4+8*(self%NSPIND*(self%LMAXD+7)+2*self%LPOT+4+2)

    ! Calculate atoms per process
    self%atoms_per_proc = self%naez / self%num_atom_procs

    call consistencyCheck01(self%IEMXD, self%LMAXD, self%NSPIND, self%SMPID)

  endsubroutine


  ! Consistency checks
  !----------------------------------------------------------------------------
  subroutine consistencyCheck01(IEMXD, LMAXD, NSPIND, SMPID)
    integer, intent(in) :: IEMXD
    integer, intent(in) :: LMAXD
    integer, intent(in) :: NSPIND
    integer, intent(in) :: SMPID

    ! -------------------------------------------------------------------------
    ! consistency checks
!    if (IEMXD < 1) then
!      write (*,*) "main2: IEMXD must be >= 1"
!      stop
!    end if

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
