!> Module that defines a datatype, which contains the parameters used to dimension important arrays
!> E.R.

! Dependencies: ConfigReader_mod

module DimParams_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: DimParams, create, destroy
  public :: createDimParams, destroyDimParams ! deprecated
  public :: writeDimParams, createDimParamsFromFile

  type DimParams
    integer :: NSYMAXD
    integer :: NAEZ
    integer :: LMAXD
    integer :: IRID
    integer :: BCPD
    integer :: IRMD
    integer :: IEMXD
    integer :: IGUESSD
    integer :: ISHLD
    integer :: IRNSD
    integer :: KPOIBZ
    integer :: NMAXD
    integer :: NSPIND
    integer :: NXIJD
    integer :: LLY
    integer :: EKMD
    integer :: XDIM
    integer :: YDIM
    integer :: ZDIM
    integer :: NATBLD
    integer :: ITDBRYD
    integer :: LMMAXD
    integer :: LMAXD1
    integer :: MMAXD
    integer :: LMXSPD
    integer :: LMPOTD
    integer :: IRMIND
    integer :: LRECRES2
    integer :: LPOT
    integer :: NGUESSD
    integer :: SMPID
    integer :: EMPID
    integer :: NTHRDS
    integer :: MAXMSHD

    integer :: num_atom_procs !< atom-parallelisation number of MPI-procs
    integer :: atoms_per_proc

  endtype ! DimParams

  
  interface create
    module procedure createDimParams
  endinterface
  
  interface destroy
    module procedure destroyDimParams
  endinterface
  
  contains

  !-----------------------------------------------------------------------------
  !> Constructs a DimParams object from FORMATTED global.conf file
  !> @param[in,out] self    The DimParams object to construct.
  subroutine createDimParamsFromFile(self, filename)
    use ConfigReader_mod, only: ConfigReader, getValue, getUnreadVariable, parseFile
    use ConfigReader_mod, only: createConfigReader, destroyConfigReader ! deprecated
    
    type(DimParams), intent(inout) :: self
    character(len=*), intent(in) :: filename ! usually 'global.conf'

    type(ConfigReader) :: conf
    character(len=40) :: variable
    integer :: next_ptr

    ! these two parameters have to be determined later
    self%IEMXD = 0
    self%EKMD = 0

    call createConfigReader(conf)
    if (parseFile(conf, filename) /= 0) die_here("parsing file"+filename+"failed!")

    if (getValue(conf, "LMAXD",   self%LMAXD) /= 0)     die_here("did not find in LMAXD file"+filename)
    if (getValue(conf, "NSPIND",  self%NSPIND) /= 0)    die_here("did not find in NSPIND file"+filename)
    if (getValue(conf, "NAEZD",   self%NAEZ) /= 0)      die_here("did not find in NAEZD file"+filename)
    if (getValue(conf, "IRNSD",   self%IRNSD) /= 0)     die_here("did not find in IRNSD file"+filename)
    if (getValue(conf, "IRMD",    self%IRMD) /= 0)      die_here("did not find in IRMD file"+filename)
    if (getValue(conf, "IRID",    self%IRID) /= 0)      die_here("did not find in IRID file"+filename)
    if (getValue(conf, "NXIJD",   self%NXIJD) /= 0)     die_here("did not find in NXIJD file"+filename)
    if (getValue(conf, "KPOIBZ",  self%KPOIBZ) /= 0)    die_here("did not find in KPOIBZ file"+filename)
    if (getValue(conf, "IGUESSD", self%IGUESSD) /= 0)   die_here("did not find in IGUESSD file"+filename)
    if (getValue(conf, "BCPD",    self%BCPD) /= 0)      die_here("did not find in BCPD file"+filename)
    if (getValue(conf, "NMAXD",   self%NMAXD) /= 0)     die_here("did not find in NMAXD file"+filename)
    if (getValue(conf, "ISHLD",   self%ISHLD) /= 0)     die_here("did not find in ISHLD file"+filename)
    if (getValue(conf, "LLY",     self%LLY) /= 0)       die_here("did not find in LLY file"+filename)
    if (getValue(conf, "SMPID",   self%SMPID) /= 0)     die_here("did not find in SMPID file"+filename)
    if (getValue(conf, "EMPID",   self%EMPID) /= 0)     die_here("did not find in EMPID file"+filename)
    if (getValue(conf, "NTHRDS",  self%NTHRDS) /= 0)    die_here("did not find in NTHRDS file"+filename)
    if (getValue(conf, "XDIM",    self%XDIM) /= 0)      die_here("did not find in XDIM file"+filename)
    if (getValue(conf, "YDIM",    self%YDIM) /= 0)      die_here("did not find in YDIM file"+filename)
    if (getValue(conf, "ZDIM",    self%ZDIM) /= 0)      die_here("did not find in ZDIM file"+filename)
    if (getValue(conf, "NATBLD",  self%NATBLD) /= 0)    die_here("did not find in NATBLD file"+filename)
    if (getValue(conf, "ITDBRYD", self%ITDBRYD) /= 0)   die_here("did not find in ITDBRYD file"+filename)
    if (getValue(conf, "num_atom_procs", self%num_atom_procs) /= 0) then
      write(*,*) "WARNING: num_atom_procs not specified, using default = NAEZD."
      self%num_atom_procs = self%NAEZ ! default value
    endif

    write(*,*) "The following variables have not been read from global.conf:"
    next_ptr = 1
    do
      if (getUnreadVariable(conf, variable, next_ptr) /= 0) exit
      write (*,*) variable
    enddo

    call destroyConfigReader(conf)
    call calculateDerivedParameters(self) ! deal with derived parameters

  endsubroutine ! create

  
  !-----------------------------------------------------------------------------
  !> Constructs a DimParams object from UNFORMATTED inp0.unf file
  !> @param[in,out] self    The DimParams object to construct.
  subroutine createDimParams(self, filename)
    type(DimParams), intent(inout) :: self
    character(len=*), intent(in) :: filename ! usually 'inp0.unf' binary file

    integer, parameter :: fu = 67

    open (fu, FILE=filename, FORM='unformatted', action='read', status='old')

    read(fu) self%LMAXD
    read(fu) self%NSPIND
    read(fu) self%NAEZ
    read(fu) self%IRNSD
    read(fu) self%IRMD
    read(fu) self%IRID
    read(fu) self%NXIJD
    read(fu) self%KPOIBZ
    read(fu) self%IGUESSD
    read(fu) self%BCPD
    read(fu) self%NMAXD
    read(fu) self%ISHLD
    read(fu) self%LLY
    read(fu) self%SMPID
    read(fu) self%EMPID
    read(fu) self%NTHRDS
    read(fu) self%XDIM
    read(fu) self%YDIM
    read(fu) self%ZDIM
    read(fu) self%NATBLD
    read(fu) self%ITDBRYD
    read(fu) self%IEMXD
    read(fu) self%EKMD
    read(fu) self%num_atom_procs

    close(fu)

    ! deal with derived parameters
    call calculateDerivedParameters(self)

  endsubroutine ! create
  
  !-----------------------------------------------------------------------------
  !> Write DimParams object to inp0.unf file
  !> @param[in,out] self    The DimParams object to construct.
  subroutine writeDimParams(self, filename)
    type(DimParams), intent(in) :: self
    character(len=*), intent(in) :: filename ! usually 'inp0.unf'

    integer, parameter :: fu = 67

    open (fu, FILE=filename, FORM='unformatted', action='write')

    write(fu) self%LMAXD
    write(fu) self%NSPIND
    write(fu) self%NAEZ
    write(fu) self%IRNSD
    write(fu) self%IRMD
    write(fu) self%IRID
    write(fu) self%NXIJD
    write(fu) self%KPOIBZ
    write(fu) self%IGUESSD
    write(fu) self%BCPD
    write(fu) self%NMAXD
    write(fu) self%ISHLD
    write(fu) self%LLY
    write(fu) self%SMPID
    write(fu) self%EMPID
    write(fu) self%NTHRDS
    write(fu) self%XDIM
    write(fu) self%YDIM
    write(fu) self%ZDIM
    write(fu) self%NATBLD
    write(fu) self%ITDBRYD
    write(fu) self%IEMXD
    write(fu) self%EKMD
    write(fu) self%num_atom_procs

    close(fu)

  endsubroutine ! write

  !-----------------------------------------------------------------------------
  !> Destroys a DimParams object.
  !> @param[in,out] self    The DimParams object to destroy.
  subroutine destroyDimParams(self)
    type(DimParams), intent(inout) :: self

!   integer :: memory_stat
    ! Nothing to do.
  endsubroutine ! destroy

!------------------------------------------------------------------------------
!  HELPER ROUTINES
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !> Helper routine to initialise some parameters derived from others
  !> correctly.
  subroutine calculateDerivedParameters(self)
    type(DimParams), intent(inout) :: self

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

  endsubroutine ! calc


  ! Consistency checks
  !----------------------------------------------------------------------------
  subroutine consistencyCheck01(IEMXD, LMAXD, NSPIND, SMPID) ! todo: remove IEMXD
    integer, intent(in) :: IEMXD, LMAXD, NSPIND, SMPID

!   if (IEMXD < 1) stop "main2: IEMXD must be >= 1"

    if (LMAXD < 0) stop "main2: LMAXD must be >= 0"

    if (SMPID /= 1 .and. SMPID /=2) stop "main2: SMPID must be 1 or 2"

    if (NSPIND /= 1 .and. NSPIND /=2) stop "main2: NSPIND must be 1 or 2"

    if ((SMPID == 2) .and. (NSPIND /= 2)) stop "main2: Spin parallelism is only possible if NSPIND=2, set SMPID=1"
    
  endsubroutine ! check

endmodule DimParams_mod
