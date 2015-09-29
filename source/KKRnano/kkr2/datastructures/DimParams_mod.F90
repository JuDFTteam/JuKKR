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
    integer :: nsymaxd
    integer :: naez
    integer :: lmaxd
    integer :: irid
    integer :: bcpd
    integer :: irmd
    integer :: iemxd
    integer :: iguessd
    integer :: irnsd
    integer :: kpoibz
    integer :: nspind
    integer :: nxijd
    integer :: lly
    integer :: ekmd
    integer :: xdim
    integer :: ydim
    integer :: zdim
    integer :: natbld
    integer :: itdbryd
    integer :: lmmaxd
    integer :: lmaxd1
    integer :: mmaxd
    integer :: lmxspd
    integer :: lmpotd
    integer :: irmind
    integer :: lrecres2
    integer :: lpot
    integer :: nguessd
    integer :: smpid
    integer :: empid
    integer :: nthrds
    integer :: maxmshd

    integer :: num_atom_procs !< atom-parallelisation number of mpi-procs

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
    use ConfigReader_mod, only: createConfigReader, destroyConfigReader ! deprecated, new: destroy, create
    
    type(DimParams), intent(inout) :: self
    character(len=*), intent(in) :: filename ! usually 'global.conf'

    type(ConfigReader) :: cr
    character(len=40) :: variable
    integer :: next_ptr

    ! these two parameters have to be determined later
    self%IEMXD = 0
    self%EKMD = 0

    call createConfigReader(cr)
    if (parseFile(cr, filename) /= 0) die_here("parsing file"+filename+"failed!")

    ! these parameters are not optional (so far)
    if (getValue(cr, "NAEZD",   self%naez) /= 0)      die_here("did not find in NAEZD file"+filename)
    if (getValue(cr, "IRNSD",   self%irnsd) /= 0)     die_here("did not find in IRNSD file"+filename)
    if (getValue(cr, "IRMD",    self%irmd) /= 0)      die_here("did not find in IRMD file"+filename)
    if (getValue(cr, "IRID",    self%irid) /= 0)      die_here("did not find in IRID file"+filename)
    ! optionals
    if (getValue(cr, "ITDBRYD", self%itdbryd, def=30) > 0)  die_here("did not find in ITDBRYD file"+filename) ! this can be moved to input params
    if (getValue(cr, "IGUESSD", self%iguessd, def=0) > 0)   die_here("did not find in IGUESSD file"+filename) ! this var has acutally become a switch (0 or 1)
    
    if (getValue(cr, "KPOIBZ",  self%kpoibz, def=1024) > 0) die_here("did not find in KPOIBZ file"+filename) ! todo: switch to dynamic allocation for kpoints
    if (getValue(cr, "NXIJD",   self%nxijd, def=1) > 0)     die_here("did not find in NXIJD file"+filename) ! todo: switch to dynamic allocation in Jij module
    if (getValue(cr, "LMAXD",   self%lmaxd, def=3) > 0)     die_here("did not find in LMAXD file"+filename)
    if (getValue(cr, "NSPIND",  self%nspind, def=1) > 0)    die_here("did not find in NSPIND file"+filename)
    if (getValue(cr, "LLY",     self%lly, def=0) > 0)       die_here("did not find in LLY file"+filename)
    if (getValue(cr, "SMPID",   self%smpid, def=1) > 0)     die_here("did not find in SMPID file"+filename)
    if (getValue(cr, "EMPID",   self%empid, def=1) > 0)     die_here("did not find in EMPID file"+filename)
    if (getValue(cr, "NTHRDS",  self%nthrds, def=0) > 0)    die_here("did not find in NTHRDS file"+filename) ! 0:automatically use environment variable OMP_NUM_THREADS
    if (getValue(cr, "BCPD",    self%bcpd, def=0) > 0)      die_here("did not find in BCPD file"+filename)
    if (getValue(cr, "NATBLD",  self%natbld, def=4) > 0)    die_here("did not find in NATBLD file"+filename)
    if (getValue(cr, "XDIM",    self%xdim, def=1) > 0)      die_here("did not find in XDIM file"+filename)
    if (getValue(cr, "YDIM",    self%ydim, def=1) > 0)      die_here("did not find in YDIM file"+filename)
    if (getValue(cr, "ZDIM",    self%zdim, def=1) > 0)      die_here("did not find in ZDIM file"+filename)

    ! new default 0: automatically adopt to the number of currently running MPI processes
    if (getValue(cr, "num_atom_procs", self%num_atom_procs, def=0) > 0) die_here("num_atom_procs could not be parsed in file"+filename) ! 0:auto

    write(*,'(9a)') " The following variables have not been read from ",trim(filename),":"
    next_ptr = 1
    do while (getUnreadVariable(cr, variable, next_ptr) == 0)
      write(*,*) trim(variable)
    enddo ! while
    
    call destroyConfigReader(cr)
    call calculateDerivedParameters(self) ! deal with derived parameters

  endsubroutine ! create

  
  !-----------------------------------------------------------------------------
  !> Constructs a DimParams object from UNFORMATTED inp0.unf file
  !> @param[in,out] self    The DimParams object to construct.
  subroutine createDimParams(self, filename)
    type(DimParams), intent(inout) :: self
    character(len=*), intent(in) :: filename ! usually 'inp0.unf' binary file

    integer, parameter :: fu = 67

    open (fu, file=filename, FORM='unformatted', action='read', status='old')

    read (fu) self%lmaxd
    read (fu) self%nspind
    read (fu) self%naez
    read (fu) self%irnsd
    read (fu) self%irmd
    read (fu) self%irid
    read (fu) self%nxijd
    read (fu) self%kpoibz
    read (fu) self%iguessd
    read (fu) self%bcpd
    read (fu) self%lly
    read (fu) self%smpid
    read (fu) self%empid
    read (fu) self%nthrds
    read (fu) self%xdim
    read (fu) self%ydim
    read (fu) self%zdim
    read (fu) self%natbld
    read (fu) self%itdbryd
    read (fu) self%iemxd
    read (fu) self%ekmd
    read (fu) self%num_atom_procs

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

    open (fu, file=filename, form='unformatted', action='write')

    write(fu) self%lmaxd
    write(fu) self%nspind
    write(fu) self%naez
    write(fu) self%irnsd
    write(fu) self%irmd
    write(fu) self%irid
    write(fu) self%nxijd
    write(fu) self%kpoibz
    write(fu) self%iguessd
    write(fu) self%bcpd
    write(fu) self%lly
    write(fu) self%smpid
    write(fu) self%empid
    write(fu) self%nthrds
    write(fu) self%xdim
    write(fu) self%ydim
    write(fu) self%zdim
    write(fu) self%natbld
    write(fu) self%itdbryd
    write(fu) self%iemxd
    write(fu) self%ekmd
    write(fu) self%num_atom_procs

    close(fu)

  endsubroutine ! write

  !-----------------------------------------------------------------------------
  !> Destroys a DimParams object.
  !> @param[in,out] self    The DimParams object to destroy.
  subroutine destroyDimParams(self)
    type(DimParams), intent(inout) :: self
    ! Nothing to do.
  endsubroutine ! destroy


  !----------------------------------------------------------------------------
  !> Helper routine to initialise some parameters derived from others correctly.
  subroutine calculateDerivedParameters(self)
    type(DimParams), intent(inout) :: self

    ! derived parameters
    self%maxmshd = 8
    self%nsymaxd = 48
    self%lpot = 2*self%lmaxd

    ! derived dimension parameters
    self%lmmaxd = (self%lmaxd+1)**2

    self%lmaxd1 = self%lmaxd+1
    self%mmaxd  = 2*self%lmaxd + 1
    self%lmxspd = (2*self%lpot+1)**2

    self%lmpotd = (self%lpot+1)**2
    !self%ntird = (self%irmd+(self%irnsd+1)*(self%lmpotd-1))*self%nspind
    self%irmind = self%irmd - self%irnsd
    self%nguessd = 1 + self%iguessd * ( self%naez * (self%lmaxd+1)**2 - 1 )

    ! record lengths
    self%lrecres2 = 4+8*(self%nspind*(self%lmaxd+7)+2*self%lpot+4+2)

    call consistencyCheck01(self%iemxd, self%lmaxd, self%nspind, self%smpid)

  endsubroutine ! calc


  ! Consistency checks
  !----------------------------------------------------------------------------
  subroutine consistencyCheck01(iemxd, lmaxd, nspind, smpid) ! todo: remove iemxd
    integer, intent(in) :: iemxd, lmaxd, nspind, smpid

!   if (iemxd < 1) stop "main2: IEMXD must be >= 1"

    if (lmaxd < 0) stop "main2: LMAXD must be >= 0"

    if (smpid /= 1 .and. smpid /=2) stop "main2: SMPID must be 1 or 2"

    if (nspind /= 1 .and. nspind /=2) stop "main2: NSPIND must be 1 or 2"

    if ((smpid == 2) .and. (nspind /= 2)) stop "main2: Spin parallelism is only possible if NSPIND=2, set SMPID=1"
    
  endsubroutine ! check

endmodule DimParams_mod
