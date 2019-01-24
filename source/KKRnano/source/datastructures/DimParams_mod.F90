!> Module that defines a datatype, which contains the parameters used to dimension important arrays
!> E.R.

! Dependencies: ConfigReader_mod

module DimParams_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: DimParams, parse, load, store, destroy

  type DimParams
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
    integer :: korbit !NOCO
    integer :: lmmaxd_noco !NOCO, lmmaxd_noco = (korbit+1)*lmmaxd
  endtype ! DimParams

  
  interface load
    module procedure createDimParamsFromDisk
  endinterface

  interface parse
    module procedure createDimParamsFromFile
  endinterface

  interface store
    module procedure writeDimParams
  endinterface
  
  interface destroy
    module procedure destroyDimParams
  endinterface
  
  contains

  !-----------------------------------------------------------------------------
  !> Constructs a DimParams object from FORMATTED global.conf file
  !> @param[in,out] self    The DimParams object to construct.
  subroutine createDimParamsFromFile(self, filename, altfile)
    use ConfigReader_mod, only: ConfigReader, create, destroy
    use ConfigReader_mod, only: getValue, getUnreadVariable, parseFile
    
    type(DimParams), intent(inout) :: self
    character(len=*), intent(in) :: filename          ! usually 'global.conf' and 'input.conf' as alternative
    character(len=*), intent(in), optional :: altfile ! usually 'global.conf' and 'input.conf' as alternative

    type(ConfigReader) :: cr
    character(len=40) :: variable
    integer :: next_ptr
    integer, parameter :: AUTO = 0 ! will be derived from other quantities

    ! these two parameters have to be determined later
    self%IEMXD = 0
    self%EKMD = 0

    call create(cr)
    
    if (parseFile(cr, filename) /= 0) then
      if (present(altfile)) then
        warn(6, "parsing file"+filename+"failed, try to parse"+altfile-"!")
        if (parseFile(cr, altfile) /= 0) then
          warn(6, "parsing file"+altfile+"failed as well, assume default values!")
        endif
      else  ! present altfile
        warn(6, "parsing file"+filename+"failed, assume default values!")
      endif ! present altfile
    endif

    ! all parameters are optional, getValue will return positive if an error occured
    if (getValue(cr, "NAEZD",   self%naez, def=AUTO) > 0)   die_here("unable to read NAEZD in file"+filename)
    if (getValue(cr, "IRNSD",   self%irnsd, def=208) > 0)   die_here("unable to read IRNSD in file"+filename)
    if (getValue(cr, "IRMD",    self%irmd, def=484) > 0)    die_here("unable to read IRMD in file"+filename)
    if (getValue(cr, "IRID",    self%irid, def=135) > 0)    die_here("unable to read IRID in file"+filename)

    if (getValue(cr, "ITDBRYD", self%itdbryd, def=30) > 0)  die_here("unable to read ITDBRYD in file"+filename) ! this can be moved to input params
    if (getValue(cr, "IGUESSD", self%iguessd, def=0) > 0)   die_here("unable to read IGUESSD in file"+filename) ! this var has acutally become a switch (0 or 1)
    
    if (getValue(cr, "KPOIBZ",  self%kpoibz, def=1024) > 0) die_here("unable to read KPOIBZ in file"+filename) ! todo: switch to dynamic allocation for kpoints
    if (getValue(cr, "NXIJD",   self%nxijd, def=1) > 0)     die_here("unable to read NXIJD in file"+filename) ! todo: switch to dynamic allocation in Jij module
    if (getValue(cr, "LMAXD",   self%lmaxd, def=3) > 0)     die_here("unable to read LMAXD in file"+filename)
    if (getValue(cr, "NSPIND",  self%nspind, def=1) > 0)    die_here("unable to read NSPIND in file"+filename)
    if (getValue(cr, "LLY",     self%lly, def=0) > 0)       die_here("unable to read LLY in file"+filename)
    if (getValue(cr, "SMPID",   self%smpid, def=1) > 0)     die_here("unable to read SMPID in file"+filename)
    if (getValue(cr, "EMPID",   self%empid, def=1) > 0)     die_here("unable to read EMPID in file"+filename)
    if (getValue(cr, "NTHRDS",  self%nthrds, def=AUTO) > 0) die_here("unable to read NTHRDS in file"+filename) ! AUTO: use environment variable OMP_NUM_THREADS
    if (getValue(cr, "BCPD",    self%bcpd, def=0) > 0)      die_here("unable to read BCPD in file"+filename)
    if (getValue(cr, "NATBLD",  self%natbld, def=4) > 0)    die_here("unable to read NATBLD in file"+filename)
    if (getValue(cr, "XDIM",    self%xdim, def=1) > 0)      die_here("unable to read XDIM in file"+filename)
    if (getValue(cr, "YDIM",    self%ydim, def=1) > 0)      die_here("unable to read YDIM in file"+filename)
    if (getValue(cr, "ZDIM",    self%zdim, def=1) > 0)      die_here("unable to read ZDIM in file"+filename)

    if (getValue(cr, "KORBIT",  self%korbit, def=0) > 0)    die_here("unable to read KORBIT in file"+filename)
    ! new default 0: automatically adopt to the number of currently running MPI processes
    if (getValue(cr, "num_atom_procs", self%num_atom_procs, def=AUTO) > 0) die_here("num_atom_procs could not be parsed in file"+filename)

    write(*,'(9a)') " The following variables have not been read from ",trim(filename),":"
    next_ptr = 1
    do while (getUnreadVariable(cr, variable, next_ptr) == 0)
      write(*,*) trim(variable)
    enddo ! while
    
    call destroy(cr)
    call calculateDerivedParameters(self) ! deal with derived parameters

  endsubroutine ! create

  
  !-----------------------------------------------------------------------------
  !> Constructs a DimParams object from UNFORMATTED inp0.unf file
  !> @param[in,out] self    The DimParams object to construct.
  subroutine createDimParamsFromDisk(self, filename)
    type(DimParams), intent(inout) :: self
    character(len=*), intent(in) :: filename ! usually 'inp0.unf' binary file

    integer, parameter :: fu = 67

    open(fu, file=filename, form='unformatted', action='read', status='old')
    read(fu) self
    close(fu)

    call calculateDerivedParameters(self) ! deal with derived parameters

  endsubroutine ! create
  
  !-----------------------------------------------------------------------------
  !> Write DimParams object to inp0.unf file
  !> @param[in,out] self    The DimParams object to construct.
  subroutine writeDimParams(self, filename)
    type(DimParams), intent(in) :: self
    character(len=*), intent(in) :: filename ! usually 'inp0.unf'

    integer, parameter :: fu = 67

    open(fu, file=filename, form='unformatted', action='write')
    write(fu) self
    close(fu)

  endsubroutine ! write

  !-----------------------------------------------------------------------------
  !> Destroys a DimParams object.
  !> @param[in,out] self    The DimParams object to destroy.
  elemental subroutine destroyDimParams(self)
    type(DimParams), intent(inout) :: self
    ! Nothing to do.
  endsubroutine ! destroy


  !----------------------------------------------------------------------------
  !> Helper routine to initialise some parameters derived from others correctly.
  subroutine calculateDerivedParameters(self)
    type(DimParams), intent(inout) :: self

    ! derived parameters
    self%maxmshd = 8
    self%lpot = 2*self%lmaxd

    ! derived dimension parameters
    self%lmmaxd = (self%lmaxd+1)**2
    self%lmmaxd_noco = (1+self%korbit)*(self%lmaxd+1)**2 !NOCO, matrix sizes msut be doubled if korbit==1   

    self%lmaxd1 = self%lmaxd+1
    self%mmaxd  = 2*self%lmaxd + 1
    self%lmxspd = (2*self%lpot+1)**2

    self%lmpotd = (self%lpot+1)**2
    !self%ntird = (self%irmd+(self%irnsd+1)*(self%lmpotd-1))*self%nspind
    self%irmind = self%irmd - self%irnsd
    self%nguessd = 1 + self%iguessd * ( self%naez * (self%lmaxd+1)**2 - 1 )

    ! record lengths
    self%lrecres2 = 4+8*(self%nspind*(self%lmaxd+7)+2*self%lpot+4+2)

    call consistencyCheck01(self%lmaxd, self%nspind, self%smpid)

  endsubroutine ! calc


  ! Consistency checks
  !----------------------------------------------------------------------------
  subroutine consistencyCheck01(lmaxd, nspind, smpid)
    integer, intent(in) :: lmaxd, nspind, smpid

    if (lmaxd < 0) stop "main2: LMAXD must be >= 0"

    if (smpid /= 1 .and. smpid /= 2) stop "main2: SMPID must be 1 or 2"
    if (nspind /= 1 .and. nspind /= 2) stop "main2: NSPIND must be 1 or 2"
    if (smpid == 2 .and. nspind /= 2) stop "main2: Spin parallelism is only possible if NSPIND=2, set SMPID=1"
    
  endsubroutine ! check

endmodule ! DimParams_mod
