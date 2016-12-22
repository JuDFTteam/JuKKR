
!------------------------------------------------------------------------------
! Automatically generated source file. Do not edit manually!
! To add/remove/modify input parameters:
! Edit InputParamsNew.txt and run 
! 'inputgenerator.py InputParams InputParamsNew.txt > InputParams_mod.F90'
! to generate source code.
!------------------------------------------------------------------------------


module InputParams_mod
  implicit none
  private
  public :: InputParams, load, store, getValues


  type InputParams
    double precision :: alat
    double precision :: bravais_a (3)
    double precision :: bravais_b (3)
    double precision :: bravais_c (3)
    logical :: cartesian
    integer :: bzdivide (3)
    double precision :: rclust
    double precision :: emin
    double precision :: emax
    integer :: npol
    integer :: npnt1
    integer :: npnt2
    integer :: npnt3
    double precision :: tempr
    double precision :: ebotsemi
    double precision :: emusemi
    double precision :: fsemicore
    integer :: npntsemi (3)
    integer :: solver
    integer :: scfsteps
    integer :: imix
    double precision :: mixing
    double precision :: fcm
    double precision :: target_rms
    double precision :: gmax
    double precision :: rmax
    integer :: kxc
    double precision :: qmrbound
    integer :: volterra
    integer :: icst
    integer :: kpre
    integer :: kforce
    logical :: ldau
    integer :: nsra
    integer :: kte
    logical :: jij
    double precision :: rcutjij
    double precision :: rclust_voronoi
    integer :: nmin_panel
    integer :: num_MT_points
    double precision :: MT_scale
    double precision :: RMT_ref_scale
    double precision :: cutoff_radius
    double precision :: lcutoff_radii (9)
    integer :: near_field
    character(len=96) :: elementdatabasepath
    integer :: write_shapes
    double precision :: mt_zero_shift
    integer :: DEBUG_morgan_electrostatics
    logical :: fullbz
    double precision :: vref
    logical :: soc
    double precision :: socscale
    integer :: npan_log
    integer :: npan_eq
    integer :: ncheb
    double precision :: r_log
  endtype ! InputParams


  interface load
    module procedure readInputParamsFromFile
  endinterface

  interface store
    module procedure writeInputParamsToFile
  endinterface

  contains
!-------------------------------------------------------------------------------
integer function getValues(filename, self) result(ierror)
  use ConfigReader_mod, only: ConfigReader, create, destroy
  use ConfigReader_mod, only: not_found => CONFIG_READER_ERR_VAR_NOT_FOUND
  use ConfigReader_mod, only: use_default => CONFIG_READER_USE_DEFAULT_VALUE
  use ConfigReader_mod, only: getValue, parseFile

  character(len=*), intent(in) :: filename

  type(InputParams), intent(inout) :: self
  type(ConfigReader) :: cr

  ierror = 0
  write(*,*) "Reading information from input.conf..."
  call create(cr) ! createConfigReader
#define destroy_and_return   call destroy(cr) ; return
  ierror = parseFile(cr, filename)
  if (ierror /= 0) then
    write(*,*) "Error reading configfile ", trim(filename)
    destroy_and_return
  endif


  ierror = getValue(cr, "alat", self%alat , def=1.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for alat. Set alat to 1.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for alat."
    destroy_and_return
  endif

  ierror = getValue(cr, "bravais_a", self%bravais_a)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_a."
    destroy_and_return
  endif

  ierror = getValue(cr, "bravais_b", self%bravais_b)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_b."
    destroy_and_return
  endif

  ierror = getValue(cr, "bravais_c", self%bravais_c)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_c."
    destroy_and_return
  endif

  ierror = getValue(cr, "cartesian", self%cartesian , def=.TRUE.)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for cartesian. Set cartesian to .TRUE."
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for cartesian."
    destroy_and_return
  endif

  ierror = getValue(cr, "bzdivide", self%bzdivide , def=8)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for bzdivide. Set bzdivide to 8"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for bzdivide."
    destroy_and_return
  endif

  ierror = getValue(cr, "rclust", self%rclust , def=1.5)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for rclust. Set rclust to 1.5"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for rclust."
    destroy_and_return
  endif

  ierror = getValue(cr, "emin", self%emin)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for emin."
    destroy_and_return
  endif

  ierror = getValue(cr, "emax", self%emax)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for emax."
    destroy_and_return
  endif

  ierror = getValue(cr, "npol", self%npol)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npol."
    destroy_and_return
  endif

  ierror = getValue(cr, "npnt1", self%npnt1)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt1."
    destroy_and_return
  endif

  ierror = getValue(cr, "npnt2", self%npnt2)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt2."
    destroy_and_return
  endif

  ierror = getValue(cr, "npnt3", self%npnt3)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt3."
    destroy_and_return
  endif

  ierror = getValue(cr, "tempr", self%tempr , def=800.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for tempr. Set tempr to 800.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for tempr."
    destroy_and_return
  endif

  ierror = getValue(cr, "ebotsemi", self%ebotsemi , def=0.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for ebotsemi. Set ebotsemi to 0.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for ebotsemi."
    destroy_and_return
  endif

  ierror = getValue(cr, "emusemi", self%emusemi , def=0.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for emusemi. Set emusemi to 0.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for emusemi."
    destroy_and_return
  endif

  ierror = getValue(cr, "fsemicore", self%fsemicore , def=1.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for fsemicore. Set fsemicore to 1.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for fsemicore."
    destroy_and_return
  endif

  ierror = getValue(cr, "npntsemi", self%npntsemi , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for npntsemi. Set npntsemi to 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for npntsemi."
    destroy_and_return
  endif

  ierror = getValue(cr, "solver", self%solver , def=3)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for solver. Set solver to 3"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for solver."
    destroy_and_return
  endif

  ierror = getValue(cr, "scfsteps", self%scfsteps)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for scfsteps."
    destroy_and_return
  endif

  ierror = getValue(cr, "imix", self%imix , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for imix. Set imix to 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for imix."
    destroy_and_return
  endif

  ierror = getValue(cr, "mixing", self%mixing)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for mixing."
    destroy_and_return
  endif

  ierror = getValue(cr, "fcm", self%fcm , def=20.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for fcm. Set fcm to 20.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for fcm."
    destroy_and_return
  endif

  ierror = getValue(cr, "target_rms", self%target_rms , def=1.0D-8)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for target_rms. Set target_rms to 1.0D-8"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for target_rms."
    destroy_and_return
  endif

  ierror = getValue(cr, "gmax", self%gmax)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for gmax."
    destroy_and_return
  endif

  ierror = getValue(cr, "rmax", self%rmax)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rmax."
    destroy_and_return
  endif

  ierror = getValue(cr, "kxc", self%kxc , def=2)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for kxc. Set kxc to 2"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for kxc."
    destroy_and_return
  endif

  ierror = getValue(cr, "qmrbound", self%qmrbound , def=1.0D-9)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for qmrbound. Set qmrbound to 1.0D-9"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for qmrbound."
    destroy_and_return
  endif

  ierror = getValue(cr, "volterra", self%volterra , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for volterra. Set volterra to 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for volterra."
    destroy_and_return
  endif

  ierror = getValue(cr, "icst", self%icst , def=4)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for icst. Set icst to 4"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for icst."
    destroy_and_return
  endif

  ierror = getValue(cr, "kpre", self%kpre , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for kpre. Set kpre to 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for kpre."
    destroy_and_return
  endif

  ierror = getValue(cr, "kforce", self%kforce , def=1)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for kforce. Set kforce to 1"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for kforce."
    destroy_and_return
  endif

  ierror = getValue(cr, "ldau", self%ldau , def=.FALSE.)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for ldau. Set ldau to .FALSE."
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for ldau."
    destroy_and_return
  endif

  ierror = getValue(cr, "nsra", self%nsra , def=2)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for nsra. Set nsra to 2"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for nsra."
    destroy_and_return
  endif

  ierror = getValue(cr, "kte", self%kte , def=1)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for kte. Set kte to 1"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for kte."
    destroy_and_return
  endif

  ierror = getValue(cr, "jij", self%jij , def=.FALSE.)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for jij. Set jij to .FALSE."
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for jij."
    destroy_and_return
  endif

  ierror = getValue(cr, "rcutjij", self%rcutjij , def=2.30)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for rcutjij. Set rcutjij to 2.30"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for rcutjij."
    destroy_and_return
  endif

  ierror = getValue(cr, "rclust_voronoi", self%rclust_voronoi)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rclust_voronoi."
    destroy_and_return
  endif

  ierror = getValue(cr, "nmin_panel", self%nmin_panel , def=7)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for nmin_panel. Set nmin_panel to 7"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for nmin_panel."
    destroy_and_return
  endif

  ierror = getValue(cr, "num_MT_points", self%num_MT_points , def=10)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for num_MT_points. Set num_MT_points to 10"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for num_MT_points."
    destroy_and_return
  endif

  ierror = getValue(cr, "MT_scale", self%MT_scale , def=0.98)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for MT_scale. Set MT_scale to 0.98"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for MT_scale."
    destroy_and_return
  endif

  ierror = getValue(cr, "RMT_ref_scale", self%RMT_ref_scale , def=0.995)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for RMT_ref_scale. Set RMT_ref_scale to 0.995"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for RMT_ref_scale."
    destroy_and_return
  endif

  ierror = getValue(cr, "cutoff_radius", self%cutoff_radius , def=-1.d0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for cutoff_radius. Set cutoff_radius to -1.d0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for cutoff_radius."
    destroy_and_return
  endif

  ierror = getValue(cr, "lcutoff_radii", self%lcutoff_radii , def=0.d0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for lcutoff_radii. Set lcutoff_radii to 0.d0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for lcutoff_radii."
    destroy_and_return
  endif

  ierror = getValue(cr, "near_field", self%near_field , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for near_field. Set near_field to 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for near_field."
    destroy_and_return
  endif

  ierror = getValue(cr, "elementdatabasepath", self%elementdatabasepath , def='~/KKR/ElementDataBase')
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for elementdatabasepath. Set elementdatabasepath to '~/KKR/ElementDataBase'"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for elementdatabasepath."
    destroy_and_return
  endif

  ierror = getValue(cr, "write_shapes", self%write_shapes , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for write_shapes. Set write_shapes to 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for write_shapes."
    destroy_and_return
  endif

  ierror = getValue(cr, "mt_zero_shift", self%mt_zero_shift , def=0.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for mt_zero_shift. Set mt_zero_shift to 0.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for mt_zero_shift."
    destroy_and_return
  endif

  ierror = getValue(cr, "DEBUG_morgan_electrostatics", self%DEBUG_morgan_electrostatics , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for DEBUG_morgan_electrostatics. Set DEBUG_morgan_electrostatics to 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for DEBUG_morgan_electrostatics."
    destroy_and_return
  endif

  ierror = getValue(cr, "fullbz", self%fullbz , def=.FALSE.)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for fullbz. Set fullbz to .FALSE."
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for fullbz."
    destroy_and_return
  endif

  ierror = getValue(cr, "vref", self%vref , def=8.d0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for vref. Set vref to 8.d0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for vref."
    destroy_and_return
  endif

  ierror = getValue(cr, "soc", self%soc , def=.FALSE.)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for soc. Set soc to .FALSE."
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for soc."
    destroy_and_return
  endif

  ierror = getValue(cr, "socscale", self%socscale , def=1.0D0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for socscale. Set socscale to 1.0D0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for socscale."
    destroy_and_return
  endif

  ierror = getValue(cr, "npan_log", self%npan_log , def=30)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for npan_log. Set npan_log to 30"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for npan_log."
    destroy_and_return
  endif

  ierror = getValue(cr, "npan_eq", self%npan_eq , def=30)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for npan_eq. Set npan_eq to 30"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for npan_eq."
    destroy_and_return
  endif

  ierror = getValue(cr, "ncheb", self%ncheb , def=10)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for ncheb. Set ncheb to 10"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for ncheb."
    destroy_and_return
  endif

  ierror = getValue(cr, "r_log", self%r_log , def=0.1D0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for r_log. Set r_log to 0.1D0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for r_log."
    destroy_and_return
  endif

  write(*,*) "Finished reading information from input.conf"
  destroy_and_return
#undef destroy_and_return
endfunction ! get

!-------------------------------------------------------------------------------
subroutine readInputParamsFromFile(self, filename, ios)
  type(InputParams), intent(inout) :: self
  character(len=*), intent(in) :: filename
  integer, intent(out) :: ios

  integer, parameter :: fu = 67
  integer :: ioc

  open(fu, file=filename, form="unformatted", action="read", status="old", iostat=ios)
  if (ios /= 0) return
  read(fu, iostat=ios) self
  close(fu, iostat=ioc)
endsubroutine ! load

!-------------------------------------------------------------------------------
subroutine writeInputParamsToFile(self, filename)
  type(InputParams), intent(inout) :: self
  character(len=*), intent(in) :: filename

  integer, parameter :: fu = 67

  open(fu, file=filename, form="unformatted", action="write")
  write(fu) self
  close(fu)
endsubroutine ! store

endmodule ! InputParams
