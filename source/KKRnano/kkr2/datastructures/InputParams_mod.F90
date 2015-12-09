
!------------------------------------------------------------------------------
! Automatically generated source file. Do not edit by hand.
! To add/remove/modify input parameters:
! Edit InputParamsNew.txt and run 
! 'inputgenerator.py InputParams InputParamsNew.txt > source.F90'
! to generate source code.
!------------------------------------------------------------------------------


module InputParams_mod
  implicit none
  private
  public :: InputParams
  public :: getInputParamsValues
  public :: readInputParamsFromFile
  public :: writeInputParamsToFile

  type InputParams
    integer :: icst
    integer :: kpre
    double precision :: bravais_a (3)
    double precision :: bravais_b (3)
    double precision :: bravais_c (3)
    double precision :: gmax
    double precision :: rmax
    double precision :: qmrbound
    integer :: nsra
    integer :: kxc
    double precision :: tempr
    double precision :: alat
    double precision :: rclust
    double precision :: mixing
    double precision :: fcm
    integer :: scfsteps
    integer :: kforce
    double precision :: emin
    double precision :: emax
    logical :: jij
    logical :: ldau
    logical :: cartesian
    double precision :: rcutjij
    integer :: kte
    integer :: imix
    integer :: bzdivide (3)
    integer :: npol
    integer :: npnt1
    integer :: npnt2
    integer :: npnt3
    double precision :: rclust_voronoi
    integer :: nmin_panel
    integer :: num_MT_points
    double precision :: MT_scale
    double precision :: RMT_ref_scale
    integer :: use_semicore
    double precision :: ebotsemi
    double precision :: emusemi
    integer :: n1semi
    integer :: n2semi
    integer :: n3semi
    double precision :: fsemicore
    integer :: volterra
    double precision :: target_rms
    integer :: near_field
    integer :: write_shapes
    double precision :: mt_zero_shift
    integer :: DEBUG_morgan_electrostatics
  endtype ! InputParams

  contains
!-------------------------------------------------------------------------------
integer function getInputParamsValues(filename, values) result(ierror)
  use ConfigReader_mod, only: ConfigReader, createConfigReader, destroy
  use ConfigReader_mod, only: not_found => CONFIG_READER_ERR_VAR_NOT_FOUND
  use ConfigReader_mod, only: use_default => CONFIG_READER_USE_DEFAULT_VALUE
  use ConfigReader_mod, only: getValue, parseFile

  character(len=*), intent(in) :: filename

  type(InputParams), intent(inout) :: values
  type(ConfigReader) :: cr

  ierror = 0
  write(*,*) "Reading information from input.conf..."
  call createConfigReader(cr)
#define destroy_and_return   call destroy(cr) ; return
  ierror = parseFile(cr, filename)
  if (ierror /= 0) then
    write(*,*) "Error reading configfile ", trim(filename)
    destroy_and_return
  endif


  ierror = getValue(cr, "icst", values%icst , def=4)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for icst. Set to icst = 4"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for icst."
    destroy_and_return
  endif

  ierror = getValue(cr, "kpre", values%kpre , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for kpre. Set to kpre = 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for kpre."
    destroy_and_return
  endif

  ierror = getValue(cr, "bravais_a", values%bravais_a)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_a."
    destroy_and_return
  endif

  ierror = getValue(cr, "bravais_b", values%bravais_b)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_b."
    destroy_and_return
  endif

  ierror = getValue(cr, "bravais_c", values%bravais_c)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_c."
    destroy_and_return
  endif

  ierror = getValue(cr, "gmax", values%gmax)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for gmax."
    destroy_and_return
  endif

  ierror = getValue(cr, "rmax", values%rmax)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rmax."
    destroy_and_return
  endif

  ierror = getValue(cr, "qmrbound", values%qmrbound , def=1.0D-6)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for qmrbound. Set to qmrbound = 1.0D-6"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for qmrbound."
    destroy_and_return
  endif

  ierror = getValue(cr, "nsra", values%nsra , def=2)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for nsra. Set to nsra = 2"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for nsra."
    destroy_and_return
  endif

  ierror = getValue(cr, "kxc", values%kxc , def=2)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for kxc. Set to kxc = 2"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for kxc."
    destroy_and_return
  endif

  ierror = getValue(cr, "tempr", values%tempr , def=800.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for tempr. Set to tempr = 800.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for tempr."
    destroy_and_return
  endif

  ierror = getValue(cr, "alat", values%alat , def=1.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for alat. Set to alat = 1.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for alat."
    destroy_and_return
  endif

  ierror = getValue(cr, "rclust", values%rclust , def=1.5)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for rclust. Set to rclust = 1.5"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for rclust."
    destroy_and_return
  endif

  ierror = getValue(cr, "mixing", values%mixing)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for mixing."
    destroy_and_return
  endif

  ierror = getValue(cr, "fcm", values%fcm , def=20.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for fcm. Set to fcm = 20.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for fcm."
    destroy_and_return
  endif

  ierror = getValue(cr, "scfsteps", values%scfsteps)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for scfsteps."
    destroy_and_return
  endif

  ierror = getValue(cr, "kforce", values%kforce , def=1)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for kforce. Set to kforce = 1"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for kforce."
    destroy_and_return
  endif

  ierror = getValue(cr, "emin", values%emin)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for emin."
    destroy_and_return
  endif

  ierror = getValue(cr, "emax", values%emax)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for emax."
    destroy_and_return
  endif

  ierror = getValue(cr, "jij", values%jij , def=.FALSE.)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for jij. Set to jij = .FALSE."
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for jij."
    destroy_and_return
  endif

  ierror = getValue(cr, "ldau", values%ldau , def=.FALSE.)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for ldau. Set to ldau = .FALSE."
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for ldau."
    destroy_and_return
  endif

  ierror = getValue(cr, "cartesian", values%cartesian , def=.TRUE.)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for cartesian. Set to cartesian = .TRUE."
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for cartesian."
    destroy_and_return
  endif

  ierror = getValue(cr, "rcutjij", values%rcutjij , def=2.30)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for rcutjij. Set to rcutjij = 2.30"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for rcutjij."
    destroy_and_return
  endif

  ierror = getValue(cr, "kte", values%kte , def=1)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for kte. Set to kte = 1"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for kte."
    destroy_and_return
  endif

  ierror = getValue(cr, "imix", values%imix , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for imix. Set to imix = 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for imix."
    destroy_and_return
  endif

  ierror = getValue(cr, "bzdivide", values%bzdivide)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bzdivide."
    destroy_and_return
  endif

  ierror = getValue(cr, "npol", values%npol)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npol."
    destroy_and_return
  endif

  ierror = getValue(cr, "npnt1", values%npnt1)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt1."
    destroy_and_return
  endif

  ierror = getValue(cr, "npnt2", values%npnt2)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt2."
    destroy_and_return
  endif

  ierror = getValue(cr, "npnt3", values%npnt3)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt3."
    destroy_and_return
  endif

  ierror = getValue(cr, "rclust_voronoi", values%rclust_voronoi)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rclust_voronoi."
    destroy_and_return
  endif

  ierror = getValue(cr, "nmin_panel", values%nmin_panel , def=7)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for nmin_panel. Set to nmin_panel = 7"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for nmin_panel."
    destroy_and_return
  endif

  ierror = getValue(cr, "num_MT_points", values%num_MT_points , def=10)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for num_MT_points. Set to num_MT_points = 10"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for num_MT_points."
    destroy_and_return
  endif

  ierror = getValue(cr, "MT_scale", values%MT_scale , def=0.98)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for MT_scale. Set to MT_scale = 0.98"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for MT_scale."
    destroy_and_return
  endif

  ierror = getValue(cr, "RMT_ref_scale", values%RMT_ref_scale , def=0.995)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for RMT_ref_scale. Set to RMT_ref_scale = 0.995"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for RMT_ref_scale."
    destroy_and_return
  endif

  ierror = getValue(cr, "use_semicore", values%use_semicore , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for use_semicore. Set to use_semicore = 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for use_semicore."
    destroy_and_return
  endif

  ierror = getValue(cr, "ebotsemi", values%ebotsemi , def=0.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for ebotsemi. Set to ebotsemi = 0.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for ebotsemi."
    destroy_and_return
  endif

  ierror = getValue(cr, "emusemi", values%emusemi , def=0.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for emusemi. Set to emusemi = 0.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for emusemi."
    destroy_and_return
  endif

  ierror = getValue(cr, "n1semi", values%n1semi , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for n1semi. Set to n1semi = 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for n1semi."
    destroy_and_return
  endif

  ierror = getValue(cr, "n2semi", values%n2semi , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for n2semi. Set to n2semi = 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for n2semi."
    destroy_and_return
  endif

  ierror = getValue(cr, "n3semi", values%n3semi , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for n3semi. Set to n3semi = 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for n3semi."
    destroy_and_return
  endif

  ierror = getValue(cr, "fsemicore", values%fsemicore , def=1.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for fsemicore. Set to fsemicore = 1.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for fsemicore."
    destroy_and_return
  endif

  ierror = getValue(cr, "volterra", values%volterra , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for volterra. Set to volterra = 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for volterra."
    destroy_and_return
  endif

  ierror = getValue(cr, "target_rms", values%target_rms , def=1.0D-8)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for target_rms. Set to target_rms = 1.0D-8"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for target_rms."
    destroy_and_return
  endif

  ierror = getValue(cr, "near_field", values%near_field , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for near_field. Set to near_field = 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for near_field."
    destroy_and_return
  endif

  ierror = getValue(cr, "write_shapes", values%write_shapes , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for write_shapes. Set to write_shapes = 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for write_shapes."
    destroy_and_return
  endif

  ierror = getValue(cr, "mt_zero_shift", values%mt_zero_shift , def=0.0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for mt_zero_shift. Set to mt_zero_shift = 0.0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for mt_zero_shift."
    destroy_and_return
  endif

  ierror = getValue(cr, "DEBUG_morgan_electrostatics", values%DEBUG_morgan_electrostatics , def=0)
  if (ierror == use_default) then
    write(*,*) "WARNING: Bad/no value given for DEBUG_morgan_electrostatics. Set to DEBUG_morgan_electrostatics = 0"
    ierror = 0 ! ok, no error
  elseif (ierror /= 0) then
    write(*,*) "Bad/no value given for DEBUG_morgan_electrostatics."
    destroy_and_return
  endif

  write(*,*) "Finished reading information from input.conf"
  destroy_and_return
#undef destroy_and_return
endfunction !

!-------------------------------------------------------------------------------
integer function readInputParamsFromFile(values, filename) result(ierror)
  type(InputParams), intent(inout) :: values
  character(len=*), intent(in) :: filename

  integer, parameter :: fu = 67

  ierror = 0
  open(fu, file=filename, form="unformatted", action="read", status="old")
  read(fu) values%icst
  read(fu) values%kpre
  read(fu) values%bravais_a
  read(fu) values%bravais_b
  read(fu) values%bravais_c
  read(fu) values%gmax
  read(fu) values%rmax
  read(fu) values%qmrbound
  read(fu) values%nsra
  read(fu) values%kxc
  read(fu) values%tempr
  read(fu) values%alat
  read(fu) values%rclust
  read(fu) values%mixing
  read(fu) values%fcm
  read(fu) values%scfsteps
  read(fu) values%kforce
  read(fu) values%emin
  read(fu) values%emax
  read(fu) values%jij
  read(fu) values%ldau
  read(fu) values%cartesian
  read(fu) values%rcutjij
  read(fu) values%kte
  read(fu) values%imix
  read(fu) values%bzdivide
  read(fu) values%npol
  read(fu) values%npnt1
  read(fu) values%npnt2
  read(fu) values%npnt3
  read(fu) values%rclust_voronoi
  read(fu) values%nmin_panel
  read(fu) values%num_MT_points
  read(fu) values%MT_scale
  read(fu) values%RMT_ref_scale
  read(fu) values%use_semicore
  read(fu) values%ebotsemi
  read(fu) values%emusemi
  read(fu) values%n1semi
  read(fu) values%n2semi
  read(fu) values%n3semi
  read(fu) values%fsemicore
  read(fu) values%volterra
  read(fu) values%target_rms
  read(fu) values%near_field
  read(fu) values%write_shapes
  read(fu) values%mt_zero_shift
  read(fu) values%DEBUG_morgan_electrostatics
  close(fu)
endfunction ! readFromFile

!-------------------------------------------------------------------------------
integer function writeInputParamsToFile(values, filename) result(ierror)
  type(InputParams), intent(inout) :: values
  character(len=*), intent(in) :: filename

  integer, parameter :: fu = 67

  ierror = 0
  open(fu, file=filename, form="unformatted", action="write")
  write(fu) values%icst
  write(fu) values%kpre
  write(fu) values%bravais_a
  write(fu) values%bravais_b
  write(fu) values%bravais_c
  write(fu) values%gmax
  write(fu) values%rmax
  write(fu) values%qmrbound
  write(fu) values%nsra
  write(fu) values%kxc
  write(fu) values%tempr
  write(fu) values%alat
  write(fu) values%rclust
  write(fu) values%mixing
  write(fu) values%fcm
  write(fu) values%scfsteps
  write(fu) values%kforce
  write(fu) values%emin
  write(fu) values%emax
  write(fu) values%jij
  write(fu) values%ldau
  write(fu) values%cartesian
  write(fu) values%rcutjij
  write(fu) values%kte
  write(fu) values%imix
  write(fu) values%bzdivide
  write(fu) values%npol
  write(fu) values%npnt1
  write(fu) values%npnt2
  write(fu) values%npnt3
  write(fu) values%rclust_voronoi
  write(fu) values%nmin_panel
  write(fu) values%num_MT_points
  write(fu) values%MT_scale
  write(fu) values%RMT_ref_scale
  write(fu) values%use_semicore
  write(fu) values%ebotsemi
  write(fu) values%emusemi
  write(fu) values%n1semi
  write(fu) values%n2semi
  write(fu) values%n3semi
  write(fu) values%fsemicore
  write(fu) values%volterra
  write(fu) values%target_rms
  write(fu) values%near_field
  write(fu) values%write_shapes
  write(fu) values%mt_zero_shift
  write(fu) values%DEBUG_morgan_electrostatics
  close(fu)
endfunction ! writeToFile

endmodule !
