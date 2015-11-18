
!------------------------------------------------------------------------------
! Automatically generated source file. Do not edit by hand.
! To add/remove/modify input parameters:
! Edit InputParamsNew.txt and run 
! 'inputgenerator.py InputParams InputParamsNew.txt > source.f90'
! to generate source code.
!------------------------------------------------------------------------------


module InputParams_mod
type InputParams
  integer :: fred
  integer :: icst
  integer :: kpre
  double precision :: bravais_b (3)
  double precision :: bravais_c (3)
  double precision :: gmax
  double precision :: qmrbound
  integer :: nsra
  integer :: kxc
  double precision :: rmax
  double precision :: tempr
  double precision :: alat
  double precision :: rclust
  double precision :: mixing
  double precision :: fcm
  integer :: scfsteps
  integer :: kforce
  double precision :: emin
  double precision :: bravais_a (3)
  logical :: jij
  logical :: ldau
  double precision :: rcutjij
  double precision :: emax
  integer :: kte
  integer :: imix
  integer :: bzdivide (3)
  logical :: cartesian
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
  double precision :: target_rms
  integer :: near_field
  integer :: write_shapes
  double precision :: mt_zero_shift
  integer :: DEBUG_morgan_electrostatics
end type InputParams

CONTAINS
!-------------------------------------------------------------------------------
integer function getInputParamsValues(filename, confvalues) result(ierror)
  use Config_Reader
  implicit none

  character(len=*) :: filename
  type (ConfigReader) :: conf
  type (InputParams), intent(inout) :: confvalues

  ierror = 0
  write(*,*) "Reading information from input.conf..."
  call createConfigReader(conf)
  call parseFile(conf, filename, ierror)
  if (ierror /= 0) then
    write(*,*) "Error reading configfile", filename
    call destroyConfigReader(conf)
    return
  end if


  call getValueInteger(conf, "fred", confvalues%fred, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%fred = 1
    write(*,*) "WARNING: Bad/no value given for fred. Set to fred = 1"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for fred."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "icst", confvalues%icst, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%icst = 4
    write(*,*) "WARNING: Bad/no value given for icst. Set to icst = 4"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for icst."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "kpre", confvalues%kpre, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%kpre = 0
    write(*,*) "WARNING: Bad/no value given for kpre. Set to kpre = 0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for kpre."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDoubleVector(conf, "bravais_b", confvalues%bravais_b, 3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_b."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDoubleVector(conf, "bravais_c", confvalues%bravais_c, 3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_c."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "gmax", confvalues%gmax, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for gmax."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "qmrbound", confvalues%qmrbound, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%qmrbound = 1.0D-6
    write(*,*) "WARNING: Bad/no value given for qmrbound. Set to qmrbound = 1.0D-6"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for qmrbound."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "nsra", confvalues%nsra, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for nsra."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "kxc", confvalues%kxc, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%kxc = 2
    write(*,*) "WARNING: Bad/no value given for kxc. Set to kxc = 2"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for kxc."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "rmax", confvalues%rmax, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rmax."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "tempr", confvalues%tempr, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%tempr = 800.0
    write(*,*) "WARNING: Bad/no value given for tempr. Set to tempr = 800.0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for tempr."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "alat", confvalues%alat, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for alat."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "rclust", confvalues%rclust, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%rclust = 1.5
    write(*,*) "WARNING: Bad/no value given for rclust. Set to rclust = 1.5"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rclust."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "mixing", confvalues%mixing, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for mixing."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "fcm", confvalues%fcm, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%fcm = 20.0
    write(*,*) "WARNING: Bad/no value given for fcm. Set to fcm = 20.0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for fcm."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "scfsteps", confvalues%scfsteps, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for scfsteps."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "kforce", confvalues%kforce, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%kforce = 1
    write(*,*) "WARNING: Bad/no value given for kforce. Set to kforce = 1"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for kforce."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "emin", confvalues%emin, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for emin."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDoubleVector(conf, "bravais_a", confvalues%bravais_a, 3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bravais_a."
    call destroyConfigReader(conf)
    return
  end if
  call getValueLogical(conf, "jij", confvalues%jij, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%jij = .FALSE.
    write(*,*) "WARNING: Bad/no value given for jij. Set to jij = .FALSE."
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for jij."
    call destroyConfigReader(conf)
    return
  end if
  call getValueLogical(conf, "ldau", confvalues%ldau, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%ldau = .FALSE.
    write(*,*) "WARNING: Bad/no value given for ldau. Set to ldau = .FALSE."
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for ldau."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "rcutjij", confvalues%rcutjij, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%rcutjij = 2.30
    write(*,*) "WARNING: Bad/no value given for rcutjij. Set to rcutjij = 2.30"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rcutjij."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "emax", confvalues%emax, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for emax."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "kte", confvalues%kte, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for kte."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "imix", confvalues%imix, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for imix."
    call destroyConfigReader(conf)
    return
  end if
  call getValueIntVector(conf, "bzdivide", confvalues%bzdivide, 3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for bzdivide."
    call destroyConfigReader(conf)
    return
  end if
  call getValueLogical(conf, "cartesian", confvalues%cartesian, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%cartesian = .TRUE.
    write(*,*) "WARNING: Bad/no value given for cartesian. Set to cartesian = .TRUE."
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for cartesian."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "npol", confvalues%npol, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npol."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "npnt1", confvalues%npnt1, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt1."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "npnt2", confvalues%npnt2, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt2."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "npnt3", confvalues%npnt3, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for npnt3."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "rclust_voronoi", confvalues%rclust_voronoi, ierror)
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for rclust_voronoi."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "nmin_panel", confvalues%nmin_panel, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%nmin_panel = 7
    write(*,*) "WARNING: Bad/no value given for nmin_panel. Set to nmin_panel = 7"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for nmin_panel."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "num_MT_points", confvalues%num_MT_points, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%num_MT_points = 10
    write(*,*) "WARNING: Bad/no value given for num_MT_points. Set to num_MT_points = 10"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for num_MT_points."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "MT_scale", confvalues%MT_scale, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%MT_scale = 0.98
    write(*,*) "WARNING: Bad/no value given for MT_scale. Set to MT_scale = 0.98"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for MT_scale."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "RMT_ref_scale", confvalues%RMT_ref_scale, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%RMT_ref_scale = 0.995
    write(*,*) "WARNING: Bad/no value given for RMT_ref_scale. Set to RMT_ref_scale = 0.995"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for RMT_ref_scale."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "use_semicore", confvalues%use_semicore, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%use_semicore = 0
    write(*,*) "WARNING: Bad/no value given for use_semicore. Set to use_semicore = 0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for use_semicore."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "ebotsemi", confvalues%ebotsemi, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%ebotsemi = 0.0
    write(*,*) "WARNING: Bad/no value given for ebotsemi. Set to ebotsemi = 0.0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for ebotsemi."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "emusemi", confvalues%emusemi, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%emusemi = 0.0
    write(*,*) "WARNING: Bad/no value given for emusemi. Set to emusemi = 0.0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for emusemi."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "n1semi", confvalues%n1semi, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%n1semi = 0
    write(*,*) "WARNING: Bad/no value given for n1semi. Set to n1semi = 0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for n1semi."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "n2semi", confvalues%n2semi, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%n2semi = 0
    write(*,*) "WARNING: Bad/no value given for n2semi. Set to n2semi = 0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for n2semi."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "n3semi", confvalues%n3semi, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%n3semi = 0
    write(*,*) "WARNING: Bad/no value given for n3semi. Set to n3semi = 0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for n3semi."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "fsemicore", confvalues%fsemicore, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%fsemicore = 1.0
    write(*,*) "WARNING: Bad/no value given for fsemicore. Set to fsemicore = 1.0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for fsemicore."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "target_rms", confvalues%target_rms, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%target_rms = 1.0D-8
    write(*,*) "WARNING: Bad/no value given for target_rms. Set to target_rms = 1.0D-8"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for target_rms."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "near_field", confvalues%near_field, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%near_field = 0
    write(*,*) "WARNING: Bad/no value given for near_field. Set to near_field = 0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for near_field."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "write_shapes", confvalues%write_shapes, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%write_shapes = 0
    write(*,*) "WARNING: Bad/no value given for write_shapes. Set to write_shapes = 0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for write_shapes."
    call destroyConfigReader(conf)
    return
  end if
  call getValueDouble(conf, "mt_zero_shift", confvalues%mt_zero_shift, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%mt_zero_shift = 0.0
    write(*,*) "WARNING: Bad/no value given for mt_zero_shift. Set to mt_zero_shift = 0.0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for mt_zero_shift."
    call destroyConfigReader(conf)
    return
  end if
  call getValueInteger(conf, "DEBUG_morgan_electrostatics", confvalues%DEBUG_morgan_electrostatics, ierror)
  if (ierror == CONFIG_READER_ERR_VAR_NOT_FOUND) then
    confvalues%DEBUG_morgan_electrostatics = 0
    write(*,*) "WARNING: Bad/no value given for DEBUG_morgan_electrostatics. Set to DEBUG_morgan_electrostatics = 0"
    ierror = 0
  end if
  if (ierror /= 0) then
    write(*,*) "Bad/no value given for DEBUG_morgan_electrostatics."
    call destroyConfigReader(conf)
    return
  end if
  call destroyConfigReader(conf)
write(*,*) "Finished reading information from input.conf"
end function

!-------------------------------------------------------------------------------
integer function readInputParamsFromFile(filename, confvalues) result(ierror)
  implicit none
  character(len=*), intent(in) :: filename
  type (InputParams), intent(inout) :: confvalues

  integer, parameter :: FILEHANDLE = 67

  ierror = 0
  open(FILEHANDLE, file=filename, form="unformatted")
  read(FILEHANDLE) confvalues%fred
  read(FILEHANDLE) confvalues%icst
  read(FILEHANDLE) confvalues%kpre
  read(FILEHANDLE) confvalues%bravais_b
  read(FILEHANDLE) confvalues%bravais_c
  read(FILEHANDLE) confvalues%gmax
  read(FILEHANDLE) confvalues%qmrbound
  read(FILEHANDLE) confvalues%nsra
  read(FILEHANDLE) confvalues%kxc
  read(FILEHANDLE) confvalues%rmax
  read(FILEHANDLE) confvalues%tempr
  read(FILEHANDLE) confvalues%alat
  read(FILEHANDLE) confvalues%rclust
  read(FILEHANDLE) confvalues%mixing
  read(FILEHANDLE) confvalues%fcm
  read(FILEHANDLE) confvalues%scfsteps
  read(FILEHANDLE) confvalues%kforce
  read(FILEHANDLE) confvalues%emin
  read(FILEHANDLE) confvalues%bravais_a
  read(FILEHANDLE) confvalues%jij
  read(FILEHANDLE) confvalues%ldau
  read(FILEHANDLE) confvalues%rcutjij
  read(FILEHANDLE) confvalues%emax
  read(FILEHANDLE) confvalues%kte
  read(FILEHANDLE) confvalues%imix
  read(FILEHANDLE) confvalues%bzdivide
  read(FILEHANDLE) confvalues%cartesian
  read(FILEHANDLE) confvalues%npol
  read(FILEHANDLE) confvalues%npnt1
  read(FILEHANDLE) confvalues%npnt2
  read(FILEHANDLE) confvalues%npnt3
  read(FILEHANDLE) confvalues%rclust_voronoi
  read(FILEHANDLE) confvalues%nmin_panel
  read(FILEHANDLE) confvalues%num_MT_points
  read(FILEHANDLE) confvalues%MT_scale
  read(FILEHANDLE) confvalues%RMT_ref_scale
  read(FILEHANDLE) confvalues%use_semicore
  read(FILEHANDLE) confvalues%ebotsemi
  read(FILEHANDLE) confvalues%emusemi
  read(FILEHANDLE) confvalues%n1semi
  read(FILEHANDLE) confvalues%n2semi
  read(FILEHANDLE) confvalues%n3semi
  read(FILEHANDLE) confvalues%fsemicore
  read(FILEHANDLE) confvalues%target_rms
  read(FILEHANDLE) confvalues%near_field
  read(FILEHANDLE) confvalues%write_shapes
  read(FILEHANDLE) confvalues%mt_zero_shift
  read(FILEHANDLE) confvalues%DEBUG_morgan_electrostatics
  close(FILEHANDLE)
end function

!-------------------------------------------------------------------------------
integer function writeInputParamsToFile(filename, confvalues) result(ierror)
  implicit none
  character(len=*), intent(in) :: filename
  type (InputParams), intent(inout) :: confvalues

  integer, parameter :: FILEHANDLE = 67

  ierror = 0
  open(FILEHANDLE, file=filename, form="unformatted")
  write(FILEHANDLE) confvalues%fred
  write(FILEHANDLE) confvalues%icst
  write(FILEHANDLE) confvalues%kpre
  write(FILEHANDLE) confvalues%bravais_b
  write(FILEHANDLE) confvalues%bravais_c
  write(FILEHANDLE) confvalues%gmax
  write(FILEHANDLE) confvalues%qmrbound
  write(FILEHANDLE) confvalues%nsra
  write(FILEHANDLE) confvalues%kxc
  write(FILEHANDLE) confvalues%rmax
  write(FILEHANDLE) confvalues%tempr
  write(FILEHANDLE) confvalues%alat
  write(FILEHANDLE) confvalues%rclust
  write(FILEHANDLE) confvalues%mixing
  write(FILEHANDLE) confvalues%fcm
  write(FILEHANDLE) confvalues%scfsteps
  write(FILEHANDLE) confvalues%kforce
  write(FILEHANDLE) confvalues%emin
  write(FILEHANDLE) confvalues%bravais_a
  write(FILEHANDLE) confvalues%jij
  write(FILEHANDLE) confvalues%ldau
  write(FILEHANDLE) confvalues%rcutjij
  write(FILEHANDLE) confvalues%emax
  write(FILEHANDLE) confvalues%kte
  write(FILEHANDLE) confvalues%imix
  write(FILEHANDLE) confvalues%bzdivide
  write(FILEHANDLE) confvalues%cartesian
  write(FILEHANDLE) confvalues%npol
  write(FILEHANDLE) confvalues%npnt1
  write(FILEHANDLE) confvalues%npnt2
  write(FILEHANDLE) confvalues%npnt3
  write(FILEHANDLE) confvalues%rclust_voronoi
  write(FILEHANDLE) confvalues%nmin_panel
  write(FILEHANDLE) confvalues%num_MT_points
  write(FILEHANDLE) confvalues%MT_scale
  write(FILEHANDLE) confvalues%RMT_ref_scale
  write(FILEHANDLE) confvalues%use_semicore
  write(FILEHANDLE) confvalues%ebotsemi
  write(FILEHANDLE) confvalues%emusemi
  write(FILEHANDLE) confvalues%n1semi
  write(FILEHANDLE) confvalues%n2semi
  write(FILEHANDLE) confvalues%n3semi
  write(FILEHANDLE) confvalues%fsemicore
  write(FILEHANDLE) confvalues%target_rms
  write(FILEHANDLE) confvalues%near_field
  write(FILEHANDLE) confvalues%write_shapes
  write(FILEHANDLE) confvalues%mt_zero_shift
  write(FILEHANDLE) confvalues%DEBUG_morgan_electrostatics
  close(FILEHANDLE)
end function

end module
