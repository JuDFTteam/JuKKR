module mod_config
  use nrtype
! integer :: INS
! integer :: IRMTD
! integer :: LMAXD
! integer :: NCORED
! integer :: lmmaxd
! integer,protected,save                                   ::     ICST     = 4   ! number of born iterations
! integer,protected,save                                   ::     INS      = 1    ! =1 full pot, =0 ASA
! integer,protected,save                                   ::     kvrel    = 1    ! =1 full pot, =0 ASA
! integer,protected,save                                   ::     nsra     = 2    ! =1 full pot, =0 ASA
! integer,protected,save                                   ::     nspin    = 2    ! =1 full pot, =0 ASA
! integer, parameter,private                               ::     dim_testflag=20 ! dimension of testflag array
! character(len=20),dimension(dim_testflag),private,save   ::     testflag =''    ! testflag array
! character(len=20),dimension(dim_testflag),private,save   ::     runflag =''     ! runflag array
! ! integer,protected,save                                   ::     naclsd=100      ! max number of atoms in a cluster
! ! integer,protected,save                                   ::     nclsd=-1        ! max number of clusters
!                                                                                 ! -1 means nclsd=natom
! real(kind=DP),protected,save                             ::     rcut=2.0        ! cluster radius
! real(kind=DP),protected,save                             ::     strmix=2.0        ! cluster radius
! real(kind=DP),protected,save                             ::     fcm=2.0        ! cluster radius

contains

subroutine config_read(config)
  use type_config
  use mod_log, only: log_write

  implicit none
!interface
type(config_type),intent(out) ::  config                       ! config type in which all
                                                               ! keywords are stored
integer,save               ::  first = 1                       ! make sure routine is called just once
character(len=*),parameter ::  cfilename_config = 'config.cfg' ! name of the config file
integer                    ::  ifile_config                    
character(len=200)         ::  string1,string2
character(len=200)         ::  keyword, keyword1, keyword2, ctemp
integer                    ::  ios,ios2,iline,itemp

! ********************************************************** 
! make sure the routine is just called once
! ********************************************************** 
ifile_config     = 1000000
if (first/=1) stop '[config_params] trying to change the config file twice is not permitted'
first=0

! ********************************************************** 
! setting default values
! ********************************************************** 


testflag='x'
runflag='x'











open(unit=ifile_config, file=cfilename_config, status='old', iostat=ios)
if (ios/=0) then
   write(*,*) '[read_config] config file does not exist'
   stop
end if

iline =0
do while (ios/=-1)
   iline=iline+1
   read(unit=ifile_config,fmt='(A)', iostat=ios) string1

   if (ios==-1) cycle
   if (string1(1:1) == "!" .or. string1(1:1) == "#") cycle
   if (len(trim(string1)) == 0) cycle

   read(string1,*) keyword2
   keyword1=trim(keyword2)
   select case (keyword1)
      case ('NSPIN=')
         read(string1,*,iostat=ios2) keyword1, config%nspin
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('INS=')
         read(string1,*,iostat=ios2) keyword1, config%INS
         config%kshape=config%INS
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('WAVEFUNC_RECALC_THRESHHOLD=')
         read(string1,*,iostat=ios2) keyword1, config%wavefunc_recalc_threshhold
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('XC=')
         read(string1,*,iostat=ios2) keyword1, config%modeexcorr
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('ICST=')
         read(string1,*,iostat=ios2) keyword1, config%icst
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('HFIELD=')
         read(string1,*,iostat=ios2) keyword1, config%hfield, config%hfield_apply_niter
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('HFIELD2=')
         read(string1,*,iostat=ios2) keyword1, config%hfield2(1),config%hfield2(2), config%hfield_apply_niter2
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if

      case ('KVREL=')
         read(string1,*,iostat=ios2) keyword1, config%kvrel
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
         if (config%kvrel==0) then 
           config%nsra=1
         elseif (config%kvrel==1) then 
           config%nsra=2
         elseif (config%kvrel==2) then 
           config%nsra=3
         elseif (config%kvrel==3) then 
           config%nsra=4
         else
           stop'[config] error KVREL=?'
         end if
      case ('MIXFAC=')
         read(string1,*,iostat=ios2) keyword1, config%mixfac
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('CALCFORCE=')
         read(string1,*,iostat=ios2) keyword1, config%calcforce
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('CALCJIJMAT=')
         read(string1,*,iostat=ios2) keyword1, config%calcJijmat
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('FCM=')
         read(string1,*,iostat=ios2) keyword1, config%fcm
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('QBOUND=')
         read(string1,*,iostat=ios2) keyword1, config%qbound
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('IMIX=')
         read(string1,*,iostat=ios2) keyword1, config%IMIX
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('NSIMPLEMIXFIRST=')
         read(string1,*,iostat=ios2) keyword1, config%NSIMPLEMIXFIRST
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('IMIXSPIN=')
         read(string1,*,iostat=ios2) keyword1, config%IMIXSPIN
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('SPINMIXFAC=')
         read(string1,*,iostat=ios2) keyword1, config%SPINMIXFAC
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('SPINMIXBOUND=')
         read(string1,*,iostat=ios2) keyword1, config%spinmixbound
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('CALCORBITALMOMENT=')
         read(string1,*,iostat=ios2) keyword1, config%calcorbitalmoment
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('ITDBRY=')
         read(string1,*,iostat=ios2) keyword1, config%ITDBRY
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('LATTICE_RELAX=')
         read(string1,*,iostat=ios2) keyword1, config%LATTICE_RELAX
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('SCFSTEPS=')
         read(string1,*,iostat=ios2) keyword1, config%SCFSTEPS
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('SPINORBIT=')
         read(string1,*,iostat=ios2) keyword1, config%kspinorbit
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('NCOLL=')
         read(string1,*,iostat=ios2) keyword1, config%ncoll
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('NPAN_LOG=')
         read(string1,*,iostat=ios2) keyword1, config%NPAN_LOG
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('NPAN_EQ=')
         read(string1,*,iostat=ios2) keyword1, config%NPAN_EQ
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('NCHEB=')
         read(string1,*,iostat=ios2) keyword1, config%NCHEB
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('NPAN_LOGPANELFAC=')
         read(string1,*,iostat=ios2) keyword1, config%npan_logfac
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('RADIUS_LOGPANELS=')
         read(string1,*,iostat=ios2) keyword1, config%RLOGPAN
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('RADIUS_MIN=')
         read(string1,*,iostat=ios2) keyword1, config%RMIN
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('TESTFLAG=')
         string2=' x x x x x x x x x x x x x x x x x x x x'
         string1=trim(string1)//string2
         read(string1,*,iostat=ios2) keyword,testflag
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
      case ('RUNFLAG=')
         string2=' x x x x x x x x x x x x x x x x x x x x'
         string1=trim(string1)//string2
         read(string1,*,iostat=ios2) keyword,runflag
         if (ios2 /= 0) then
            write(*,*) 'error in config ', keyword1
            stop
         end if
   end select
end do
! ********************************************************** 
! set some default values                          --
! ********************************************************** 

if ( .not. config_testflag('tmatnew') .and. config%kspinorbit==1 ) then
  stop '[config] spinorbit only works with the new solver'
end if

if ( .not. config_testflag('tmatnew') .and. config%ncoll==1 ) then
  stop '[config] non-collinear magnetism only works with the new solver'
end if

if (config%kvrel==2 .and. config%nspin==1) then 
  stop '[config] nspin=1 and kvrel=2 conflicts'
end if

if (config%kspinorbit==1 .and. config%ncoll==0) then 
  stop '[config] ncoll=0 does not work if spinorbit=1'
end if

if (config%kspinorbit==1 .and. config%NSRA==3) then 
  stop '[config] config%kspinorbit==1 does not work if KVREL=2'
end if

if (config%calcjijmat==1 .and. config%ncoll==0) then 
  stop '[config] config%calcjijmat==1 does not work if NCOLL=0'
end if



! none jet ;-)

!--------------------------------------------------------
!--  write out config                                  --
!--------------------------------------------------------

write(1337,*) '########################################'
write(1337,*) '#######  Config Parameter   ############'
write(1337,*) '########################################'
write(1337,*)          'ICST     = ',config%ICST
write(1337,*)          'INS      = ',config%INS
write(1337,*)          'KVREL    = ',config%KVREL
write(1337,*)          'NSRA     = ',config%NSRA
write(1337,*)          'NSPIN    = ',config%NSPIN
write(1337,*)          'IMIX      = ',config%IMIX
write(1337,'(A,F17.3)') 'MIXFAC   = ',config%MIXFAC
write(1337,'(A,F17.3)') 'FCM      = ',config%FCM
write(1337,'(A,E22.3)') 'QBOUND   = ',config%QBOUND
write(1337,*)          'SCFSTEPS =  ',config%SCFSTEPS
! write(1337,*) 'RCUT   =  ',config%RCUT
! write out test- and running flag information
write(1337,*) '########################################'
write(1337,*) ' TESTFLAG=         '
write(1337,*) '########################################'

do itemp = 1,dim_flags
  if (trim(testflag(itemp))/='x') then
    write(1337,'(A)',advance='no') testflag(itemp)
  end if 
end do

write(1337,'(A)') ''
write(1337,*) '########################################'
write(1337,*) ' RUNFLAGS         '
write(1337,*) '########################################'
do itemp = 1,dim_flags
  if (trim(runflag(itemp))/='x') then
    write(1337,'(A)',advance='no') runflag(itemp)
  end if 
end do
write(1337,*) ''
write(1337,*) '########################################'

end subroutine config_read


function config_testflag(ctestflag) result(istestflag)
  use type_config, only: testflag,dim_flags
  implicit none
!interface variables
  character(len=*),intent(in)       :: ctestflag
  logical                           :: istestflag
!local variables
  integer                           :: itemp

istestflag=.false.
do itemp=1,dim_flags
    if (trim(ctestflag)==trim(testflag(itemp))) istestflag=.true.
end do 
end function config_testflag

function config_runflag(crunflag) result(istestflag)
  use type_config, only: runflag,dim_flags
  implicit none
!interface variables
  character(len=*),intent(in)      :: crunflag
  logical                          :: istestflag
!local variables
  integer                :: itemp

istestflag=.false.
do itemp=1,dim_flags
  if (trim(crunflag)==trim(runflag(itemp))) istestflag=.true.
end do 
end function config_runflag

end module mod_config