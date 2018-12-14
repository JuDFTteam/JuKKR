!------------------------------------------------------------------------------------
!> Summary: KKRimp program
!> Author: 
!> 
!------------------------------------------------------------------------------------
program kkrflex


#ifdef CPP_MPI
  use mpi
#endif
! modules
  use nrtype
  use mod_read_potential
  use mod_config, only: config_read, config_testflag,config_runflag
  use mod_read_atominfo
  use mod_gauntharmonics, only: gauntharmonics_set
!   use mod_gauntharmonics_test, only: gauntharmonics_set_test
  use mod_calctmat
  use mod_rhocore
  use mod_gauntshape, only: gen_gauntshape
  use arrayparams, only: arrayparams_set
  use mod_energyloop
  use mod_rhototb
  use mod_preconditioning
  use mod_mpienergy, only: mpienergy_distribute
  use mod_vinters2010
  use mod_vintras
  use mod_convol
  use MOD_MIXSTR
  use mod_mixbroyden
  use mod_vxcdrv, only: vxcdrv
  use mod_rites
  use mod_epotinb
  use mod_ESPCB
  use mod_ECOUB
  use mod_etotb1
  use mod_wrmoms
  use mod_calcforce
  use mod_utrafo
  use mod_initldau     ! lda+u
  use mod_calcwldau    ! lda+u
  use mod_averagewldau ! lda+u

! type definitions
  use type_cell
  use type_shapefun
  use type_gauntshape
  use type_corestate
  use type_density
  use type_config
  use type_energyparts
  use type_gmat
  use type_gmatonsite
  use type_tmat
  use type_gmatbulk
  use type_ldau               ! lda+u
  use mod_jmtrx !test
  use mod_timing
  use mod_log, only: log_write
  use mod_calctmat_bauernew !test
  use mod_change_nrmin

  use global_variables, only: ipand
#ifdef CPP_MPI
  use mod_mympi, only: myrank, master
#endif


  use mod_mathtools
  use mod_version_info
  implicit none 

!***********************************
! main variables
!***********************************
  integer                               :: natom                                 ! number of impurity atoms
  integer                               :: nspin                                 ! number of spin directions
                                                                                 ! 1= paramagnetic
                                                                                 ! 2= collinear,non collinear, SOC, ...
  integer                               :: lmaxd                                 ! maximum lmax of all atoms 
                                                                                 ! meaning max(lmaxatom(iatom),iatom=1,natom)
  integer,allocatable                   :: lmaxatom(:)                           ! lmax value for each atom
  integer,allocatable                   :: killatom(:)                           ! host sites which are removed from the cluster
  integer,allocatable                   :: isvatom(:)                            ! host site which is a calculated as a
                                                                                 ! virtual atom (meaning delta t = 0)
  real(kind=dp),allocatable             :: zatom(:)                              ! nucleus charge
  real(kind=dp),allocatable             :: rimpatom_host(:,:)                    ! real space impurity positions given by the host
  real(kind=dp),allocatable             :: rimpatom(:,:)                         ! real space impurity positions shifted by rimpshift
                                                                                 ! rimpatom = rimpatom_host + rimpshift
  real(kind=dp),allocatable             :: rimpshift(:,:)                        ! shift of the atomic positions by the U-trafo
  real(kind=dp),allocatable             :: vpot(:,:,:,:)                         ! input potential array
  real(kind=dp)                         :: vpottemp
  real(kind=dp),allocatable             :: vpot_out(:,:,:,:)                     ! output potential array
  real(kind=dp),allocatable             :: cmom(:,:)                             ! charge moments
  real(kind=dp),allocatable             :: cmom_interst(:,:)                             ! charge moments
  real(kind=dp)                         :: alat                                  ! bulk lattice constant
  integer                               :: itscf                                 ! selfconsistency running variable
  real(kind=dp)                         :: vmtzero(2)                               ! muffin tin zero from the bulk code
  integer                               :: iatom, ispin, idummy, ilm             ! running variables

!***********************************
! type variables
!***********************************
  type(cell_type),allocatable           :: cell(:)                               ! cell properties
  type(energyparts_type)                :: energyparts                           ! total energy stuff
  type(shapefun_type),allocatable       :: shapefun(:)!(natom)                   ! information on the shape functions
  type(gauntshape_type),allocatable    :: gauntshape(:)                          ! gaunt coefficients used for convoluting
                                                                                 ! shape functions
  type(corestate_type),allocatable      :: corestate(:)!(natom)                  ! derived type containing all shape information 
  type(density_type),allocatable        :: density(:)!(natom)                    ! type containing information on the density
                                                                                 ! of each atom
  type(config_type)                     :: config                                ! type containing all config variables
  type(gmat_type)                       :: gmat                                  ! type containing information on the real
                                                                                 ! space Greens function Gnn'
  type(ldau_type),allocatable           :: ldau(:)                               ! lda+u variables, intended dimension: (NATOM)  ! lda+u
                                                                                 ! space Greens function Gnn'
!***********************************
! energy variables 
! needs maybe to be combined to an derived type
!***********************************
  double complex,allocatable            :: ez(:)                                 ! energy mesh from the host
  double complex,allocatable            :: wez(:)                                ! integration weights
  integer                               :: ielast                                ! number of energy points in ez, wez
  double complex                        :: llyfac                                ! renormalization factor for the simplest form of lloyds formula (opton LLYsimple)
!   integer,allocatable                   :: mpi_ielast(:)
  real(kind=dp),allocatable             :: intercell_ach(:,:)                    ! intercell potential of the host
                                                       !(lmpotd,ntotatom)        ! without the contribution of the 
                                                                                 ! impurity atoms
!***********************************
! mpi stuff
!***********************************
  integer                               :: my_rank,mpi_size,ierror
!***********************************
! constants
!***********************************
  real(kind=8),parameter                  :: rfpi=3.5449077018110318d0 !sqrt(4*pi)
!***********************************
! temp variables
!***********************************
  character(len=20)                        :: ctemp
!***********************************
! Bulk Properties
!***********************************
 type(gmatbulk_type)                      :: gmatbulk
!***********************************
! mixing stuff
!***********************************
  real(kind=8)                          :: rmsavq, rmsavm
  real(kind=8)                          :: mixldau ! lda+u
! !   real(kind=8)                          :: mixing
  real(kind=8)                          :: sum,rv
  integer                               :: irmin1,irc1
  integer                               :: ir        ! running variable for the radial mesh
  integer                              :: lmpotin


  type(gmatonsite_type),allocatable   :: gmatonsite(:,:)
  type(tmat_type),allocatable    :: tmatll(:,:) !(lmmaxd,lmmaxd)
  integer                              :: istop_selfcons





! print *, rotvector(0.0D0,0.0D0,(/ -0.001D0,0.001D0,-1.0D0 /),1.000000D0)

! stop





! **************************************************************************
! **************************************************************************
!                               KKR FLEX IMPURITY CODE
! **************************************************************************
! **************************************************************************

! mpi stuff
my_rank=0
mpi_size=1

#ifdef CPP_MPI
      CALL MPI_INIT(ierror)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierror)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierror)

      myrank = my_rank
      master = 0
#endif
! find serial number that is printed to files
call construct_serialnr()

call timing_init(my_rank)
call timing_start('Total running time')
call timing_start('time until scf starts')

if (my_rank==0) then
  write(*,*) ' **************************************************************************'
  write(*,*) ' **************************************************************************'
  write(*,*) '                               KKR FLEX IMPURITY CODE'
  write(*,'(2A)') '                          Version: ',serialnr
  write(*,*) ' **************************************************************************'
  write(*,*) ' **************************************************************************'
end if

#ifdef CPP_MPI
if (my_rank==0) then
  write(*,*) ' ###############################################'
  write(*,*) ' ##########   MPI Initialization    ############'
  write(*,*) ' ###############################################'
  write(*,*) ' ###    using ',mpi_size,' processors'
  write(*,*) ' ###############################################'
end if 
#endif

write(*,*) 'check all matrix inversions. There might be an error due: Hermitian'

! ********************************************************** 
! open the log file for each processor 
! file is called out_log.xxx.txt where xxx is 
! the processor id (my_rank)
! ********************************************************** 
write(ctemp,'(I03.3)') my_rank
open(unit=1337, file='out_log.'//trim(ctemp)//'.txt')
call version_print_header(1337)
! ********************************************************** 
! first all parameters are read in from the config
! file and stored into the config type
! ********************************************************** 
call log_write('>>>>>>>>>>>>>>>>>>>>> read_config >>>>>>>>>>>>>>>>>>>>>')
call config_read(config)
call log_write('<<<<<<<<<<<<<<<<<<< end read_config <<<<<<<<<<<<<<<<<<<')
nspin=config%nspin

! Check compatibility of config flags   ! lda+u
if (config_runflag('LDA+U').and..not.config_testflag('tmatnew')) stop &
 'testflag tmatnew should be applied if runflag LDA+U is used.' 


! ********************************************************** 
! the gaunt coefficients are calculated using the Juelich Muenichen
! routines. They were developed for a fixes LMAX cutoff. In the 
! KKRFLEX Code an LMAX can be defined for each atom. Therefore
! an array is created containing the data for all LMAX components
! which are used.
! ********************************************************** 

call log_write('>>>>>>>>>>>>>>>>>>>>> set_gaunt_harmonics >>>>>>>>>>>>>>>>>>>>>')
! allocates gaunt coefficients from l=2,..4
call gauntharmonics_set()
call log_write('<<<<<<<<<<<<<<<<<<< end set_gaunt_harmonics <<<<<<<<<<<<<<<<<<<')



! call log_write('>>>>>>>>>>>>>>>>>>>>> set_gaunt_harmonics test>>>>>>>>>>>>>>>>>>>>>')
! call gauntharmonics_set_test()
! call log_write('<<<<<<<<<<<<<<<<<<< end set_gaunt_harmonics test<<<<<<<<<<<<<<<<<<<')
! stop

!test
! allocate(test1(16,16))
! call  JMTRX(0.000001D0,0.000001D0,0.0000001D0,(1.5D0,0.0D0),3,test1,test2)
! do iatom=1,16
!   write(123456,'(50000F)') test1(iatom,:)
! end do
! stop


! ********************************************************** 
! The routine does :
! -  a dyson step (with -t_i, i=all atoms to be
!    removed ) to kill all atoms which are substituted within the
!    impurity calculation
! -  cut off all LM components in the greensfunction which are
!    not used  
! -  reads the intercell potential and substracts all contributions
!    of all atoms which are removed or substituted by other
!    atoms
! ********************************************************** 
! - ez, wez, ielast are read from the xxxxxx file
! - intercell_ach is the intercell potential (without the contribution
!   of all impurity atoms)
! - alat is the lattice constant read from the intercell file
! ********************************************************** 
call log_write('>>>>>>>>>>>>>>>>>>>>> PRECONDITIONING_start >>>>>>>>>>>>>>>>>>>>>')
call timing_start('PRECONDITIONING_start')
call PRECONDITIONING_start(my_rank,mpi_size,ez, wez, ielast, intercell_ach,alat,vmtzero,config%lattice_relax,gmatbulk)
call timing_stop('PRECONDITIONING_start')
call log_write('<<<<<<<<<<<<<<<<<<< end PRECONDITIONING_start <<<<<<<<<<<<<<<<<<<')



! LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple
if ( config_runflag('LLYsimple') ) then
  call log_write('>>>>>>>>>>>>>>>>>>>>> NEWWEIGHTS >>>>>>>>>>>>>>>>>>>>>')
  if (my_rank==0) then
    write(*,'(A)') 'Found option LLYsimple: reading in renormalization factor from file kkrflex_llyfac'
    open(192837, file='kkrflex_llyfac', form='formatted', iostat=ierror)
    if(ierror/=0) stop 'Error: File kkrflex_llyfac not found, needed for LLYsimple option'
    read(192837, *) llyfac
    close(192837)
    write(*,*) 'Renormalize weights with factor:',llyfac
  end if
#ifdef CPP_MPI
  call MPI_Bcast(llyfac, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
  if(ierror/=0) stop 'Error in MPI_Bcast for llyfac'
#endif
  !renormalize weights on every rank
  wez(:) = wez(:)*llyfac
  do idummy=1,ielast
    write(1337, '(A,I5,A,F20.14,A,F20.14)') 'IE: ',idummy,' new weight: ',real(wez(idummy)), ' ', imag(wez(idummy))
  end do
  call log_write('<<<<<<<<<<<<<<<<<<< end NEWWEIGHTS <<<<<<<<<<<<<<<<<<<')
end if
! LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple LLYsimple



! ********************************************************** 
! Reads the atomic positions, LMAX cutoff per Atom and other properties
! of the impurity atoms
! ********************************************************** 
call log_write('>>>>>>>>>>>>>>>>>>>>> read_atominfo >>>>>>>>>>>>>>>>>>>>>')
call read_atominfo('imp','kkrflex_atominfo',natom,idummy,rimpatom_host,zatom, &
!                      <           <            >     >       >       >  
                     lmaxd,lmaxatom,killatom,isvatom)
!                      >       >        >       >
call log_write('<<<<<<<<<<<<<<<<<<< end read_atominfo <<<<<<<<<<<<<<<<<<<')


! ********************************************************** 
! sets up some array parameters (will probably be removed)
! ********************************************************** 
call log_write('>>>>>>>>>>>>>>>>>>>>> set_array_params >>>>>>>>>>>>>>>>>>>>>')
call arrayparams_set(lmaxd)
call log_write('<<<<<<<<<<<<<<<<<<< end set_array_params <<<<<<<<<<<<<<<<<<<')



! ********************************************************** 
!  read in the initial potentials of the impurity atoms
!  and does some consistency checks 
! ********************************************************** 
call timing_start('read potential')
call log_write('>>>>>>>>>>>>>>>>>>>>> read_potential >>>>>>>>>>>>>>>>>>>>>')
call read_potential('potential','shapefun',natom,lmaxd,zatom,lmaxatom, &
!                          <             <           <    <     <      <
                        cell,shapefun, corestate, vpot,alat,nspin,config%ins)
!                        >       >         >       >    <     <       <
call log_write('<<<<<<<<<<<<<<<<<<< end read_potential <<<<<<<<<<<<<<<<<<<')
call timing_stop('read potential')


! ********************************************************** 
! Set up exchange correlation mode
! ********************************************************** 

if (my_rank==0) print *,'****************************************************************'
if ((trim(config%modeexcorr)=='LDA').OR.(trim(config%modeexcorr)=='LDA-VWN')) then
  do iatom=1,natom
    cell(iatom)%kxc=2
  end do
  if (my_rank==0) print *,'using LDA (Vosko-Wilk-Nusair)'
elseif (trim(config%modeexcorr)=='LDA-vBH') then
  do iatom=1,natom
    cell(iatom)%kxc=1
  end do
  if (my_rank==0) print *,'using LDA (von Bath-Hedin)'
elseif (trim(config%modeexcorr)=='LDA-MJW') then
  do iatom=1,natom
    cell(iatom)%kxc=0
  end do
  if (my_rank==0) print *,'using LDA (Moruzzi-Janak-Williams)'
elseif (trim(config%modeexcorr)=='GGA') then
  do iatom=1,natom
    cell(iatom)%kxc=3
  end do
  if (my_rank==0) print *,'using GGA (PW91)'
elseif (trim(config%modeexcorr)=='File') then
  open(unit=23437243,file='kkrflex_xc')
  if (my_rank==0)  then
    print *,' Using different exchange correlation functionals for each atom '
    print *,' kxc=2 (LDA), kxc=3 (GGA)'
  end if
  do iatom=1,natom
    read(23437243,*) cell(iatom)%kxc
    if (my_rank==0) print *,iatom,cell(iatom)%kxc
    if (cell(iatom)%kxc/=2 .and. cell(iatom)%kxc/=3) stop 'indiv. ex.corr. functional not known'
  end do
  close(23437243)
else
  stop 'ex.corr. functional not known'
end if
if (my_rank==0) print *,'****************************************************************'

if ( config_runflag('SIMULASA') ) then
  if (my_rank==0) then
    write(*,*) 'Run mode SIMULASA'
    write(*,*) 'Cutting of all spherical components'
  end if
  vpot(:,2:,:,:)=0.0D0
end if

! **********************************************************                         !lda+u
! Start LDA+U initialization                                                         !lda+u
! In particular, U-matrix ULDAU and basis PHI are created.                           !lda+u
! Allocate arrays that depend on NATOM. (Arrays dependent on                         !lda+u
! NATLDAU are allocated in routines INITLDAU,RWLDAUPOT.)                             !lda+u
      ALLOCATE( LDAU(NATOM) )                                                        !lda+u
      LDAU(:)%LOPT = -1                                                              !lda+u
      IF ( CONFIG_RUNFLAG('LDA+U') ) THEN                                            !lda+u
         WRITE(*,*) 'LDA+U calculation'                                              !lda+u
         CALL INITLDAU(LMAXD,NATOM,NSPIN,VPOT,ZATOM,1,CELL,LDAU)                     !lda+u
         DO IATOM = 1,NATOM                                                          !lda+u
            IF (LDAU(IATOM)%LOPT.GE.0) WRITE(*,FMT='(A12,I4,a3,I3,3(A6,F6.3))') &    !lda+u
            'LDA+U: Atom',IATOM,' l=',LDAU(IATOM)%LOPT,' UEFF=',LDAU(IATOM)%UEFF, &  !lda+u
            ' JEFF=',LDAU(IATOM)%JEFF,' EREF=',LDAU(IATOM)%EREFLDAU                  !lda+u
            IF (LDAU(IATOM)%LOPT.GT.LMAXATOM(IATOM)) THEN                            !lda+u
               WRITE(*,*) 'Atom:',IATOM,' LDA+U orbital=',LDAU(IATOM)%LOPT,  &       !lda+u
                    ' but lmax=',LMAXATOM(IATOM)                                     !lda+u
               STOP 'LDA+U control: lmax'                                            !lda+u
            ENDIF                                                                    !lda+u
            LDAU(IATOM)%IELDAUSTART = 1                                              !lda+u
            LDAU(IATOM)%IELDAUEND = IELAST                                           !lda+u
         ENDDO                                                                       !lda+u
! CALL AVERAGEWLDAU   ! average read-in interaction over m's                         !lda+u
         CALL AVERAGEWLDAU(NATOM,NSPIN,LDAU)                                         !lda+u
      ENDIF                                                                          !lda+u
! **********************************************************                         !lda+u


if ( config_testflag('step_potential') ) then
  do iatom=1,natom
    do ispin=1,natom
      vpottemp=1.0D0
      vpot(1,:,ispin,iatom)=vpottemp  
      do ir=2,cell(iatom)%nrmax
        if (cell(iatom)%rmesh(ir)==cell(iatom)%rmesh(ir-1)) vpottemp=vpottemp+1.0D0
        do ilm=1,(2*lmaxatom(iatom)+1)**2
          vpot(ir,ilm,ispin,iatom)=vpottemp
        end do
      end do
    end do
  end do
  do iatom=1,natom
    do ispin=1,natom
      do ilm=1,(2*lmaxatom(iatom)+1)**2
        write(424,'(50000E)') vpot(:,ilm,ispin,iatom)
      end do
    end do
  end do
end if


if ( config_testflag('change_nrmin') ) then
  call change_nrmin(natom,cell)

  call rites(nspin,natom,zatom,alat,cell(1)%nrmaxd,lmaxd,config%ins, &
            config%qbound,dreal(ez(ielast)),vmtzero,cell,corestate,lmaxatom,vpot)

  stop
end if

! read(5339,*) vpot
!     write(5338,'(50000F)') vpot
! vpot(:cell(1)%NRMIN_NS,2:,:,:)=0.0D0
! vpot(:210,2:,:,:)=0.0D0


if ( config_testflag('write_vpotin') ) then
  open(unit=43623223,file='test_vpotin')
  do iatom=1,natom
    do ispin=1,nspin
      write(43623223,*) '# atom ',iatom,'ispin ',ispin
      write(43623223,'(50000g24.16)') vpot(:,:,ispin,iatom)
    end do !ispin
  end do !iatom
end if ! config_testflag('write_gmatonsite')


! ********************************************************** 
! generates the gaunt coefficients used to convolute the potential
! etc with the shape functions
! ********************************************************** 
! if (config%ins/=0) then
  call log_write('>>>>>>>>>>>>>>>>>>>>> gen_gauntshape >>>>>>>>>>>>>>>>>>>>>')
  call gen_gauntshape(lmaxd,natom,lmaxatom,gauntshape,shapefun,config%ins)
!    if (ins==0) allocate thetas(1,1)
  !                     <     <       <       >          <
  call log_write('<<<<<<<<<<<<<<<<<<< end gen_gauntshape <<<<<<<<<<<<<<<<<<<')
! else 
! allocate( gauntshape%gsh(1),gauntshape%imaxsh(0:(2*lmaxd+1)**2) )
! gauntshape%imaxsh=0
! gauntshape%gsh(1)=0
! end if

if (config_testflag('writeout_shapefun')) then
  open(unit=324,file='test_shapefun')
  do iatom=1,natom
    write(324,*) 'iatom',iatom
    write(324,*) 'shapefun(iatom)%nrshaped'
    write(324,*) shapefun(iatom)%nrshaped
    write(324,*) 'shapefun(iatom)%nlmshaped'
    write(324,*) shapefun(iatom)%nlmshaped
    write(324,*) 'shapefun(iatom)%nlmshape'
    write(324,*) shapefun(iatom)%nlmshape
    write(324,*) 'shapefun(iatom)%index2lm'
    write(324,*) shapefun(iatom)%index2lm
    write(324,*) 'shapefun(iatom)%lmused'
    write(324,*) shapefun(iatom)%lmused
    write(324,*) 'shapefun(iatom)%lm2index'
    write(324,*) shapefun(iatom)%lm2index
    write(324,*) 'shapefun(iatom)%thetas'
    write(324,*) shapefun(iatom)%thetas
  end do
  close(324)
end if ! testflag
call timing_stop('time until scf starts')


! ********************************************************** 
! if relaxations are used the Greens function of the
! host is beed transformed to a greensfunction of atomic
! positions which are shifted by a vector s
! This is done with the precision of the bulk lmmax value.
! Afterwards the Greensfunction is cut back to the lmmax
! value specified in the kkrflex_atominfo file
! ********************************************************** 


! for now the impurity atoms are not shifted
! need to write an input file for the shift
if (.not. allocated(rimpshift)) allocate(rimpshift(3,natom))
rimpshift=0.0D0

if ( config%lattice_relax==1 ) then
  call log_write('>>>>>>>>>>>>>>>>>>>>>  U-transformation >>>>>>>>>>>>>>>>>>>>> ')
    if ( my_rank==0 ) then
      call utrafo(ielast,ez,natom,lmaxatom,gmatbulk,rimpshift)
    end if
  call log_write('<<<<<<<<<<<<<<<<<<<  end U-transformation <<<<<<<<<<<<<<<<<<< ')
end if

! calculate the new impurity positions by considering
! the shift by rimpshift
if (.not. allocated(rimpatom)) allocate( rimpatom(3,natom) )

if ( config%lattice_relax==1 ) then
  rimpatom = rimpatom_host + rimpshift
else
  rimpatom = rimpatom_host
end if

! ********************************************************** 
! start of the selfconsistency cycle
! ********************************************************** 

!          if ( config_testflag('tmatnew') ) then


! call calctmat_bauernew(vpot(:,:,1,1),cell(1),lmaxatom(1),ez(1),zatom(1))
 
!         end if


! ********************************************************** 
! start of the selfconsistency cycle
! ********************************************************** 
istop_selfcons=0    ! if set to 1 the scf cycle will stop 
if (config%scfsteps<1)  stop 'scfsteps is too small'

call log_write('***********************************************************')
call log_write('**************  starting selfconsistency ******************')
call log_write('***********************************************************')

do itscf=1,config%scfsteps
  write(ctemp,'(I03.3)') itscf
  call timing_start('Iteration number '//ctemp)
  write(1337,*) ' Iteration Number ',itscf
  if (my_rank==0) write(   *,*) ' Iteration Number ',itscf

  if (itscf==1) then
    allocate(energyparts%espc(0:3,nspin,natom))
    allocate(energyparts%EPOTIN(natom))
    allocate(energyparts%ECOU(0:2*LMAXD,NATOM))
    allocate(energyparts%ESPV(0:lmaxd+1,NSPIN,NATOM))
  end if
  energyparts%espc=0.0D0
  energyparts%EPOTIN=0.0D0
  energyparts%ECOU=0.0D0
  energyparts%ESPV=0.0D0


if (itscf<=config%hfield_apply_niter) then
  do iatom=1,natom
    do ispin=1,nspin
      do ir = 1,cell(iatom)%NRCORE !irmin1-1 
        vpot(ir,1,ispin,iatom) = vpot(ir,1,ispin,iatom) - DBLE(2*ISPIN-3)*config%HFIELD
      end do
    end do !ispin
  end do !natom
end if

if (itscf<=config%hfield_apply_niter2) then
  do iatom=1,19
    if(dabs(config%HFIELD2(ispin))>0.0d0) then
      do ispin=1,nspin
        write(6,*) 'atom',iatom,'spin',ispin,'shifted by',config%HFIELD2(ispin)
          do ir = 1,cell(iatom)%NRCORE !irmin1-1 
            vpot(ir,1,ispin,iatom) = vpot(ir,1,ispin,iatom) + config%HFIELD2(ispin)
          end do
      end do !ispin
    end if
  end do !natom
end if



  call timing_start('energy loop')
  call log_write('>>>>>>>>>>>>>>>>>>>>>  energyloop >>>>>>>>>>>>>>>>>>>>> ')
  ! ********************************************************** 
  ! for all energies IE in the energy loop:
  !  -  t-matrix is calculated
  !  -  dyson equation is solved
  !  -  density is calculated
  ! ********************************************************** 
  call energyloop(my_rank,mpi_size,itscf,cell, vpot, shapefun,zatom,natom,nspin,lmaxatom, &
  !                   <    <     <        <      <     <     <       <
                    lmaxd,density,ielast,ez,wez,config,gmat,gmatonsite,tmatll,energyparts, &
                    ldau)                                                                ! lda+u

  !                   <    >>>      <    <   <     <    <>

  call log_write('<<<<<<<<<<<<<<<<<<<  end energyloop <<<<<<<<<<<<<<<<<<< ')
  call timing_stop('energy loop')

if (itscf<=config%hfield_apply_niter) then
  do iatom=1,natom
    if(dabs(config%HFIELD2(ispin))>0.0d0) then
      do ispin=1,nspin
        write(6,*) 'atom',iatom,'spin',ispin,'shifted back by',DBLE(2*ISPIN-3)*config%HFIELD
        do ir = 1,cell(iatom)%NRCORE !irmin1-1 
          vpot(ir,1,ispin,iatom) = vpot(ir,1,ispin,iatom) + DBLE(2*ISPIN-3)*config%HFIELD
        end do
      end do !ispin
    end if
  end do !natom
end if

if (itscf<=config%hfield_apply_niter2) then
  do iatom=1,19
    do ispin=1,nspin
    write(6,*) 'atom',iatom,'spin',ispin,'shifted back by',-config%HFIELD2(ispin)
      do ir = 1,cell(iatom)%NRCORE !irmin1-1 
        vpot(ir,1,ispin,iatom) = vpot(ir,1,ispin,iatom) - config%HFIELD2(ispin)
      end do
    end do !ispin
  end do !natom
end if


  if (my_rank==0) then 
    call timing_start('time after energy loop')
    call log_write('>>>>>>>>>>>>>>>>>>>>>  WRMOMS >>>>>>>>>>>>>>>>>>>>> ')
    call wrmoms(natom,nspin, density,lmaxd,lmaxd+1,lmaxatom)
    call log_write('<<<<<<<<<<<<<<<<<<<  end WRMOMS <<<<<<<<<<<<<<<<<<< ')



    call log_write('>>>>>>>>>>>>>>>>>>>>> rhocore >>>>>>>>>>>>>>>>>>>>>')
    ! ********************************************************** 
    ! the core charge is calculated 
    ! ********************************************************** 
    do iatom=1,natom
      do ispin=1,nspin
        call  rhocore(ispin,nspin,iatom,natom,cell(iatom),vpot(:,1,ispin,iatom),zatom(iatom), &
                      corestate(iatom),density(iatom),config,itscf)
      end do !ispin=1,nspin
    end do !iatom=1,natom
    call log_write('<<<<<<<<<<<<<<<<<<< end rhocore <<<<<<<<<<<<<<<<<<<')



    ! ********************************************************** 
    ! combine rhoval and rhocore
    ! ********************************************************** 
    call log_write('>>>>>>>>>>>>>>>>>>>>> rhototb >>>>>>>>>>>>>>>>>>>>>')
    call rhototb(itscf,natom,nspin,lmaxatom,density, &
                zatom,cell,shapefun,config%ins)
    call log_write('<<<<<<<<<<<<<<<<<<< end rhototb <<<<<<<<<<<<<<<<<<<')





  
  !       DO ISPIN = 1,NSPIN
  ! C ----------------------------------------------------------------------
  !          IF (KTE.EQ.1) THEN
  !             DO IATOM = 1,NATOM
  ! !                IPOT = (I1-1)*NSPIN + ISPIN
  !                ESPV(0,ISPIN,IATOM) = ESPV(0,ISPIN,IATOM) -
  !      +                        DREAL(EZ(IELAST))*CHRGNT/DBLE(NSPIN*NATOM)
  !             END DO
  !          END IF                 ! (KTE.EQ.1)
  ! C ----------------------------------
  !      end do !ispin
  
    if ( config_testflag('write_totden') ) then
      open(unit=34329823,file='test_totrho2ns')
      if (nspin==1) then
        do iatom=1,natom
          write(34329823,*) '# atom ',iatom
          write(34329823,'(50000g24.16)') density(iatom)%rho2ns(:,:,1)
        end do !iatom
      elseif (nspin==2) then
        do iatom=1,natom
          write(34329823,*) '# atom ',iatom,' charge density'
          write(34329823,'(50000g24.16)') density(iatom)%rho2ns(:,:,1)
          write(34329823,*) '# atom ',iatom,' spin density'
          write(34329823,'(50000g24.16)') density(iatom)%rho2ns(:,:,2)
        end do !iatom
      else
        stop '[energyloop] nspin error'
      end if
    end if ! config_testflag('write_gmatonsite')

    call log_write('>>>>>>>>>>>>>>>>>>>>> vintras >>>>>>>>>>>>>>>>>>>>>')
    ! ********************************************************** 
    ! calculate from the charge density:
    !  -  the intracell potential -> vpot_out
    !  -  the charge moments      -> cmom
    ! ********************************************************** 
    if (itscf==1) then 
      allocate(vpot_out(cell(1)%nrmaxd,(2*lmaxd+1)**2,nspin,natom), &
               cmom_interst( (2*lmaxd+1)**2,natom ), &
               cmom        ( (2*lmaxd+1)**2,natom ) )
      vpot_out=0.0D0;cmom_interst=0.0D0;cmom=0.0D0
    end if
    call vintras(natom, nspin, cell(1)%nrmaxd,lmaxd, lmaxatom,cell, vpot_out, &
    !              <      <          <          <        <     <       >
                shapefun, gauntshape, density,cmom,cmom_interst,config%ins)
    !                <         <          <      >       <
    call log_write('<<<<<<<<<<<<<<<<<<< end vintras <<<<<<<<<<<<<<<<<<<')
    


    if ( config_testflag('write_vintrascmom') ) then
      if (itscf==1) open(unit=7836785,file='test_vintrascmom')
      do iatom=1,natom
        do ispin=1,nspin
          write(7836785,*) '# atom ',iatom,'ispin ',ispin
          write(7836785,'(50000g24.16)') cmom(:,iatom)
        end do !ispin
      end do !iatom
    end if ! config_testflag('write_gmatonsite')
    
    
    if ( config_testflag('write_vintraspot') ) then
      if (itscf==1) open(unit=786785,file='test_vintraspot')
      do iatom=1,natom
        do ispin=1,nspin
          write(786785,*) '# atom ',iatom,'ispin ',ispin
          write(786785,'(50000g24.16)') vpot_out(:,:,ispin,iatom)
        end do !ispin
      end do !iatom
    end if ! config_testflag('write_gmatonsite')
    
    if ( config_testflag('vinters=0') ) then
      intercell_ach=0.0d0
      cmom=0.0D0
    end if
    call log_write('>>>>>>>>>>>>>>>>>>>>> vinters2010 >>>>>>>>>>>>>>>>>>>>> ')
    ! ********************************************************** 
    ! calculate the intercell potential and adds it to VPOT_OUT
    ! ********************************************************** 
    call vinters2010(natom,nspin,lmaxd,lmaxatom,cell,rimpatom,cmom,alat,&
    !                    <     <      <      <      <     <      <   <
                      config%ins,cell(1)%nrmaxd,vpot_out,intercell_ach,zatom,config%lattice_relax,rimpshift)
    !                       <             <          >          <
    call log_write('<<<<<<<<<<<<<<<<<<< end vinters2010 <<<<<<<<<<<<<<<<<<<')
    
    if ( config_testflag('write_vpotout') ) then
      if (itscf==1) open(unit=54633563,file='test_vpotinters')
      do iatom=1,natom
        do ispin=1,nspin
          write(54633563,*) '# atom ',iatom,'ispin ',ispin
          write(54633563,'(50000g24.16)') vpot_out(:,:,ispin,iatom)
        end do !ispin
      end do !iatom
    end if ! config_testflag('write_gmatonsite')
  
!-------------------------------------------------------------------                       ! lda+u
! Calculate output interaction potential for lda+u and mix.                                ! lda+u
    mixldau = config%mixfac                                                                ! lda+u
    if ( config_testflag('stepmixldau').or.config_runflag('stepmixldau') ) then            ! lda+u
      mixldau = 0.d0                                                                       ! lda+u
      if (mod(itscf,config%ITDBRY).eq.0) mixldau = config%mixfac                           ! lda+u
    endif                                                                                  ! lda+u
    if ( config_testflag('freezeldau').or.config_runflag('freezeldau') ) mixldau = 0.d0    ! lda+u
    call calcwldau(nspin,natom,lmaxd,cell(1)%nrmaxd,lmaxatom,density,mixldau,ldau)         ! lda+u
!-------------------------------------------------------------------                       ! lda+u


    if (config%calcforce==1) then
      call log_write('>>>>>>>>>>>>>>>>>>>>> force part 1 >>>>>>>>>>>>>>>>>>>>>')
      call calcforce('part1',cmom,cmom_interst,lmaxatom,lmaxd,nspin,natom,density,VPOT_OUT, &
                       cell,config%ins, zatom,cell(1)%nrmaxd,alat)
      call log_write('<<<<<<<<<<<<<<<<<<< force part 1 <<<<<<<<<<<<<<<<<<<')
    end if


    call log_write('>>>>>>>>>>>>>>>>>>>>> total energy >>>>>>>>>>>>>>>>>>>>>')
    ! ********************************************************** 
    ! calculate parts of the single particle energy
    ! ********************************************************** 
    ! single particle core energies
    call espcb(energyparts%espc,nspin,natom,corestate)
    ! input potential
    ipand = cell(1)%npand
    call epotinb(energyparts%epotin,nspin,natom,vpot,config%ins, &
                            lmaxatom,zatom,cell,density, &
                            ipand, cell(1)%nrmaxd,2*lmaxd)
    ! the electrostatic potential-energies
    call ecoub(cmom,cmom_interst,energyparts%ecou,cell,density,shapefun,gauntshape,lmaxatom,lmaxd,nspin,natom,vpot_out,zatom, &
                          config%ins,cell(1)%nrmaxd,(2*lmaxd+1)**2,2*lmaxd)
    
    call log_write('>>>>>>>>>>>>>>>>>>>>> total energy >>>>>>>>>>>>>>>>>>>>>')
  
  
    call log_write('>>>>>>>>>>>>>>>>>>>>> vxcdrv >>>>>>>>>>>>>>>>>>>>>')
    ! ********************************************************** 
    ! calculate the exchange-correlation potential and add it 
    !  to the potential
    ! ********************************************************** 
    if (itscf==1) allocate( energyparts%exc(0:2*lmaxd,natom) )
    call vxcdrv(energyparts%exc,config%kte,nspin,natom,density, & 
                  vpot_out, cell,config%kshape,gauntshape, shapefun,lmaxd, & 
                  (2*lmaxd), (2*lmaxd+1)**2, (4*lmaxd+1)**2, cell(1)%nrmaxd, lmaxatom,config%ins)
    call log_write('<<<<<<<<<<<<<<<<<<< end vxcdrv <<<<<<<<<<<<<<<<<<<')
    
    if ( config_testflag('write_vpotout') ) then
      if (itscf==1) open(unit=54644563,file='test_vpot_xc')
      do iatom=1,natom
        do ispin=1,nspin
          write(54644563,*) '# atom ',iatom,'ispin ',ispin
          write(54644563,'(50000g24.16)') vpot_out(:,:,ispin,iatom)
        end do !ispin
      end do !iatom
    end if ! config_testflag('write_gmatonsite')


    if (config%calcforce==1) then
      call log_write('>>>>>>>>>>>>>>>>>>>>> force part 2 >>>>>>>>>>>>>>>>>>>>>')
      call calcforce('part2',cmom,cmom_interst,lmaxatom,lmaxd,nspin,natom,density,VPOT_OUT, &
                     cell,config%ins, zatom,cell(1)%nrmaxd,alat)
      call log_write('<<<<<<<<<<<<<<<<<<< force part 2 <<<<<<<<<<<<<<<<<<<')
    end if

    call log_write('>>>>>>>>>>>>>>>>>>>>> shift >>>>>>>>>>>>>>>>>>>>>')
    ! ********************************************************** 
    ! shift the potential using the muffin thin zero
    ! ********************************************************** 
    ! muffin tin zero of the bulk needs to be implemented
    ! vmtzero=0.0D0
    do ispin = 1,nspin
      do iatom = 1,natom
        do ir = 1,cell(iatom)%nrmax !irceq(iatyp)
          vpot_out(ir,1,ispin,iatom) = vpot_out (ir,1,ispin,iatom) + vmtzero(ispin)*rfpi
  !       attention mtzero is not jet read in by the code
        end do
      end do
    end do
    call log_write('<<<<<<<<<<<<<<<<<<< end shift <<<<<<<<<<<<<<<<<<<')
  
    if ( config_testflag('write_vpotout') ) then
      if (itscf==1) open(unit=126456563,file='test_vpot_shift')
      do iatom=1,natom
        do ispin=1,nspin
          write(126456563,*) '# atom ',iatom,'ispin ',ispin
          write(126456563,'(50000g24.16)') vpot_out(:,:,ispin,iatom)
        end do !ispin
      end do !iatom
    end if ! config_testflag('write_gmatonsite')
  
    if (config%ins.ne.0) then
      call log_write('>>>>>>>>>>>>>>>>>>>>> convoldrv >>>>>>>>>>>>>>>>>>>>>')
      ! ********************************************************** 
      ! convolute the potential
      ! ********************************************************** 
      call convoldrv(cell(1)%nrmaxd,lmaxd,nspin,natom,lmaxatom,cell,gauntshape, &
      !                        <          <     <     <       <      <      < 
                    shapefun,zatom, vpot_out)
      !                 <      <      >>>
      call log_write('<<<<<<<<<<<<<<<<<<< end convoldrv <<<<<<<<<<<<<<<<<<<')
    end if

    if ( config_runflag('SIMULASA') ) then
      vpot_out(:,2:,:,:)=0.0D0
    end if
  
    if ( config_testflag('write_vpotout') ) then
      if (itscf==1) open(unit=546456563,file='test_vpot_conv')
      do iatom=1,natom
        do ispin=1,nspin
          write(546456563,*) '# atom ',iatom,'ispin ',ispin
          write(546456563,'(50000g24.16)') vpot_out(:,:,ispin,iatom)
        end do !ispin
      end do !iatom
    end if ! config_testflag('write_gmatonsite')
    
    
    
    if ( config_testflag('write_vpotout') ) then
      if (itscf==1) open(unit=54633563,file='test_vpotout')
      do iatom=1,natom
        do ispin=1,nspin
          write(54633563,*) '# atom ',iatom,'ispin ',ispin
          write(54633563,'(50000g24.16)') vpot_out(:,:,ispin,iatom)
        end do !ispin
      end do !iatom
    end if ! config_testflag('write_gmatonsite')
    
    
    call log_write('>>>>>>>>>>>>>>>>>>>>> mixstr >>>>>>>>>>>>>>>>>>>>>')
    ! ********************************************************** 
    ! straight mixing of the potential
    ! ********************************************************** 
    lmpotin=1
!     if (config%ins==1) lmpotin=(2*lmaxd+1)**2
    lmpotin=(2*lmaxd+1)**2
    call mixstr(rmsavq,rmsavm,config%ins,natom,lmaxatom,lmaxd, &
                  nspin,itscf, config%mixfac,config%fcm,vpot, vpot_out, &
                  cell,cell(1)%nrmaxd,lmpotin )
    call log_write('<<<<<<<<<<<<<<<<<<< end MIXSTR <<<<<<<<<<<<<<<<<<<')
  
    if ( config_testflag('write_vpotin') ) then
      open(unit=43623223,file='test_vpotin2')
      do iatom=1,natom
        do ispin=1,nspin
          write(43623223,*) '# atom ',iatom,'ispin ',ispin
          write(43623223,'(50000g24.16)') vpot(:,:,ispin,iatom)
        end do !ispin
      end do !iatom
    end if ! config_testflag('write_gmatonsite')

  
    if (config%imix>=3 .and. itscf>config%nsimplemixfirst) then
    call log_write('>>>>>>>>>>>>>>>>>>>>> mix broy >>>>>>>>>>>>>>>>>>>>>')
      call mixbroyden(vpot,vpot_out,config%ins,config%mixfac,nspin,cell,lmaxatom, &
                            natom,config%ITDBRY,config%imix,(2*lmaxd+1)**2,cell(1)%nrmaxd)
    call log_write('<<<<<<<<<<<<<<<<<<< end mix broy <<<<<<<<<<<<<<<<<<<')
    end if
  

    if ( config_testflag('write_vpotmixed') ) then
      if (itscf==1) open(unit=54633563,file='test_vpotmixed')
      do iatom=1,natom
        do ispin=1,nspin
          write(54633563,*) '# atom ',iatom,'ispin ',ispin
          write(54633563,'(50000g24.16)') vpot_out(:,:,ispin,iatom)
        end do !ispin
      end do !iatom
    end if ! config_testflag('write_gmatonsite')




    call log_write('>>>>>>>>>>>>>>>>>>>>> copy pot >>>>>>>>>>>>>>>>>>>>>')
    ! ********************************************************** 
    ! this does :
    !  - copy VPOT -> VPOT_OUT
    !  - remove the non-spherical part if rms of (potential-0) is below
    !    the scf cut-off value (because of numerical reasons)
    ! ********************************************************** 
            do ispin = 1,nspin
              do iatom = 1,natom
                
                if (vpot_out(1,1,ispin,iatom)>1.0e3) then 
                  write(*,*) '[WARNING] potential of atom', iatom, 'spin', ispin,'is too high'
                end if
                irc1 = cell(iatom)%nrmax !irc(it)
    !            call dcopy(irc1,vpot_out(1,1,ispin,iatom),1,vpot(1,1,ispin,iatom),1)
                if (config%ins==1) then
                  vpot    (:,:,ispin,iatom) = vpot_out(:,:,ispin,iatom)
                elseif ( config%ins==0 ) then
                  vpot    (:,1,ispin,iatom) = vpot_out(:,1,ispin,iatom)
                end if
                vpot_out(:,:,ispin,iatom) = 0.0d0
                if ( ( config%ins.ne.0 ) ) then
                  irmin1 = cell(iatom)%nrmin_ns !irmin(it)
                  do ilm = 2,(2*lmaxatom(iatom)+1)**2!lmpot
                    sum = 0.0d0
                    do ir = irmin1,irc1
                      rv = vpot(ir,ilm,ispin,iatom)*cell(iatom)%rmesh(ir)
!                       write(*,*) 'rv',vpot(ir,ilm,ispin,iatom),cell(iatom)%rmesh(ir)
                      sum = sum + rv*rv*cell(iatom)%drmeshdi(ir)
                    end do
                    if ( sqrt(sum).lt.config%qbound ) then
                      write(1337,*) 'cutting pot. ','ispin',ispin,'iatom',iatom,'ilm',ilm,'rms',sqrt(sum)
                      do ir = 1,irc1 ! the 1 was irmin1 before!
                          vpot(ir,ilm,ispin,iatom) = 0.0d0
                      end do
                    end if
                    if ( .not. config_testflag('nocut_sphpart') ) then
                      do ir = 1,irmin1-1 
                          vpot(ir,ilm,ispin,iatom) = 0.0d0
                      end do
                    end if
                  end do
                end if
              end do !iatom
            end do !ispin
    call log_write('<<<<<<<<<<<<<<<<<<< end copy pot <<<<<<<<<<<<<<<<<<<')
  
    do ispin = 1,nspin
      do iatom = 1,natom
        do ilm = 2,(2*lmaxatom(iatom)+1)**2!lmpot
          sum = 0.0d0
          irmin1 = cell(iatom)%nrmin_ns !irmin(it)
          do ir = 1,irmin1-1
            rv = vpot(ir,ilm,ispin,iatom)*cell(iatom)%rmesh(ir)
            sum = sum + rv*rv*cell(iatom)%drmeshdi(ir)
          end do
          if ( sqrt(sum)>1d-3 ) then
            write(*,'(A,I3,A)') '[WARNING] the non-sph. potential component',ilm,' is non-zero:'
            write(*,*) '          RMS is ',sqrt(sum)
          end if
        end do
      end do
    end do

!     write(5339,'(50000F)') vpot

    call log_write('>>>>>>>>>>>>>>>>>>>>> rites >>>>>>>>>>>>>>>>>>>>>')
    call rites(nspin,natom,zatom,alat,cell(1)%nrmaxd,lmaxd,config%ins, &
              config%qbound,dreal(ez(ielast)),vmtzero,cell,corestate,lmaxatom,vpot)
    call log_write('<<<<<<<<<<<<<<<<<<< end rites <<<<<<<<<<<<<<<<<<<')
  
  
    call log_write('>>>>>>>>>>>>>>>>>>>>> ETOTB1 >>>>>>>>>>>>>>>>>>>>>')
    if (config%kte.eq.1) then !.and. icc.eq.0)      why icc==0 ???????
      call etotb1(dreal(ez(ielast)),lmaxatom,energyparts,corestate, &
                  nspin,natom,2*lmaxd)
    end if
    call log_write('<<<<<<<<<<<<<<<<<<< end ETOTB1 <<<<<<<<<<<<<<<<<<<')



  ! ********************************************************** 
  ! Section for MYRANK==0 ends here 
  ! ********************************************************** 
  end if ! my_rank==0

  if (my_rank==0) then
    IF (MAX(RMSAVQ,RMSAVM).LT.config%QBOUND) THEN
      istop_selfcons=1
      WRITE(*,'(17X,A)') '++++++ SCF ITERATION CONVERGED ++++++'
      WRITE(*,'(79(1H*))')
    END IF

    if ( config_runflag('force_angles') ) then
      if (itscf==density(1)%nangleconfigur) then
        write(*,*) 'All magnetic configurations have been have been calculated'
        write(*,*) 'Selfconsistency will be stopped now!'
        istop_selfcons=1
      else
        istop_selfcons=0
      end if
    end if




  end if

#ifdef CPP_MPI
    if (my_rank==0 .and. istop_selfcons==1) call log_write('send stop signal')
    call mpi_bcast( istop_selfcons,1,MPI_INTEGER,0,MPI_COMM_WORLD, ierror)

  call mpi_bcast( vpot,cell(1)%NRMAXD*(2*LMAXD+1)**2*NSPIN*NATOM,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierror)
  do iatom = 1,natom                                                                                       ! lda+u
  if (ldau(iatom)%lopt.ge.0) then                                                                          ! lda+u
    call mpi_bcast( ldau(iatom)%wldau,(2*LMAXD+1)**2*NSPIN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, ierror)  ! lda+u
  endif                                                                                                    ! lda+u
  enddo                                                                                                    ! lda+u
#endif


  write(ctemp,'(I03.3)') itscf
  if (my_rank==0)   call timing_stop('time after energy loop')
  call timing_stop('Iteration number '//ctemp)

  if (istop_selfcons==1) then 
     call log_write('received stop signal: exiting scf cycle')
     exit !exit loop over selfconsitancy
  end if
end do ! selfconsistency

! ********************************************************** 
! delete the kkrflex_greennew file
if (my_rank==0) then
  write(*,*) 'Deleting file kkrflex_greennew, processor',my_rank
  close(3434560,status='delete')
end if
! ********************************************************** 

#ifdef CPP_MPI
call log_write('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
call log_write('>>>>>>>>>>>>>>>>>>>>> mpi_finalize >>>>>>>>>>>>>>>>>>>>>')
call log_write('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
call mpi_finalize(ierror)
#endif


call log_write('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
call log_write('>>>>>>>>>>>>>>>>>>>>> end program >>>>>>>>>>>>>>>>>>>>>')
call log_write('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

call timing_stop('Total running time')


end program kkrflex
