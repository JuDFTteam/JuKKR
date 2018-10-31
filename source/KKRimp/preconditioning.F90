module mod_preconditioning
  contains
   subroutine preconditioning_start(my_rank,mpi_size,ez, wez, ielast, intercell_ach, alat,vmtzero,lattice_relax,gmatbulk) 
#ifdef CPP_MPI
  use mpi
#endif
     use nrtype
     use mod_dysonvirtatom, only: dysonvirtatom
     use mod_read_atominfo
     use mod_config, only: config_testflag,config_runflag
     use mod_log, only: log_write
     use mod_mpienergy, only: mpienergy_distribute
     use type_gmatbulk

     implicit none
!interface
     double complex,allocatable,intent(out)  :: ez(:)
     double complex,allocatable,intent(out)  :: wez(:)
     integer,intent(out)                     :: ielast
     real(kind=DP),allocatable,intent(out)   :: intercell_ach(:,:)   ! intercell potential
                                               !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
     real(kind=DP),intent(out)               :: alat
     real(kind=dp)                           :: vmtzero(2)  
     integer                                 :: lattice_relax
     type(gmatbulk_type)                     :: gmatbulk
!local

     integer                     :: natom                ! number of impurity atoms
     integer                     :: ntotatom             ! number of imp atoms+killatoms
     integer,allocatable         :: lmaxatom(:)          ! lmax for all atoms
     integer,allocatable         :: isvatom(:)           ! 1=vatom (delta t=0), 0=no vatom
     integer,allocatable         :: killatom(:)          ! 1=atom will be removed in a dyson step
     integer                     :: KGREFSOC             ! =1 if SOC for the GREF 
     integer                     :: NSOC                 ! =2 if SOC for the GREF else: =1 

     integer                     :: ie
     integer                     :: natomimpd,lmsizehost
     double complex,allocatable  :: gmathost(:,:)
     double complex,allocatable  :: gmathostnew(:,:)
     integer                     :: recl1,recl2,nspin,ispin,iatom
     double complex,allocatable  ::  tmat(:,:,:,:,:)
     integer                     :: nlmhostnew
     real(kind=dp),allocatable   :: RIMPATOM(:,:),zatom(:)
     integer                      :: lmaxd
!mpi
      integer,allocatable                      :: mpi_iebounds(:,:)
      integer                                  :: my_rank
      integer                                  :: mpi_size,ierror

call log_write('>>>>>>>>>>>> preconditioning read_atominfo >>>>>>>>>>>>>>>>>>>>>')
call read_atominfo('total','kkrflex_atominfo',natom,ntotatom,RIMPATOM,&
                   zatom,lmaxd,lmaxatom,killatom,isvatom)
call log_write('<<<<<<<<<<<< preconditioning end read_atominfo <<<<<<<<<<<<<<<<<<<')


call log_write('>>>>>>>>>>>>>>>>>>>>> preconditioning_readtmatinfo >>>>>>>>>>>>>>>>>>>>>')
! ###########################################
! read out some information of the bulk 
! ###########################################
call preconditioning_readtmatinfo(NTOTATOM,NSPIN,IELAST,lmsizehost,KGREFSOC)
!                                     in    out    out      out 
call log_write('<<<<<<<<<<<<<<<<<<< end preconditioning_readtmatinfo <<<<<<<<<<<<<<<<<<<')

allocate(tmat(lmsizehost,lmsizehost,Ntotatom,IELAST,NSPIN-KGREFSOC))

! ###########################################
! read in the t-matricies
! ###########################################
call log_write('>>>>>>>>>>>>>>>>>>>>> preconditioning_readtmatinfo >>>>>>>>>>>>>>>>>>>>>')

if( .not. config_testflag('tmat_mpicomm')) then
! ###########################################
! every process reads the file itself
! ###########################################
  call preconditioning_readtmat(IELAST,lmsizehost,NTOTATOM,NSPIN,TMAT,KGREFSOC)
else
! ###########################################
! my_rank=0 reads and sends the information to the other processes
! ###########################################
  if (my_rank==0) then 
    write(*,*) 'my_rank=1 reads tmat and communicates to other processes'
    call preconditioning_readtmat(IELAST,lmsizehost,NTOTATOM,NSPIN,TMAT,KGREFSOC)
  end if
#ifdef CPP_MPI
   call mpi_bcast( TMAT, lmsizehost*lmsizehost*Ntotatom*IELAST*(NSPIN-KGREFSOC),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, ierror)
#endif
end if

call log_write('<<<<<<<<<<<<<<<<<<< end preconditioning_readtmatinfo <<<<<<<<<<<<<<<<<<<')

if (KGREFSOC==1) then
  NSOC=2
else
  NSOC=1
end if

call log_write('>>>>>>>>>>>>>>>>>>>>> preconditioning_readenergy >>>>>>>>>>>>>>>>>>>>>')
! ###########################################
! read in the energy mesh
! ###########################################
call preconditioning_readenergy(IELAST,NSPIN,EZ,WEZ,NATOMIMPD,NTOTATOM,lmsizehost,kgrefsoc)
call log_write('<<<<<<<<<<<<<<<<<< < end preconditioning_readenergy <<<<<<<<<<<<<<<<<<<')

if (NATOMIMPD/=NTOTATOM) stop '[preconditioning] NATOMIMPD/=NTOTATOM'

gmatbulk%natom=natom

if (KGREFSOC==0) then
  gmatbulk%lmax=int(sqrt(real(lmsizehost  ))-1)
  if ( (gmatbulk%lmax+1)**2 /= lmsizehost) then 
    print *,'KGREFSOC',KGREFSOC
    print *,'gmatbulk%lmax',gmatbulk%lmax
    print *,'lmsizehost',lmsizehost
    print *, (gmatbulk%lmax+1)**2, lmsizehost
    stop 'lmmax error'
  end if
else
  gmatbulk%lmax=int(sqrt(real(lmsizehost/2))-1)
  if ( 2* (gmatbulk%lmax+1)**2 /= lmsizehost) then
    print *,'KGREFSOC',KGREFSOC
    print *,'gmatbulk%lmax',gmatbulk%lmax
    print *,'lmsizehost',lmsizehost
    print *, 2*(gmatbulk%lmax+1)**2, lmsizehost
    stop 'lmmax error'
  end if
end if

gmatbulk%lmmax=lmsizehost
gmatbulk%nspin=nspin
gmatbulk%hostdim=natom*lmsizehost

! ###########################################
! check for inconsitencies in the host greens function file
! ###########################################
!    do iatom=1,ntotatom
!      if ( (lmaxatom(iatom)+1)**2>lmsizehost) stop 'preconditioning error'
!    end do

! ###########################################
! check for inconsitencies in the host greens function file and
! determine the size of the new matrix (sum of all lm-components)
! ###########################################
nlmhostnew=0
do iatom=1,ntotatom
  if (killatom(iatom)/=1) then
    if (NSOC*(lmaxatom(iatom)+1)**2>lmsizehost) then
      write(*,*) 'lmax value of atom ',iatom,' is greater than the lmax value of the host'
      stop
    end if
    nlmhostnew=nlmhostnew+NSOC*(lmaxatom(iatom)+1)**2  
  end if
end do

if (lattice_relax==1) nlmhostnew=natom*lmsizehost
write(1337,*) 'Number of lm-components for the ghost matrix:',nlmhostnew
if (lattice_relax==1) then
  write(1337,*) 'The lm-components are not cut of because they are needed'
  write(1337,*) 'by the U-transformation'
end if

allocate(gmathostnew(nlmhostnew,nlmhostnew))


! ###########################################
! now start the energy loop which reads in the host greens function of a given
! energy and then kills all the atoms which are to be removed using a Dyson step
! with -t_kill. Then all lm-components and atoms are removed which are not used
! in the calculation
! ###########################################

recl1=wlength*2*natomimpd*lmsizehost*natomimpd*lmsizehost
! print *, 'natomimpd',natomimpd,'lmsizehost',lmsizehost

! Obsolete after introducinf wlength
!!#ifdef GFORT
!! the gfortran compiler assumes a rec length of 1 byte whereas the 
!! ifc compiler assumes a rec length of 4byte by default.
!! one should change the bulk program in the future!!!
!!recl1=4*recl1
!!#endif

recl2=wlength*4*nlmhostnew*nlmhostnew
open (88,access='direct',recl=recl1,file='kkrflex_green',form='unformatted')
if (lattice_relax==0) then
  open (89,access='direct',recl=recl2,file='kkrflex_greennew',form='unformatted')
else
  open (89,access='direct',recl=recl2,file='kkrflex_greenvoid',form='unformatted')
end if

if(config_testflag('gmat_plain')) open (8989,file='kkrflex_greennew.txt')

call log_write('**********************************************************')
call log_write('preconditioning dyson steps')
call log_write('**********************************************************')
! ###################################
! # Distribute energy points IE to processors
! ###################################
call mpienergy_distribute(my_rank,mpi_size,ielast,mpi_iebounds)


if (my_rank==0) then
  write(*,*) ' ###############################################'
  write(*,*) ' ###    Starting pre dyson energy loop'
  write(*,*) ' ###############################################'
end if
#ifdef CPP_MPI
call mpi_barrier(MPI_COMM_WORLD,ierror)
! print *, 'ierror=',ierror
#endif
if (my_rank==0) then
  write(*,*) ' ###############################################'
end if
if (mpi_size/=1) write(*,*) 'Proc = ',my_rank,'ie= ',mpi_iebounds(1,my_rank), ' .. ',mpi_iebounds(2,my_rank)

do ie=mpi_iebounds(1,my_rank),mpi_iebounds(2,my_rank)
  do ispin=1,nspin-KGREFSOC
    write(1337,*) 'proc = ',my_rank,' IE = ',ie,' ispin= ',ispin

    call preconditioning_readgreenfn(ie,ispin,ielast,lmsizehost,ntotatom,gmathost,'singleprecision')
!                                      in    in         in        in     out
    if(config_testflag('gtest')) write(10000+my_rank,'(832E)') gmathost

    call dysonvirtatom(natom,ntotatom,lmsizehost,gmathost,tmat(:,:,:,ie,ispin),killatom, &
                       isvatom,lmaxatom,gmathostnew,nlmhostnew,lattice_relax,NSOC)

    call preconditioning_writegreenfn(ie,ispin,ielast,gmathostnew,nlmhostnew)

       if(config_testflag('gtest')) write(20000+my_rank,'(832E)') gmathostnew


    write(1337,*) 'proc = ',my_rank,' IE = ',ie,'ispin= ',ispin,'done...'
  end do
end do

!  write header
if (my_rank==0) write(89,rec=1) natom,nlmhostnew,nspin,kgrefsoc

close(88)
close(89)
if(config_testflag('gmat_plain')) close(8989)
deallocate(gmathostnew)

#ifdef CPP_MPI
call mpi_barrier(MPI_COMM_WORLD,ierror)
#endif


call log_write('**********************************************************')
call log_write(' starting intercell preparation..')
call log_write('**********************************************************')
call preconditioning_intercell(natom,ntotatom,lmsizehost,rimpatom,killatom,isvatom,intercell_ach,alat,vmtzero,nsoc)
call log_write('**********************************************************')
call log_write('done preconditioning')
call log_write('**********************************************************')

if (vmtzero(1)>100.0D0 .or. vmtzero(2)>100.0D0) then
  stop '[preconditioning_start] vmtzero too high'
end if

! if (my_rank==0) then
  if(config_testflag('writeout_intercell')) then
    write(1337,*) 'Intercell potential after substraction'
    do iatom=1,natom
      write(1337,'(5000g24.16)') intercell_ach(:,iatom)
      write(1337,*) ' '
    end do
  end if
! end if !my_rank

end subroutine



!  *******************************************************************
!  *******************************************************************
!                                  SUBROUTINES
!  *******************************************************************
!  *******************************************************************



subroutine preconditioning_intercell(natom,ntotatom,lmsizehost,ratom,killatom,isvatom,intercell_ach,alat,vmtzero,nsoc)
use nrtype
use mod_amn2010
use mod_gauntharmonics, only: gauntcoeff
use mod_log, only: log_write
use mod_config, only: config_testflag
use mod_shftvout

implicit none
!interface
integer                    :: natom
integer                    :: ntotatom
integer                    :: lmsizehost
real(kind=DP)              :: ratom(3,ntotatom)
integer                    :: killatom(ntotatom)
integer                    :: isvatom(ntotatom)
real(kind=DP),allocatable  :: intercell_ach(:,:) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
real(kind=DP)              :: alat
real(kind=dp)              :: vmtzero(2)
integer                    :: nsoc

!local
integer                    :: lmaxhost
integer                    :: iatom,jatom,lval,mval,lm,lm2
real(kind=DP),allocatable  :: cmom(:,:) !(lmpotd,ntotatom)
real(kind=DP),allocatable  :: ach(:,:) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
real(kind=DP),allocatable  :: achnew(:,:) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
real(kind=DP),allocatable  :: zatref(:) !(lmpotd,ntotatom),achnew(lmpotd,ntotatom)
real(kind=DP)              :: ac
integer                    :: lmpothost
integer                    :: lpothost
real(kind=DP),allocatable  :: amat(:,:,:)
real(kind=DP),allocatable  :: bmat(:,:)

integer                    :: substract_cmoms(ntotatom,ntotatom)
character(len=*),parameter :: correct_mode='before_substractcmom'

lmaxhost=int(sqrt(dble(lmsizehost/nsoc))-1)

if (nsoc*(lmaxhost+1)**2/=lmsizehost) then
  print *, nsoc,lmaxhost,nsoc*(lmaxhost+1)**2,lmsizehost
  stop 'lmaxhost conversion error'
end if

lpothost=lmaxhost*2
lmpothost=(lpothost+1)**2

if (.not. allocated(gauntcoeff)) stop '[preconditioning_intercell] gauntcoeff not allocated'

allocate( amat(ntotatom,lmpothost,lmpothost), &
          bmat(ntotatom,lmpothost ) ) 


allocate( cmom(lmpothost,ntotatom),   ach(lmpothost,ntotatom), &
          achnew(lmpothost,ntotatom), zatref(ntotatom) )

! read the charge moments of each cell
call log_write('call preconditioning_readcmoms')
call preconditioning_readcmoms(ntotatom,lmpothost,cmom,zatref)


! read the intercell potential of the bulk expanded in lm components: ach
call log_write('call preconditioning_readintercell')
call preconditioning_readintercell(ntotatom,lmpothost,ach,alat,vmtzero)

call log_write('substracting bulk cmoms')

substract_cmoms(:,:)=1
if (correct_mode=='before_substractcmom') then
  call correct_virtualintercell(ntotatom,lmaxhost,lmpothost,ach,ratom,isvatom,substract_cmoms)
elseif (correct_mode=='after_substractcmom') then
else
  stop '[preconditioning] error in correct_mode' 
end if

do iatom=1,ntotatom
  call amn2010(iatom,amat,bmat,alat,gauntcoeff(lmaxhost),lpothost,ntotatom,ratom)
  do lval = 0,lpothost
    do mval = -lval,lval
      lm = lval*lval + lval + mval + 1
      ac = 0.0d0
      if (.not. config_testflag('hostintercell')) then
        if (ntotatom/=1) then
          do jatom = 1,ntotatom

            if (substract_cmoms(iatom,jatom)==1) then
              do lm2 = 1,lmpothost
                ac = ac - amat(jatom,lm,lm2)*cmom(lm2,jatom)
              end do !lm2
              ac = ac - bmat(jatom,lm)*zatref(jatom)
            end if
          end do
        end if  
      end if !(.not. config_testflag('hostintercell'))
      ac = ac + ach(lm,iatom)
      achnew(lm,iatom)=ac
    end do    ! m
  end do       ! l
end do !iatom=1,ntotatom

if (correct_mode=='after_substractcmom') then
  call correct_virtualintercell(ntotatom,lmaxhost,lmpothost,achnew,ratom,isvatom,substract_cmoms)
end if








jatom=0
allocate(intercell_ach(lmpothost,ntotatom)) ! Phivos 2.6.14: before this was: allocate(intercell_ach(lmpothost,natom)), too small in case of "killed" atoms
do iatom=1,ntotatom
!   write(*,'(A,5000F)') 'test',achnew(:,iatom)
  if (killatom(iatom)==0) then
    jatom=jatom+1
    intercell_ach(:,jatom)=achnew(:,iatom)
  end if
end do !iatom

call log_write('substracting bulk cmoms done')
end subroutine preconditioning_intercell !(natotatom)

subroutine preconditioning_readcmoms(ntotatom,lmpothost,cmom,zatref)
use nrtype
use mod_version_info
implicit none
real(kind=DP)        :: cmom(lmpothost,ntotatom)
real(kind=DP)        :: zatref(ntotatom)
character(len=10000)   :: cline
integer              :: lmpothost,ntotatom
integer              :: ios,iatom
open(unit=88, file='kkrflex_intercell_cmoms')
call version_check_header(88)
do iatom=1,ntotatom
  cline=this_readline(88,ios)
  if (ios/=0) stop 'error in preconditioning_readcmoms'
  read(cline,*) zatref(iatom),cmom(:,iatom) 
end do !natom
close(88)
end subroutine preconditioning_readcmoms

subroutine preconditioning_readintercell(ntotatom,lmpothost,ach,alat,vmtzero)
use nrtype
use mod_version_info
implicit none
integer              :: lmpothost,ntotatom
real(kind=DP)        :: ach(lmpothost,ntotatom)
real(kind=DP)        :: alat
! real(kind=DP)        :: zatref(ntotatom)
character(len=10000)   :: cline
integer              :: ios,iatom
integer              :: ntotatomtemp,lmpothosttemp
real(kind=DP)              :: vmtzero(2)

open(unit=88, file='kkrflex_intercell_ref')
call version_check_header(88)
 cline=this_readline(88,ios)

read(cline,*) ntotatomtemp,lmpothosttemp,alat,vmtzero(1),vmtzero(2)
if (abs(vmtzero(2))<10e-10) then
  write(1337,*) 'WARNING : VMTZERO(2) is zero, I will set VMTZERO(2)=VMTZERO(1)'
  VMTZERO(2)=VMTZERO(1)
end if
if (abs(vmtzero(2)-vmtzero(1))>10e-10) stop '[preconditioning_readintercell] vmtzero do not match'
if (ntotatomtemp/=ntotatom) stop '[preconditioning_readintercell] ntotatom error'
if (lmpothosttemp/=lmpothost) stop '[preconditioning_readintercell] lmpothost error'

do iatom=1,ntotatom
  cline=this_readline(88,ios)
  if (ios/=0) stop 'error in preconditioning_readcmoms'
  read(cline,*) ach(:,iatom) 
end do !natom
close(88)
end subroutine preconditioning_readintercell



function this_readline(ifile,ios)
!--------------------------------------------------------
!--  reads the next line in file unit ifile            --
!--------------------------------------------------------
!--  files starting with a dash (#) are treated as     --
!--  a comment !!!                                     --
!--  OUTPUT: next line which is not commented out      --
!--          error variable IOS (should be zero)       --
!--------------------------------------------------------
! input variables
  implicit none
integer,intent(in)               :: ifile
integer,intent(out)              :: ios
! local variables
character(len=10000)  ::this_readline
do
  read(unit=ifile,fmt='(A)', iostat=ios) this_readline
  if (ios/=0 ) exit
  if (this_readline(1:1)/='#' .and. this_readline(2:2)/='#') exit
end do
end function this_readline


   subroutine preconditioning_readenergy(IELAST,NSPIN,EZ,WEZ,NATOMIMPD,NTOTATOM,lmsizehost,kgrefsoc)
     use mod_config, only: config_testflag,config_runflag
     use nrtype, only: wlength
     implicit none
     integer          ::   ielast,nspin
     integer          ::   ielasttemp,nspintemp
     DOUBLE COMPLEX,allocatable   ::   EZ(:),WEZ(:)
!      DOUBLE COMPLEX               ::   EZtemp(1000),WEZtemp(1000)
     integer                      ::   natomimpd,NATOMIMP,ie,recl1,ntotatom,lmsizehost
     integer :: kgrefsoc,kgrefsoctemp,lmsizehosttemp
!  *******************************************************************
!  first open the header of green function file to get information of
!  dimensions
!  *******************************************************************
   OPEN (88,ACCESS='direct',RECL=wlength*2*16,FILE='kkrflex_green',FORM='unformatted')
     read(88,rec=1) ielasttemp,nspintemp,natomimpd,natomimp,lmsizehosttemp,kgrefsoctemp

      if (kgrefsoctemp/=kgrefsoc) then
        print *, '[warning] the value for kgrefsoc does not agree between'
        print *, '          tmat and green function file'
        print *, 'Do you use an old JM code? Maybe you should update your JM code'
        print *, 'Results are probably still correct'
        if (.not. config_runflag('oldJMcode')) then
          print *,'please use runflag "oldJMcode" '
          stop
        end if
      end if

   close(88)
   write(1337,*) '[read_energy] number of energy points ',ielast
   write(1337,*) '[read_energy] number of spins ',nspintemp
   write(1337,*) '[read_energy] number of host atoms ',natomimp
   write(1337,*) '[read_energy] maximum number of host atoms ',natomimpd
   write(1337,*) '[read_energy] host lmsizehost',lmsizehosttemp
   write(1337,*) '[read_energy] Spin orbit coupling used? (1=yes,0=no)',kgrefsoc

  if (ntotatom/=natomimp) then
    write(*,*) ntotatom, natomimp
    stop '[preconditioning_readenergy] ntotatom error'
  endif
  if (natomimpd/=natomimp) then
    write(*,*) natomimpd, natomimp
    stop '[preconditioning_readenergy] ntotatom error'
  endif
  if (ielast/=ielasttemp) then
    write(*,*) ielast, ielasttemp
    stop '[preconditioning_readenergy] ielast error'
  endif
  if (nspin/=nspintemp) then
    write(*,*) nspin, nspintemp
    stop '[preconditioning_readenergy] nspin error'
  endif
  if (lmsizehost/=lmsizehosttemp) then
    write(*,*) lmsizehost, lmsizehosttemp
    stop '[preconditioning_readenergy] lmsizehost error'
  endif

!  *******************************************************************
!  now reopen the file to read the energy weights
!  *******************************************************************
   allocate( ez(ielast), wez(ielast) )
   RECL1=wlength*2*NATOMIMPD*lmsizehost*NATOMIMPD*lmsizehost
   OPEN (88,ACCESS='direct',RECL=RECL1,FILE='kkrflex_green',FORM='unformatted')
   if ( .not. config_runflag('oldJMcode') ) then
     read(88,rec=1) ielasttemp,nspintemp,natomimpd,NATOMIMP,lmsizehosttemp,kgrefsoctemp,(ez(ie),ie=1,ielast),(wez(ie),ie=1,ielast)
   else
     ! old version does not include kgrefsoc value
     read(88,rec=1) ielasttemp,nspintemp,natomimpd,NATOMIMP,lmsizehosttemp,(ez(ie),ie=1,ielast),(wez(ie),ie=1,ielast)
   end if
   write(1337,*) '[read_energy] energies and weights are:'
   write(1337,*) '[read_energy]   energie           weights'
   write(1337,*) '[read_energy]  real   imag      real    imag'
   write(1337,*) '--------------------------------------------'
!    write(*,*) '[read_energy] energie weights are:',wez
   do ie=1,ielast
     write(1337,'(A,I3,2F,A,2F)') '       ',ie,ez(ie),' ',wez(ie)
   end do
!    read(88,rec=1) ielast,nspin,natomimp,lmmaxd,(ez(ie),ie=1,ielast),(wez(ie),ie=1,ielast)
   close(88)
   end subroutine !read_energy!(IE,LMMAXD,NATOMIMP,GMATHOST)




   subroutine preconditioning_readgreenfn(IE,ISPIN,IELAST,lmsizehost,NATOMIMP,GCLUST,CMODE)
     implicit none
     double complex,allocatable   ::  gclust(:,:)
     complex,allocatable   ::  gclustsingle(:,:)
     integer                      ::  natomimp,ie,lmsizehost,ispin,ielast
     integer                      ::  ierror,ngclus,irec
     character(len=*)             ::  cmode
   ngclus=natomimp*lmsizehost
!    write(*,*) natomimp,lmsizehost

   irec = ielast*(ispin-1)+ ie+1
!     write(*,*) 'irec',irec
   if (cmode=='singleprecision') then 
     allocate (gclust(ngclus,ngclus),stat=ierror)
     allocate (gclustsingle(ngclus,ngclus),stat=ierror)
     read(88,rec=irec) gclustsingle

!      write(*,'(50000F)') gclustsingle

     gclust=DCMPLX(gclustsingle)
     deallocate (gclustsingle)
   else if (cmode=='doubleprecision') then
     allocate (gclust(ngclus,ngclus),stat=ierror)
     read(88,rec=irec) gclust
   else
     stop '[preconditioning_readgreenfn] wrong mode'
   end if
   end subroutine !precontitioning_start

   subroutine preconditioning_writegreenfn(IE,ISPIN,IELAST,GCLUST,nlmhostnew)
     use mod_config, only: config_testflag
     implicit none
     double complex             ::  gclust(nlmhostnew,nlmhostnew)
!      complex,allocatable   ::  gclustsingle(:,:)
     integer                      ::  ie,ispin,ielast,nlmhostnew
     integer                      ::  irec
   irec = ielast*(ispin-1)+ ie+1
   write(89,rec=irec) gclust
   if(config_testflag('gmat_plain')) write(8989,'(65000f)') gclust
   end subroutine !precontitioning_start


     subroutine preconditioning_readtmatinfo(NTOTATOM,NSPIN,IELAST,lmsizehost,KGREFSOC)
     use mod_version_info
     implicit none
     integer          ::  NTOTATOM,NATOMTEMP,NSPIN,IELAST,lmsizehost,KGREFSOC
     character(len=5) ::  CHAR1 
     character(len=100) ::  CHAR2 
     character(len=103) ::  CHAR3 
     integer :: ios

     open(unit=6699, file='kkrflex_tmat', status='old', iostat=ios)
     call version_check_header(6699)

     read(unit=6699,fmt='(A)', iostat=ios) CHAR2
     CHAR3=CHAR2 !//' 0' ! older version of the JM do not have the KGREFSOC flag. We, therefore, add a zero 
                    !  to the string to prevent a conversion error
     read(CHAR3,*) CHAR1,NATOMTEMP,NSPIN,IELAST,lmsizehost,KGREFSOC


!      end if
     close(6699)
     write(1337,*) '[preconditioning_readtmatinfo] Information from the T-matrix file'
     write(1337,*) '[preconditioning_readtmatinfo] NTOTATOMIMP',NATOMTEMP
     write(1337,*) '[preconditioning_readtmatinfo] NSPIN',NSPIN
     write(1337,*) '[preconditioning_readtmatinfo] IELAST',IELAST
     write(1337,*) '[preconditioning_readtmatinfo] lmsizehost',lmsizehost
     write(1337,*) '[preconditioning_readtmatinfo] KGREFSOC',KGREFSOC
     if (ntotatom/=natomtemp) stop 'preconditioning_readtmatinfo] natom error in tmat file'

   end subroutine !precontitioning_start

   subroutine preconditioning_readtmat(IELAST,lmsizehost,NATOMIMP,NSPIN,TMAT,KGREFSOC)
     use mod_version_info
     implicit none
     double complex   ::  tmat(lmsizehost,lmsizehost,NATOMIMP,IELAST,NSPIN-KGREFSOC)
     integer          ::  ie,ielast,lmsizehost,natomimp
     integer          ::  NATOMIMPtemp,NSPINtemp,IELASTtemp,lmsizehosttemp
     integer          ::  IATOMtemp,ISPINtemp,IEtemp,temp1
     integer          ::  IATOM,ISPIN,NSPIN,KGREFSOC
     character(len=50) ::  CHAR1 
     OPEN (88,FILE='kkrflex_tmat',STATUS='unknown')
     call version_check_header(88)
     read(88,*)     CHAR1,NATOMIMPtemp,NSPINtemp,IELASTtemp,lmsizehosttemp
     if (natomimp/=natomimptemp)    stop '[preconditioning_readtmat] error'
     if (NSPIN/=NSPINtemp)          stop '[preconditioning_readtmat] error'
     if (IELAST/=IELASTtemp)        stop '[preconditioning_readtmat] error'
     if (lmsizehost/=lmsizehosttemp)  stop '[preconditioning_readtmat] error'
     DO IATOM = 1,NATOMIMP
       DO ISPIN=1,NSPIN-KGREFSOC
         DO IE=1,IELAST
!            write(*,*) iatom,ispin,ie
!            read(88,*)  IATOMtemp,ISPINtemp,IEtemp,temp1,TMAT(:,:,iatom,ie)
           read(88,'(4I12,(50000E25.16))')  IATOMtemp,ISPINtemp,IEtemp,temp1,TMAT(:,:,iatom,ie,ispin)
           if (iatomtemp/= iatom) stop 'preconditioning_readtmat iatom error'
           if (ispin    /= ispin) stop 'preconditioning_readtmat ispin error'
           if (ietemp   /= ie)    stop 'preconditioning_readtmat iettemp error'
         END DO !IE=1,IELAST
       END DO !ISPIN=1,NSPIN
     END DO !IATOM = 1,NATOMIMP
     close(88)
   end subroutine !precontitioning_start

subroutine correct_virtualintercell(ntotatom,lmaxhost,lmpothost,ach,ratom,isvatom,substract_cmoms)
use mod_gauntharmonics, only: gauntcoeff
use mod_shftvout
implicit none
integer               :: ntotatom
integer               :: lmaxhost
integer               :: lmpothost
double precision      :: ach(lmpothost,ntotatom)
double precision      :: ratom(3,ntotatom)
integer               :: isvatom(ntotatom)
! type(gauntcoeff_type) :: gauntcoeff(lmaxhost)
integer               :: substract_cmoms(ntotatom,ntotatom)
!local
integer               :: iatom,jatom
double precision      :: deltar(3)


do iatom=1,ntotatom
  if (isvatom(iatom)==1) then
    do jatom=1,ntotatom
      if (.not. isvatom(jatom)==1) then
        deltar=ratom(:,iatom)-ratom(:,jatom)
        write(*,*) deltar
        if ( sqrt(deltar(1)**2+deltar(2)**2+deltar(3)**2 ) <0.3D0) then
          write(*   ,'(A,I3,A   )') '###################################################################'
          write(*   ,'(A,I3,A   )') 'intercell potential of virt. atom',iatom,' is replaced'
          write(*   ,'(A,I3,A   )') 'by the intercell potential of the real atom',jatom,' and shifted'
          write(*   ,'(A        )') 'by the U-transformation'
          write(*   ,'(A,I3,A   )') '###################################################################'
          write(1337,'(A,I3,A   )') '###################################################################'
          write(1337,'(A,I3,A   )') 'intercell potential of virt. atom',iatom,' is replaced'
          write(1337,'(A,I3,A   )') 'by the intercell potential of the real atom',jatom,' and shifted'
          write(1337,'(A        )') 'by the U-transformation'
          write(1337,'(A,I3,A   )') '###################################################################'
        !  copy intercell_ach(lmpothost,jatom) -> intercell_ach(lmpothost,iatom)
        !  intercell_temp=intercell_ach(lmpothost,iatom)
        !  shift intercell_ach(lmpothost,iatom) by 
        
           if ( substract_cmoms(iatom,jatom)/=0) then
             substract_cmoms(iatom,jatom)=0
           else
             stop '[correct_virtualintercell] some atoms seem to be too close together'
           end if

!          deltar=0
!   write(*,*) ach(:,jatom),ach(:,iatom),deltar, &
!                             2*lmaxhost,gauntcoeff(lmaxhost)%WG,gauntcoeff(lmaxhost)%YRG,(lmaxhost+1)**2, &
!                             4*lmaxhost,2*lmaxhost,(4*lmaxhost+1)**2,gauntcoeff(lmaxhost)%NCLEB
            call SHFTVOUT(ach(:,jatom),ach(:,iatom),deltar, &
                            2*lmaxhost,gauntcoeff(lmaxhost)%WG,gauntcoeff(lmaxhost)%YRG,(lmaxhost+1)**2, &
                            4*lmaxhost,2*lmaxhost,(4*lmaxhost+1)**2,gauntcoeff(lmaxhost)%NCLEB)
        
        ! ach(:,iatom)=ach(:,jatom)
        
        !  intercell_temp=intercell_ach(lmpothost,iatom)
        end if
      end if
    end do
  end if
end do

end subroutine





end module mod_preconditioning
