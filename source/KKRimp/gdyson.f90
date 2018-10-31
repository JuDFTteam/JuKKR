module mod_gdyson
complex(8),allocatable,private       ::  Gbulk_storage(:,:,:,:)
contains

subroutine gdyson(igmatnewfile,ie,ispin,nspin,natom,lmaxatom,tmatll,use_fullgmat,gmat,gmatonsite,ielast,mpi_iebounds,ITSCF,saveGmat)
use mod_mathtools, only: linearsolve_dc
use type_tmat
use type_gmat
use type_gmatonsite
use mod_config, only: config_testflag

implicit none
!interface
integer                                :: ie
integer                                :: ispin
integer                                :: nspin
integer                                :: natom
integer                                :: lmaxatom(natom)
type(TMAT_TYPE),allocatable            :: TMATLL(:,:) !(LMMAXD,LMMAXD)
integer                                :: use_fullgmat
type(gmat_type)                        :: gmat
type(gmatonsite_type),allocatable       :: gmatonsite(:,:)
integer                                :: ielast
integer                                  :: mpi_iebounds(2)
integer                                :: ITSCF
integer                                :: saveGmat

!local
integer                                :: iatom,ilm
double complex,parameter               :: cone=(1.0D0,0.0D0)
double complex,parameter               :: czero=(0.0D0,0.0D0)
double complex,allocatable             ::  gref1(:,:)

double complex,allocatable                         ::  tmat_big(:,:)
double complex,allocatable                         ::  gref(:,:)
double complex,allocatable                         ::  gtmat(:,:)
! double complex                         ::  tmat_big(gmat%gmatdim,gmat%gmatdim)
! double complex                         ::  gref(gmat%gmatdim,gmat%gmatdim)
! double complex                         ::  gtmat(gmat%gmatdim,gmat%gmatdim) 

integer                                ::   nlmindex1,nlmindex2
integer                                ::  tmatubound,tmatlbound
integer                                :: igmatnewfile


allocate (tmat_big(gmat%gmatdim,gmat%gmatdim),gref(gmat%gmatdim,gmat%gmatdim),gtmat(gmat%gmatdim,gmat%gmatdim)) 

allocate (gref1(gmat%gmathostdim,gmat%gmathostdim))

call gdyson_readgmat(igmatnewfile,use_fullgmat,ielast,ie,Gref,gref1,mpi_iebounds,ITSCF,gmat%gmathostdim,gmat,ispin)


! open(32492345,file='test_gdyson_gref')
! write(32492345,'(5000E)') Gref

! #############################
! set up T-matrix
! #############################
if (.not. config_testflag('nodyson')) then

  tmat_big=CZERO
  do iatom=1,natom
    nlmindex1=gmat%iatom2nlmindex(1,iatom)
    nlmindex2=gmat%iatom2nlmindex(3,iatom)
    tmatubound=ubound(tmatll(iatom,ispin)%tmat,1)
    tmatlbound=lbound(tmatll(iatom,ispin)%tmat,1)
    if (nlmindex2-nlmindex1/=tmatubound-tmatlbound) then
      write(*,*) '[gdyson] error in tmat dimension'
      stop
    end if
    tmat_big(nlmindex1:nlmindex2,nlmindex1:nlmindex2) = -tmatll(iatom,ispin)%tmat
  end do !iatom

  if (config_testflag('write_gdyson')) then
    open(243234,file='test_gdyson1')
    write(243234,'(50000E)') tmat_big
    open(243235,file='test_gdyson2')
    write(243235,'(50000E)') Gref
  end if

  ! #############################
  ! calculate Gref*T
  ! #############################
  
  call matmat_zmzm(Gref,tmat_big,GTMAT)
  ! CALL ZGEMM('N','N',gmat%gmatdim,gmat%gmatdim,gmat%gmatdim,-CONE,GRef,gmat%gmatdim,tmat_big, &
  !            gmat%gmatdim,CZERO,GTMAT,gmat%gmatdim)
  if (config_testflag('write_gdyson')) then
    open(243236,file='test_gdyson3')
    write(243236,'(50000E)') GTMAT
  end if
  ! #############################
  ! calculate 1- Gref*T
  ! #############################
  do ilm=1,gmat%gmatdim
    GTMAT(ilm,ilm)=GTMAT(ilm,ilm)+CONE
  end do
  
  if (config_testflag('write_gdyson')) then
    open(243237,file='test_gdyson4')
    write(243237,'(50000E)') GTMAT
  end if
  ! #############################
  ! calculate G=(1-Gref*T)**-1*Gref
  ! #############################
  call linearsolve_dc(GTMAT,Gref)
  ! the Greensfunction is now stored in Gref
  if (config_testflag('write_gdyson')) then
    open(243238,file='test_gdyson5')
    write(243238,'(50000E)') Gref
  end if
end if !(.not. config_testopt('nodyson')) then

! #############################
! save onsite Greensfunction
! #############################
do iatom=1,natom
  nlmindex1=gmat%iatom2nlmindex(1,iatom)
  nlmindex2=gmat%iatom2nlmindex(3,iatom)
!   write(*,*) nlmindex1,nlmindex2
  gmatonsite(iatom,ispin)%gmat=Gref(nlmindex1:nlmindex2,nlmindex1:nlmindex2)
end do !iatom

if (saveGmat==1 .or. config_testflag('gtest')) then 
  if (.not. allocated (gmat%gmat)) then
    allocate( gmat%gmat(gmat%gmatdim,gmat%gmatdim) )
  end if
  gmat%gmat=Gref!(nlmindex1:nlmindex2,nlmindex1:nlmindex2)
end if


end subroutine

! subroutine read_gmat()
! 
! end subroutine !read_gmat()

      subroutine matmat_zmzm(mat1,mat2,matout)
      implicit none
      complex(8), intent(in) :: mat1(:,:),mat2(:,:)
      complex(8)             :: matout(size(mat1,1),size(mat2,2))
      integer                :: n1,n,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,2)
      if(size(mat2,1).ne.n) stop 'matmat_zmzm: dimensions of matrices are inconsistent.'
      call zgemm('N','N',n1,n2,n,(1d0,0d0),mat1,n1,mat2,n,(0d0,0d0),matout,n1)
      end subroutine matmat_zmzm

! call gdyson_readgmat      (use_fullgmat,ielast,ie,Gref,gref1,mpi_iebounds,ITSCF)

subroutine gdyson_readgmat(igmatnewfile,use_fullgmat,ielast,ie,Gref,greftemp,mpi_iebounds,ITSCF,nlmhost,gmat,ispin)
use mod_config, only: config_runflag
use type_gmat
  implicit none
integer :: igmatnewfile
integer :: use_fullgmat
integer :: ielast
integer :: ie
type(gmat_type)                        :: gmat
! integer                                :: gmatdim
integer :: ispin
complex(8) :: Gref(gmat%gmatdim,gmat%gmatdim)
complex(8) :: greftemp(nlmhost,nlmhost)
integer                                  :: mpi_iebounds(2)
integer                                :: ITSCF
integer                                :: nlmhost
!local
integer :: irec
integer :: temp1,temp2
integer,save :: nspinhost
integer,save :: first=1
integer,save :: kgrefsoc
integer :: iatom1,iatom2,natom

if (first==1) then
  read(igmatnewfile,rec=1) temp1,temp2,nspinhost,kgrefsoc
  if (kgrefsoc/=1 .and. kgrefsoc/=0) then
    stop '[gdyson_readgmat] Gref SOC value not read correctly. Old Juelich-Muenchen code?'
  end if
  write(1337,*) 'Read Gmat bulk header'
  write(1337,*) 'NSPINBULK ==',nspinhost
  write(1337,*) 'kgrefsoc ==',kgrefsoc
end if

if (config_runflag('GBULKtomemory') .and. first==1) then
  allocate(Gbulk_storage(nlmhost,nlmhost,nspinhost-kgrefsoc,mpi_iebounds(1):mpi_iebounds(2)))
end if


Gref = (0.0D0,0.0D0)
if (use_fullgmat==1) then
  if (config_runflag('GBULKtomemory')) then
    stop 'GBULKtomemory not implemeted for use_fullgmat==1'
  else
    if      (kgrefsoc == 1) then 
      irec =                   ie+1
      read(igmatnewfile,rec=irec) Gref
    else if (kgrefsoc == 0) then

      irec = ielast*(1-1)+ ie+1
      read(igmatnewfile,rec=irec) greftemp
      natom = ubound(gmat%iatom2nlmindex,2)
  !      write(*,*) natom
      do iatom1=1,natom
        do iatom2=1,natom
          Gref(gmat%iatom2nlmindex(1,iatom1):gmat%iatom2nlmindex(2,iatom1),    &
              gmat%iatom2nlmindex(1,iatom2):gmat%iatom2nlmindex(2,iatom2)) =  &
          greftemp(gmat%iatom2nlmindexhost(1,iatom1):gmat%iatom2nlmindexhost(2,iatom1), &
              gmat%iatom2nlmindexhost(1,iatom2):gmat%iatom2nlmindexhost(2,iatom2))
        end do !iatom2=1,natom
      end do !iatom1=1,natom
      if (nspinhost==2) then
        irec = ielast*(2-1)+ ie+1
        read(igmatnewfile,rec=irec) greftemp
      end if
      do iatom1=1,natom
        do iatom2=1,natom
          Gref(gmat%iatom2nlmindex(2,iatom1)+1:gmat%iatom2nlmindex(3,iatom1),    &
              gmat%iatom2nlmindex(2,iatom2)+1:gmat%iatom2nlmindex(3,iatom2)) =  &
          greftemp(gmat%iatom2nlmindexhost(1,iatom1):gmat%iatom2nlmindexhost(2,iatom1), &
              gmat%iatom2nlmindexhost(1,iatom2):gmat%iatom2nlmindexhost(2,iatom2))
        end do !iatom2=1,natom
      end do !iatom1=1,natom

    end if

  end if
else
  if      (kgrefsoc == 1) then 
    stop 'noSOC calculation with a SOC reference Green function not possible'
  end if

  if (config_runflag('GBULKtomemory')) then
    if (ITSCF==1) then
      if (nspinhost==2) then
        irec = ielast*(ispin-1)+ ie+1
        read(igmatnewfile,rec=irec) Gref
        Gbulk_storage(:,:,ispin,IE)=Gref
      else
        irec =                   ie+1
        read(igmatnewfile,rec=irec) Gref
        Gbulk_storage(:,:,1,IE)=Gref
      end if
    else
      if (nspinhost==2) then
        Gref=Gbulk_storage(:,:,ispin,IE)
      else
        Gref=Gbulk_storage(:,:,1,IE)
      end if
    end if
  else
    ! read gmat for spin=ispin
    if (nspinhost==2) then 
      irec = ielast*(ispin-1)+ ie+1
!       write(*,*) irec 
      read(igmatnewfile,rec=irec) Gref
    else
      irec =                   ie+1
!       write(*,*) irec 
      read(igmatnewfile,rec=irec) Gref
    end if
  end if
end if 
first=0
end subroutine gdyson_readgmat


subroutine gdyson_read_kgrefsoc(kgrefsoc)
use nrtype, only: wlength
implicit none
 integer :: temp1,temp2,nspinhost,kgrefsoc
  open (13234,access='direct',recl=wlength*4,file='kkrflex_greennew',form='unformatted')
  read(13234,rec=1) temp1,temp2,nspinhost,kgrefsoc
  if (kgrefsoc/=1 .and. kgrefsoc/=0) then
    stop '[gdyson_read_kgrefsoc] Gref SOC value not read correctly. Old Juelich-Muenchen code?'
  end if
  write(1337,*) '[gdyson_read_kgrefsoc] kgrefsoc ==',kgrefsoc
  close(13234)

end subroutine gdyson_read_kgrefsoc



end module
