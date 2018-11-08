module mod_calcsph

contains

subroutine calcsph(nsra,cellnew,zatom,use_fullgmat,nspin,ispin,lmax,eryd, &
                   jlk_index,hlk,jlk,hlk2,jlk2,GMATPREFACTOR,gauntcoeff,tmat,idotime )
!interface
use type_cellnew
use mod_vllmat
use mod_vllmatsra
use mod_rllsll
use type_gauntcoeff
use mod_config, only: config_testflag
use mod_basistransform

implicit none

integer                                   :: nsra
type(cell_typenew)                        :: cellnew
double precision                          :: zatom
integer                                   :: use_fullgmat
integer                                   :: nspin
integer                                   :: ispin
integer                                   :: lmax
double complex                            :: eryd
integer                                   :: jlk_index(:)
double complex,allocatable                :: hlk(:,:),jlk(:,:),hlk2(:,:),jlk2(:,:)
double complex                            :: gmatprefactor
!local
integer                                   :: lmsize,lmsize2,lmpot,nspintemp,jspin,jspin2
integer                                   :: lmaxatom,lval,ir,lshift_sra,lshift_spin
double complex,allocatable                :: vpotll(:,:,:)
type(gauntcoeff_type)                     :: gauntcoeff
double complex,allocatable                :: hlktemp(:,:),jlktemp(:,:),hlk2temp(:,:),jlk2temp(:,:)
integer,allocatable                       :: jlk_indextemp(:)
double complex,allocatable                :: RLLtemp(:,:,:,:),SLLtemp(:,:,:,:)
double complex,allocatable                :: JLKnew(:,:),HLKnew(:,:)
integer                                   :: nvec,idotime
integer                                   :: ivec,ispinfullgmat,l1,m1,lm1
double complex,allocatable                :: tmattemp(:,:)
double complex,allocatable                :: tmat(:)

!#######################################################
! Setting up dimension variables. Since in this routine
! just the sperical part of the potential is used lmsize
! is eather 1 or 2 (for SRA).
!#######################################################
write(1337,*)  'Starting spherical calculation'

lmsize =1
if (nsra==2 .or. nsra==4) then
  lmsize2=2
  nvec=2
elseif (nsra==1) then
  lmsize2=1
  nvec=1
else
  stop '[calcsph] nsra error'
end if

lmaxatom=0  ! lmax used for the expansion of the green function
lmpot=ubound(cellnew%Vpotnew(:,:,:),2)

if (use_fullgmat==1) then 
  nspintemp=2
else
  nspintemp=1
end if

!#######################################################
! allocation of local arrays and setting up pointers
!#######################################################
allocate(RLLtemp(lmsize2,lmsize,cellnew%nrmaxnew,1), &
         SLLtemp(lmsize2,lmsize,cellnew%nrmaxnew,1))
allocate(hlktemp(nvec,cellnew%nrmaxnew),jlktemp(nvec,cellnew%nrmaxnew),&
         hlk2temp(nvec,cellnew%nrmaxnew),jlk2temp(nvec,cellnew%nrmaxnew))
allocate(jlk_indextemp(lmsize2))
allocate(tmattemp(lmsize,lmsize))

if (.not. allocated(tmat))   allocate( tmat( nvec*(lmax+1) ) )

do ivec=1,nvec
  jlk_indextemp(ivec)=ivec
end do

allocate(JLKnew(nvec*(use_fullgmat+1)*(lmax+1),cellnew%nrmaxnew),&
          HLKnew(nvec*(use_fullgmat+1)*(lmax+1),cellnew%nrmaxnew) )

!#######################################################
! depending on use_fullgmat eather one or two spin components
! are calculated
!#######################################################
do jspin=1,nspintemp
  if (use_fullgmat==1) then
    jspin2=jspin
  else
    jspin2=ispin
  end if
  lshift_spin=(lmax+1)*(jspin-1)
  lshift_sra=(lmax+1)*nspintemp

  !#######################################################
  ! Now, for each value of l the Lippmann-Schwinger equation 
  ! is solved using the free-potential wavefunctions and potentials
  ! corresponding to l-value.
  !#######################################################
  do lval=0,lmax

    allocate(Vpotll(lmsize,lmsize,cellnew%nrmaxnew))

    call VLLMAT(Vpotll,cellnew%Vpotnew(:,:,:),LMAXATOM,(LMAXATOM+1)**2,LMPOT,1, &
                cellnew%nrmaxnew,gauntcoeff,zatom,cellnew%rmeshnew,lmsize,0 &! changed from use_fullgmat 
                ,NSPIN,JSPIN2,'SPH' )

    if (nsra==2 .or. nsra==4 ) call vllmatsra(Vpotll,cellnew%rmeshnew,eryd,lmaxatom,lval,'Ref=0')



    jlktemp (1,:)=jlk (lval+1,:)
    hlktemp (1,:)=hlk (lval+1,:)
    jlk2temp(1,:)=jlk2(lval+1,:)
    hlk2temp(1,:)=hlk2(lval+1,:)

    if (nsra==2 .or. nsra==4) then
      jlktemp (2,:)=jlk (lmax+lval+2,:)
      hlktemp (2,:)=hlk (lmax+lval+2,:)
      jlk2temp(2,:)=jlk2(lmax+lval+2,:)
      hlk2temp(2,:)=hlk2(lmax+lval+2,:)
    end if
    if (nsra==4) then
      stop 'does not work up here'
    end if



    call RLLSLL(cellnew%rpan_intervall,cellnew%rmeshnew,VPOTLL,&
                RLLtemp(:,:,:,1),SLLtemp(:,:,:,1),tmattemp, &
                cellnew%ncheb,cellnew%npan_tot,lmsize,lmsize2,cellnew%nrmaxnew,nvec, &
                jlk_indextemp,hlktemp,jlktemp,hlk2temp,jlk2temp,GMATPREFACTOR,'1','1','0',idotime)

    do ir=1,cellnew%nrmaxnew
      JLKnew(lshift_spin+lval+1,           ir)=RLLtemp(1,1,ir,1)/cellnew%rmeshnew(ir)
      HLKnew(lshift_spin+lval+1,           ir)=SLLtemp(1,1,ir,1)/cellnew%rmeshnew(ir)
    end do!ir

    if (nsra==2) then
      do ir=1,cellnew%nrmaxnew
        JLKnew(lshift_spin+lshift_sra+lval+1,ir)=RLLtemp(2,1,ir,1)/cellnew%rmeshnew(ir)
        HLKnew(lshift_spin+lshift_sra+lval+1,ir)=SLLtemp(2,1,ir,1)/cellnew%rmeshnew(ir)
      end do!ir
    end if

    tmat(lshift_spin+lval+1)=tmattemp(1,1)
    deallocate(Vpotll)

  end do !lval

end do !jspin

deallocate(jlk,hlk,jlk2,hlk2)

allocate(JLK (nvec*(use_fullgmat+1)*(lmax+1),cellnew%nrmaxnew),&
           HLK (nvec*(use_fullgmat+1)*(lmax+1),cellnew%nrmaxnew),&
           JLK2(nvec*(use_fullgmat+1)*(lmax+1),cellnew%nrmaxnew),&
           HLK2(nvec*(use_fullgmat+1)*(lmax+1),cellnew%nrmaxnew) )
jlk =(0.0D0,0.0D0)
hlk =(0.0D0,0.0D0)
jlk2=(0.0D0,0.0D0)
hlk2=(0.0D0,0.0D0)


lm1 = 1
do ivec=1,nvec
  do ispinfullgmat=0,use_fullgmat
    do l1 = 0,lmax
      do m1 = -l1,l1
        jlk_index(lm1) = l1+(ivec-1)*nspintemp*(lmax+1)+ispinfullgmat*(lmax+1)+1
        lm1 = lm1 + 1
      end do   
    end do  
  end do!ispinorbit=0,use_fullgmat
end do !nvec

do ir=1,cellnew%nrmaxnew
  do l1 = 1,nvec*(lmax+1)*nspintemp
    jlk(l1,ir) =jlknew(l1,ir)
    hlk(l1,ir) =hlknew(l1,ir)
  end do
end do

if (nsra==2) then
  do ir=1,cellnew%nrmaxnew
    do l1 = 1,(lmax+1)*nspintemp
      jlk2(l1,ir) = -jlknew(l1+lmax+1,ir)
      hlk2(l1,ir) = -hlknew(l1+lmax+1,ir)
    end do
    do l1 = nspintemp*(lmax+1)+1,2*(lmax+1)*nspintemp
      jlk2(l1,ir) = jlknew(l1-(lmax+1)*nspintemp,ir)
      hlk2(l1,ir) = hlknew(l1-(lmax+1)*nspintemp,ir)
    end do
  end do !ir
elseif (nsra==1 .or. nsra==4) then
  do ir=1,cellnew%nrmaxnew
    do l1 = 1,nvec*(lmax+1)*nspintemp
      jlk2(l1,ir) =jlknew(l1,ir)
      hlk2(l1,ir) =hlknew(l1,ir)
    end do
  end do
end if

end subroutine calcsph

end module mod_calcsph
