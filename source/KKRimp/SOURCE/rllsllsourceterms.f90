module mod_rllsllsourceterms

contains

subroutine rllsllsourceterms(nsra,nvec,eryd,rmesh,nrmax,lmax,lmsize,use_fullgmat,jlk_index,hlk,jlk,hlk2,jlk2,GMATPREFACTOR)
use mod_physic_params, only: cvlight
use mod_timing
use mod_beshank
use mod_chebint
use mod_config, only: config_testflag
use mod_rllslltools
use mod_physic_params,only: cvlight
use sourceterms
implicit none
! ************************************************************************
! calculates the source terms J,H and the left solution J2, H2 for:
! - non-relativistic
! - scalar-relativistic
! - full-relativistic
! calculations
! ************************************************************************
double complex,parameter   :: ci=(0.0d0,1.0d0)
integer                    :: nsra
integer                    :: nvec
double complex             :: eryd
double precision           :: rmesh(nrmax)
integer,allocatable        :: jlk_index(:)
integer                    :: l1,lm1,m1,ivec,ispinfullgmat,ir
integer                    :: use_fullgmat
integer                    :: lmsize

double complex             :: ek,ek2,gmatprefactor
double complex,allocatable :: hlk(:,:),  jlk(:,:)
double complex,allocatable :: hlk2(:,:), jlk2(:,:)
integer                    :: lmax
integer                    :: nrmax

if (nsra==2 .or. nsra==3) then 
  nvec=2
elseif (nsra==1) then 
  nvec=1
elseif (nsra==4) then 
  nvec=1
else 
  stop '[rllsllsourceterms] error'
end if

allocate ( jlk_index(nvec*lmsize) )
if (nsra<=2 .or. nsra==4) then
  allocate ( hlk(1:(lmax+1)*nvec,nrmax ) )
  allocate ( jlk(1:(lmax+1)*nvec,nrmax ) )
  allocate ( hlk2(1:(lmax+1)*nvec,nrmax ) )
  allocate ( jlk2(1:(lmax+1)*nvec,nrmax ) )
elseif (nsra==3) then
  allocate(  hlk (2*lmsize,nrmax),&
             jlk (2*lmsize,nrmax),&
             hlk2(2*lmsize,nrmax),&
             jlk2(2*lmsize,nrmax))
else 
        stop '[rllsll] nsra not known'
end if


if (nsra<=2 .or. nsra==4) then 
  lm1 = 1
  do ivec=1,nvec
    do ispinfullgmat=0,use_fullgmat
      do l1 = 0,lmax
        do m1 = -l1,l1
          jlk_index(lm1) = l1+(ivec-1)*(lmax+1)+1
          lm1 = lm1 + 1
        end do   
      end do  
    end do!ispinorbit=0,use_fullgmat
  end do !nvec
else if (nsra==3) then
  do lm1=1,lmsize*nvec
    jlk_index(lm1) = lm1
  end do !nvec
end if

if (nsra==1) then 
  ek = sqrt(eryd)
  ek2 = sqrt(eryd)
elseif (nsra==2) then
  ek = sqrt(eryd+(eryd/cvlight)**2)
  ek2 = sqrt(eryd+(eryd/cvlight)**2) *(1.0d0+eryd/cvlight**2)
elseif (nsra==4) then
  ek = sqrt(eryd+(eryd/cvlight)**2)
  ek2 = sqrt(eryd+(eryd/cvlight)**2) *(1.0d0+eryd/cvlight**2)
elseif (nsra==3) then
  ek =  sqrt(eryd+(eryd/cvlight)**2)
  ek2 = sqrt(eryd+(eryd/cvlight)**2)
else
  stop'[rllsll] wrong value for nvec'
end if
write(1337,*) '**********************************************'
write(1337,*) '  rllsllsourceterms'
write(1337,*) '**********************************************'
write(1337,*) 'ek',ek
write(1337,*) 'ek2',ek2
write(1337,*) '**********************************************'


do ir = 1,nrmax
  if (nsra<=2 .or. nsra==4) then ! SRA or non-relatistic

    call beshank(hlk(:,ir),jlk(:,ir),ek*rmesh(ir),lmax)

    if (nsra==2) then
      call beshank_smallcomp(hlk(:,ir),jlk(:,ir),&
                        ek*rmesh(ir),rmesh(ir),eryd,lmax)
    end if

    do l1 = 1,nvec*(lmax+1)
      hlk(l1,ir) = -ci*hlk(l1,ir)
    end do

    if (nsra==1 .or. nsra==4) then
      do l1 = 1,nvec*(lmax+1)
        jlk2(l1,ir) = jlk(l1,ir)
        hlk2(l1,ir) = hlk(l1,ir)
      end do
    else if (nsra==2) then
      do l1 = 1,lmax+1
        jlk2(l1,ir) = jlk(l1,ir)
        hlk2(l1,ir) = hlk(l1,ir)
      end do
      do l1 = lmax+2,2*(lmax+1)
        jlk2(l1,ir) = -jlk(l1,ir)
        hlk2(l1,ir) = -hlk(l1,ir)
      end do
    end if

  else if (nsra==3) then ! Dirac

  call sourcetermsupervector(lmax,eryd,rmesh(ir),jlk(:,ir), &
        hlk(:,ir),jlk2(:,ir),hlk2(:,ir))
  end if
end do

gmatprefactor=ek2 ! prefacor for the Green function: 
                  ! non-relativistic  = kappa
                  ! sra = M_0 * kappa (check Bauer, PhD, section 4.3)

end subroutine rllsllsourceterms

end module mod_rllsllsourceterms
