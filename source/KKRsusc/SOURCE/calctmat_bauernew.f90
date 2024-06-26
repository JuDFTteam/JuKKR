module mod_calctmat_bauernew
contains
subroutine calctmat_bauernew(cell,tmat,lmaxatom,eryd_in,ZATOM,cellnew,wavefunction, &
                             ispin,nspin,kspinorbit,use_fullgmat,theta,phi,ncoll,nsra,config,idotime, &
                             ie,ldau,iatom,inpsusc)        ! lda+u

use mod_gauntharmonics, only: gauntcoeff
use mod_timing
use type_tmat
use type_cell
use type_cellnew
use type_wavefunction
use type_config
use type_ldau                                   ! lda+u
use mod_interpolpot
use mod_vllmat
use mod_vllmatsra
use mod_rllsll
use mod_calccouplingconstants, only: calccouplingdeltat
use mod_config, only: config_testflag
use mod_spinorbit
use mod_rllsllsourceterms
use mod_physic_params, only: cvlight
use Potential
use mod_calcsph
use mod_checknan
use mod_basistransform
use mod_mathtools
use mod_calctmat_bauernew_testtools
use mod_wronskian

! Connection to KKRSUSC                                ! susc |||| added by Juba (2015)
use type_inpsusc                                       ! susc |||| added by Juba (2015)
use mod_config, only: config_testflag, config_runflag  ! susc |||| config_runflag added by Juba (2015)
use mod_cheb2oldgridc                                  ! susc |||| added by Juba (2015)
 
implicit none
!interface
type(cell_type),intent(in)                :: cell
type(tmat_type)                           :: tmat
integer                                   :: lmaxatom
double complex                            :: eryd_in
double precision                          :: zatom
type(cell_typenew)                        :: cellnew
type(wavefunction_type)                   :: wavefunction
integer                                   :: ispin,nspin
integer                                   :: kspinorbit,use_fullgmat
double precision                          :: theta,phi
integer                                   :: ncoll
integer                                   :: nsra
type(config_type)                         :: config
integer                                   :: ie                         ! lda+u
type(ldau_type)                           :: ldau                       ! lda+u variables
!local
double complex,allocatable                ::  vll2ddr(:,:)
double complex,allocatable                ::  vll2ddr2(:,:)
double complex,allocatable                ::  vpotll(:,:,:)
double complex,allocatable                ::  vpotll2(:,:,:)
integer                                   :: ir,lm1,lm2,lmpot,lmmax,ipan,ipan2,lmsize,lmsize2
integer                                   :: irmin,irmax,irminnew,irmaxnew
double precision                          :: rmin,rmax,rval
double precision,allocatable              :: c1(:,:)
double complex,allocatable                :: srllp(:,:), &
                                             ull(:,:,:),sll(:,:,:),&
                                             rll(:,:,:)
double complex                            :: gmatprefactor
double complex,allocatable                :: hlk(:,:),jlk(:,:),hlk2(:,:),jlk2(:,:)
integer,allocatable                       :: jlk_index(:)

integer                                   :: ierror,idotime

complex(kind=dpc),allocatable             :: tmattemp(:,:)
complex(kind=dpc),allocatable             :: tmatsph(:)
double complex                            :: eryd
integer                                   :: idoldau,lmlo,lmhi,mmax,m1,imt1  ! lda+u 
! kkrsusc
integer,intent(in)                        :: iatom                           ! susc |||| added by Juba (2015) 
type(inpsusc_type),intent(inout)          :: inpsusc                         ! susc |||| added by Juba (2015)       
integer                                   :: nrmaxd,ilm,ilm2,isn,i           ! susc |||| added by Juba (2015) 
integer,allocatable                       :: i2lm(:,:)
! use_spinorbit = 0
complex(kind=dpc),allocatable             :: work1_in(:,:), work2_in(:,:)
complex(kind=dpc),allocatable             :: work1_out(:,:), work2_out(:,:)
! use_spinorbit = 1
complex(kind=dpc),allocatable             :: work1r_in(:,:), work2r_in(:,:), work3r_in(:,:), work4r_in(:,:)
complex(kind=dpc),allocatable             :: work1l_in(:,:), work2l_in(:,:), work3l_in(:,:), work4l_in(:,:) 
complex(kind=dpc),allocatable             :: work1r_out(:,:), work2r_out(:,:), work3r_out(:,:), work4r_out(:,:)
complex(kind=dpc),allocatable             :: work1l_out(:,:), work2l_out(:,:), work3l_out(:,:), work4l_out(:,:)

write(1337,*) 'starting calctmatnew'

eryd=(eryd_in)

lmmax=(lmaxatom+1)**2
lmpot=(2*lmaxatom+1)**2

!#######################################################
! set the energy in case of a non-rel or SRA calculation
!#######################################################
if (nsra==1) then                   ! non-relativistc calculation
  wavefunction%nvec=1               !   wave function just has one component
elseif (nsra==2 .or. nsra==3) then  ! sra or full-relativistic
  wavefunction%nvec=2               !   spinor with 2 components
elseif (nsra==4) then               ! test option (might be deleted)
  wavefunction%nvec=1               
end if


!#######################################################
! set the size of the t-matrix and the wavefunctions
! in case of a spin-orbit calculation the tmatrix is twice as big!
!#######################################################
if (use_fullgmat==1) then           ! use_fullgmat means we treat a GF with spin up/down components.
  lmsize=2*lmmax                    ! Thus, we need to multiply by 2
else
  lmsize=lmmax
end if
lmsize2=wavefunction%nvec*lmsize    ! lmsize2 is a combined index of (nvec, lmsize). A factor of 2
wavefunction%lmsize=lmsize          ! is included in case of sra or full-relativistic calculation
wavefunction%lmsize2=lmsize2

allocate(tmattemp(lmsize,lmsize))

if (ubound(tmat%tmat,1)/=lmsize) stop 'calctmat: error in tmat dim'

if ( config_testflag('tmatdebug') ) then
  do ir=1,cellnew%nrmaxnew
    write(1001,'(50000E)') cellnew%rmeshnew(ir),(cellnew%vpotnew(ir,lm1,ispin),lm1=1,lmpot)
  end do
  write(5000,'(50000F)') cellnew%rmeshnew
  write(5001,'(50000F)') cell%rmesh
end if

!#######################################################
! set up the calculation of the VLL matrix
!  according to formula 4.12 of Bauer,PhD Thesis
!#######################################################
IF (NSRA<=2) then ! non-relativistic case
  allocate(Vpotll(lmsize,lmsize,cellnew%nrmaxnew))
  if ( .not. config_testflag('sph') .or. nsra==5 ) then ! set up VLL by omitting the spherical
                                                        ! potential
    call VLLMAT(Vpotll,cellnew%Vpotnew(:,:,:),LMAXATOM,(LMAXATOM+1)**2,LMPOT,1, &
                cellnew%nrmaxnew,gauntcoeff(lmaxatom),zatom,cellnew%rmeshnew,lmsize,use_fullgmat,NSPIN,ISPIN,'NS')
  else                                                  ! use all components of Vpot

    call VLLMAT(Vpotll,cellnew%Vpotnew(:,:,:),LMAXATOM,(LMAXATOM+1)**2,LMPOT,1, &
                cellnew%nrmaxnew,gauntcoeff(lmaxatom),zatom,cellnew%rmeshnew,lmsize,use_fullgmat,NSPIN,ISPIN,'NOSPH')

  end if
else if (nsra==3) then ! use the Dirac solver by Pascal Kordt, PhD thesis
  allocate(Vpotll(2*lmsize,2*lmsize,cellnew%nrmaxnew))
  call PotentialMatrixArray(lmaxatom,lmaxatom,zatom,cellnew%rmeshnew,cellnew%nrmaxnew,eryd,cellnew%vpotnew,Vpotll) ! cellnew%vpotnew(ir,lm1,ispin)
end if

!----------------------------------------------------------------------------------- ! lda+u
! Add wldau to vpotll                                                                ! lda+u
imt1 = cellnew%ipan_intervall(cellnew%npan_log+cellnew%npan_eq) + 1                  ! lda+u
if (ldau%lopt.ge.0.and.ie.ge.ldau%ieldaustart.and.ie.le.ldau%ieldauend) then         ! lda+u
   if (nsra==3) stop 'lda+u not implemented for nsra=3'                              ! lda+u

   lmlo = ldau%lopt**2 + 1                                                           ! lda+u
   lmhi = (ldau%lopt + 1)**2                                                         ! lda+u

 
   if (use_fullgmat.eq.1) then    ! 2x2 in spin space                                ! lda+u
      
      do ir = 1,imt1                                                                 ! lda+u
         vpotll(lmlo:lmhi,lmlo:lmhi,ir) = vpotll(lmlo:lmhi,lmlo:lmhi,ir) +     &     ! lda+u
                                          ldau%wldau(1:mmax,1:mmax,1)                ! lda+u
      enddo                                                                          ! lda+u
      lmlo = lmlo + lmmax                                                            ! lda+u
      lmhi = lmhi + lmmax                                                            ! lda+u
      do ir = 1,imt1                                                                 ! lda+u
         vpotll(lmlo:lmhi,lmlo:lmhi,ir) = vpotll(lmlo:lmhi,lmlo:lmhi,ir) +     &     ! lda+u
                                                  ldau%wldau(1:mmax,1:mmax,2)        ! lda+u
      enddo                                                                          ! lda+u

   else       ! 1x1 in spin space                                                    ! lda+u 

      do ir = 1,imt1
         vpotll(lmlo:lmhi,lmlo:lmhi,ir) = vpotll(lmlo:lmhi,lmlo:lmhi,ir) +     &     ! lda+u
                                          ldau%wldau(1:mmax,1:mmax,ispin)            ! lda+u
      enddo                                                                          ! lda+u
   endif  

endif                                                                                ! lda+u
!----------------------------------------------------------------------------------- ! lda+u


if ( config_testflag('vlldebug') ) then
  do ir=1,cellnew%nrmaxnew
      write(1233,'(50000E)') Vpotll(:,:,ir)
      if (kspinorbit==1) write(11233,'(50000E)') Vpotll2(:,:,ir)
  end do
end if


!#######################################################
! add the spin-orbit Hamiltonian to the VLL matrix
!#######################################################
! in case of spin-orbit coupling V_LL ist not any more symmetric in L-space. Thus,
! the left- and right solutions need to be calculated explicitly. We transpose the
! potential in L-space in order to calculate the the left solution.
if (kspinorbit==1) then

  allocate(Vpotll2(lmsize,lmsize,cellnew%nrmaxnew))
  Vpotll2=Vpotll

  if (cellnew%use_spinorbit==1) then

    call spinorbit(lmaxatom,zatom,       eryd,cellnew,cellnew%nrmaxnew,nspin,Vpotll,theta,phi,ncoll,'1')

    call spinorbit(lmaxatom,zatom,       eryd,cellnew,cellnew%nrmaxnew,nspin,Vpotll2,theta,phi,ncoll,'transpose')

  end if

end if

if ( config_testflag('vlldebug') ) then
  do ir=1,cellnew%nrmaxnew
      write(1234,'(50000E)') Vpotll(:,:,ir)
      if (kspinorbit==1) write(11234,'(50000E)') Vpotll2(:,:,ir)
  end do
end if

!#######################################################
! Extend matrix for the SRA treatment
! according to formula 4.107 of Bauer, PhD thesis
!
!V = ( 1/2M = 1/2M_0 l(l+1)/r**2 + V_LL;     0     )
!    (             0                       2M-2M_0 )
!#######################################################

!Write(*,*)'NSRA = ', nsra

if (nsra==2) then

  if ( .not. config_testflag('sph') .or. nsra==5 ) then
    call vllmatsra(Vpotll,cellnew%rmeshnew,eryd,lmaxatom,0,'Ref=0')
  else
    call vllmatsra(Vpotll,cellnew%rmeshnew,eryd,lmaxatom,0,'Ref=Vsph')
  end if

  if (kspinorbit==1) then    ! do the same with the potential matrix
                             ! used for the left solution
    if ( .not. config_testflag('sph') .or. nsra==5 ) then
      call vllmatsra(Vpotll2,cellnew%rmeshnew,eryd,lmaxatom,0,'Ref=0')
    else
      call vllmatsra(Vpotll2,cellnew%rmeshnew,eryd,lmaxatom,0,'Ref=Vsph')
    end if

  end if

end if


if ( config_testflag('vlldebug') ) then
  do ir=1,cellnew%nrmaxnew
    write(1235,'(50000E)') Vpotll(:,:,ir)
    if (kspinorbit==1) write(11235,'(50000E)') Vpotll2(:,:,ir)
  end do
end if

! might be deleted in the future
if ( config_testflag('kappamutest')) then
  call VLL_TRANSFORM(Vpotll,lmaxatom)
end if

if ( config_testflag('vlldebug') ) then
  do ir=1,cellnew%nrmaxnew
    write(12351,'(50000E)') Vpotll(:,:,ir)
    if (kspinorbit==1) write(112351,'(50000E)') Vpotll2(:,:,ir)
  end do
end if


!#######################################################
! calculate the source terms in the Lippmann-Schwinger equation
!#######################################################

!#######################################################
! these are in priciple spherical hankel and bessel functions
! which are extended with derivates of j and h's for a SR calculation
! for details, check chapter 4 of Bauer, PhD thesis
!#######################################################

! calculate the Bessel and Hankel functions
call rllsllsourceterms( nsra,wavefunction%nvec,eryd,cellnew%rmeshnew,cellnew%nrmaxnew,lmaxatom,lmsize,use_fullgmat,jlk_index,hlk,jlk,hlk2,jlk2,GMATPREFACTOR)

! might be deleted in the future
if ( config_testflag('kappamutest')) then
  call JLK_TRANSFORM(jlk,lmaxatom,jlk_index)
  call JLK_TRANSFORM(hlk,lmaxatom,jlk_index)
  call JLK_TRANSFORM(jlk2,lmaxatom,jlk_index)
  call JLK_TRANSFORM(hlk2,lmaxatom,jlk_index)
        DO LM1=1,2*LMSIZE
          jlk_index(LM1)=LM1
        END DO
end if

if ( config_testflag('writesourceterms')) then
  do lm1=1,ubound(jlk,1)
    write(1661,'(50000E)') jlk(lm1,:)
    write(1662,'(50000E)') hlk(lm1,:)
    write(1663,'(50000E)') jlk2(lm1,:)
    write(1664,'(50000E)') hlk2(lm1,:)
  end do
end if

if ( config_testflag('conjgtest')) then
  call conjugate2(jlk )
  call conjugate2(jlk2)
  call conjugate2(hlk )
  call conjugate2(hlk2)
  GMATPREFACTOR=conjg(GMATPREFACTOR)
  call conjugate3(VPOTLL)
end if

!#######################################################
! if the option 'sph' is set then the wave functions of the 
! spherical part of the potential is used as a reference
! function. The solutions of a sperical potential are calculated by using
! bessel and hanel functions and are stored in the same array:
!#######################################################

if ( config_testflag('sph') .and. nsra/=5 ) then
      call calcsph(nsra,cellnew,zatom,use_fullgmat,nspin,ispin,lmaxatom,eryd, &
                   jlk_index,hlk,jlk,hlk2,jlk2,GMATPREFACTOR,gauntcoeff(lmaxatom) ,tmatsph ,idotime)
end if

if ( config_testflag('writesourceterms')) then
  do lm1=1,ubound(jlk,1)
    write(1671,'(50000E)') jlk(lm1,:)
    write(1672,'(50000E)') hlk(lm1,:)
    write(1673,'(50000E)') jlk2(lm1,:)
    write(1674,'(50000E)') hlk2(lm1,:)
  end do
end if


!#######################################################
! The Lippmann-Schwinger equation for the full-potential
! is solved using using the Chebyshev integration method
! see Chapter 5 of Bauer, PhD thesis
!#######################################################

if (.not. allocated(wavefunction%SLL)) then
  allocate (wavefunction%SLL(lmsize2,lmsize,cellnew%nrmaxnew,1),&
            wavefunction%RLL(lmsize2,lmsize,cellnew%nrmaxnew,1))
end if

wavefunction%rll=(0.0D0,0.0D0)
wavefunction%sll=(0.0D0,0.0D0)

! might be deleted in the future
if ( config_testflag('sw') ) then
  write(*,*) 'switch source terms'
  call switch_jlk(jlk)
  call switch_jlk(jlk2)
  call switch_jlk(hlk)
  call switch_jlk(hlk2)
  call switch_vll(VPOTLL)
end if

call timing_start('---rll call---')

if (nsra==4) then
  hlk = hlk  / sqrt( (1.0D0,0.0D0)+(eryd)/cvlight**2)
  jlk = jlk  / sqrt( (1.0D0,0.0D0)+(eryd)/cvlight**2)
  hlk2= hlk2 / sqrt( (1.0D0,0.0D0)+(eryd)/cvlight**2)
  jlk2= jlk2 / sqrt( (1.0D0,0.0D0)+(eryd)/cvlight**2)
end if


!#######################################################! 
! calculate the right-hand side solution of the single-site wave functions
!#######################################################! 
call RLLSLL(cellnew%rpan_intervall,cellnew%rmeshnew,VPOTLL,&
            wavefunction%RLL(:,:,:,1),wavefunction%SLL(:,:,:,1),tmat%tmat, &
            cellnew%ncheb,cellnew%npan_tot,lmsize,lmsize2,cellnew%nrmaxnew,wavefunction%nvec,jlk_index,hlk,jlk,hlk2,jlk2,GMATPREFACTOR,'1','1','0',idotime)

if (nsra==2) then ! for nummerical reasons a factor of cvlight has been added to the equations
                  ! which needs to be removed now
  wavefunction%RLL(lmsize+1:,:,:,1)=wavefunction%RLL(lmsize+1:,:,:,1)/cvlight
  wavefunction%SLL(lmsize+1:,:,:,1)=wavefunction%SLL(lmsize+1:,:,:,1)/cvlight
end if

!#######################################################! 
! Full-relativistic: Transformation from kappa-mu to Lms basis
!#######################################################! 
if (nsra==3) then
! Pascal
  write(*,*) "Transforming wavefunctions to real spherical harmonics basis."
  do ir=1,cellnew%nrmaxnew
      call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,wavefunction%RLL(1:lmsize,1:lmsize,ir,1),'REL>RLM')
      call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,wavefunction%RLL(lmsize+1:2*lmsize,1:lmsize,ir,1),'REL>RLM')
      call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,wavefunction%SLL(1:lmsize,1:lmsize,ir,1),'REL>RLM')
      call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,wavefunction%SLL(lmsize+1:2*lmsize,1:lmsize,ir,1),'REL>RLM')
  end do
  call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,tmat%tmat,'REL>RLM')
  write(*,*) "done."
end if

if ( config_testflag('conjgtest')) then
  call  conjugate4(wavefunction%RLL)
  call  conjugate4(wavefunction%SLL)
end if

call timing_stop('---rll call---')

if ( config_testflag('kappamutest')) then
  call RLL_TRANSFORM(wavefunction%RLL(:,:,:,1),lmaxatom,'REL>RLM')
  call RLL_TRANSFORM(wavefunction%SLL(:,:,:,1),lmaxatom,'REL>RLM')
end if


!#######################################################
! In case the option 'sph' is set. The output t-matrix
! just contains the non-sph part of the t-matrix. Thus,
! the sperical needs to be added
!#######################################################
if ( config_testflag('sph') .or. nsra==5 ) then
  do lm1=1,lmsize
    tmat%tmat(lm1,lm1)=tmat%tmat(lm1,lm1)+tmatsph(jlk_index(lm1))
  end do
end if

!#######################################################
! If spin-orbit coupling is used the left solution of the
! Hamiltonian is non-trivial and needs to be calculated explicitly
!#######################################################
if (kspinorbit==1) then

  if (.not. allocated(wavefunction%SLLleft)) then
    allocate (wavefunction%SLLleft(lmsize2,lmsize,cellnew%nrmaxnew,1),&
              wavefunction%RLLleft(lmsize2,lmsize,cellnew%nrmaxnew,1))
  end if

  wavefunction%SLLleft=(0.0D0,0.0D0)
  wavefunction%RLLleft=(0.0D0,0.0D0)

  deallocate(jlk_index,hlk,jlk,hlk2,jlk2)
  call rllsllsourceterms(nsra,wavefunction%nvec,(eryd),cellnew%rmeshnew,cellnew%nrmaxnew,lmaxatom,lmsize,use_fullgmat,jlk_index,hlk,jlk,hlk2,jlk2,GMATPREFACTOR)
  
  if ( config_testflag('sph') .and. nsra/=5 ) then
      call calcsph(nsra,cellnew,zatom,use_fullgmat,nspin,ispin,lmaxatom,(eryd), &
                   jlk_index,hlk2,jlk2,hlk,jlk,GMATPREFACTOR,gauntcoeff(lmaxatom) ,tmatsph ,idotime)
  end if

  call RLLSLL(cellnew%rpan_intervall,cellnew%rmeshnew,VPOTLL2,&
              wavefunction%RLLleft(:,:,:,1),wavefunction%SLLleft(:,:,:,1),tmattemp, &
!                                                                ------------>    watch out here changed the order for left and right solution <-----------
              cellnew%ncheb,cellnew%npan_tot,lmsize,lmsize2,cellnew%nrmaxnew,wavefunction%nvec,jlk_index,hlk2,jlk2,hlk,jlk,GMATPREFACTOR,'1','1','0',idotime)
  if (nsra==2) then
    wavefunction%RLLleft(lmsize+1:,:,:,1)=wavefunction%RLLleft(lmsize+1:,:,:,1)/(cvlight)
    wavefunction%SLLleft(lmsize+1:,:,:,1)=wavefunction%SLLleft(lmsize+1:,:,:,1)/(cvlight)
  end if

  if ( config_testflag('tmatdebug') ) then
    do lm1=1,lmsize
      do lm2=1,lmsize
        write(4100,'(50000E)') wavefunction%rllleft(lm2,lm1,:,1)
        write(4101,'(50000E)') wavefunction%sllleft(lm2,lm1,:,1)
      end do
    end do
  end if
  
  if (wavefunction%nvec==2) then
    if ( config_testflag('tmatdebug') ) then
      do lm1=1,lmsize
        do lm2=lmsize+1,2*lmsize
          write(4110,'(50000E)') wavefunction%rllleft(lm2,lm1,:,1)
          write(4111,'(50000E)') wavefunction%sllleft(lm2,lm1,:,1)
        end do
      end do
    end if
  end if

end if !(kspinorbit==1)

if ( config_testflag('tmatdebug') ) then
  do lm1=1,lmsize
    do lm2=1,lmsize
      write(4000,'(50000E)') wavefunction%rll(lm2,lm1,:,1)
      write(4001,'(50000E)') wavefunction%sll(lm2,lm1,:,1)
    end do
  end do
end if

if (wavefunction%nvec==2) then
  if ( config_testflag('tmatdebug') ) then
    do lm1=1,lmsize
      do lm2=lmsize+1,2*lmsize
        write(4010,'(50000E)') wavefunction%rll(lm2,lm1,:,1)
        write(4011,'(50000E)') wavefunction%sll(lm2,lm1,:,1)
      end do
    end do
  end if
end if

! kkrsusc
if ( config_runflag('kkrsusc') ) then

  ! Passing t-matrix to KKRSUSC
   inpsusc%tmatll(:,:,iatom) = tmat%tmat                                                                                   ! susc |||| added by Juba (2016)

! ! NOSOC (use_spinorbit = 0) or SOC (use_spinorbit = 1) with the new solver           
! if (cellnew%use_spinorbit==0) then
!   ! kkrsusc saving the wave functions for the projections                                                                 ! susc |||| added by Juba (2015) 
!   ! Allocate the regular and irregular solutions defined on the old grid                                                  ! susc |||| added by Juba (2015)
!   if (.not. allocated(wavefunction%SLL_new_grid)) then                                                                    ! susc |||| added by Juba (2015)
!     allocate (wavefunction%SLL_new_grid(0:2*lmaxatom+1,cellnew%nrmaxnew),&                                                ! susc |||| added by Juba (2015)
!             wavefunction%RLL_new_grid(0:2*lmaxatom+1,cellnew%nrmaxnew))                                                   ! susc |||| added by Juba (2015)
!   endif

!   ! Pointers from Full potential to ASA L --> (l,m)
!   if (.not.allocated(i2lm)) then
!     allocate (i2lm(2,lmmax))                       
!     allocate (work1_in(1:cellnew%nrmaxnew,0:lmaxatom),work2_in(1:cellnew%nrmaxnew,0:lmaxatom)&
!              ,work1_out(1:cell%nrmax,0:lmaxatom),work2_out(1:cell%nrmax,0:lmaxatom) )                                        ! susc |||| added by Juba (2015)
!   endif

!   ! Define the pointer
!   i = 0
!   do ilm = 0, lmaxatom
!     do ilm2 = -ilm, ilm
!       i = i + 1
!       i2lm(:,i) = (/ilm2,ilm/)
!     end do
!   end do

!   ! Save the diagonal parts in from the wave functions      
!   do isn = 1, nspin                                                             
!     ! Initialize Arrays for cheb2oldgride
!     work1_out = 0.d0
!     work2_out = 0.d0

!     do ilm = lmmax*(isn-1) + 1, lmmax*isn
!       wavefunction%RLL_new_grid(i2lm(2,ilm-lmmax*(isn-1)),:) = wavefunction%RLL(ilm,ilm,:,1)
!       wavefunction%SLL_new_grid(i2lm(2,ilm-lmmax*(isn-1)),:) = wavefunction%SLL(ilm,ilm,:,1)
!       write(99999999,*)'Pointers again ------>'
!       write(99999999,*) 'l component', i2lm(2,ilm-lmmax*(isn-1)), 'spin', isn
!     enddo

!     ! Transposing the wave functions for cheb2oldgrid
!     do ilm = 0, lmaxatom
!       write(99999999,'(1I12)') ilm
!       do i = 1, cellnew%nrmaxnew
!         work1_in(i,ilm) = wavefunction%RLL_new_grid(ilm,i) 
!         work2_in(i,ilm) = wavefunction%SLL_new_grid(ilm,i)
!         write(99999999,'(3F12.6)') cellnew%rmeshnew(i), work1_in(i,ilm)
!       enddo
!     enddo

!     ! Convert regular and irregular to the old grid                                                                ! susc |||| added by Juba (2015) 
!     call cheb2oldgridc(cell,cellnew,cellnew%ncheb,lmaxatom+1,work1_in,work1_out)                                   ! susc |||| added by Juba (2015)
!     call cheb2oldgridc(cell,cellnew,cellnew%ncheb,lmaxatom+1,work2_in,work2_out)                                   ! susc |||| added by Juba (2015) 
!    
!     ! Fill pzsusc array, etc. with radial wave functions for every energy:                                         ! susc |||| added by Juba (2015)
!     ! Only the big component is stored 
!     if ( config_runflag('kkrsusc') ) then                                                                          ! susc |||| added by Juba (2015)
!       inpsusc%pzsusc(1:cell%nrmax,0:lmaxatom,iatom,isn) = work1_out                                                ! susc |||| added by Juba (2015)
!       inpsusc%qzsusc(1:cell%nrmax,0:lmaxatom,iatom,isn) = work2_out                                                ! susc |||| added by Juba (2015)
!       inpsusc%szsusc(1:cell%nrmax,0:lmaxatom,iatom,isn) = 0.d0                                                     ! susc |||| added by Juba (2015)
!       inpsusc%fzsusc(1:cell%nrmax,0:lmaxatom,iatom,isn) = 0.d0                                                     ! susc |||| added by Juba (2015)
!     end if                                                                                                         ! susc |||| added by Juba (2015)

!   enddo ! spin

! else ! use_spinorbit 


  ! Allocate arrays for interpolation non-spherical solutions 
  if (.not.allocated(work1r_in)) then
  ! Big component
    allocate  (work1r_in(cellnew%nrmaxnew,1:lmsize),work2r_in(cellnew%nrmaxnew,1:lmsize)&
              ,work1l_in(cellnew%nrmaxnew,1:lmsize),work2l_in(cellnew%nrmaxnew,1:lmsize)&              
              ,work1r_out(cell%nrmax,1:lmsize),work2r_out(cell%nrmax,1:lmsize)&
              ,work1l_out(cell%nrmax,1:lmsize),work2l_out(cell%nrmax,1:lmsize))
    ! Small component
    if (nsra == 2) then 
      allocate (work3r_in(cellnew%nrmaxnew,1:lmsize),work4r_in(cellnew%nrmaxnew,1:lmsize)&
               ,work3l_in(cellnew%nrmaxnew,1:lmsize),work4l_in(cellnew%nrmaxnew,1:lmsize)&
               ,work3r_out(cell%nrmax,1:lmsize),work4r_out(cell%nrmax,1:lmsize)&
               ,work3l_out(cell%nrmax,1:lmsize),work4l_out(cell%nrmax,1:lmsize)) 
    endif ! SRA

  endif

  do ilm = 1, lmsize 
    ! Initialize Arrays for cheb2oldgrid 
    ! Big component
    work1r_out = 0.d0
    work1l_out = 0.d0
    work2r_out = 0.d0
    work2l_out = 0.d0
    ! Small component
    if (nsra == 2) then
      work3r_out = 0.d0
      work3l_out = 0.d0
      work4r_out = 0.d0
      work4l_out = 0.d0
    endif ! SRA
     
    ! Transposing the wave functions for cheb2oldgrid
    do ilm2 = 1, lmsize
      do i = 1, cellnew%nrmaxnew      
        ! Big components
        work1r_in(i,ilm2) = wavefunction%RLL(ilm2,ilm,i,1)
        work1l_in(i,ilm2) = wavefunction%RLLleft(ilm2,ilm,i,1)
        work2r_in(i,ilm2) = wavefunction%SLL(ilm2,ilm,i,1)
        work2l_in(i,ilm2) = wavefunction%SLLleft(ilm2,ilm,i,1)
        ! Small components
        if (nsra == 2) then        
          work3r_in(i,ilm2) = wavefunction%RLL(ilm2+lmsize,ilm,i,1)
          work3l_in(i,ilm2) = wavefunction%RLLleft(ilm2+lmsize,ilm,i,1)
          work4r_in(i,ilm2) = wavefunction%SLL(ilm2+lmsize,ilm,i,1)
          work4l_in(i,ilm2) = wavefunction%SLLleft(ilm2+lmsize,ilm,i,1)
        endif ! SRA
      enddo
    enddo

    ! Convert regular and irregular solutions to the old grid
    ! Single column array --> See cheb2olgrid input                                                                                                               
    ! Big components
    call cheb2oldgridc(cell,cellnew,cellnew%ncheb,lmsize,work1r_in,work1r_out)                 ! susc |||| added by Juba (2015)
    call cheb2oldgridc(cell,cellnew,cellnew%ncheb,lmsize,work1l_in,work1l_out)                 ! susc |||| added by Juba (2015) 
    call cheb2oldgridc(cell,cellnew,cellnew%ncheb,lmsize,work2r_in,work2r_out)                 ! susc |||| added by Juba (2015)
    call cheb2oldgridc(cell,cellnew,cellnew%ncheb,lmsize,work2l_in,work2l_out)                 ! susc |||| added by Juba (2015) 
    ! Small components 
    if (nsra == 2) then
      call cheb2oldgridc(cell,cellnew,cellnew%ncheb,lmsize,work3r_in,work3r_out)               ! susc |||| added by Juba (2015)
      call cheb2oldgridc(cell,cellnew,cellnew%ncheb,lmsize,work3l_in,work3l_out)               ! susc |||| added by Juba (2015) 
      call cheb2oldgridc(cell,cellnew,cellnew%ncheb,lmsize,work4r_in,work4r_out)               ! susc |||| added by Juba (2015)
      call cheb2oldgridc(cell,cellnew,cellnew%ncheb,lmsize,work4l_in,work4l_out)               ! susc |||| added by Juba (2015) 
    endif ! SRA

    ! Fill pzsusc array, etc. with radial wave functions for every energy:                                         ! susc |||| added by Juba (2015)
    ! Only the big component is stored for now 
    inpsusc%pzsusc_soc(1:cell%nrmax,1:lmsize,ilm,iatom)     = work1r_out   
    inpsusc%pzsusc_socleft(1:cell%nrmax,1:lmsize,ilm,iatom) = work1l_out
    inpsusc%qzsusc_soc(1:cell%nrmax,1:lmsize,ilm,iatom)     = work2r_out
    inpsusc%qzsusc_socleft(1:cell%nrmax,1:lmsize,ilm,iatom) = work2l_out
    ! Small components 
    if (nsra == 2) then
      inpsusc%fzsusc_soc(1:cell%nrmax,1:lmsize,ilm,iatom)      = work3r_out                                  ! susc |||| added by Juba (2015)
      inpsusc%fzsusc_socleft(1:cell%nrmax,1:lmsize,ilm,iatom)  = work3l_out                                  ! susc |||| added by Juba (2015)
      inpsusc%szsusc_soc(1:cell%nrmax,1:lmsize,ilm,iatom)      = work4r_out
      inpsusc%szsusc_socleft(1:cell%nrmax,1:lmsize,ilm,iatom)  = work4l_out
    else 
      inpsusc%fzsusc_soc(1:cell%nrmax,1:lmsize,ilm,iatom)      = 0.d0                                        ! susc |||| added by Juba (2015)
      inpsusc%fzsusc_socleft(1:cell%nrmax,1:lmsize,ilm,iatom)  = 0.d0                                        ! susc |||| added by Juba (2015)
      inpsusc%szsusc_soc(1:cell%nrmax,1:lmsize,ilm,iatom)      = 0.d0
      inpsusc%szsusc_socleft(1:cell%nrmax,1:lmsize,ilm,iatom)  = 0.d0
    endif ! SRA
  
  enddo ! ilm

!  endif ! use_spinorbit

endif ! kkrsusc


!#######################################################
! calculation of Jij's by a Lichtenstein-like approach
! check section 6.3.3 Bauer, PhD
!#######################################################

if (.not. allocated (tmat%deltaT_Jij)) then
  allocate(tmat%deltaT_Jij(lmsize,lmsize,3) )
end if

if (config%calcJijmat==1) then
call  calccouplingdeltat(wavefunction,tmat%deltaT_Jij,cellnew,gauntcoeff(lmaxatom),theta,phi,lmmax,lmsize,lmaxatom,lmpot,cellnew%nrmaxnew)
end if

!#######################################################
! calculation of the Wronskian. Just for nummerical checks
!#######################################################

if ( config_testflag('wronskian') ) then
  if (kspinorbit==1) then
    call calcwronskian(wavefunction%rll(:,:,:,1),wavefunction%sll(:,:,:,1), &
                       wavefunction%rllleft(:,:,:,1),wavefunction%sllleft(:,:,:,1), &
                       cellnew)
  else
    call calcwronskian(wavefunction%rll(:,:,:,1),wavefunction%sll(:,:,:,1), &
                       wavefunction%rll(:,:,:,1),wavefunction%sll(:,:,:,1), &
                       cellnew)
  end if
end if

if ( config_testflag('checknan') ) then

      call checknan( wavefunction%rll ,ierror)
      if (ierror==1) stop '[calctmat_bauernew] wavefunction%rll nan'
      call checknan( wavefunction%sll ,ierror)
      if (ierror==1) stop '[calctmat_bauernew] wavefunction%sll nan'

end if


end subroutine calctmat_bauernew

end module mod_calctmat_bauernew
