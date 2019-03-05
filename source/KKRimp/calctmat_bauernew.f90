module mod_calctmat_bauernew
!-------------------------------------------------------------------------------
!> Summary: Calculate t-matrices including non-collinear magnetism and
!> spin-orbit coupling
!> Author: David Bauer
!> Category: KKRimp, single-site, spin-orbit-coupling, dirac, potential   
!>           
!-------------------------------------------------------------------------------

contains
!-------------------------------------------------------------------------------
!> Summary: Calculate t-matrices including non-collinear magnetism and
!> spin-orbit coupling
!> Author: David Bauer
!> Category: KKRimp, single-site, spin-orbit-coupling, dirac, potential   
!>           
!-------------------------------------------------------------------------------
subroutine calctmat_bauernew(cell,tmat,lmaxatom,eryd_in,ZATOM,cellnew,wavefunction, &
  ispin,nspin,kspinorbit,use_fullgmat,theta,phi,ncoll,nsra,config,idotime, &
  ie,ldau,iatom,cellorbit,calcleft)        ! lda+u

use mod_datatypes, only: dp
use mod_constants, only: czero
use global_variables, only: korbit
use mod_gauntharmonics, only: gauntcoeff
use mod_timing, only: timing_start, timing_stop
use type_tmat, only: tmat_type
use type_cell, only: cell_type
use type_cellnew, only: cell_typenew
use type_cellorbit, only: cell_typeorbit
use type_wavefunction, only: wavefunction_type
use type_config, only: config_type
use type_ldau, only: ldau_type                ! lda+u
use mod_interpolpot, only: interpolpot
use mod_vllmat, only: vllmat
use mod_vllmatsra, only: vllmatsra
use mod_rllsll, only: rllsll
use mod_calccouplingconstants, only: calccouplingdeltat
use mod_config, only: config_testflag
use mod_spinorbit_ham, only: spinorbit_ham
use mod_rllsllsourceterms, only: rllsllsourceterms
use mod_physic_params, only: cvlight
use Potential, only:  PotentialMatrixArray
use mod_calcsph, only: calcsph
use mod_checknan, only: checknan
use mod_basistransform, only: rll_transform, vll_transform, jlk_transform, single_transform
use mod_mathtools, only: conjugate2, conjugate3, conjugate4
use mod_calctmat_bauernew_testtools, only: switch_jlk, switch_vll
use mod_wronskian, only: calcwronskian
implicit none
!interface
type(cell_type),intent(in)                :: cell
type(tmat_type)                           :: tmat
integer                                   :: lmaxatom
double complex                            :: eryd_in
double precision                          :: zatom
type(cell_typenew)                        :: cellnew
type(cell_typeorbit)                      :: cellorbit
integer                                   :: iatom
type(wavefunction_type)                   :: wavefunction
integer                                   :: ispin,nspin
integer                                   :: kspinorbit,use_fullgmat
double precision                          :: theta,phi
integer                                   :: ncoll
integer                                   :: nsra
type(config_type)                         :: config
integer                                   :: ie                         ! lda+u
type(ldau_type)                           :: ldau                       ! lda+u variables
logical                                   :: calcleft
!local
real (kind=dp),allocatable                :: vins(:,:,:)
double complex,allocatable                :: vpotll(:,:,:)
double complex,allocatable                :: vpotll2(:,:,:)
double complex,allocatable                :: vnspll0(:,:,:)
double complex,allocatable                :: vnspll1(:,:,:)
double complex,allocatable                :: vnspll2(:,:,:)
integer                                   :: ir,lm1,lm2,lmpot,lmmax,lmsize,lmsize2, use_sratrick, istat
double complex                            :: gmatprefactor
double complex,allocatable                :: hlk(:,:),jlk(:,:),hlk2(:,:),jlk2(:,:)
integer,allocatable                       :: jlk_index(:)

integer                                   :: ierror,idotime

complex(kind=dp),allocatable              :: tmattemp(:,:)
complex(kind=dp),allocatable              :: tmatsph(:)
double complex                            :: eryd
integer                                   :: lmlo,lmhi,mmax,imt1  ! lda+u 


    character (len=100) :: filename

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
tmattemp = czero

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
IF (NSRA<=2) then ! non/scalar-relativistic case (+SOC)

  allocate(Vpotll(2*lmsize,2*lmsize,cellnew%nrmaxnew), stat=istat)
  if (istat/=0) stop 'Error allocating Vpotll in calctmat_bauernew'
  Vpotll = czero

  if ( config_testflag('nosph') ) then
    use_sratrick = 0
  else
    use_sratrick = 1
  endif

  allocate(vins(cellnew%nrmaxnew, lmpot, nspin), vnspll0(lmsize, lmsize, cellnew%nrmaxnew), stat=istat)
  if (istat/=0) stop 'Error allocating vins in calctmat_bauernew'
  vins = czero
  vnspll0 = czero
  vins(1:cellnew%nrmaxnew, 1:lmpot, 1) = cellnew%Vpotnew(:,:,1)
  vins(1:cellnew%nrmaxnew, 1:lmpot, nspin) = cellnew%Vpotnew(:,:,nspin)
  call vllmat(1, cellnew%nrmaxnew, cellnew%nrmaxnew, lmmax, lmsize, vnspll0, vins, lmpot, gauntcoeff(lmaxatom)%cleb, gauntcoeff(lmaxatom)%icleb, &
    gauntcoeff(lmaxatom)%iend, nspin, zatom, cellnew%rmeshnew, use_sratrick, gauntcoeff(lmaxatom)%ncleb)

else if (nsra==3) then ! use the full Dirac solver by Pascal Kordt, PhD thesis

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
         vnspll0(lmlo:lmhi,lmlo:lmhi,ir) = vnspll0(lmlo:lmhi,lmlo:lmhi,ir) + ldau%wldau(1:mmax,1:mmax,1)                ! lda+u
      enddo                                                                          ! lda+u
      lmlo = lmlo + lmmax                                                            ! lda+u
      lmhi = lmhi + lmmax                                                            ! lda+u
      do ir = 1,imt1                                                                 ! lda+u
         vnspll0(lmlo:lmhi,lmlo:lmhi,ir) = vnspll0(lmlo:lmhi,lmlo:lmhi,ir) + ldau%wldau(1:mmax,1:mmax,2)        ! lda+u
      enddo                                                                          ! lda+u

   else       ! 1x1 in spin space                                                    ! lda+u 

      do ir = 1,imt1
         vnspll0(lmlo:lmhi,lmlo:lmhi,ir) = vnspll0(lmlo:lmhi,lmlo:lmhi,ir) + ldau%wldau(1:mmax,1:mmax,ispin)            ! lda+u
      enddo                                                                          ! lda+u
   endif  

endif                                                                                ! lda+u
!----------------------------------------------------------------------------------- ! lda+u


if ( config_testflag('vlldebug') ) then
  do ir=1,cellnew%nrmaxnew
      write(1233,'(50000E)') Vnspll0(:,:,ir)
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

  allocate(Vpotll2(2*lmsize,2*lmsize,cellnew%nrmaxnew))
  allocate(vnspll1(lmsize,lmsize,cellnew%nrmaxnew))
  allocate(vnspll2(lmsize,lmsize,cellnew%nrmaxnew))
  vpotll2 = czero
  vnspll1 = czero
  vnspll2 = czero

  write(1337,*) 'spinorbit index','','atom',iatom,cellorbit%use_spinorbit(iatom)
  if (cellorbit%use_spinorbit(iatom)==1) then

    ! for right solution
    call spinorbit_ham(lmaxatom, lmmax, vins, cellnew%rmeshnew, eryd, zatom, cvlight, 1.0_dp, nspin, &
      lmpot, theta, phi, cellnew%ipan_intervall, cellnew%rpan_intervall, cellnew%npan_tot, &
      cellnew%ncheb, cellnew%nrmaxnew, cellnew%nrmaxnew, vnspll0, vnspll1, '1')

    ! for left solution
    call spinorbit_ham(lmaxatom, lmmax, vins, cellnew%rmeshnew, eryd, zatom, cvlight, 1.0_dp, nspin, &
      lmpot, theta, phi, cellnew%ipan_intervall, cellnew%rpan_intervall, cellnew%npan_tot, &
      cellnew%ncheb, cellnew%nrmaxnew, cellnew%nrmaxnew, vnspll0, vnspll2, 'transpose')

  else

    vnspll1 = vnspll0
    vnspll2 = vnspll0

  end if

end if

if ( config_testflag('vlldebug') ) then
  do ir=1,cellnew%nrmaxnew
      write(1234,'(50000E)') vnspll1(:,:,ir)
      if (kspinorbit==1) write(11234,'(50000E)') vnspll2(:,:,ir)
  end do
  open (7352834, file='vnspll_SOC.txt', form='formatted')
  write (7352834, '(A,3I9)') '# LMMAXSO,LMMAXSO,IRMDNEW=', lmsize, lmsize, cellnew%nrmaxnew
  write (7352834, '(2F25.14)') vnspll1(:, :, :)
  close (7352834)
end if


!#######################################################
! Extend matrix for the SRA treatment
! according to formula 4.107 of Bauer, PhD thesis
!
!V = ( 1/2M = 1/2M_0 l(l+1)/r**2 + V_LL;     0     )
!    (             0                       2M-2M_0 )
!#######################################################

if (nsra==2) then

  if ( config_testflag('nosph') .or. nsra==5 ) then
    call vllmatsra(vnspll1, Vpotll, cellnew%rmeshnew, lmsize, cellnew%nrmaxnew, cellnew%nrmaxnew, &
      eryd, lmaxatom, 0, 'Ref=0')
  else
    call vllmatsra(vnspll1, Vpotll, cellnew%rmeshnew, lmsize, cellnew%nrmaxnew, cellnew%nrmaxnew, &
      eryd, lmaxatom, 0, 'Ref=Vsph')
  end if

  if (kspinorbit==1) then    ! do the same with the potential matrix
                             ! used for the left solution
    if ( config_testflag('nosph') .or. nsra==5 ) then
      call vllmatsra(vnspll2, Vpotll2, cellnew%rmeshnew, lmsize, cellnew%nrmaxnew, cellnew%nrmaxnew, &
        eryd, lmaxatom, 0, 'Ref=0')
    else
      call vllmatsra(vnspll2, Vpotll2, cellnew%rmeshnew, lmsize, cellnew%nrmaxnew, cellnew%nrmaxnew, &
        eryd, lmaxatom, 0, 'Ref=Vsph')
    end if

  end if

end if

if ( config_testflag('vlldebug') ) then
  do ir=1,cellnew%nrmaxnew
    write(1235,'(50000E)') Vpotll(:,:,ir)
    if (kspinorbit==1) write(11235,'(50000E)') Vpotll2(:,:,ir)
  end do
  open (7352834, file='vnspll_sra.txt', form='formatted')
  if (nsra==2) then
    write (7352834, '(A,3I9)') '# 2*LMMAXSO,2*LMMAXSO,IRMDNEW=', 2*lmsize, 2*lmsize, cellnew%nrmaxnew
  else
    write (7352834, '(A,3I9)') '# LMMAXSO,LMMAXSO,IRMDNEW=', lmsize, lmsize, cellnew%nrmaxnew
  end if
  write (7352834, '(2F25.14)') Vpotll(:, :, :)
  close (7352834)
end if

! might be deleted in the future
if ( config_testflag('kappamutest')) then
  call VLL_TRANSFORM(Vpotll,lmaxatom)
end if

if ( config_testflag('vlldebug') ) then
  do ir=1,cellnew%nrmaxnew
    write(12351,'(50000E)') Vpotll(1:lmsize,1:lmsize,ir)
    if (kspinorbit==1) write(112351,'(50000E)') Vpotll2(1:lmsize,1:lmsize,ir)
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
allocate(jlk_index(wavefunction%nvec*lmsize), &
  hlk(1:nsra*(1+kspinorbit)*(lmaxatom+1),cellnew%nrmaxnew), &
  jlk(1:nsra*(1+kspinorbit)*(lmaxatom+1),cellnew%nrmaxnew), &
  hlk2(1:nsra*(1+kspinorbit)*(lmaxatom+1),cellnew%nrmaxnew), &
  jlk2(1:nsra*(1+kspinorbit)*(lmaxatom+1),cellnew%nrmaxnew), &
  stat=istat)
if(istat/=0) stop 'Error allocating jlk_index etc. in calctmat_bauernew'
jlk = czero
hlk = czero
jlk2 = czero
hlk2 = czero
! this is needed since rllsllsourceterms uses this from global variables
korbit = kspinorbit
call rllsllsourceterms(nsra, wavefunction%nvec, eryd, cellnew%rmeshnew, cellnew%nrmaxnew, cellnew%nrmaxnew, &
  lmaxatom, lmsize, use_fullgmat, jlk_index, hlk, jlk, hlk2, jlk2, GMATPREFACTOR)

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
  write (filename, '(A,I0.3,A,I0.3,A)') 'rll_source_jlk_atom_', 1, '_energ_', 1, '.dat'
  open (888888, file=trim(filename), form='formatted')
  write (888888, '(2ES21.9)') jlk(:, :)
  close (888888)
  write (filename, '(A,I0.3,A,I0.3,A)') 'rll_source_hlk_atom_', 1, '_energ_', 1, '.dat'
  open (888888, file=trim(filename), form='formatted')
  write (888888, '(2ES21.9)') hlk(:, :)
  close (888888)
  write (filename, '(A,I0.3,A,I0.3,A)') 'rll_source_jlk2_atom_', 1, '_energ_', 1, '.dat'
  open (888888, file=trim(filename), form='formatted')
  write (888888, '(2ES21.9)') jlk2(:, :)
  close (888888)
  write (filename, '(A,I0.3,A,I0.3,A)') 'rll_source_hlk2_atom_', 1, '_energ_', 1, '.dat'
  open (888888, file=trim(filename), form='formatted')
  write (888888, '(2ES21.9)') hlk2(:, :)
  close (888888)
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
! if the option 'nosph' is not set then the wave functions of the 
! spherical part of the potential is used as a reference
! function. The solutions of a sperical potential are calculated by using
! bessel and hanel functions and are stored in the same array:
!#######################################################

if ( .not. config_testflag('nosph') .and. nsra/=5 ) then
      allocate(tmatsph(nspin*(lmaxatom+1)), stat=istat)
      tmatsph = czero
      if(istat/=0) stop 'Error allocating tmatsph in calctmat_bauernew'
      call calcsph(nsra, cellnew%nrmaxnew, cellnew%nrmaxnew, lmaxatom, nspin/(2-kspinorbit), zatom, eryd, &
        lmpot, lmsize, cellnew%rmeshnew, vins, cellnew%ncheb, cellnew%npan_tot, cellnew%rpan_intervall, jlk_index, &
        hlk, jlk, hlk2, jlk2, gmatprefactor, tmatsph, tmattemp, use_sratrick)
end if

if ( config_testflag('writesourceterms')) then
  do lm1=1,ubound(jlk,1)
    write(1671,'(50000E)') jlk(lm1,:)
    write(1672,'(50000E)') hlk(lm1,:)
    write(1673,'(50000E)') jlk2(lm1,:)
    write(1674,'(50000E)') hlk2(lm1,:)
  end do
  write (filename, '(A,I0.3,A,I0.3,A)') 'tmatsph_atom_', 1, '_energ_', 1, '.dat'
  open (888888, file=trim(filename), form='formatted')
  write (888888, '(2ES21.9)') tmatsph(:)
  close (888888)
  write (filename, '(A,I0.3,A,I0.3,A)') 'rll_sph_jlk_atom_', 1, '_energ_', 1, '.dat'
  open (888888, file=trim(filename), form='formatted')
  write (888888, '(2ES21.9)') jlk(:, :)
  close (888888)
  write (filename, '(A,I0.3,A,I0.3,A)') 'rll_sph_hlk_atom_', 1, '_energ_', 1, '.dat'
  open (888888, file=trim(filename), form='formatted')
  write (888888, '(2ES21.9)') hlk(:, :)
  close (888888)
  write (filename, '(A,I0.3,A,I0.3,A)') 'rll_sph_jlk2_atom_', 1, '_energ_', 1, '.dat'
  open (888888, file=trim(filename), form='formatted')
  write (888888, '(2ES21.9)') jlk2(:, :)
  close (888888)
  write (filename, '(A,I0.3,A,I0.3,A)') 'rll_sph_hlk2_atom_', 1, '_energ_', 1, '.dat'
  open (888888, file=trim(filename), form='formatted')
  write (888888, '(2ES21.9)') hlk2(:, :)
  close (888888)
end if


!#######################################################
! The Lippmann-Schwinger equation for the full-potential
! is solved using using the Chebyshev integration method
! see Chapter 5 of Bauer, PhD thesis
!#######################################################

if (.not. allocated(wavefunction%SLL)) then
  allocate (wavefunction%SLL(lmsize2,lmsize,cellnew%nrmaxnew,1),&
            wavefunction%RLL(lmsize2,lmsize,cellnew%nrmaxnew,1))
  wavefunction%SLL = czero
  wavefunction%RLL = czero
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
tmat%tmat = czero
call rllsll(cellnew%rpan_intervall, cellnew%rmeshnew, Vpotll, wavefunction%RLL(:,:,:,1), wavefunction%SLL(:,:,:,1), &
  tmat%tmat, cellnew%ncheb, cellnew%npan_tot, lmsize, lmsize2, nsra*(1+kspinorbit)*(lmaxatom+1), cellnew%nrmaxnew, nsra, &
  jlk_index, hlk, jlk, hlk2, jlk2, GMATPREFACTOR, '1', '1', '0', use_sratrick, tmattemp)

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
! In case the option 'nosph' is not set. The output t-matrix
! just contains the non-sph part of the t-matrix. Thus,
! the sperical needs to be added
!#######################################################
if ( .not. config_testflag('nosph') .or. nsra==5 ) then
  do lm1=1,lmsize
    tmat%tmat(lm1,lm1)=tmat%tmat(lm1,lm1)+tmatsph(jlk_index(lm1))
  end do
end if

!#######################################################
! If spin-orbit coupling is used the left solution of the
! Hamiltonian is non-trivial and needs to be calculated explicitly
!#######################################################
if ((kspinorbit==1).and.calcleft) then

  if (.not. allocated(wavefunction%SLLleft)) then
    allocate (wavefunction%SLLleft(lmsize2,lmsize,cellnew%nrmaxnew,1),&
              wavefunction%RLLleft(lmsize2,lmsize,cellnew%nrmaxnew,1))
    wavefunction%SLLleft = czero
    wavefunction%RLLleft = czero
  end if

  wavefunction%SLLleft=(0.0D0,0.0D0)
  wavefunction%RLLleft=(0.0D0,0.0D0)

  jlk_index = czero
  hlk = czero
  jlk = czero
  jlk2 = czero
  hlk2 = czero
  call rllsllsourceterms(nsra, wavefunction%nvec, eryd, cellnew%rmeshnew, cellnew%nrmaxnew, cellnew%nrmaxnew, &
    lmaxatom, lmsize, use_fullgmat, jlk_index, hlk, jlk, hlk2, jlk2, GMATPREFACTOR)
   
  
  if ( .not. config_testflag('nosph') .and. nsra/=5 ) then
      call calcsph(nsra, cellnew%nrmaxnew, cellnew%nrmaxnew, lmaxatom, nspin/(2-kspinorbit), zatom, eryd, &
        lmpot, lmsize, cellnew%rmeshnew, vins, cellnew%ncheb, cellnew%npan_tot, cellnew%rpan_intervall, jlk_index, &
        hlk2, jlk2, hlk, jlk, gmatprefactor, tmatsph, tmattemp, use_sratrick)
  end if

  call rllsll(cellnew%rpan_intervall, cellnew%rmeshnew, Vpotll2, wavefunction%RLLleft(:,:,:,1), wavefunction%SLLleft(:,:,:,1), &
    tmattemp, cellnew%ncheb, cellnew%npan_tot, lmsize, lmsize2, nsra*(1+kspinorbit)*(lmaxatom+1), cellnew%nrmaxnew, nsra, &
    ! ------------>    watch out here changed the order for left and right solution <-----------
    jlk_index, hlk2, jlk2, hlk, jlk, GMATPREFACTOR, '1', '1', '0', use_sratrick, tmattemp)
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
      cellnew%ncheb, cellnew%npan_tot, cellnew%ipan_intervall, cellnew%rpan_intervall)
      !cellnew)
  else
    call calcwronskian(wavefunction%rll(:,:,:,1),wavefunction%sll(:,:,:,1), &
      wavefunction%rll(:,:,:,1),wavefunction%sll(:,:,:,1), &
      cellnew%ncheb, cellnew%npan_tot, cellnew%ipan_intervall, cellnew%rpan_intervall)
      !cellnew)
  end if
end if

if ( config_testflag('checknan') ) then

      call checknan( wavefunction%rll ,ierror)
      if (ierror==1) stop '[calctmat_bauernew] wavefunction%rll nan'
      call checknan( wavefunction%sll ,ierror)
      if (ierror==1) stop '[calctmat_bauernew] wavefunction%sll nan'

end if

! clean up allocations
deallocate(Vpotll, tmattemp, vins, jlk_index, hlk, jlk, hlk2, jlk2, stat=istat)
if (istat/=0) stop 'Error deallocating arrays in calctmat_bauernew'
if ( .not. config_testflag('nosph') .and. nsra/=5 ) then
  deallocate(tmatsph, stat=istat)
  if (istat/=0) stop 'Error deallocating arrays in calctmat_bauernew'
end if
  


end subroutine calctmat_bauernew

end module mod_calctmat_bauernew
