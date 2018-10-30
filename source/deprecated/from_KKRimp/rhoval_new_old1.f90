module mod_rhoval_new

contains

subroutine rhoval_new(eryd,ie,wez,cellnew,wavefunction,cell,gmatll,iatom,ispin,nspin,SHAPEFUN,GAUNTCOEFF, ZATOM, DENSITY, &
                         LMAXATOM,LMMAXATOM,config,lmaxd,energyparts,kspinorbit,use_fullgmat,nspinden)
! use mod_gauntharmonics, only: gauntcoeff
use type_cell
use type_cellnew
use type_shapefun
use type_gauntcoeff
use type_density
use type_config
use type_wavefunction
use type_energyparts
use mod_cheb2oldgridc
! use mod_interpolpot
use mod_vllmat
use mod_checknan
! use mod_gauntharmonics, only: gauntcoeff
use mod_rllsll
use mod_rhooutnew
use mod_config, only: config_testflag
use mod_physic_params, only: cvlight
use mod_spinorbit
implicit none
!interface
complex(kind=dpc)                         :: eryd
integer                                   :: ie
complex(kind=dpc)                         :: wez
type(cell_typenew)                        :: cellnew
type(wavefunction_type)                   :: wavefunction
type(cell_type),intent(in)                :: cell
complex(kind=dpc),intent(in)              :: gmatll(:,:)
integer                                   :: iatom
integer                                   :: ispin
integer                                   :: nspin
type(shapefun_type),intent(in)            :: shapefun
type(gauntcoeff_type),intent(in)          :: gauntcoeff
real(kind=dp),intent(in)                  :: zatom
type(density_type)                        :: density
integer                                   :: lmaxatom
integer,intent(in)                        :: lmmaxatom
type(config_type)                         :: config
integer                                   :: lmaxd
type(energyparts_type)                    :: energyparts
integer                                   :: kspinorbit
integer                                   :: use_fullgmat
integer                                   :: nspinden

!local
double complex,allocatable                :: rho2ns_complex(:,:,:)
!       real(kind=dp)                         ::   espv(0:lmaxatom+1,2)
integer                                   :: ierror
complex(kind=dpc)                         ::   df,ek
double complex                            :: rho2ns_integrated(4)
double complex,allocatable                ::  vpotll(:,:,:)!(1:cellnew%nrmax,(2*lmaxatom+1)**2)
integer                                   :: lval,lmsize
double complex                            :: tmat( (lmaxatom+1)**2, (lmaxatom+1)**2)

double precision,allocatable              :: c1(:,:)
double complex,allocatable                :: srllp(:,:), &
                                             ull(:,:,:,:)
double complex                            :: cden(nspinden,cellnew%nrmaxnew,0:lmaxd),cdenns(cellnew%nrmaxnew),efac(lmmaxatom) &
                                            ,cdenlm(nspinden,cellnew%nrmaxnew,lmmaxatom)  ! lm-dos

integer                                   :: lm1,lm2
integer                                   :: vlllmmax,jspin
integer                                   :: nspinstart, nspinstop

if (use_fullgmat==1) then
  lmsize=2*lmmaxatom
  nspinden=4
  nspinstart=1
  nspinstop=nspinden
else
  nspinden=nspin
  lmsize=lmmaxatom
  nspinstart=ispin
  nspinstop=ispin
end if

allocate(rho2ns_complex(cellnew%nrmaxnew,(2*lmaxatom+1)**2,nspinden))

if (ubound(gmatll,1)/=lmsize ) then
  stop 'error in rhovalnew'
end if

rho2ns_complex=(0.0d0,0.0d0)

df = wez/dble(nspin)
if (config%nsra.eq.1) ek = sqrt(eryd)
if (config%nsra.eq.2) ek = sqrt(eryd+eryd*eryd/ (cvlight*cvlight))


!#######################################################
! if the wavefunctions are stored in memory, then there
! is no need to recalculate them
!#######################################################

if ( .not. allocated(wavefunction%rll)) then

  allocate(vpotll(lmsize,lmsize,cellnew%nrmaxnew))

  !#######################################################
  ! calculate the VLL matrix using the expansion of the Potential
  ! in spherical harmonics
  !#######################################################
  call vllmat(vpotll,cellnew%vpotnew(:,:,:),lmaxatom,(lmaxatom+1)**2,(2*lmaxatom+1)**2,1, &
              cellnew%nrmaxnew,gauntcoeff,zatom,cellnew%rmeshnew,lmsize,use_fullgmat,nspin,ispin )

  !#######################################################
  ! add the spin-orbit Hamiltonian to the VLL matrix
  !#######################################################
  if (kspinorbit==1) then
    call spinorbit(lmaxatom,zatom,eryd,cellnew,cellnew%nrmaxnew,nspin,vpotll,density%theta,density%phi,config%ncoll)
  end if

  allocate( c1(0:cellnew%ncheb,0:cellnew%ncheb))

  allocate(srllp(lmmaxatom,lmmaxatom), &
           ull(lmmaxatom,lmmaxatom,cellnew%nrmaxnew,2),wavefunction%sll(lmmaxatom,lmmaxatom,cellnew%nrmaxnew,2),&
           wavefunction%rll(lmmaxatom,lmmaxatom,cellnew%nrmaxnew,2))!,hll(lmmaxatom,lmmaxatom,cellnew%nrmaxnew,2))

stop 'check lmsize'
! watch out the wave functions need to be multiplied by a factor of 1/sqrt(e) to
! be consistant with pns qns
  call rllsll(cellnew%rpan_intervall,cellnew%rmeshnew,vpotll,eryd,&
              c1,srllp,wavefunction%rll(:,:,:,1),ull(:,:,:,1),wavefunction%sll(:,:,:,1),tmat, &
              cellnew%ncheb,cellnew%npan_tot,lmaxatom,lmaxatom+1,lmsize,lmsize,(cellnew%ncheb+1)*cellnew%npan_tot,use_fullgmat,1,1)
  wavefunction%rll=wavefunction%rll/ek
  ull=ull/ek
  wavefunction%sll=wavefunction%sll/-ek
end if

if ( config_testflag('tmatnewdebug') ) then

write(8000,'(50000F)') cellnew%rmeshnew
do lm1=1,(use_fullgmat+1)*(lmaxatom+1)**2
  do lm2=1,(use_fullgmat+1)*(lmaxatom+1)**2
    write(8001,'(50000F)') wavefunction%rll(lm2,lm1,:,1)
  end do
end do

do lm1=1,(use_fullgmat+1)*(lmaxatom+1)**2
  do lm2=1,(use_fullgmat+1)*(lmaxatom+1)**2
    write(8002,'(50000E)') wavefunction%sll(lm2,lm1,:,1)
  end do
end do

end if

call rhooutnew(gauntcoeff,df,gmatll,ek,cellnew,wavefunction,rho2ns_complex(:,:,:), &
               config%nsra, &
               lmaxd,lmmaxatom,lmsize,(2*lmaxd+1)**2,cellnew%nrmaxnew,cellnew%npan_sph*cellnew%ncheb+1, &
               ispin,nspinden)


do jspin=nspinstart,nspinstop
    call intcheb_cell(rho2ns_complex(:,1,jspin),cellnew,jspin,rho2ns_integrated)
    density%rho2ns_integrated(jspin)=density%rho2ns_integrated(jspin)+rho2ns_integrated(jspin)
    call cheb2oldgridc(cell,cellnew,cellnew%ncheb,(2*lmaxatom+1)**2,rho2ns_complex(:,:,jspin),density%rho2ns_complex(:,:,jspin))
end do


! c
! c---> calculate complex density of states
! c
!       DO 30 L = 0,LMAXD
! c
! c---> call integration subroutine
! c
!         CALL CSIMPK(CDEN(1,L),DEN(L),IPAN,IRCUT,DRDI)
!    30 CONTINUE

!       DO 40 L = 1,LMMAXD  ! lm-dos
!         CALL CSIMPK(CDENLM(1,L),DENLM2(L),IPAN,IRCUT,DRDI)  ! lm-dos
!         DENLM(L)=DENLM2(L)
!    40 CONTINUE  ! lm-dos

! ! c Energy depends on EK and NSRA:                            ! lm-dos
! c     IF (NSRA.EQ.1) EK = SQRT(E)                           ! lm-dos
! c     IF (NSRA.EQ.2) EK = SQRT(E+E*E/ (CVLIGHT*CVLIGHT))    ! lm-dos
! c     CVLIGHT=274.0720442D0                                 ! lm-dos 
! c Therefore the following is a good approximation           ! lm-dos
! c for energies of a few Ryd:                                ! lm-dos
!       ENERG = EK**2                                         ! lm-dos

!       WRITE(30,9000) DREAL(ENERG),(-DIMAG(DENLM(LM))/PI,LM=1,LMMAXD)
!  9000 FORMAT(30E12.4)


!       IF (IPAN.GT.1) THEN
!         CALL CSIMPK(CDENNS,DENNS,IPAN,IRCUT,DRDI)
!         DEN(LMAXD+1) = DENNS
!       END IF




end subroutine rhoval_new

end module mod_rhoval_new
