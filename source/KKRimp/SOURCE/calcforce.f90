module mod_calcforce
      DOUBLE PRECISION,allocatable :: FLM(:,:),FLMC(:,:)
contains

subroutine calcforce(mode,cmom,cmom_interst,lmaxatom,lmaxd,nspin,natom,density,VPOT_OUT, &
                     cell,ins, zatom,nrmaxd,alat)
use type_density
use type_cell
use mod_forceh
use mod_force
use mod_forcxc
implicit none
!interface
character(len=*) :: mode
double precision              :: cmom((2*lmaxd+1)**2,natom) !(lmpotd,ntotatom)
double precision          :: cmom_interst((2*lmaxd+1)**2,natom) !(lmpotd,ntotatom)
integer             :: lmaxatom(natom)
integer             :: lmaxd
integer             :: nspin
integer             :: natom
type(density_type)  :: density(natom)
real*8              :: vpot_out(nrmaxd,(2*lmaxd+1)**2,nspin,natom) !thetas(iri,nfund,*),
type(cell_type)     :: cell(natom)
integer             :: ins
double precision    :: zatom(natom)
integer             :: nrmaxd
real*8              :: alat
!local
double precision          :: cmom_mt((2*lmaxd+1)**2,natom) !(lmpotd,ntotatom)

! allocate arrays
if (.not. allocated(flm)) allocate (flm(-1:1,natom),flmc(-1:1,natom))

! by the kkrflex definition cmom is the total charge moment
! for the force calculation the mt-part is needed. that's
! why the interstitial part is substracted
 cmom_mt=cmom-cmom_interst

! the routines are called in two part. in the first part the madelung
! potential is needed. in the second part the full kohn-sham potential
! is needed
if (mode=='part1') then
  flm=0.0d0
  flmc=0.0d0
  call forceh(cmom_mt,flm,2*lmaxd,(2*lmaxd+1)**2,nspin,natom,cell,density, &
              vpot_out,cell(1)%nrmaxd,zatom,ins)

  call force(flm,flmc,2*lmaxd,nspin,natom,vpot_out, &
             density,cell,cell(1)%nrmaxd, (2*lmaxd+1)**2,ins)

else if (mode=='part2') then
  call forcxc(flm,flmc,2*lmaxd,nspin,natom,vpot_out,density,cell,alat,(2*lmaxd+1)**2, &
              cell(1)%nrmaxd,ins)
else 
  stop '[calcforce] mode not known'
end if

end subroutine

end module mod_calcforce