module mod_vllmatsra2

contains

subroutine vllmatsra2(vll,vpot,cellnew,eryd,zatom)
!************************************************************************************
! The perubation matrix for the SRA-equations are set up
!************************************************************************************
use mod_physic_params, only: cvlight
use mod_mathtools, only: inverse
use type_cellnew
use mod_chebyshev, only : diff_cell_complex,diff_cell,diff2_cell
implicit none
!interface
  DOUBLE COMPLEX             :: VLL(:,:,:)
  DOUBLE COMPLEX,allocatable :: VLL2DDR(:,:)
  type(cell_typenew)          :: cellnew
!   double precision,allocatable :: rmesh(:)
  double complex              :: eryd
  double precision            :: zatom
  integer                      :: lmsize
!local

  double precision             :: vpot(:)
  double precision,allocatable :: dvpotdr(:),d2vpotdr2(:)
  double complex,allocatable  :: vpottemp(:),dvlldr(:,:)
  integer                     :: nrmax
  integer                     :: ival,ival2,ir
  double complex,parameter    :: cone=(1.0D0,0.0D0)
  double complex,parameter    :: czero=(0.0D0,0.0D0)
!  double precision,parameter :: cvlight = 274.0720442d0
  double complex               :: Mass,Mass0,dMassdr,d2Massdr2,vpotsra

lmsize=ubound(vll,1)
nrmax=ubound(vll,3)
! write(*,*) lmsize,nrmax
! allocate( vll2ddr(lmsize,nrmax) )
! allocate( dvlldr(nrmax,lmsize) )
! allocate( vpottemp(nrmax) )
! allocate( Mass(lmsize),Mass0(lmsize),dMassdr(lmsize) )
allocate( dvpotdr(nrmax),d2vpotdr2(nrmax) )
! vpot=0.0D0
dvpotdr=0.0D0
d2vpotdr2=0.0D0

! do ival=1, lmsize
!   write(8584,'(5000E)') vll(ival,ival,:)
! end do !ival

! do ival=1, lmsize
!   vpottemp=vll(ival,ival,:)
!   write(8585,'(5000E)') vpottemp

   do ir=1,nrmax
     vpot(ir)=cellnew%rmeshnew(ir)**2+1.0D0
   end do !ir


  write(8585,'(5000E)') vpot
  call  diff_cell(vpot,cellnew,dvpotdr)
  write(8586,'(5000E)') dvpotdr
  call  diff2_cell(vpot,cellnew,d2vpotdr2)
  write(8587,'(5000E)') d2vpotdr2
!   write(8586,'(5000E)') dvlldr(:,ival)
! end do !ival
! stop
! do ival=1, lmsize
!   do ir=1,nrmax
!     dvlldr(ir,ival)=2.0D0*26/cellnew%rmeshnew(ir)**2
!   end do !ir
! end do !ival
! -2.0D0*ZATOM/rmesh(ir)

! do ival=1, lmsize
!   write(8585,'(5000E)') dvlldr(ival,ival,:)
! end do !ival


do ir=1,nrmax

!   do ival=1, lmsize  
     Mass      = cone + (eryd+2*zatom/cellnew%rmeshnew(ir)    - vpot(ir)      ) /cvlight**2
     dMassdr   =        (    -2*zatom/cellnew%rmeshnew(ir)**2 - dvpotdr(ir)   ) /cvlight**2
     d2Massdr2 =        (     4*zatom/cellnew%rmeshnew(ir)**3 - d2vpotdr2(ir) ) /cvlight**2
!     dMassdr(ival) =-dvlldr(ir,ival)/cvlight**2
     Mass0     =cone+eryd/cvlight**2
!   end do

  do ival=1, lmsize  
    do ival2=1,lmsize
      vll(ival2,ival,ir)=Mass*vll(ival2,ival,ir)
    end do !ival2
  end do !ival


!    vpotsra=0.75*(-2*zatom/cellnew%rmeshnew(ir)**2/cvlight**2/Mass)**2 !-0.5*d2Massdr2/Mass-dMassdr/Mass/cellnew%rmeshnew(ir)
  vpotsra=             0.75*(dMassdr/Mass)**2 -0.5*d2Massdr2/Mass -dMassdr/Mass/cellnew%rmeshnew(ir)



  write(3244,'(50000E)') vpotsra,0.75*(dMassdr/Mass)**2-0.5*d2Massdr2/Mass-dMassdr/Mass/cellnew%rmeshnew(ir)
  write(3243,'(50000E)') 0.75*(dMassdr/Mass)**2,-0.5*d2Massdr2/Mass,-dMassdr/Mass/cellnew%rmeshnew(ir)


  do ival=1, lmsize  
    vll(ival,ival,ir)=vll(ival,ival,ir)+vpotsra+(Mass0-Mass) * eryd
  end do !ival

!   do ival=1, lmsize  
!     vll2ddr(ival,ir)=dMassdr(ival)/Mass(ival)
!     vll(ival,ival,ir)=vll(ival,ival,ir)+(Mass0(ival)-Mass(ival))*eryd-vll2ddr(ival,ir)/cellnew%rmeshnew(ir)
!   end do
end do

! do ival=1, lmsize
!   write(8587,'(5000E)') vll2ddr(ival,:)
! end do !ival

! writE(*,*) 'dealloc'
! deallocate(vll2ddr)
end subroutine vllmatsra2

end module mod_vllmatsra2
