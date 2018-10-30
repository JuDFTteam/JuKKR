module mod_vllmatsra2

contains

subroutine vllmatsra2(vll,vll2ddr,cellnew,eryd)
!************************************************************************************
! The perubation matrix for the SRA-equations are set up
!************************************************************************************
use mod_physic_params, only: cvlight
use mod_mathtools, only: inverse
use type_cellnew
use mod_chebyshev, only : diff_cell_complex
implicit none
!interface
  DOUBLE COMPLEX             :: VLL(:,:,:)
  DOUBLE COMPLEX,allocatable :: VLL2DDR(:,:)
  type(cell_typenew)          :: cellnew
!   double precision,allocatable :: rmesh(:)
  double complex              :: eryd
  integer                      :: lmsize
!local
  double complex,allocatable :: vpottemp(:),dvlldr(:,:)
  integer                     :: nrmax
  integer                     :: ival,ival2,ir
  double complex,parameter    :: cone=(1.0D0,0.0D0)
  double complex,parameter    :: czero=(0.0D0,0.0D0)
!  double precision,parameter :: cvlight = 274.0720442d0
  double complex,allocatable   :: Mass(:),Mass0(:),dMassdr(:)

lmsize=ubound(vll,1)
nrmax=ubound(vll,3)
! write(*,*) lmsize,nrmax
allocate( vll2ddr(lmsize,nrmax) )
allocate( dvlldr(nrmax,lmsize) )
allocate( vpottemp(nrmax) )
allocate( Mass(lmsize),Mass0(lmsize),dMassdr(lmsize) )

! do ival=1, lmsize
!   write(8584,'(5000E)') vll(ival,ival,:)
! end do !ival

! do ival=1, lmsize
!   vpottemp=vll(ival,ival,:)
!   write(8585,'(5000E)') vpottemp
!   call  diff_cell_complex(vpottemp,cellnew,dvlldr(:,ival))
!   write(8586,'(5000E)') dvlldr(:,ival)
! end do !ival

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
     Mass(ival) =cone + (eryd+2*zatom/cellnew%rmeshnew(ir) + vll(ival,ival,ir) ) /cvlight**2
     dMassdr(ival) =cone + (eryd+2*zatom/cellnew%rmeshnew(ir) + vll(ival,ival,ir) ) /cvlight**2
!     dMassdr(ival) =-dvlldr(ir,ival)/cvlight**2
!     Mass0(ival)=cone+eryd/cvlight**2
!   end do

  do ival=1, lmsize  
    do ival2=1,lmsize
      vll(ival2,ival,ir)=Mass(ival2)*vll(ival2,ival,ir)
    end do !ival2
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
