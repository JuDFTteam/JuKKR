module mod_checkinterpolation
!-------------------------------------------------------------------------------
!> Summary: Checks the validity of the interpolation between 'old' mesh
!> and Chebyshev mesh
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Summary: Checks the validity of the interpolation between 'old' mesh
!> and Chebyshev mesh
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine checkinterpolation(cell,cellnew,vpotin,config)!,lmpot)
use type_cell, only: cell_type
use type_cellnew, only: cell_typenew
use type_config, only: config_type
use mod_interpolatecell, only: interpolatecell
use mod_cheb2oldgrid, only: cheb2oldgrid

implicit none
type(cell_type)    :: cell
type(cell_typenew) :: cellnew
double precision   :: vpotin(cell%nrmax,1)
type(config_type)  :: config
!local
double precision,allocatable   :: testpot(:,:)
double complex,allocatable   :: testpot2(:,:)
double complex               :: testpot_out(cell%nrmax)
double precision               :: testpot_out2(cell%nrmax)
double precision   :: rms
integer            :: ir

testpot_out=0.0D0
testpot_out2=0.0D0
print *,'Interpol old->new'
call interpolatecell(1,VPOTIN,1,CELL,1,cellnew,1,1,config,'test',testpot)

print *,'allocate ',ubound(testpot,1),ubound(testpot,2)
allocate (testpot2(ubound(testpot,1),ubound(testpot,2)))
testpot2=testpot

print *,'Interpol new->old'
call cheb2oldgrid(cell%nrmax, cellnew%nrmaxnew, 1, cell%rmesh, cellnew%ncheb, cellnew%npan_tot, cellnew%rpan_intervall, cellnew%ipan_intervall, testpot2, testpot_out, cell%nrmax)

testpot_out2=testpot_out

print *,'Calculate rms'
rms=0
print *, 'vpotin        vpotout'
do ir=1,cell%nrmax
  write(*,'(I5,2F25.14,E25.14)'),ir,vpotin(ir,1),testpot_out2(ir),(vpotin(ir,1)-testpot_out2(ir))
  write(112,'(5I5,2F25.14,E25.14)'),ir,vpotin(ir,1),testpot_out2(ir),(vpotin(ir,1)-testpot_out2(ir))
  rms = rms + (vpotin(ir,1)-testpot_out2(ir))**2
end do
rms = sqrt(rms)

write(*,*) 'Interpolation rms error is ',rms/cell%nrmax
write(*,*) 'done'
stop

end subroutine checkinterpolation

end module mod_checkinterpolation
