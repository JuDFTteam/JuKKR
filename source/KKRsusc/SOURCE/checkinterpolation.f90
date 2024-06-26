module mod_checkinterpolation

contains

subroutine checkinterpolation(cell,cellnew,vpotin,config)!,lmpot)
use type_cell
use type_cellnew
use type_config
use mod_interpolatecell
use mod_cheb2oldgridc

implicit none
type(cell_type)    :: cell
type(cell_typenew) :: cellnew
double precision   :: vpotin(cell%nrmax,1)
type(config_type)  :: config
! integer             :: lmpot
!local
double precision,allocatable   :: testpot(:,:)
double complex,allocatable   :: testpot2(:,:)
double complex               :: testpot_out(cell%nrmax)
double precision               :: testpot_out2(cell%nrmax)
double precision   :: rms
integer            :: ir
! do ilm = 1 , 
testpot_out=0.0D0
testpot_out2=0.0D0
print *,'Interpol old->new'
call interpolatecell(1,VPOTIN,1,CELL,1,cellnew,1,1,config,'test',testpot)

print *,'allocate ',ubound(testpot,1),ubound(testpot,2)
allocate (testpot2(ubound(testpot,1),ubound(testpot,2)))
testpot2=testpot

! allocate ( testpot_out(cellnew%nrmaxnew) )
! allocate ( testpot_out2(cellnew%nrmaxnew) )

print *,'Interpol new->old'
call cheb2oldgridc(cell,cellnew,cellnew%ncheb,1,testpot2,testpot_out)

testpot_out2=testpot_out

print *,'Calculate rms'
rms=0
print *, 'vpotin        vpotout'
do ir=1,cell%nrmax
  write(*,'(I,2F,E)'),ir,vpotin(ir,1),testpot_out2(ir),(vpotin(ir,1)-testpot_out2(ir))
  write(112,'(I,2F,E)'),ir,vpotin(ir,1),testpot_out2(ir),(vpotin(ir,1)-testpot_out2(ir))
  rms = rms + (vpotin(ir,1)-testpot_out2(ir))**2
end do
rms = sqrt(rms)

write(*,*) 'Interpolation rms error is ',rms/cell%nrmax
write(*,*) 'done'
stop

end subroutine checkinterpolation

end module mod_checkinterpolation
