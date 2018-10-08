 include 'dsort.f90'
program test2
use MOD_DSORT
real*8 :: array(6)
integer :: index1(6),pos

 array=(/1.0,5.0,2.0,3.0,10.0,-1.0/)
 call dsort(array,index1,6,pos)

 write(*,*) array
 write(*,*) index1






end program test2
