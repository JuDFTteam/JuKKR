module mod_umatrix

contains

subroutine umatrix(natom,lmaxhost,umat,dimgmathost,svec,eryd)
use mod_jmtrx
implicit none
integer :: natom
integer :: lmaxhost
integer :: dimgmathost
double complex :: umat(dimgmathost,dimgmathost)
double precision :: svec(3,natom)
double complex :: eryd
!local
integer       :: iatom,ilm,lmmax
logical       :: lcall
double complex :: umatrix_small( (lmaxhost+1)**2,(lmaxhost+1)**2 )

lcall=.False.
lmmax=(lmaxhost+1)**2
umat=(0.0D0,0.0D0)
do iatom=1,natom
 umatrix_small=(0.0D0,0.0D0)
 call  JMTRX(svec(1,iatom),svec(2,iatom),svec(3,iatom),eryd,lmaxhost,umatrix_small,lcall)
 do ilm=1,lmmax
   umat((iatom-1)*lmmax+1:(iatom)*lmmax,(iatom-1)*lmmax+ilm)=umatrix_small(:,ilm)
 end do

end do

end subroutine umatrix





end module mod_umatrix