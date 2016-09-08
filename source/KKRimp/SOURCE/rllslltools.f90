module mod_rllslltools
contains

subroutine rllslltools_solvesra(nmat,ns,mat,s1,s2)
! use mod_wrapper
use mod_mathtools
use mod_timing
implicit none
double complex :: mat(nmat,nmat),s1(nmat,ns),s2(nmat,ns)
double complex :: sout(nmat,ns)
integer :: nmat,nmat2,ns
double complex,allocatable :: mattemp(:,:),mattemp2(:,:)

nmat=ubound(mat,1)
nmat2=nmat/2
allocate(mattemp(nmat2,nmat2),mattemp2(nmat2,nmat2))

!set 1x1 block
mattemp=mat(1:nmat/2,1:nmat/2)

call timing_start('inverse')
call inverse(nmat2,mattemp)
call timing_stop('inverse')

mat(1:nmat2,1:nmat2)=mattemp

!set 2x2 block
mattemp2=-mat(1+nmat2:nmat,1:nmat2)
mattemp2=matmat_zmzm(mattemp2,mattemp)
mat(1+nmat/2:nmat,1:nmat/2)=mattemp2

!calc mat**-1 *s1
sout=matmat_zmzm(mat,s1)
s1=sout

!calc mat**-1 *s1
sout=matmat_zmzm(mat,s2)
s2=sout

end subroutine rllslltools_solvesra !(mat,s1,s2)

      function matmat_zmzm(mat1,mat2)
      implicit none
      complex(8), intent(in) :: mat1(:,:),mat2(:,:)
      complex(8)             :: matmat_zmzm(size(mat1,1),size(mat2,2))
      integer                :: n1,n,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,2)
      if(size(mat2,1).ne.n) stop 'matmat_zmzm: dimensions of matrices are inconsistent.'
      call zgemm('N','N',n1,n2,n,(1d0,0d0),mat1,n1,mat2,n,(0d0,0d0),matmat_zmzm,n1)
      end function matmat_zmzm

subroutine inverse(nmat,mat,work,ipiv)
!interface
integer        :: nmat
double complex :: mat(nmat,nmat)
double complex :: work(nmat,nmat)
!local
integer        :: IPIV(nmat)
integer        :: info

call ZGETRF( nmat, nmat, mat, nmat, IPIV, INFO )
if (info/=0) stop '[inverse] error INFO' 
call ZGETRI( nmat, mat, nmat, IPIV, WORK, nmat*nmat, INFO )
if (info/=0) stop '[inverse] error INFO' 
end subroutine inverse


subroutine rllslltools(dim1,dimblock,nblock,mat1,mat2,mat3)
use mod_mathtools
implicit none

integer          :: dim1
integer          :: dimblock
integer          :: nblock
double complex :: mat1(dim1,dim1)
double complex :: mat2(dim1,nblock)
double complex :: mat3(dim1,nblock)

!local
double complex,allocatable :: matsmall(:,:)
double complex,allocatable :: mat2small(:,:)
double complex,allocatable :: mat3small(:,:)
integer          :: ival1, ival2,ival3, ival4
integer          :: iblock
! nblock = dim1/dimblock
if (nblock*dimblock/=dim1) stop 'error in rllslltools'

allocate ( matsmall(dimblock,dimblock), &
           mat2small(dimblock,nblock), &
           mat3small(dimblock,nblock) )

! mat2=(0.0D0,0.0D0)
! mat3=(0.0D0,0.0D0)
! print *,dim1
! print *,dimblock
! print *,nblock

do iblock = 1, nblock
!          print *,nblock
        ival1=1+(iblock-1)*dimblock
        ival2=   iblock   *dimblock

        ival3=1+(iblock-1)*dimblock
        ival4=   iblock   *dimblock


!         print *,ival1,ival2
        matsmall = mat1(ival1:ival2 , ival1:ival2)
!         print *,ival1,ival2
!         CALL ZGETRF(dimblock,dimblock,matsmall,dimblock,IPIV,INFO)


        mat2small            = mat2(ival3:ival4 , :)
        mat3small            = mat3(ival1:ival2 , :)
        call linearsolve2_dc(matsmall,mat2small,mat3small)
!         CALL ZGETRS('N',dimblock,nblock,matsmall,dimblock,IPIV,mat2small,dimblock,INFO)
        mat2(ival3:ival4 , :)= mat2small 
        mat3(ival3:ival4 , :)= mat3small 

!         call linearsolve_dc(matsmall,mat3small)
!         CALL ZGETRS('N',dimblock,nblock,matsmall,dimblock,IPIV,mat3small,dimblock,INFO)
end do !nblock
! print *,'end'
!         CALL ZGETRF(NPLM,NPLM,SLV,NPLM,IPIV,INFO)
!         CALL ZGETRS('N',NPLM,LMMAX,SLV,NPLM,IPIV,YRLL,NPLM,INFO)
!         CALL ZGETRS('N',NPLM,LMMAX,SLV,NPLM,IPIV,ZRLL,NPLM,INFO)
end subroutine rllslltools

end module mod_rllslltools
