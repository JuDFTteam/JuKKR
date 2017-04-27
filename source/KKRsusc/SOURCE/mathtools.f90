module mod_mathtools
contains
subroutine inverse(nmat,mat)
! inverts a matrix mat of size nmat x nmat
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

subroutine transpose(mat)
! transpose a matrix mat of size nmat x nmat
 implicit none
  double complex mat(:,:)
  double complex temp
  integer         :: ndim1,ndim2
  integer         :: idim1,idim2
ndim1=ubound(mat,1)
ndim2=ubound(mat,2)
if (ndim1/=ndim2) stop '[transpose] matrix is not a square matrix'
do idim1=1,ndim1
  do idim2=1,idim1-1
    temp             =  mat(idim1,idim2)
    mat(idim1,idim2) =  mat(idim2,idim1)
    mat(idim2,idim1) =  temp
  end do
end do
end subroutine transpose

subroutine transpose_dm(mat)
! inverts a matrix mat of size nmat x nmat - double precision
 implicit none
  double precision mat(:,:)
  double precision temp
  integer         :: ndim1,ndim2
  integer         :: idim1,idim2
ndim1=ubound(mat,1)
ndim2=ubound(mat,2)
if (ndim1/=ndim2) stop '[transpose] matrix is not a square matrix'
do idim1=1,ndim1
  do idim2=1,idim1-1
    temp             =  mat(idim1,idim2)
    mat(idim1,idim2) =  mat(idim2,idim1)
    mat(idim2,idim1) =  temp
  end do
end do
end subroutine transpose_dm

subroutine conjugate2(mat)
 implicit none
  double complex mat(:,:)
  integer         :: ndim1,ndim2
  integer         :: idim1,idim2
ndim1=ubound(mat,1)
ndim2=ubound(mat,2)

do idim2=1,ndim2
  do idim1=1,ndim1
    mat(idim1,idim2) =  conjg( mat(idim1,idim2) )
  end do
end do

end subroutine conjugate2

subroutine conjugate3(mat)
 implicit none
  double complex mat(:,:,:)
  integer         :: ndim1,ndim2,ndim3
  integer         :: idim1,idim2,idim3

ndim1=ubound(mat,1)
ndim2=ubound(mat,2)
ndim3=ubound(mat,3)

do idim3=1,ndim3
  do idim2=1,ndim2
    do idim1=1,ndim1
      mat(idim1,idim2,idim3) =  conjg( mat(idim1,idim2,idim3) )
    end do
  end do
end do

end subroutine conjugate3

subroutine conjugate4(mat)
 implicit none
  double complex mat(:,:,:,:)
  integer         :: ndim1,ndim2,ndim3,ndim4
  integer         :: idim1,idim2,idim3,idim4

ndim1=ubound(mat,1)
ndim2=ubound(mat,2)
ndim3=ubound(mat,3)
ndim4=ubound(mat,4)

do idim4=1,ndim4
  do idim3=1,ndim3
    do idim2=1,ndim2
      do idim1=1,ndim1
        mat(idim1,idim2,idim3,idim4) =  conjg( mat(idim1,idim2,idim3,idim4) )
      end do
    end do
  end do
end do

end subroutine conjugate4


subroutine linearsolve_dc(Amat,bmat)!,ndim,lmgf0d,ngd)
  use nrtype
  implicit none
!interface variables
      complex(kind=dpc),intent(in)        :: Amat(:,:)
      complex(kind=dpc),intent(inout)     :: bmat(:,:)
!local variables
      integer                             :: nrow,ncol
      integer                             :: info
      integer                             :: temp(size(bmat,1))

nrow=size(bmat,1)
ncol=size(bmat,2)

if (nrow/=size(Amat,1) .or. nrow/=size(Amat,2)) then
   stop '[linearsolve] dimension error while solving Ax=b' 
end if
! allocate (temp(nrow))

!--------------------------------------------------
!-- solve the system of linear equations         --
!--------------------------------------------------
call zgetrf(nrow,nrow,Amat,nrow,temp,info)
      if(info.ne.0) then 
        WRITE(*,*) 'info',info
        stop '[linearsolve_dc] zgetrf failed' 
      endif
call zgetrs('n',nrow,ncol,Amat,nrow,temp,bmat,nrow,info)
      if(info.ne.0) then 
        WRITE(*,*) 'info',info
        stop '[linearsolve_dc] zgetrf failed' 
      endif
! deallocate (temp)

end subroutine

subroutine linearsolve2_dc(Amat,bmat,cmat)!,ndim,lmgf0d,ngd)
  use nrtype
  implicit none
!interface variables
      complex(kind=dpc),intent(in)        :: Amat(:,:)
      complex(kind=dpc),intent(inout)     :: bmat(:,:)
      complex(kind=dpc),intent(inout)     :: cmat(:,:)
!local variables
      integer                             :: nrow,ncol
      integer                             :: info
      integer                             :: temp(size(bmat,1))

nrow=size(bmat,1)
ncol=size(bmat,2)

if (nrow/=size(Amat,1) .or. nrow/=size(Amat,2)) then
   stop '[linearsolve] dimension error while solving Ax=b' 
end if
if (nrow/=size(cmat,1) .or. ncol/=size(cmat,2)) then
   stop '[linearsolve] dimension error while solving Ax=b' 
end if

! allocate (temp(nrow))

!--------------------------------------------------
!-- solve the system of linear equations         --
!--------------------------------------------------
call zgetrf(nrow,nrow,Amat,nrow,temp,info)
      if(info.ne.0) then 
        WRITE(*,*) 'info',info
        stop '[linearsolve_dc] zgetrf failed' 
      endif
call zgetrs('n',nrow,ncol,Amat,nrow,temp,bmat,nrow,info)
      if(info.ne.0) then 
        WRITE(*,*) 'info',info
        stop '[linearsolve_dc] zgetrf failed' 
      endif
call zgetrs('n',nrow,ncol,Amat,nrow,temp,cmat,nrow,info)
      if(info.ne.0) then 
        WRITE(*,*) 'info',info
        stop '[linearsolve_dc] zgetrf failed' 
      endif
! deallocate (temp)

end subroutine



      function matmat(mat1,mat2)
      implicit none
      double complex, intent(in) :: mat1(:,:),mat2(:,:)
      double complex             :: matmat(size(mat1,1),size(mat2,2))
      integer                :: n1,n,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,2)
      if(size(mat2,1).ne.n) stop 'matmat: dimensions of matrices are inconsistent.'
      call zgemm('N','N',n1,n2,n,(1d0,0d0),mat1,n1,mat2,n,(0d0,0d0),matmat,n1)
      end function matmat

      function matmat1T(mat1,mat2)
      implicit none
      complex(8), intent(in) :: mat1(:,:),mat2(:,:)
      complex(8)             :: matmat1T(size(mat1,1),size(mat2,1))
      integer                :: n11,n12,n21,n22
      n11 = size(mat1,1)
      n12  = size(mat1,2)
      n21 = size(mat2,1)
      n22 = size(mat2,2)
      if(n12.ne.n22) stop 'matmat1T: dimensions of matrices are inconsistent.'
      call zgemm('N','T',n11,n21,n12,(1d0,0d0),mat1,n11,mat2,n11,(0d0,0d0),matmat1T,n21)
      end function matmat1T

      function matmatT1(mat1,mat2)
      implicit none
      double complex, intent(in) :: mat1(:,:),mat2(:,:)
      double complex             :: matmatT1(size(mat1,2),size(mat2,2))
      integer                    :: n11,n12,n21,n22
      n11 = size(mat1,1)
      n12 = size(mat1,2)
      n21 = size(mat2,1)
      n22 = size(mat2,2)
      if(n11.ne.n21) stop 'matmatT1: dimensions of matrices are inconsistent.'
      call zgemm('T','N',n12,n22,n11,(1d0,0d0),mat1,n11,mat2,n21,(0d0,0d0),matmatT1,n22)
      end function matmatT1

      function matmat_dmdm(mat1,mat2)
      implicit none
      double precision, intent(in) :: mat1(:,:),mat2(:,:)
      double precision             :: matmat_dmdm(size(mat1,1),size(mat2,2))
      integer                :: n1,n,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,2)
      if(size(mat2,1).ne.n) stop 'matmat: dimensions of matrices are inconsistent.'
      call dgemm('N','N',n1,n2,n,(1d0,0d0),mat1,n1,mat2,n,(0d0,0d0),matmat_dmdm,n1)
      end function matmat_dmdm

      function matmat1T_dmdm(mat1,mat2)
      implicit none
      double precision, intent(in) :: mat1(:,:),mat2(:,:)
      double precision             :: matmat1T_dmdm(size(mat1,1),size(mat2,1))
      integer                :: n11,n12,n21,n22
      n11 = size(mat1,1)
      n12  = size(mat1,2)
      n21 = size(mat2,1)
      n22 = size(mat2,2)
      if(n12.ne.n22) stop 'matmat1T: dimensions of matrices are inconsistent.'
      call dgemm('N','T',n11,n21,n12,(1d0,0d0),mat1,n11,mat2,n11,(0d0,0d0),matmat1T_dmdm,n21)
      end function matmat1T_dmdm

      function matmatT1_dmdm(mat1,mat2)
      implicit none
      double precision, intent(in) :: mat1(:,:),mat2(:,:)
      double precision             :: matmatT1_dmdm(size(mat1,2),size(mat2,2))
      integer                    :: n11,n12,n21,n22
      n11 = size(mat1,1)
      n12 = size(mat1,2)
      n21 = size(mat2,1)
      n22 = size(mat2,2)
      if(n11.ne.n21) stop 'matmatT1: dimensions of matrices are inconsistent.'
      call dgemm('T','N',n12,n22,n11,(1d0,0d0),mat1,n11,mat2,n21,(0d0,0d0),matmatT1_dmdm,n22)
      end function matmatT1_dmdm


      function matvec_dmdm(mat1,vec1)
      implicit none
      real(8), intent(in) :: mat1(:,:),vec1(:)
      real(8)             :: matvec_dmdm(size(mat1,1))
      integer             :: n,m
      m = size(mat1,1)
      n = size(mat1,2)
      if(size(vec1,1).ne.n) stop 'matvec_dmdm: dimensions of first input array differ.'
!       if(size(mat2,1).ne.n) stop 'matmat_dmdm: second input array has wrong dimensions.'
!       if(size(mat2,2).ne.n) stop 'matmat_dmdm: dimensions of second input array differ.'
      call DGEMV('N',M,N,1.0D0,mat1,M,vec1,1,0.0D0,matvec_dmdm,1)
      end function matvec_dmdm


      function matvec_dzdz(mat1,vec1)
      implicit none
      double complex, intent(in) :: mat1(:,:),vec1(:)
      double complex             :: matvec_dzdz(size(mat1,1))
      integer             :: n,m
      m = size(mat1,1)
      n = size(mat1,2)
      if(size(vec1,1).ne.n) stop 'matvec_dzdz: dimensions of first input array differ.'
!       if(size(mat2,1).ne.n) stop 'matmat_dmdm: second input array has wrong dimensions.'
!       if(size(mat2,2).ne.n) stop 'matmat_dmdm: dimensions of second input array differ.'
      call ZGEMV('N',M,N,(1.0D0,0.0D0),mat1,M,vec1,1,(0.0D0,0.0D0),matvec_dzdz,1)
      end function matvec_dzdz




function cross(a,b) result(axb)
 
implicit none
integer,parameter :: wp=selected_real_kind(15, 307) !double precision
real(wp),dimension(3) :: axb
real(wp),dimension(3),intent(in) :: a
real(wp),dimension(3),intent(in) :: b 
 
axb(1) = a(2)*b(3) - a(3)*b(2)
axb(2) = a(3)*b(1) - a(1)*b(3)
axb(3) = a(1)*b(2) - a(2)*b(1)
 
end function cross




function  rotvector(theta1,phi1,m2,fac1)
implicit none
real(kind=8) theta1,phi1
real(kind=8) m2(3)
real(kind=8) fac1
real(kind=8) rotvector(3)

real(kind=8) :: absm2,e1(3),m1(3),nvec(3),absnvec,alpha
real(kind=8) :: cosa,sina,n1,n2,n3,rotmat(3,3)



absm2=dsqrt(m2(1)**2+m2(2)**2+m2(3)**2)
print *,'absm2',absm2
e1=(/dcos(phi1)*dsin(theta1),dsin(phi1)*dsin(theta1),dcos(theta1) /)
print *,'e1',e1
m1=absm2*e1
print *,'m1',m1

nvec=cross(m1,m2)/absm2**2
print *,'nvec',nvec
absnvec=dsqrt ( nvec(1)**2+nvec(2)**2+nvec(3)**2)
print *,'absnvec',absnvec

! This gives errors for theta>180deg
!alpha = dasin( absnvec ) 

! better do this
alpha=dacos(DOT_PRODUCT(m1,m2)/absm2**2)

print *,'alpha',alpha, 'and in deg',alpha/3.1515926*180

alpha=alpha*fac1
print *,'new alpha',alpha

 cosa=dcos(alpha)
 sina=dsin(alpha)
 n1=nvec(1)/absnvec
 n2=nvec(2)/absnvec
 n3=nvec(3)/absnvec

rotmat(1,:) =  (/ n1*n1*(1-cosa) +    cosa, n1*n2*(1-cosa) - n3*sina, n1*n3*(1-cosa)+n2*sina /)
rotmat(2,:) =  (/ n2*n1*(1-cosa) + n3*sina, n2*n2*(1-cosa) +    cosa, n2*n3*(1-cosa)-n1*sina /)
rotmat(3,:) =  (/ n3*n1*(1-cosa) - n2*sina, n3*n2*(1-cosa) + n1*sina, n3*n3*(1-cosa)+cosa    /) 


rotvector=matvec_dmdm(rotmat,m1)

print *,m2
print *,rotvector

! stop'rotvector'
end function rotvector

end module mod_mathtools





