!------------------------------------------------------------------------------------
!> Summary: Math tools collection
!> Author: 
!> 
!------------------------------------------------------------------------------------
module mod_mathtools

contains

  !-------------------------------------------------------------------------------
  !> Summary: Inverts a double complex matrix mat of size nmat x nmat
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  subroutine inverse(nmat,mat)
    ! inverts a matrix mat of size nmat x nmat
    use mod_datatypes, only: dp
    implicit none
    !interface
    integer        :: nmat
    complex (kind=dp) :: mat(nmat,nmat)
    complex (kind=dp) :: work(nmat,nmat)
    !local
    integer        :: IPIV(nmat)
    integer        :: info
   
    call ZGETRF( nmat, nmat, mat, nmat, IPIV, INFO )
    if (info/=0) stop '[inverse] error INFO' 
    call ZGETRI( nmat, mat, nmat, IPIV, WORK, nmat*nmat, INFO )
    if (info/=0) stop '[inverse] error INFO' 
  end subroutine inverse

  !-------------------------------------------------------------------------------
  !> Summary: Transposes a matrix mat of size nmat x nmat
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  subroutine transpose(mat)
    use mod_datatypes, only: dp
    implicit none
    complex (kind=dp) mat(:,:)
    complex (kind=dp) temp
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

  !-------------------------------------------------------------------------------
  !> Summary: Transposes a matrix mat of size nmat x nmat - real (kind=dp)
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  subroutine transpose_dm(mat)
    use mod_datatypes, only: dp
    implicit none
    real (kind=dp) mat(:,:)
    real (kind=dp) temp
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

  !-------------------------------------------------------------------------------
  !> Summary: Conjugates a 2-dim matrix mat
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  subroutine conjugate2(mat)
    use mod_datatypes, only: dp
    implicit none
    complex (kind=dp) mat(:,:)
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

  !-------------------------------------------------------------------------------
  !> Summary: Conjugates a 3-dim matrix mat
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  subroutine conjugate3(mat)
    use mod_datatypes, only: dp
    implicit none
    complex (kind=dp) mat(:,:,:)
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

  !-------------------------------------------------------------------------------
  !> Summary: Conjugates a 4-dim matrix mat
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  subroutine conjugate4(mat)
    use mod_datatypes, only: dp
    implicit none
    complex (kind=dp) mat(:,:,:,:)
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

  !-------------------------------------------------------------------------------
  !> Summary: Solves the system of linear equations
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  subroutine linearsolve_dc(Amat,bmat)!,ndim,lmgf0d,ngd)
    use mod_datatypes, only: dp
    implicit none
    !interface variables
    complex(kind=dp),intent(in)        :: Amat(:,:)
    complex(kind=dp),intent(inout)     :: bmat(:,:)
    !local variables
    integer                             :: nrow,ncol
    integer                             :: info
    integer                             :: temp(size(bmat,1))
   
    nrow=size(bmat,1)
    ncol=size(bmat,2)
   
    if (nrow/=size(Amat,1) .or. nrow/=size(Amat,2)) then
      stop '[linearsolve] dimension error while solving Ax=b' 
    end if
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
  end subroutine linearsolve_dc

  !-------------------------------------------------------------------------------
  !> Summary: Solves the system of linear equations
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  subroutine linearsolve2_dc(Amat,bmat,cmat)!,ndim,lmgf0d,ngd)
    use mod_datatypes, only: dp
    implicit none
    !interface variables
    complex(kind=dp),intent(in)        :: Amat(:,:)
    complex(kind=dp),intent(inout)     :: bmat(:,:)
    complex(kind=dp),intent(inout)     :: cmat(:,:)
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
  end subroutine linearsolve2_dc

  !-------------------------------------------------------------------------------
  !> Summary: Multiplies matrices mat1 and mat2 with zgemm: N N
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  function matmat(mat1,mat2)
    use mod_datatypes, only: dp
    use mod_constants, only: cone, czero
    implicit none
    complex (kind=dp), intent(in) :: mat1(:,:),mat2(:,:)
    complex (kind=dp)             :: matmat(size(mat1,1),size(mat2,2))
    integer                :: n1,n,n2
    n1 = size(mat1,1)
    n  = size(mat1,2)
    n2 = size(mat2,2)
    if(size(mat2,1).ne.n) stop 'matmat: dimensions of matrices are inconsistent.'
    call zgemm('N','N',n1,n2,n,cone,mat1,n1,mat2,n,czero,matmat,n1)
  end function matmat

  !-------------------------------------------------------------------------------
  !> Summary: Multiplies matrices mat1 and mat2 with zgemm: N T
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  function matmat1T(mat1,mat2)
    use mod_datatypes, only: dp
    use mod_constants, only: cone, czero
    implicit none
    complex(kind=dp), intent(in) :: mat1(:,:),mat2(:,:)
    complex(kind=dp)             :: matmat1T(size(mat1,1),size(mat2,1))
    integer                :: n11,n12,n21,n22
    n11 = size(mat1,1)
    n12  = size(mat1,2)
    n21 = size(mat2,1)
    n22 = size(mat2,2)
    if(n12.ne.n22) stop 'matmat1T: dimensions of matrices are inconsistent.'
    call zgemm('N','T',n11,n21,n12,cone,mat1,n11,mat2,n11,czero,matmat1T,n21)
  end function matmat1T

  !-------------------------------------------------------------------------------
  !> Summary: Multiplies matrices mat1 and mat2 with zgemm: T N
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  function matmatT1(mat1,mat2)
    use mod_datatypes, only: dp
    use mod_constants, only: cone, czero
    implicit none
    complex (kind=dp), intent(in) :: mat1(:,:),mat2(:,:)
    complex (kind=dp)             :: matmatT1(size(mat1,2),size(mat2,2))
    integer                    :: n11,n12,n21,n22
    n11 = size(mat1,1)
    n12 = size(mat1,2)
    n21 = size(mat2,1)
    n22 = size(mat2,2)
    if(n11.ne.n21) stop 'matmatT1: dimensions of matrices are inconsistent.'
    call zgemm('T','N',n12,n22,n11,cone,mat1,n11,mat2,n21,czero,matmatT1,n22)
  end function matmatT1

  !-------------------------------------------------------------------------------
  !> Summary: Multiplies matrices mat1 and mat2 with dgemm: N N
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  function matmat_dmdm(mat1,mat2)
    use mod_datatypes, only: dp
    use mod_constants, only: cone, czero
    implicit none
    real (kind=dp), intent(in) :: mat1(:,:),mat2(:,:)
    real (kind=dp)             :: matmat_dmdm(size(mat1,1),size(mat2,2))
    integer                :: n1,n,n2
    n1 = size(mat1,1)
    n  = size(mat1,2)
    n2 = size(mat2,2)
    if(size(mat2,1).ne.n) stop 'matmat: dimensions of matrices are inconsistent.'
    call dgemm('N','N',n1,n2,n,cone,mat1,n1,mat2,n,czero,matmat_dmdm,n1)
  end function matmat_dmdm

  !-------------------------------------------------------------------------------
  !> Summary: Multiplies matrices mat1 and mat2 with dgemm: N T
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  function matmat1T_dmdm(mat1,mat2)
    use mod_datatypes, only: dp
    use mod_constants, only: cone, czero
    implicit none
    real (kind=dp), intent(in) :: mat1(:,:),mat2(:,:)
    real (kind=dp)             :: matmat1T_dmdm(size(mat1,1),size(mat2,1))
    integer                :: n11,n12,n21,n22
    n11 = size(mat1,1)
    n12  = size(mat1,2)
    n21 = size(mat2,1)
    n22 = size(mat2,2)
    if(n12.ne.n22) stop 'matmat1T: dimensions of matrices are inconsistent.'
    call dgemm('N','T',n11,n21,n12,cone,mat1,n11,mat2,n11,czero,matmat1T_dmdm,n21)
  end function matmat1T_dmdm

  !-------------------------------------------------------------------------------
  !> Summary: Multiplies matrices mat1 and mat2 with dgemm: T N
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  function matmatT1_dmdm(mat1,mat2)
    use mod_datatypes, only: dp
    use mod_constants, only: cone, czero
    implicit none
    real (kind=dp), intent(in) :: mat1(:,:),mat2(:,:)
    real (kind=dp)             :: matmatT1_dmdm(size(mat1,2),size(mat2,2))
    integer                    :: n11,n12,n21,n22
    n11 = size(mat1,1)
    n12 = size(mat1,2)
    n21 = size(mat2,1)
    n22 = size(mat2,2)
    if(n11.ne.n21) stop 'matmatT1: dimensions of matrices are inconsistent.'
    call dgemm('T','N',n12,n22,n11,cone,mat1,n11,mat2,n21,czero,matmatT1_dmdm,n22)
  end function matmatT1_dmdm

  !-------------------------------------------------------------------------------
  !> Summary: Multiplies matrix mat1 and vector vec1 with dgemv: N
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  function matvec_dmdm(mat1,vec1)
    use mod_datatypes, only: dp
    implicit none
    real (kind=dp), intent(in) :: mat1(:,:),vec1(:)
    real (kind=dp)             :: matvec_dmdm(size(mat1,1))
    integer             :: n,m
    m = size(mat1,1)
    n = size(mat1,2)
    if(size(vec1,1).ne.n) stop 'matvec_dmdm: dimensions of first input array differ.'
    ! if(size(mat2,1).ne.n) stop 'matmat_dmdm: second input array has wrong dimensions.'
    ! if(size(mat2,2).ne.n) stop 'matmat_dmdm: dimensions of second input array differ.'
    call DGEMV('N',M,N,1.0D0,mat1,M,vec1,1,0.0D0,matvec_dmdm,1)
  end function matvec_dmdm

  !-------------------------------------------------------------------------------
  !> Summary: Multiplies matrix mat1 and vector vec1 with zgemv: N
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  function matvec_dzdz(mat1,vec1)
    use mod_datatypes, only: dp
    use mod_constants, only: cone, czero
    implicit none
    complex (kind=dp), intent(in) :: mat1(:,:),vec1(:)
    complex (kind=dp)             :: matvec_dzdz(size(mat1,1))
    integer             :: n,m
    m = size(mat1,1)
    n = size(mat1,2)
    if(size(vec1,1).ne.n) stop 'matvec_dzdz: dimensions of first input array differ.'
    ! if(size(mat2,1).ne.n) stop 'matmat_dmdm: second input array has wrong dimensions.'
    ! if(size(mat2,2).ne.n) stop 'matmat_dmdm: dimensions of second input array differ.'
    call ZGEMV('N',M,N,cone,mat1,M,vec1,1,czero,matvec_dzdz,1)
  end function matvec_dzdz

  !-------------------------------------------------------------------------------
  !> Summary: Cross product of 2 vectors
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  function cross(a,b) result(axb)
    use mod_datatypes, only: dp
    implicit none
    real(kind=dp),dimension(3) :: axb
    real(kind=dp),dimension(3),intent(in) :: a
    real(kind=dp),dimension(3),intent(in) :: b 
    
    axb(1) = a(2)*b(3) - a(3)*b(2)
    axb(2) = a(3)*b(1) - a(1)*b(3)
    axb(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

  !-------------------------------------------------------------------------------
  !> Summary: Rotates vector m2 by theta1 and phi1 angles
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !-------------------------------------------------------------------------------
  function  rotvector(theta1,phi1,m2,fac1)
    use mod_datatypes, only: dp
    implicit none
    real(kind=dp) theta1,phi1
    real(kind=dp) m2(3)
    real(kind=dp) fac1
    real(kind=dp) rotvector(3)
    
    real(kind=dp) :: absm2,e1(3),m1(3),nvec(3),absnvec,alpha
    real(kind=dp) :: cosa,sina,n1,n2,n3,rotmat(3,3)
    
    absm2=sqrt(m2(1)**2+m2(2)**2+m2(3)**2)
    print *,'absm2',absm2
    e1=(/cos(phi1)*sin(theta1),sin(phi1)*sin(theta1),cos(theta1) /)
    print *,'e1',e1
    m1=absm2*e1
    print *,'m1',m1
    
    nvec=cross(m1,m2)/absm2**2
    print *,'nvec',nvec
    absnvec=sqrt ( nvec(1)**2+nvec(2)**2+nvec(3)**2)
    print *,'absnvec',absnvec
    
    ! This gives errors for theta>180deg
    !alpha = dasin( absnvec ) 
    ! better do this
    alpha=acos(DOT_PRODUCT(m1,m2)/absm2**2)
    
    print *,'alpha',alpha, 'and in deg',alpha/3.1515926*180
    
    alpha=alpha*fac1
    print *,'new alpha',alpha
    
    cosa=cos(alpha)
    sina=sin(alpha)
    n1=nvec(1)/absnvec
    n2=nvec(2)/absnvec
    n3=nvec(3)/absnvec
    
    rotmat(1,:) =  (/ n1*n1*(1-cosa) +    cosa, n1*n2*(1-cosa) - n3*sina, n1*n3*(1-cosa)+n2*sina /)
    rotmat(2,:) =  (/ n2*n1*(1-cosa) + n3*sina, n2*n2*(1-cosa) +    cosa, n2*n3*(1-cosa)-n1*sina /)
    rotmat(3,:) =  (/ n3*n1*(1-cosa) - n2*sina, n3*n2*(1-cosa) + n1*sina, n3*n3*(1-cosa)+cosa    /) 


    rotvector=matvec_dmdm(rotmat,m1)
    
    print *,m2
    print *,rotvector
  end function rotvector

end module mod_mathtools

