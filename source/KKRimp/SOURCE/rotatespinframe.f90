!-------------------------------------------------------------------------------
!> Summary: Subroutines for rotation of vector and matrix between global and local spin frame
!> Author: 
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module mod_rotatespinframe

  character(len=*),parameter  :: spinmode='kkr'

  private
  public :: spinmode, rotatematrix, rotatevector

contains
  
  !-------------------------------------------------------------------------------
  !> Summary: Rotates a matrix in the local frame pointing in the direction of phi and theta to the global frame
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculate matrix in the global frame:
  !>
  !>  t_glob = U * t_loc * Udegga
  !-------------------------------------------------------------------------------
  subroutine rotatematrix(mat,theta,phi,lmmax,cmode)
    implicit none
    !interface
    double complex,intent(inout)    ::  mat(2*lmmax,2*lmmax)
    double precision,intent(in)     :: phi
    double precision,intent(in)     :: theta
    integer                         :: lmmax
    character(len=*)                :: cmode
    !local
    double complex   :: Umat(2*lmmax,2*lmmax)
    double complex   :: Udeggamat(2*lmmax,2*lmmax)
    double complex   :: mattemp(2*lmmax,2*lmmax)
    
    !***********************************************************************
    ! create the rotation matrix:
    !     | cos(theta/2) exp(-i/2 phi)   -sin(theta/2) exp(-i/2 phi) |
    !  U= |                                                          |
    !     | sin(theta/2) exp( i/2 phi)    cos(theta/2) exp( i/2 phi) |
    !
    !  Udegga = transpose(complex conjug ( U ) )
    !***********************************************************************
    call create_Umatrix(theta,phi,lmmax,Umat,Udeggamat)

    !***********************************************************************
    ! calculate matrix in the global frame:
    !
    !  t_glob = U * t_loc * Udegga
    !***********************************************************************
    if (cmode=='loc->glob') then
      mattemp = matmat_zmzm(mat,Udeggamat)
      mat  = matmat_zmzm(Umat,mattemp)
    elseif (cmode=='glob->loc') then
      mattemp = matmat_zmzm(mat,Umat)
      mat  = matmat_zmzm(Udeggamat,mattemp)
    else
      stop '[rotatematrix] mode not known'
    end if
  
  end subroutine rotatematrix
  
  
  !-------------------------------------------------------------------------------
  !> Summary: Rotate densities and charges from old to new local frame
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> does the rotation from the old local to the new local spin frame reference
  !> for densities and charges 
  !>  rho_loc(ir,lm)= W1 * rho_glob(ir,lm) * W2 
  !> where rho and W are matricies in spin space
  !-------------------------------------------------------------------------------
  subroutine rotatevector(rho2nsc,rho2ns,nrmaxd,lmpotd,theta,phi,theta_old,phi_old)
    implicit none
    !interface
    double complex   :: rho2nsc(nrmaxd,lmpotd,4)
    double precision :: rho2ns(nrmaxd,lmpotd,4)
    integer          :: nrmaxd,lmpotd
    double precision :: theta,phi
    double precision :: theta_old,phi_old
    !local
    integer          :: ir,ilm
    double complex   :: W1(2,2),W2(2,2)
    double complex   :: W1_11W2_11, W1_11W2_22, W1_11W2_12, W1_11W2_21
    double complex   :: W1_22W2_11, W1_22W2_22, W1_22W2_12, W1_22W2_21
    double complex   :: W1_12W2_11, W1_12W2_22, W1_12W2_12, W1_12W2_21
    double complex   :: W1_21W2_11, W1_21W2_22, W1_21W2_12, W1_21W2_21
    
    call create_Wmatrix(theta,phi,theta_old,phi_old,1,W1,W2)
    
    W1_11W2_11=W1(1,1)*W2(1,1)
    W1_11W2_22=W1(1,1)*W2(2,2)
    W1_11W2_12=W1(1,1)*W2(1,2)
    W1_11W2_21=W1(1,1)*W2(2,1)
    W1_22W2_11=W1(2,2)*W2(1,1)
    W1_22W2_22=W1(2,2)*W2(2,2)
    W1_22W2_12=W1(2,2)*W2(1,2)
    W1_22W2_21=W1(2,2)*W2(2,1)
    W1_12W2_11=W1(1,2)*W2(1,1)
    W1_12W2_22=W1(1,2)*W2(2,2)
    W1_12W2_12=W1(1,2)*W2(1,2)
    W1_12W2_21=W1(1,2)*W2(2,1)
    W1_21W2_11=W1(2,1)*W2(1,1)
    W1_21W2_22=W1(2,1)*W2(2,2)
    W1_21W2_12=W1(2,1)*W2(1,2)
    W1_21W2_21=W1(2,1)*W2(2,1)
    
    do ir=1,nrmaxd
      do ilm=1,lmpotd
        rho2ns(ir,ilm,1)= dimag(&
        +(rho2nsc(ir,ilm,1)*W1_11W2_11) &
        +(rho2nsc(ir,ilm,3)*W1_11W2_21) &
        +(rho2nsc(ir,ilm,4)*W1_12W2_11) &
        +(rho2nsc(ir,ilm,2)*W1_12W2_21)) 
    
        rho2ns(ir,ilm,2)= dimag (&
        +(rho2nsc(ir,ilm,1)*W1_21W2_12)  &
        -(rho2nsc(ir,ilm,3)*W1_21W2_22) &
        -(rho2nsc(ir,ilm,4)*W1_22W2_12) &
        +(rho2nsc(ir,ilm,2)*W1_22W2_22))
      end do
    enddo
  
  end subroutine

  
  !-------------------------------------------------------------------------------
  !> Summary: Create rotation matrices W, W^dagger for rotation between two local frames 
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !> 
  !> create the rotation matrix W:
  !>
  !>  W1= U_degga_locnew * U_locold
  !>  W2= U_degga_locold * U_locnew
  !>
  !>  the rotation matrix is created such that it rotates an operator
  !>  which is in a local frame (locold) to another local frame (locnew)
  !>  This is done by first transforming the old local frame to the
  !>  global frame using the U matrix and then transforming the global
  !>  frame to the new local frame
  !>
  !>  A_locnew = W1 * A_locold * W2
  !>
  !>  Udegga = transpose(complex conjug ( U ) )
  !-------------------------------------------------------------------------------
  subroutine create_Wmatrix(theta,phi,theta_old,phi_old,lmmax,Wmat1,Wmat2)
    implicit none
    !interface
    double precision,intent(in)     :: phi
    double precision,intent(in)     :: theta
    double precision,intent(in)     :: phi_old
    double precision,intent(in)     :: theta_old
    integer,intent(in)              :: lmmax
    double complex,intent(out)      :: Wmat1(2*lmmax,2*lmmax)
    double complex,intent(out)      :: Wmat2(2*lmmax,2*lmmax)
    !local
    double complex                  :: Umat1     (2*lmmax,2*lmmax)
    double complex                  :: Udeggamat1(2*lmmax,2*lmmax)
    double complex                  :: Umat2     (2*lmmax,2*lmmax)
    double complex                  :: Udeggamat2(2*lmmax,2*lmmax)
    
    call create_Umatrix(theta_old,phi_old,lmmax,Umat1,Udeggamat1)
    
    call create_Umatrix(theta,phi,lmmax,Umat2,Udeggamat2)
    
    Wmat1 = matmat_zmzm(Udeggamat2,Umat1)
    
    Wmat2 = matmat_zmzm(Udeggamat1,Umat2)
  
  end subroutine create_Wmatrix
  

  !-------------------------------------------------------------------------------
  !> Summary: Create rotation matrices U and U^dagger
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !> create the rotation matrix:
  !>     | cos(theta/2) exp(-i/2 phi)   -sin(theta/2) exp(-i/2 phi) |
  !>  U= |                                                          |
  !>     | sin(theta/2) exp( i/2 phi)    cos(theta/2) exp( i/2 phi) |
  !>
  !>  Udegga = transpose(complex conjug ( U ) )
  !-------------------------------------------------------------------------------
  subroutine create_Umatrix(theta,phi,lmmax,Umat,Udeggamat)
    implicit none
    !interface
    double precision,intent(in)     :: phi
    double precision,intent(in)     :: theta
    integer,intent(in)              :: lmmax
    double complex,intent(out)      :: Umat(2*lmmax,2*lmmax)
    double complex,intent(out)      :: Udeggamat(2*lmmax,2*lmmax)
    !local
    double complex                  :: Umat11,Umat12,Umat21,Umat22
    double complex                  :: Udeggamat11,Udeggamat12,Udeggamat21,Udeggamat22
    integer                         :: ival
    double complex,parameter        :: ci=(0.0D0,1.0D0)
    
    if (spinmode=='regular') then
      Umat11      =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
      Umat12      = -sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
      Umat21      =  sin(theta/2.0D0)*exp( ci/2.0D0*phi)
      Umat22      =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
    else if (spinmode=='kkr') then
      Umat11      =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
      Umat12      =  sin(theta/2.0D0)*exp( ci/2.0D0*phi)
      Umat21      = -sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
      Umat22      =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
    else 
      stop '[create_Umatrix] mode not known'
    end if
    
    
    Umat=(0.0D0,0.0D0)
    do ival=1,lmmax
      Umat(      ival,      ival) = Umat11
      Umat(      ival,lmmax+ival) = Umat12
      Umat(lmmax+ival,ival)       = Umat21
      Umat(lmmax+ival,lmmax+ival) = Umat22
    end do
    
    if (spinmode=='regular') then
      Udeggamat11 =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
      Udeggamat12 =  sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
      Udeggamat21 = -sin(theta/2.0D0)*exp( ci/2.0D0*phi)
      Udeggamat22 =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
    else if (spinmode=='kkr') then
      Udeggamat11 =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
      Udeggamat12 = -sin(theta/2.0D0)*exp( ci/2.0D0*phi)
      Udeggamat21 =  sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
      Udeggamat22 =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
    else 
      stop '[create_Umatrix] mode not known'
    end if
    

    Udeggamat=(0.0D0,0.0D0)
    do ival=1,lmmax
      Udeggamat(      ival,      ival) = Udeggamat11
      Udeggamat(      ival,lmmax+ival) = Udeggamat12
      Udeggamat(lmmax+ival,ival)       = Udeggamat21
      Udeggamat(lmmax+ival,lmmax+ival) = Udeggamat22
    end do
    
  end subroutine create_Umatrix
  
  !-------------------------------------------------------------------------------
  !> Summary: Wrapper for complex matrix-matrix multiplication using zgemm (LAPACK)
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: True ! This needs to be set to True for deprecated subroutines
  !> @note Multiple implementation of this exist... @endnote
  !-------------------------------------------------------------------------------
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


end module mod_rotatespinframe
