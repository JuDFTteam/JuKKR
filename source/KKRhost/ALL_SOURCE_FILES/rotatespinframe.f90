!-------------------------------------------------------------------------------
! SUBROUTINE: rotatematrix
!> @brief Rotates a matrix in the local frame pointing in
!> the direction of phi and theta to the global frame
!-------------------------------------------------------------------------------
subroutine rotatematrix(mat,theta,phi,lmmax,mode)
   implicit none
   !interface
   integer, intent(in)                                         :: lmmax
   integer, intent(in)                                         :: mode
   double precision,intent(in)                                 :: phi
   double precision,intent(in)                                 :: theta
   double complex, dimension(2*lmmax,2*lmmax), intent(inout)   :: mat
   !local
   double complex, dimension(2*lmmax,2*lmmax)   :: Umat
   double complex, dimension(2*lmmax,2*lmmax)   :: Udeggamat
   double complex, dimension(2*lmmax,2*lmmax)   :: mattemp

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


   if (mode==0) then ! 'loc->glob'
      call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),mat,2*lmmax,Udeggamat,2*lmmax,(0d0,0d0),mattemp,2*lmmax)
      call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),Umat,2*lmmax,mattemp,2*lmmax,(0d0,0d0),mat,2*lmmax)
   elseif (mode==1) then !'glob->loc'
      call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),mat,2*lmmax,Umat,2*lmmax,(0d0,0d0),mattemp,2*lmmax)
      call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),Udeggamat,2*lmmax,mattemp,2*lmmax,(0d0,0d0),mat,2*lmmax)
   else
      stop '[rotatematrix] mode not known'
   end if

end subroutine rotatematrix

!-------------------------------------------------------------------------------
! SUBROUTINE: rotatevector
!> @brief Does the rotation from the old local to the new local spin frame reference
!> for densities and charges
!> \f$\rho_{loc}(ir,lm)= W1 * \rho_{glob}(ir,lm) * W2\f$
!> where \f$\rho\f$ and \f$W\f$ are matricies in spin space
!-------------------------------------------------------------------------------
subroutine rotatevector(rho2nsc,rho2ns,nrmax,lmpotd,theta,phi,theta_old,phi_old,nrmaxd)

   implicit none
   !interface
   double complex   :: rho2nsc(nrmaxd,lmpotd,4)
   double precision :: rho2ns(nrmaxd,lmpotd,4)
   integer          :: nrmaxd,lmpotd,nrmax
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

   do ir=1,nrmax
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
! SUBROUTINE: create_Wmatrix
!> @brief Create the rotation matrix \f$W\f$:
!>
!>  \f$W1= U_{degga_{locnew}} * U_{locold}\f$
!>
!>  \f$W2= U_{degga_{locold}} * U_{locnew}\f$
!!>
!> @detail The rotation matrix is created such that it rotates an operator
!>  which is in a local frame (locold) to another local frame (locnew)
!>  This is done by first transforming the old local frame to the
!>  global frame using the U matrix and then transforming the global
!>  frame to the new local frame
!>
!>  \f$A_{locnew} = W1 * A_{locold} * W2\f$
!>
!>  \f$Udegga = transpose(complex conjug ( U ) )\f&
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

   call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),Udeggamat2,2*lmmax,Umat1,2*lmmax,(0d0,0d0),Wmat1,2*lmmax)
   call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),Udeggamat1,2*lmmax,Umat2,2*lmmax,(0d0,0d0),Wmat2,2*lmmax)


end subroutine create_Wmatrix

!-------------------------------------------------------------------------------
!> @brief create the rotation matrix:
!>  \f$U= \left(\begin{array}{cc} cos(\theta/2) exp(-i/2 \phi) & -sin(\theta/2) exp(-i/2 \phi) \\ sin(\theta/2) exp( i/2 \phi)  &  cos(\theta/2) exp( i/2 \phi)\end{array}\right)\f$
!>  \f$Udegga = transpose(complex conjug ( U ) )\f$
!-------------------------------------------------------------------------------
subroutine create_Umatrix(theta,phi,lmmax,Umat,Udeggamat)

   use Constants

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
   character*25               :: spinmode

   spinmode='kkr'
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

   Umat=CZERO
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

   Udeggamat=CZERO
   do ival=1,lmmax
      Udeggamat(      ival,      ival) = Udeggamat11
      Udeggamat(      ival,lmmax+ival) = Udeggamat12
      Udeggamat(lmmax+ival,ival)       = Udeggamat21
      Udeggamat(lmmax+ival,lmmax+ival) = Udeggamat22
   end do

end subroutine create_Umatrix
