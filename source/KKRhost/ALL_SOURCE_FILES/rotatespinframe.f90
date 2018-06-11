!-------------------------------------------------------------------------------
! SUBROUTINE: rotatematrix
!> @brief Rotates a matrix in the local frame pointing in
!> the direction of phi and theta to the global frame
!-------------------------------------------------------------------------------
subroutine rotatematrix(mat, theta, phi, lmmax, mode)
  implicit none
!interface
  integer, intent (in) :: lmmax
  integer, intent (in) :: mode
  double precision, intent (in) :: phi
  double precision, intent (in) :: theta
  double complex, dimension (2*lmmax, 2*lmmax), intent (inout) :: mat
!local
  double complex, dimension (2*lmmax, 2*lmmax) :: umat
  double complex, dimension (2*lmmax, 2*lmmax) :: udeggamat
  double complex, dimension (2*lmmax, 2*lmmax) :: mattemp

!***********************************************************************
! create the rotation matrix:
!     | cos(theta/2) exp(-i/2 phi)   -sin(theta/2) exp(-i/2 phi) |
!  U= |                                                          |
!     | sin(theta/2) exp( i/2 phi)    cos(theta/2) exp( i/2 phi) |
!
!  Udegga = transpose(complex conjug ( U ) )
!***********************************************************************


  call create_umatrix(theta, phi, lmmax, umat, udeggamat)
!***********************************************************************
! calculate matrix in the global frame:
!
!  t_glob = U * t_loc * Udegga
!***********************************************************************


  if (mode==0) then ! 'loc->glob'
    call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1d0,0d0), mat, 2*lmmax, &
      udeggamat, 2*lmmax, (0d0,0d0), mattemp, 2*lmmax)
    call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1d0,0d0), umat, 2*lmmax, &
      mattemp, 2*lmmax, (0d0,0d0), mat, 2*lmmax)
  else if (mode==1) then !'glob->loc'
    call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1d0,0d0), mat, 2*lmmax, &
      umat, 2*lmmax, (0d0,0d0), mattemp, 2*lmmax)
    call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1d0,0d0), udeggamat, &
      2*lmmax, mattemp, 2*lmmax, (0d0,0d0), mat, 2*lmmax)
  else
    stop '[rotatematrix] mode not known'
  end if

end subroutine

!-------------------------------------------------------------------------------
! SUBROUTINE: rotatevector
!> @brief Does the rotation from the old local to the new local spin frame reference
!> for densities and charges
!> \f$\rho_{loc}(ir,lm)= W1 * \rho_{glob}(ir,lm) * W2\f$
!> where \f$\rho\f$ and \f$W\f$ are matricies in spin space
!-------------------------------------------------------------------------------
subroutine rotatevector(rho2nsc, rho2ns, nrmax, lmpotd, theta, phi, theta_old, &
  phi_old, nrmaxd)

  implicit none
!interface
  integer :: nrmaxd, lmpotd, nrmax
  double complex :: rho2nsc(nrmaxd, lmpotd, 4)
  double precision :: rho2ns(nrmaxd, lmpotd, 4)
  double precision :: theta, phi
  double precision :: theta_old, phi_old
!local
  integer :: ir, ilm
  double complex :: w1(2, 2), w2(2, 2)
  double complex :: w1_11w2_11, w1_11w2_22, w1_11w2_12, w1_11w2_21
  double complex :: w1_22w2_11, w1_22w2_22, w1_22w2_12, w1_22w2_21
  double complex :: w1_12w2_11, w1_12w2_22, w1_12w2_12, w1_12w2_21
  double complex :: w1_21w2_11, w1_21w2_22, w1_21w2_12, w1_21w2_21

  call create_wmatrix(theta, phi, theta_old, phi_old, 1, w1, w2)

  w1_11w2_11 = w1(1, 1)*w2(1, 1)
  w1_11w2_22 = w1(1, 1)*w2(2, 2)
  w1_11w2_12 = w1(1, 1)*w2(1, 2)
  w1_11w2_21 = w1(1, 1)*w2(2, 1)
  w1_22w2_11 = w1(2, 2)*w2(1, 1)
  w1_22w2_22 = w1(2, 2)*w2(2, 2)
  w1_22w2_12 = w1(2, 2)*w2(1, 2)
  w1_22w2_21 = w1(2, 2)*w2(2, 1)
  w1_12w2_11 = w1(1, 2)*w2(1, 1)
  w1_12w2_22 = w1(1, 2)*w2(2, 2)
  w1_12w2_12 = w1(1, 2)*w2(1, 2)
  w1_12w2_21 = w1(1, 2)*w2(2, 1)
  w1_21w2_11 = w1(2, 1)*w2(1, 1)
  w1_21w2_22 = w1(2, 1)*w2(2, 2)
  w1_21w2_12 = w1(2, 1)*w2(1, 2)
  w1_21w2_21 = w1(2, 1)*w2(2, 1)

  do ir = 1, nrmax
    do ilm = 1, lmpotd
      rho2ns(ir, ilm, 1) = aimag(+(rho2nsc(ir,ilm,1)*w1_11w2_11)+(rho2nsc(ir, &
        ilm,3)*w1_11w2_21)+(rho2nsc(ir,ilm,4)*w1_12w2_11)+(rho2nsc(ir,ilm, &
        2)*w1_12w2_21))

      rho2ns(ir, ilm, 2) = aimag(+(rho2nsc(ir,ilm,1)*w1_21w2_12)-(rho2nsc(ir, &
        ilm,3)*w1_21w2_22)-(rho2nsc(ir,ilm,4)*w1_22w2_12)+(rho2nsc(ir,ilm, &
        2)*w1_22w2_22))
    end do
  end do

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
subroutine create_wmatrix(theta, phi, theta_old, phi_old, lmmax, wmat1, wmat2)
  implicit none
!interface
  double precision, intent (in) :: phi
  double precision, intent (in) :: theta
  double precision, intent (in) :: phi_old
  double precision, intent (in) :: theta_old
  integer, intent (in) :: lmmax
  double complex, intent (out) :: wmat1(2*lmmax, 2*lmmax)
  double complex, intent (out) :: wmat2(2*lmmax, 2*lmmax)
!local
  double complex :: umat1(2*lmmax, 2*lmmax)
  double complex :: udeggamat1(2*lmmax, 2*lmmax)
  double complex :: umat2(2*lmmax, 2*lmmax)
  double complex :: udeggamat2(2*lmmax, 2*lmmax)

  call create_umatrix(theta_old, phi_old, lmmax, umat1, udeggamat1)

  call create_umatrix(theta, phi, lmmax, umat2, udeggamat2)

  call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1d0,0d0), udeggamat2, &
    2*lmmax, umat1, 2*lmmax, (0d0,0d0), wmat1, 2*lmmax)
  call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1d0,0d0), udeggamat1, &
    2*lmmax, umat2, 2*lmmax, (0d0,0d0), wmat2, 2*lmmax)


end subroutine

!-------------------------------------------------------------------------------
!> @brief create the rotation matrix:
!>  \f$U= \left(\begin{array}{cc} cos(\theta/2) exp(-i/2 \phi) & -sin(\theta/2) exp(-i/2 \phi) \\ sin(\theta/2) exp( i/2 \phi)  &  cos(\theta/2) exp( i/2 \phi)\end{array}\right)\f$
!>  \f$Udegga = transpose(complex conjug ( U ) )\f$
!-------------------------------------------------------------------------------
subroutine create_umatrix(theta, phi, lmmax, umat, udeggamat)

  use :: constants

  implicit none
!interface
  double precision, intent (in) :: phi
  double precision, intent (in) :: theta
  integer, intent (in) :: lmmax
  double complex, intent (out) :: umat(2*lmmax, 2*lmmax)
  double complex, intent (out) :: udeggamat(2*lmmax, 2*lmmax)
!local
  double complex :: umat11, umat12, umat21, umat22
  double complex :: udeggamat11, udeggamat12, udeggamat21, udeggamat22
  integer :: ival
  character (len=25) :: spinmode

  spinmode = 'kkr'
  if (spinmode=='regular') then
    umat11 = cos(theta/2.0d0)*exp(-ci/2.0d0*phi)
    umat12 = -sin(theta/2.0d0)*exp(-ci/2.0d0*phi)
    umat21 = sin(theta/2.0d0)*exp(ci/2.0d0*phi)
    umat22 = cos(theta/2.0d0)*exp(ci/2.0d0*phi)
  else if (spinmode=='kkr') then
    umat11 = cos(theta/2.0d0)*exp(ci/2.0d0*phi)
    umat12 = sin(theta/2.0d0)*exp(ci/2.0d0*phi)
    umat21 = -sin(theta/2.0d0)*exp(-ci/2.0d0*phi)
    umat22 = cos(theta/2.0d0)*exp(-ci/2.0d0*phi)
  else
    stop '[create_Umatrix] mode not known'
  end if

  umat = czero
  do ival = 1, lmmax
    umat(ival, ival) = umat11
    umat(ival, lmmax+ival) = umat12
    umat(lmmax+ival, ival) = umat21
    umat(lmmax+ival, lmmax+ival) = umat22
  end do

  if (spinmode=='regular') then
    udeggamat11 = cos(theta/2.0d0)*exp(ci/2.0d0*phi)
    udeggamat12 = sin(theta/2.0d0)*exp(-ci/2.0d0*phi)
    udeggamat21 = -sin(theta/2.0d0)*exp(ci/2.0d0*phi)
    udeggamat22 = cos(theta/2.0d0)*exp(-ci/2.0d0*phi)
  else if (spinmode=='kkr') then
    udeggamat11 = cos(theta/2.0d0)*exp(-ci/2.0d0*phi)
    udeggamat12 = -sin(theta/2.0d0)*exp(ci/2.0d0*phi)
    udeggamat21 = sin(theta/2.0d0)*exp(-ci/2.0d0*phi)
    udeggamat22 = cos(theta/2.0d0)*exp(ci/2.0d0*phi)
  else
    stop '[create_Umatrix] mode not known'
  end if

  udeggamat = czero
  do ival = 1, lmmax
    udeggamat(ival, ival) = udeggamat11
    udeggamat(ival, lmmax+ival) = udeggamat12
    udeggamat(lmmax+ival, ival) = udeggamat21
    udeggamat(lmmax+ival, lmmax+ival) = udeggamat22
  end do

end subroutine
