!-------------------------------------------------------------------------------
! SUBROUTINE: rotatematrix
!> @brief Rotates a matrix in the local frame pointing in
!> the direction of phi and theta to the global frame
!-------------------------------------------------------------------------------
    Subroutine rotatematrix(mat, theta, phi, lmmax, mode)
      Use mod_datatypes, Only: dp
      Implicit None
!interface
      Integer, Intent (In) :: lmmax
      Integer, Intent (In) :: mode
      Real (Kind=dp), Intent (In) :: phi
      Real (Kind=dp), Intent (In) :: theta
      Complex (Kind=dp), Dimension (2*lmmax, 2*lmmax), Intent (Inout) :: mat
!local
      Complex (Kind=dp), Dimension (2*lmmax, 2*lmmax) :: umat
      Complex (Kind=dp), Dimension (2*lmmax, 2*lmmax) :: udeggamat
      Complex (Kind=dp), Dimension (2*lmmax, 2*lmmax) :: mattemp

!***********************************************************************
! create the rotation matrix:
!     | cos(theta/2) exp(-i/2 phi)   -sin(theta/2) exp(-i/2 phi) |
!  U= |                                                          |
!     | sin(theta/2) exp( i/2 phi)    cos(theta/2) exp( i/2 phi) |
!
!  Udegga = transpose(complex conjug ( U ) )
!***********************************************************************


      Call create_umatrix(theta, phi, lmmax, umat, udeggamat)
!***********************************************************************
! calculate matrix in the global frame:
!
!  t_glob = U * t_loc * Udegga
!***********************************************************************


      If (mode==0) Then ! 'loc->glob'
        Call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1E0_dp,0E0_dp), mat, &
          2*lmmax, udeggamat, 2*lmmax, (0E0_dp,0E0_dp), mattemp, 2*lmmax)
        Call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1E0_dp,0E0_dp), umat, &
          2*lmmax, mattemp, 2*lmmax, (0E0_dp,0E0_dp), mat, 2*lmmax)
      Else If (mode==1) Then !'glob->loc'
        Call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1E0_dp,0E0_dp), mat, &
          2*lmmax, umat, 2*lmmax, (0E0_dp,0E0_dp), mattemp, 2*lmmax)
        Call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1E0_dp,0E0_dp), &
          udeggamat, 2*lmmax, mattemp, 2*lmmax, (0E0_dp,0E0_dp), mat, 2*lmmax)
      Else
        Stop '[rotatematrix] mode not known'
      End If

    End Subroutine

!-------------------------------------------------------------------------------
! SUBROUTINE: rotatevector
!> @brief Does the rotation from the old local to the new local spin frame reference
!> for densities and charges
!> \f$\rho_{loc}(ir,lm)= W1 * \rho_{glob}(ir,lm) * W2\f$
!> where \f$\rho\f$ and \f$W\f$ are matricies in spin space
!-------------------------------------------------------------------------------
    Subroutine rotatevector(rho2nsc, rho2ns, nrmax, lmpotd, theta, phi, &
      theta_old, phi_old, nrmaxd)
      Use mod_datatypes, Only: dp

      Implicit None
!interface
      Integer :: nrmaxd, lmpotd, nrmax
      Complex (Kind=dp) :: rho2nsc(nrmaxd, lmpotd, 4)
      Real (Kind=dp) :: rho2ns(nrmaxd, lmpotd, 4)
      Real (Kind=dp) :: theta, phi
      Real (Kind=dp) :: theta_old, phi_old
!local
      Integer :: ir, ilm
      Complex (Kind=dp) :: w1(2, 2), w2(2, 2)
      Complex (Kind=dp) :: w1_11w2_11, w1_11w2_22, w1_11w2_12, w1_11w2_21
      Complex (Kind=dp) :: w1_22w2_11, w1_22w2_22, w1_22w2_12, w1_22w2_21
      Complex (Kind=dp) :: w1_12w2_11, w1_12w2_22, w1_12w2_12, w1_12w2_21
      Complex (Kind=dp) :: w1_21w2_11, w1_21w2_22, w1_21w2_12, w1_21w2_21

      Call create_wmatrix(theta, phi, theta_old, phi_old, 1, w1, w2)

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

      Do ir = 1, nrmax
        Do ilm = 1, lmpotd
          rho2ns(ir, ilm, 1) = aimag(+(rho2nsc(ir,ilm, &
            1)*w1_11w2_11)+(rho2nsc(ir,ilm,3)*w1_11w2_21)+(rho2nsc(ir,ilm, &
            4)*w1_12w2_11)+(rho2nsc(ir,ilm,2)*w1_12w2_21))

          rho2ns(ir, ilm, 2) = aimag(+(rho2nsc(ir,ilm, &
            1)*w1_21w2_12)-(rho2nsc(ir,ilm,3)*w1_21w2_22)-(rho2nsc(ir,ilm, &
            4)*w1_22w2_12)+(rho2nsc(ir,ilm,2)*w1_22w2_22))
        End Do
      End Do

    End Subroutine

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
    Subroutine create_wmatrix(theta, phi, theta_old, phi_old, lmmax, wmat1, &
      wmat2)
      Use mod_datatypes, Only: dp
      Implicit None
!interface
      Real (Kind=dp), Intent (In) :: phi
      Real (Kind=dp), Intent (In) :: theta
      Real (Kind=dp), Intent (In) :: phi_old
      Real (Kind=dp), Intent (In) :: theta_old
      Integer, Intent (In) :: lmmax
      Complex (Kind=dp), Intent (Out) :: wmat1(2*lmmax, 2*lmmax)
      Complex (Kind=dp), Intent (Out) :: wmat2(2*lmmax, 2*lmmax)
!local
      Complex (Kind=dp) :: umat1(2*lmmax, 2*lmmax)
      Complex (Kind=dp) :: udeggamat1(2*lmmax, 2*lmmax)
      Complex (Kind=dp) :: umat2(2*lmmax, 2*lmmax)
      Complex (Kind=dp) :: udeggamat2(2*lmmax, 2*lmmax)

      Call create_umatrix(theta_old, phi_old, lmmax, umat1, udeggamat1)

      Call create_umatrix(theta, phi, lmmax, umat2, udeggamat2)

      Call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1E0_dp,0E0_dp), &
        udeggamat2, 2*lmmax, umat1, 2*lmmax, (0E0_dp,0E0_dp), wmat1, 2*lmmax)
      Call zgemm('N', 'N', 2*lmmax, 2*lmmax, 2*lmmax, (1E0_dp,0E0_dp), &
        udeggamat1, 2*lmmax, umat2, 2*lmmax, (0E0_dp,0E0_dp), wmat2, 2*lmmax)


    End Subroutine

!-------------------------------------------------------------------------------
!> @brief create the rotation matrix:
!>  \f$U= \left(\begin{array}{cc} cos(\theta/2) exp(-i/2 \phi) & -sin(\theta/2) exp(-i/2 \phi) \\ sin(\theta/2) exp( i/2 \phi)  &  cos(\theta/2) exp( i/2 \phi)\end{array}\right)\f$
!>  \f$Udegga = transpose(complex conjug ( U ) )\f$
!-------------------------------------------------------------------------------
    Subroutine create_umatrix(theta, phi, lmmax, umat, udeggamat)

      Use constants
      Use mod_datatypes, Only: dp

      Implicit None
!interface
      Real (Kind=dp), Intent (In) :: phi
      Real (Kind=dp), Intent (In) :: theta
      Integer, Intent (In) :: lmmax
      Complex (Kind=dp), Intent (Out) :: umat(2*lmmax, 2*lmmax)
      Complex (Kind=dp), Intent (Out) :: udeggamat(2*lmmax, 2*lmmax)
!local
      Complex (Kind=dp) :: umat11, umat12, umat21, umat22
      Complex (Kind=dp) :: udeggamat11, udeggamat12, udeggamat21, udeggamat22
      Integer :: ival
      Character (Len=25) :: spinmode

      spinmode = 'kkr'
      If (spinmode=='regular') Then
        umat11 = cos(theta/2.0E0_dp)*exp(-ci/2.0E0_dp*phi)
        umat12 = -sin(theta/2.0E0_dp)*exp(-ci/2.0E0_dp*phi)
        umat21 = sin(theta/2.0E0_dp)*exp(ci/2.0E0_dp*phi)
        umat22 = cos(theta/2.0E0_dp)*exp(ci/2.0E0_dp*phi)
      Else If (spinmode=='kkr') Then
        umat11 = cos(theta/2.0E0_dp)*exp(ci/2.0E0_dp*phi)
        umat12 = sin(theta/2.0E0_dp)*exp(ci/2.0E0_dp*phi)
        umat21 = -sin(theta/2.0E0_dp)*exp(-ci/2.0E0_dp*phi)
        umat22 = cos(theta/2.0E0_dp)*exp(-ci/2.0E0_dp*phi)
      Else
        Stop '[create_Umatrix] mode not known'
      End If

      umat = czero
      Do ival = 1, lmmax
        umat(ival, ival) = umat11
        umat(ival, lmmax+ival) = umat12
        umat(lmmax+ival, ival) = umat21
        umat(lmmax+ival, lmmax+ival) = umat22
      End Do

      If (spinmode=='regular') Then
        udeggamat11 = cos(theta/2.0E0_dp)*exp(ci/2.0E0_dp*phi)
        udeggamat12 = sin(theta/2.0E0_dp)*exp(-ci/2.0E0_dp*phi)
        udeggamat21 = -sin(theta/2.0E0_dp)*exp(ci/2.0E0_dp*phi)
        udeggamat22 = cos(theta/2.0E0_dp)*exp(-ci/2.0E0_dp*phi)
      Else If (spinmode=='kkr') Then
        udeggamat11 = cos(theta/2.0E0_dp)*exp(-ci/2.0E0_dp*phi)
        udeggamat12 = -sin(theta/2.0E0_dp)*exp(ci/2.0E0_dp*phi)
        udeggamat21 = sin(theta/2.0E0_dp)*exp(-ci/2.0E0_dp*phi)
        udeggamat22 = cos(theta/2.0E0_dp)*exp(ci/2.0E0_dp*phi)
      Else
        Stop '[create_Umatrix] mode not known'
      End If

      udeggamat = czero
      Do ival = 1, lmmax
        udeggamat(ival, ival) = udeggamat11
        udeggamat(ival, lmmax+ival) = udeggamat12
        udeggamat(lmmax+ival, ival) = udeggamat21
        udeggamat(lmmax+ival, lmmax+ival) = udeggamat22
      End Do

    End Subroutine
