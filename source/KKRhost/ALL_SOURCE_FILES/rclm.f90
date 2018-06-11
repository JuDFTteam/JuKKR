    Subroutine rclm(key, ll, ldim, vmat)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * Transform complex matrix VMAT:                                     *
! * - From real to complex spherical harmonics basis (KEY=1)           *
! * - From complex to real spherical harmonics basis (KEY=2)           *
! *                                                                    *
! * Transformation matrix is A: Yreal = A * Ycmplx                     *
! *                                                                    *
! *       || I*i/sqrt(2)     0      J*(-i)/sqrt(2)  ||                 *
! *   A = ||      0          1            0         ||                 *
! *       || J*1/sqrt(2)     0      I*1/sqrt(2)     ||                 *
! *                                                                    *
! * (of dimension (2l+1)), where I is the l*l unit matrix,             *
! * J is the l*l antidiagonal unit matrix                              *
! *                                                                    *
! *     || 0 0 1 ||                                                    *
! * J = || 0 1 0 ||                                                    *
! *     || 1 0 0 ||                                                    *
! *                                                                    *
! * A_{mm'} = Int Y_{lm} Y^{c *}_{lm'} d\Omega                         *
! *                                                                    *
! * Transformation rule:                                               *
! * Complex --> Real (VC --> VR) :                                     *
! *                                                                    *
! * VR_{mm'} = Sum_{m1,m2} A_{m,m1} VC_{m1,m2} A^*_{m'm2}              *
! * with VC_{m1,m2} = Int Y^{c *}_{m2} V Y^{c}_{m1}                    *
! *                                                                    *
! * LDIM corresponds to the dimension of the array VMAT as             *
! * VMAT(2*LDIM+1,2*LDIM+1).                                           *
! * LL is the angular momentum l, corresponding to the part            *
! *    of VMAT used -- VMAT(2*LL+1,2*LL+1)                             *
! *                                                                    *
! * Result returns in again in VMAT (original VMAT is destroyed).      *
! *                                                                    *
! *                                ph. mavropoulos, juelich 2004       *
! **********************************************************************
      Implicit None
!..
      Complex (Kind=dp) :: cone, ci, czero
      Parameter (cone=(1E0_dp,0E0_dp), ci=(0E0_dp,1E0_dp), &
        czero=(0E0_dp,0E0_dp))
!..
      Integer :: ldim, key, ll
      Complex (Kind=dp) :: vmat(2*ldim+1, 2*ldim+1)
!.. Locals
      Complex (Kind=dp) :: vtmp(2*ldim+1, 2*ldim+1)
      Complex (Kind=dp) :: aa(2*ldim+1, 2*ldim+1), aac(2*ldim+1, 2*ldim+1)
      Real (Kind=dp) :: ovsqrtwo
      Complex (Kind=dp) :: oneovrt, ciovrt, cimovrt, blj
      Complex (Kind=dp) :: a11, a13, a31, a33
      Integer :: mdim, icall
      Integer :: m1, m2, m3, m4, mm, midl

      Data icall/0/
      Save :: icall, oneovrt, ciovrt, cimovrt
! ======================================================================
      mdim = 2*ldim + 1
      If (ll>ldim) Stop 'ERROR IN RCLM: LL>LDIM'
! ----------------------------------------------------------------------
      icall = icall + 1
      If (icall==1) Then
        ovsqrtwo = 1E0_dp/sqrt(2E0_dp)
        oneovrt = ovsqrtwo*cone
        ciovrt = ovsqrtwo*ci
        cimovrt = -ciovrt
      End If
! ----------------------------------------------------------------------

      If (key==2) Then ! matrix elements
        a11 = ciovrt
        a13 = cimovrt
        a31 = oneovrt
        a33 = oneovrt
      Else If (key==1) Then ! adjoint matrix elements
        a11 = conjg(ciovrt)
        a13 = oneovrt
        a31 = conjg(cimovrt)
        a33 = oneovrt
      Else
        Stop 'ERROR IN RCLM: KEY NOT EQUAL 1 OR 2'
      End If

! -> Construct transformation matrix AA

      Call cinit(mdim*mdim, aa)
      mm = 2*ll + 1
      midl = ll + 1
      Do m1 = 1, ll
        aa(m1, m1) = a11
        aa(m1, mm+1-m1) = a13
      End Do

      aa(midl, midl) = cone
      Do m1 = midl + 1, mm
        aa(m1, m1) = a33
        aa(m1, mm+1-m1) = a31
      End Do

! -> Construct transformation matrix AA+

      Do m2 = 1, mm
        Do m1 = 1, mm
          aac(m1, m2) = conjg(aa(m2,m1))
        End Do
      End Do

! -> Multiply from left

      Call cinit(mdim*mdim, vtmp)
      Do m2 = 1, mm
        Do m4 = 1, mm
          blj = aac(m4, m2)
          If (blj/=czero) Then
            Do m3 = 1, mm
              vtmp(m3, m2) = vtmp(m3, m2) + vmat(m3, m4)*blj
            End Do
          End If
        End Do
      End Do

      Call cinit(mdim*mdim, vmat)
      Do m2 = 1, mm
        Do m3 = 1, mm
          blj = vtmp(m3, m2)
          If (blj/=czero) Then
            Do m1 = 1, mm
              vmat(m1, m2) = vmat(m1, m2) + aa(m1, m3)*blj
            End Do
          End If
        End Do
      End Do
    End Subroutine
