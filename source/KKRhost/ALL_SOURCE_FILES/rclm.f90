subroutine rclm(key, ll, ldim, vmat)
  use :: mod_datatypes, only: dp
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
  implicit none
  real (kind=dp), parameter :: eps=1.0D-12
  ! ..
  complex (kind=dp) :: cone, ci, czero
  parameter (cone=(1e0_dp,0e0_dp), ci=(0e0_dp,1e0_dp), czero=(0e0_dp,0e0_dp))
  ! ..
  integer :: ldim, key, ll
  complex (kind=dp) :: vmat(2*ldim+1, 2*ldim+1)
  ! .. Locals
  complex (kind=dp) :: vtmp(2*ldim+1, 2*ldim+1)
  complex (kind=dp) :: aa(2*ldim+1, 2*ldim+1), aac(2*ldim+1, 2*ldim+1)
  real (kind=dp) :: ovsqrtwo
  complex (kind=dp) :: oneovrt, ciovrt, cimovrt, blj
  complex (kind=dp) :: a11, a13, a31, a33
  integer :: mdim, icall
  integer :: m1, m2, m3, m4, mm, midl

  data icall/0/
  save :: icall, oneovrt, ciovrt, cimovrt
  ! ======================================================================
  mdim = 2*ldim + 1
  if (ll>ldim) stop 'ERROR IN RCLM: LL>LDIM'
  ! ----------------------------------------------------------------------
  icall = icall + 1
  if (icall==1) then
    ovsqrtwo = 1e0_dp/sqrt(2e0_dp)
    oneovrt = ovsqrtwo*cone
    ciovrt = ovsqrtwo*ci
    cimovrt = -ciovrt
  end if
  ! ----------------------------------------------------------------------

  if (key==2) then                 ! matrix elements
    a11 = ciovrt
    a13 = cimovrt
    a31 = oneovrt
    a33 = oneovrt
  else if (key==1) then            ! adjoint matrix elements
    a11 = conjg(ciovrt)
    a13 = oneovrt
    a31 = conjg(cimovrt)
    a33 = oneovrt
  else
    stop 'ERROR IN RCLM: KEY NOT EQUAL 1 OR 2'
  end if

  ! -> Construct transformation matrix AA

  call cinit(mdim*mdim, aa)
  mm = 2*ll + 1
  midl = ll + 1
  do m1 = 1, ll
    aa(m1, m1) = a11
    aa(m1, mm+1-m1) = a13
  end do

  aa(midl, midl) = cone
  do m1 = midl + 1, mm
    aa(m1, m1) = a33
    aa(m1, mm+1-m1) = a31
  end do

  ! -> Construct transformation matrix AA+

  do m2 = 1, mm
    do m1 = 1, mm
      aac(m1, m2) = conjg(aa(m2,m1))
    end do
  end do

  ! -> Multiply from left

  call cinit(mdim*mdim, vtmp)
  do m2 = 1, mm
    do m4 = 1, mm
      blj = aac(m4, m2)
      if (abs(blj)>eps) then
        do m3 = 1, mm
          vtmp(m3, m2) = vtmp(m3, m2) + vmat(m3, m4)*blj
        end do
      end if
    end do
  end do

  call cinit(mdim*mdim, vmat)
  do m2 = 1, mm
    do m3 = 1, mm
      blj = vtmp(m3, m2)
      if (abs(blj)>eps) then
        do m1 = 1, mm
          vmat(m1, m2) = vmat(m1, m2) + aa(m1, m3)*blj
        end do
      end if
    end do
  end do
end subroutine rclm
