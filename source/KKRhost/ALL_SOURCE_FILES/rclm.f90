SUBROUTINE rclm(key,ll,ldim,vmat)
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
IMPLICIT NONE
!..
DOUBLE COMPLEX CONE,CI,CZERO
PARAMETER ( CONE=(1D0,0D0), CI=(0D0,1D0), CZERO=(0D0,0D0) )
!..
INTEGER LDIM,KEY,LL
DOUBLE COMPLEX VMAT(2*LDIM+1,2*LDIM+1)
!.. Locals
DOUBLE COMPLEX VTMP(2*LDIM+1,2*LDIM+1)
DOUBLE COMPLEX AA(2*LDIM+1,2*LDIM+1),AAC(2*LDIM+1,2*LDIM+1)
DOUBLE PRECISION OVSQRTWO
DOUBLE COMPLEX ONEOVRT,CIOVRT,CIMOVRT,BLJ
DOUBLE COMPLEX A11,A13,A31,A33
INTEGER MDIM,ICALL
INTEGER M1,M2,M3,M4,MM,MIDL

DATA ICALL /0/
SAVE ICALL, ONEOVRT,CIOVRT,CIMOVRT
! ======================================================================
mdim = 2*ldim+1
IF (ll > ldim) STOP 'ERROR IN RCLM: LL>LDIM'
! ----------------------------------------------------------------------
icall = icall + 1
IF ( icall == 1 ) THEN
  ovsqrtwo = 1D0/DSQRT(2D0)
  oneovrt = ovsqrtwo * cone
  ciovrt =  ovsqrtwo * ci
  cimovrt = -ciovrt
END IF
! ----------------------------------------------------------------------

IF (key == 2) THEN        ! matrix elements
  a11 = ciovrt
  a13 = cimovrt
  a31 = oneovrt
  a33 = oneovrt
ELSE IF (key == 1) THEN   ! adjoint matrix elements
  a11 = DCONJG(ciovrt)
  a13 = oneovrt
  a31 = DCONJG(cimovrt)
  a33 = oneovrt
ELSE
  STOP 'ERROR IN RCLM: KEY NOT EQUAL 1 OR 2'
END IF

! -> Construct transformation matrix AA

CALL cinit(mdim*mdim,aa)
mm = 2*ll + 1
midl = ll+1
DO m1 = 1,ll
  aa(m1,m1) = a11
  aa(m1,mm+1-m1) = a13
END DO

aa(midl,midl) = cone
DO m1 = midl+1,mm
  aa(m1,m1) = a33
  aa(m1,mm+1-m1) = a31
END DO

! -> Construct transformation matrix AA+

DO m2 = 1,mm
  DO m1 = 1,mm
    aac(m1,m2) = DCONJG(aa(m2,m1))
  END DO
END DO

! -> Multiply from left

CALL cinit(mdim*mdim,vtmp)
DO m2 = 1,mm
  DO m4 = 1,mm
    blj = aac(m4,m2)
    IF ( blj /= czero ) THEN
      DO m3 = 1,mm
        vtmp(m3,m2) = vtmp(m3,m2) + vmat(m3,m4)*blj
      END DO
    END IF
  END DO
END DO

CALL cinit(mdim*mdim,vmat)
DO m2 = 1,mm
  DO m3 = 1,mm
    blj = vtmp(m3,m2)
    IF ( blj /= czero ) THEN
      DO m1 = 1,mm
        vmat(m1,m2) = vmat(m1,m2) + aa(m1,m3)*blj
      END DO
    END IF
  END DO
END DO
END SUBROUTINE rclm

