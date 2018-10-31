      SUBROUTINE RCLM(KEY,LL,LDIM,VMAT)
C **********************************************************************
C *                                                                    *
C * Transform complex matrix VMAT:                                     *
C * - From real to complex spherical harmonics basis (KEY=1)           *
C * - From complex to real spherical harmonics basis (KEY=2)           *
C *                                                                    *
C * Transformation matrix is A: Yreal = A * Ycmplx                     *
C *                                                                    *
C *       || I*i/sqrt(2)     0      J*(-i)/sqrt(2)  ||                 *
C *   A = ||      0          1            0         ||                 *
C *       || J*1/sqrt(2)     0      I*1/sqrt(2)     ||                 *
C *                                                                    *
C * (of dimension (2l+1)), where I is the l*l unit matrix,             *
C * J is the l*l antidiagonal unit matrix                              *
C *                                                                    *
C *     || 0 0 1 ||                                                    *
C * J = || 0 1 0 ||                                                    *
C *     || 1 0 0 ||                                                    *
C *                                                                    *
C * A_{mm'} = Int Y_{lm} Y^{c *}_{lm'} d\Omega                         *
C *                                                                    *
C * Transformation rule:                                               *
C * Complex --> Real (VC --> VR) :                                     *
C *                                                                    *
C * VR_{mm'} = Sum_{m1,m2} A_{m,m1} VC_{m1,m2} A^*_{m'm2}              *
C * with VC_{m1,m2} = Int Y^{c *}_{m2} V Y^{c}_{m1}                    *
C *                                                                    *
C * LDIM corresponds to the dimension of the array VMAT as             *
C * VMAT(2*LDIM+1,2*LDIM+1).                                           *
C * LL is the angular momentum l, corresponding to the part            *
C *    of VMAT used -- VMAT(2*LL+1,2*LL+1)                             *
C *                                                                    *
C * Result returns in again in VMAT (original VMAT is destroyed).      *
C *                                                                    *
C *                                ph. mavropoulos, juelich 2004       *
C **********************************************************************
      IMPLICIT NONE
C     ..
      DOUBLE COMPLEX CONE,CI,CZERO
      PARAMETER ( CONE=(1D0,0D0), CI=(0D0,1D0), CZERO=(0D0,0D0) )
C     ..
      INTEGER LDIM,KEY,LL
      DOUBLE COMPLEX VMAT(2*LDIM+1,2*LDIM+1)
C     .. Locals
      DOUBLE COMPLEX VTMP(2*LDIM+1,2*LDIM+1)
      DOUBLE COMPLEX AA(2*LDIM+1,2*LDIM+1),AAC(2*LDIM+1,2*LDIM+1)
      DOUBLE PRECISION OVSQRTWO
      DOUBLE COMPLEX ONEOVRT,CIOVRT,CIMOVRT,BLJ
      DOUBLE COMPLEX A11,A13,A31,A33
      INTEGER MDIM,ICALL
      INTEGER M1,M2,M3,M4,MM,MIDL
C
      DATA ICALL /0/
      SAVE ICALL, ONEOVRT,CIOVRT,CIMOVRT
C ======================================================================
      MDIM = 2*LDIM+1
      IF (LL.GT.LDIM) STOP 'ERROR IN RCLM: LL>LDIM'
C ----------------------------------------------------------------------
      ICALL = ICALL + 1
      IF ( ICALL.EQ.1 ) THEN
         OVSQRTWO = 1D0/DSQRT(2D0)
         ONEOVRT = OVSQRTWO * CONE
         CIOVRT =  OVSQRTWO * CI
         CIMOVRT = -CIOVRT
      END IF
C ----------------------------------------------------------------------
C
      IF (KEY.EQ.2) THEN        ! matrix elements
         A11 = CIOVRT
         A13 = CIMOVRT
         A31 = ONEOVRT
         A33 = ONEOVRT
      ELSE IF (KEY.EQ.1) THEN   ! adjoint matrix elements
         A11 = DCONJG(CIOVRT)
         A13 = ONEOVRT
         A31 = DCONJG(CIMOVRT)
         A33 = ONEOVRT
      ELSE
         STOP 'ERROR IN RCLM: KEY NOT EQUAL 1 OR 2'
      END IF
C
C -> Construct transformation matrix AA
C
      CALL CINIT(MDIM*MDIM,AA)
      MM = 2*LL + 1
      MIDL = LL+1
      DO M1 = 1,LL
         AA(M1,M1) = A11
         AA(M1,MM+1-M1) = A13
      ENDDO
C
      AA(MIDL,MIDL) = CONE
      DO M1 = MIDL+1,MM
         AA(M1,M1) = A33
         AA(M1,MM+1-M1) = A31
      ENDDO
C
C -> Construct transformation matrix AA+
C
      DO M2 = 1,MM
         DO M1 = 1,MM
            AAC(M1,M2) = DCONJG(AA(M2,M1))
         END DO
      END DO
C
C -> Multiply from left
C
      CALL CINIT(MDIM*MDIM,VTMP)
      DO M2 = 1,MM
         DO M4 = 1,MM
            BLJ = AAC(M4,M2)
            IF ( BLJ.NE.CZERO ) THEN
               DO M3 = 1,MM
                  VTMP(M3,M2) = VTMP(M3,M2) + VMAT(M3,M4)*BLJ
               END DO
            END IF
         END DO
      END DO
C
      CALL CINIT(MDIM*MDIM,VMAT)
      DO M2 = 1,MM
         DO M3 = 1,MM
            BLJ = VTMP(M3,M2)
            IF ( BLJ.NE.CZERO ) THEN
               DO M1 = 1,MM
                  VMAT(M1,M2) = VMAT(M1,M2) + AA(M1,M3)*BLJ
               END DO
            END IF
         END DO
      END DO
      END

