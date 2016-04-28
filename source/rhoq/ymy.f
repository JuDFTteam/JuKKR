c **********************************************************************
      SUBROUTINE YMY(V1,V2,V3,R,YLM,LMAX)
c **********************************************************************
c    this subroutine calculates real spherical harmonics with the
c     normalization : <y|y> =1
c    returns also r = length of vector v
c
c     generate the complex spherical harmonics for the vector v
c     using a stable upward recursion in l.  (see notes
c     by m. weinert.)
c                                  m.weinert  1982
c
c     converted to real spherical harmonics .
c                                  b.drittler 1987
c-----------------------------------------------------------------------
      IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION SZERO
      PARAMETER (SZERO=1.0D-20)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION, intent(out) :: R
      DOUBLE PRECISION, intent(in) :: V1,V2,V3
      INTEGER, intent(in) :: LMAX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION, intent(out) :: YLM((2*LMAX+1)**2)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,CD,CPH,CTH,FAC,FPI,PI,RTWO,SGM,SPH,STH,T,XY,XYZ
      INTEGER I,L,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION C(0:LMAX),P(0:LMAX,0:LMAX),S(0:LMAX)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
!       
!       A   = 0.0d0
!       CD  = 0.0d0
!       CPH = 0.0d0
!       CTH = 0.0d0
!       FAC = 0.0d0
!       FPI = 0.0d0
!       PI  = 0.0d0
!       RTWO= 0.0d0
!       SGM = 0.0d0
!       SPH = 0.0d0
!       STH = 0.0d0
!       T   = 0.0d0
!       XY  = 0.0d0
!       XYZ = 0.0d0
!       I   = 0
!       L   = 0
!       M   = 0
! 
!       C(0:LMAX)= 0.0d0
!       P(0:LMAX,0:LMAX)= 0.0d0
!       S(0:LMAX)= 0.0d0
!       
!       
!       
!       
C     ..
      PI = 4.D0*ATAN(1.D0)
      FPI = 4.D0*PI
      RTWO = SQRT(2.0D0)
      
c
c--->    calculate sin and cos of theta and phi
c
      XY = V1**2 + V2**2
      XYZ = XY + V3**2
!       write(*,'(A,4ES15.7)')'ymy input', V1,V2,V3, XYZ
c
      R = SQRT(XYZ)
      IF (XYZ.LE.0.0D0) THEN
        write(*,*) xyz
        CALL RCSTOP('ylm=0   ')

      ELSE

        IF (XY.GT.SZERO*XYZ) THEN
          XY = SQRT(XY)
          XYZ = SQRT(XYZ)
          CTH = V3/XYZ
          STH = XY/XYZ
          CPH = V1/XY
          SPH = V2/XY

        ELSE

          STH = 0.0D0
          CTH = 1.0D0
          IF (V3.LT.0) CTH = -1.0D0
          CPH = 1.0D0
          SPH = 0.0D0
        END IF
c
c--->    generate associated legendre functions for m.ge.0
c        loop over m values
c
        FAC = 1.0D0
        DO 20 M = 0,LMAX - 1
          FAC = - (2*M-1)*FAC
          P(M,M) = FAC
          P(M+1,M) = (2*M+1)*CTH*FAC
c
c--->    recurse upward in l
c
          DO 10 L = M + 2,LMAX
            P(L,M) = ((2*L-1)*CTH*P(L-1,M)- (L+M-1)*P(L-2,M))/ (L-M)
   10     CONTINUE
          FAC = FAC*STH
   20   CONTINUE
        P(LMAX,LMAX) = - (2*LMAX-1)*FAC
c
c--->    determine sin and cos of phi
c
        S(0) = 0.0D0
        S(1) = SPH
        C(0) = 1.0D0
        C(1) = CPH
        DO 30 M = 2,LMAX
          S(M) = 2*CPH*S(M-1) - S(M-2)
          C(M) = 2*CPH*C(M-1) - C(M-2)
   30   CONTINUE
c
c--->    multiply in the normalization factors
c
        I = 0
        DO 50 L = 0,LMAX
          I = I + L + 1
          A = SQRT((2*L+1)/FPI)
          CD = 1
          YLM(I) = A*P(L,0)
          SGM = -RTWO
          DO 40 M = 1,L
            T = (L+1-M)* (L+M)
            CD = CD/T
            T = A*SQRT(CD)
            YLM(I+M) = SGM*T*P(L,M)*C(M)
            YLM(I-M) = SGM*T*P(L,M)*S(M)
            SGM = -SGM
   40     CONTINUE
          I = I + L
   50   CONTINUE

      END IF

      RETURN

      END
