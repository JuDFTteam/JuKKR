      SUBROUTINE cSPHER(YLM,L,X)
!
!      spherical harmonics except the facter exp(i*m*phi)
!
!      m=-l to l , for given l.
!      x=cos(theta)
!     .. Scalar Arguments ..
      DOUBLE PRECISION X
      INTEGER L
!     ..
!     .. Array Arguments ..
      DOUBLE complex YLM(*)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION FAC,OVR1,PI,QQ
      INTEGER I,II,L2,LM,LN,M,NN
!     ..
!     .. Intrinsi! Functions ..
      INTRINSIC ABS,ATAN,DBLE,SQRT
!     ..

      PI = 4.0D0*ATAN(1.0D0)

      OVR1 = ABS(X) - 1.D0
      IF (OVR1.GT.0.1D-12) THEN
        WRITE (6,FMT=9000) X
        STOP
      ELSE IF (ABS(OVR1).LT.1.D-10) THEN
        IF (X.GT.0.0D0) THEN
          FAC = 1.0D0
        ELSE
          FAC = (-1)**L
        END IF
        L2 = 2*L + 1
        DO I = 1,L2
          YLM(I) = 0.0D0
        END DO
        YLM(L+1) = SQRT(DBLE(L2)/ (4.0D0*PI))*FAC
        RETURN
      END IF

! l<0
      IF (L.LT.0) THEN
        WRITE (6,FMT=*) ' === l=',L,' < 0  : in sub.spher. ==='
        STOP '=== stop in sub.spher. (l<0) ==='
! l=0
      ELSE IF (L.EQ.0) THEN
        YLM(1) = SQRT(1.0D0/ (4.0D0*PI))
! l=1
      ELSE IF (L.EQ.1) THEN
        FAC = SQRT(3.0D0/ (4.0D0*PI))
        YLM(1) = FAC*SQRT((1.0D0-X*X)/2.0D0)
        YLM(2) = FAC*X
        YLM(3) = -YLM(1)
! l>1
      ELSE
        YLM(1) = 1.0D0
        YLM(2) = X
        DO I = 2,L
          YLM(I+1) = ((2*I-1)*X*YLM(I)- (I-1)*YLM(I-1))/I
        END DO
        FAC = 1.0D0/SQRT(1.0D0-X*X)
        DO M = 1,L
          LM = L + M
          YLM(LM+1) = FAC* (- (L-M+1)*X*YLM(LM)+ (LM-1)*YLM(L))
          IF (M.LT.L) THEN
            NN = M + 1
            DO I = NN,L
              II = L - I + NN
              YLM(II) = FAC* (- (II-M)*X*YLM(II)+ (II+M-2)*YLM(II-1))
            END DO
          END IF
        END DO
        FAC = SQRT((2*L+1)/ (4.0D0*PI))
        YLM(L+1) = FAC*YLM(L+1)
        DO M = 1,L
          FAC = -FAC/SQRT(DBLE((L+M)* (L-M+1)))
          LM = L + 1 + M
          LN = L + 1 - M
          QQ = YLM(LM)
          YLM(LM) = FAC*QQ
          YLM(LN) = ABS(FAC)*QQ
        END DO
      END IF

      RETURN
 9000 FORMAT (/,/,3x,'==invalid argument for spher; x=',D24.16,' ==')

      END SUBROUTINE cSPHER
