module BesselHankelNeumann

contains

      SUBROUTINE BESHANK(HL,JL,NL,Z,LMAX)
!-----------------------------------------------------------------------
!  calculates spherical bessel, hankel and neumann functions
!  for the orders lmin .le. l .le. lmax.
!  For |z| .lt. l+1 the taylor expansions of jl and nl are used.
!  For |z| .ge. l+1 the explicit expressions for hl(+), hl(-) are used.
!-----------------------------------------------------------------------
!     .. Parameters ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
!     ..
!     .. Scalar Arguments ..
      DOUBLE COMPLEX Z
      INTEGER LMAX
!     ..
!     .. Array Arguments ..
      DOUBLE COMPLEX HL(0:LMAX),JL(0:LMAX),NL(0:LMAX)
!     ..
!     .. Local Scalars ..
      DOUBLE COMPLEX TERMJ,TERMN,Z2,ZJ,ZN
      DOUBLE PRECISION RL,RN,RNM
      INTEGER L,M,N
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CDABS,EXP
!     ..
      ZJ = 1.D0
      ZN = 1.D0
      Z2 = Z*Z
      IF (CDABS(Z).LT.LMAX+1.D0) THEN
        DO 20 L = 0,LMAX
          RL = L + L
          TERMJ = -0.5D0/ (RL+3.D0)*Z2
          TERMN = 0.5D0/ (RL-1.D0)*Z2
          JL(L) = 1.D0
          NL(L) = 1.D0
          DO 10 N = 2,25
            JL(L) = JL(L) + TERMJ
            NL(L) = NL(L) + TERMN
            RN = N + N
            TERMJ = -TERMJ/ (RL+RN+1.D0)/RN*Z2
            TERMN = TERMN/ (RL-RN+1.D0)/RN*Z2
   10     CONTINUE
          JL(L) = JL(L)*ZJ
          NL(L) = -NL(L)*ZN/Z
          HL(L) = JL(L) + NL(L)*CI

          ZJ = ZJ*Z/ (RL+3.D0)
          ZN = ZN/Z* (RL+1.D0)
   20   CONTINUE
      END IF

      DO 40 L = 0,LMAX
        IF (CDABS(Z).GE.L+1.D0) THEN
          HL(L) = 0.D0
          NL(L) = 0.D0
          RNM = 1.D0
          DO 30 M = 0,L
            HL(L) = HL(L) + RNM/ (-CI* (Z+Z))**M
            NL(L) = NL(L) + RNM/ (CI* (Z+Z))**M
            RNM = RNM* (L*L+L-M*M-M)/ (M+1.D0)
   30     CONTINUE
          HL(L) = HL(L)* (-CI)**L*EXP(CI*Z)/ (CI*Z)
          NL(L) = NL(L)*CI**L*EXP(-CI*Z)/ (-CI*Z)
          JL(L) = (HL(L)+NL(L))*0.5D0
          NL(L) = (HL(L)-JL(L))/CI
        END IF
   40 CONTINUE

      RETURN

      END SUBROUTINE BESHANK

end module BesselHankelNeumann