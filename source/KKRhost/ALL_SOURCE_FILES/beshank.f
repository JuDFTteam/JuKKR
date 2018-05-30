      SUBROUTINE BESHANK(HL,JL,Z,LMAX)
c-----------------------------------------------------------------------
c  calculates spherical bessel, hankel and neumann functions
c  for the orders lmin .le. l .le. lmax.
c  For |z| .lt. l+1 the taylor expansions of jl and nl are used.
c  For |z| .ge. l+1 the explicit expressions for hl(+), hl(-) are used.
c-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX Z
      INTEGER LMAX
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX HL(0:LMAX),JL(0:LMAX),NL(0:LMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX TERMJ,TERMN,Z2,ZJ,ZN
      DOUBLE PRECISION RL,RN,RNM
      INTEGER L,M,N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CDABS,EXP
C     ..
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

      END SUBROUTINE

      SUBROUTINE BESHANK_SMALLCOMP(HL,JL,ZVAL,TAU,ERYD,LMAX)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c  takes the spherical bessel etc functions stored in an array up to LMAX
c  array entries from LMAX+1 to 2*LMAX are assumed to be empty
c  these values are filled with the potential-free solution of the 
c  SRA-equations
c-----------------------------------------------------------------------
      DOUBLE COMPLEX HL(0:2*(LMAX+1)-1),
     +               JL(0:2*(LMAX+1)-1),
     +               NL(0:2*(LMAX+1)-1)
      DOUBLE PRECISION CVLIGHT
      PARAMETER (CVLIGHT=274.0720442D0)
      DOUBLE COMPLEX ZVAL
      DOUBLE COMPLEX ERYD
      DOUBLE PRECISION TAU
      INTEGER LMAX

!       DOUBLE PRECISION CVLIGHT
      DOUBLE COMPLEX PREFAC
      INTEGER IL,IL2


      PREFAC = 1.0D0 / (1.0D0+eryd/cvlight**2) / tau !/cvlight  !last cvlight for small component test

      IL=0
      IL2=IL+LMAX+1
      NL(IL2)=PREFAC * (ZVAL* (-NL(IL+1)) )
      JL(IL2)=PREFAC * (ZVAL* (-JL(IL+1)) )
!       HL(IL2)=JL(IL2)+ CI*NL(IL2)
      HL(IL2)=PREFAC * (ZVAL* (-HL(IL+1)) )
!       write(*,'(5000E)') tau,HL(IL2),JL(IL2)+ (0.0D0,1.0D0)*NL(IL2)
!       write(*,'(5000E)') tau,HL(0),JL(0)+ (0.0D0,1.0D0)*NL(0)

      PREFAC = 1.0D0 / (1.0D0+eryd/cvlight**2) / tau !/cvlight !last cvlight for small component test

      DO IL=1,LMAX
        IL2=IL+LMAX+1
        NL(IL2)=PREFAC * ( ZVAL * NL(IL-1)-(IL+1)*NL(IL) )
        JL(IL2)=PREFAC * ( ZVAL * JL(IL-1)-(IL+1)*JL(IL) )
!         HL(IL2)=JL(IL2)+ CI*NL(IL2)
        HL(IL2)=PREFAC * ( ZVAL * HL(IL-1)-(IL+1)*HL(IL) )
!         HL(IL2)=PREFAC * ( ZVAL * HL(IL-1)-(IL+1)*HL(IL) )
!         write(*,'(5000E)') tau,HL(IL2),JL(IL2)+ (0.0D0,1.0D0)*NL(IL2)
      END DO

      END SUBROUTINE BESHANK_SMALLCOMP

