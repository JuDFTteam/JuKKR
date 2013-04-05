      SUBROUTINE BESSEL(JL,NL,HL,Z,LMX)
c**********************************************************************
c
c    attention : contrary to abramowitz and stegun and
c                contrary to SUBROUTINE BESHAN
c
c                the bessel functions of third kind ( hankel functions)
c                are defined as:      hl(l) = nl(l) - i * jl(l)
c
c**********************************************************************
C     .. Parameters ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.D0,1.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX Z
      INTEGER LMX
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX HL(0:LMX),JL(0:LMX),NL(0:LMX)
C     ..
C     .. Local Scalars ..
      INTEGER L
C     ..
C     .. External Subroutines ..
      EXTERNAL BESHAN
C     ..
      CALL BESHAN(HL,JL,NL,Z,LMX)
      DO L = 0,LMX
        HL(L) = -CI*HL(L)
      END DO

      END
