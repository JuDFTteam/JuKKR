  !-------------------------------------------------------------------------------
  !> Summary: Helper for PW91 correlation
  !> Author: 
  !> Category: KKRimp, xc-potential
  !> Deprecated: False
  !>
  !> Called by corlsd
  !-------------------------------------------------------------------------------
      SUBROUTINE GCOR91(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,A1,B1,B2,B3,B4,GG,GGRS,P,RS
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION P1,Q0,Q1,Q2,Q3,RS12,RS32,RSP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG,SQRT
C     ..
      P1 = P + 1.d0
      Q0 = -2.d0*A* (1.d0+A1*RS)
      RS12 = SQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.d0*A* (B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = LOG(1.d0+1.d0/Q1)
      GG = Q0*Q2
      Q3 = A* (B1/RS12+2.d0*B2+3.d0*B3*RS12+2.d0*B4*P1*RSP)
      GGRS = -2.d0*A*A1*Q2 - Q0*Q3/ (Q1**2+Q1)
      RETURN
      END SUBROUTINE
