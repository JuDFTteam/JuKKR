  !-------------------------------------------------------------------------------
  !> Summary: Exchange part of PW91 GGA
  !> Author: 
  !> Category: xc-potential, KKRimp
  !> Deprecated: False 
  !> gga91 exchange for a spin-unpolarized electronic system
  !> 
  !> input d: density
  !> s:  abs(grad d)/(2*kf*d)
  !> u:  (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
  !> v: (laplacian d)/(d*(2*kf)**2)
  !> output:  exchange energy (ex) and potential (vx) in ry.
  !> kf=cbrt(3*pai**2*d).
  !------------------------------------------------------------------------------- 
        SUBROUTINE EXCH91(D,S,U,V,EXL,EXG,VXL,VXG)
C     .. Scalar Arguments ..
      DOUBLE PRECISION D         !! Density
      DOUBLE PRECISION EXG,EXL,S,U,V,VXG,VXL
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A,A1,A2,A3,A4,AX,B1,EX,F,FAC,FS,FSS,P0,P1,P10,
     +                 P11,P2,P3,P4,P5,P6,P7,P8,P9,S2,S3,S4,THRD,THRD4,
     +                 VX
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP,LOG,SQRT
C     ..
C     .. Save statement ..
      SAVE A1,A2,A3,A4,AX,A,B1,THRD,THRD4
C     ..
C     .. Data statements ..
c.....-----------------------------------------------------------------
      DATA A1,A2,A3,A4/0.19645d0,0.27430d0,0.15084d0,100.d0/
      DATA AX,A,B1/-0.7385588d0,7.7956d0,0.004d0/
      DATA THRD,THRD4/0.333333333333d0,1.33333333333d0/
C     ..
c.....-----------------------------------------------------------------
      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.d0/SQRT(1.d0+A*A*S2)
      P1 = LOG(A*S+1.d0/P0)
      P2 = EXP(-A4*S2)
      P3 = 1.d0/ (1.d0+A1*S*P1+B1*S4)
      P4 = 1.d0 + A1*S*P1 + (A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*F
c  local exchange exl
      EXL = FAC
      EXG = EX - EXL

c  energy done. now the potential:
      P5 = B1*S2 - (A2-A3*P2)
      P6 = A1*S* (P1+A*S*P0)
      P7 = 2.d0* (A2-A3*P2) + 2.d0*A3*A4*S2*P2 - 4.d0*B1*S2*F
      FS = P3* (P3*P5*P6+P7)
      P8 = 2.d0*S* (B1-A3*A4*P2)
      P9 = A1*P1 + A*A1*S*P0* (3.d0-A*A*S2*P0*P0)
      P10 = 4.d0*A3*A4*S*P2* (2.d0-A4*S2) - 8.d0*B1*S*F - 4.d0*B1*S3*FS
      P11 = -P3*P3* (A1*P1+A*A1*S*P0+4.d0*B1*S3)
      FSS = P3*P3* (P5*P9+P6*P8) + 2.d0*P3*P5*P6*P11 + P3*P10 + P7*P11
      VX = FAC* (THRD4*F- (U-THRD4*S3)*FSS-V*FS)
c  local exchange vxl:
      VXL = FAC*THRD4
      VXG = VX - VXL

c in ry and energy density.
      EXL = EXL*2.d0*D
      EXG = EXG*2.d0*D
      VXL = VXL*2.d0
      VXG = VXG*2.d0
      RETURN
      END SUBROUTINE
