    Subroutine exch91(d, s, u, v, exl, exg, vxl, vxg)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------
!gga91 exchange for a spin-unpolarized electronic system
!-----------------------------------------------------------------
!input d: density
!      s:  abs(grad d)/(2*kf*d)
!      u:  (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
!      v: (laplacian d)/(d*(2*kf)**2)
!output:  exchange energy (ex) and potential (vx) in ry.
!  kf=cbrt(3*pai**2*d).
!-----------------------------------------------------------------
      Implicit None
!.. Scalar Arguments ..
      Real (Kind=dp) :: d, exg, exl, s, u, v, vxg, vxl
!..
!.. Local Scalars ..
      Real (Kind=dp) :: a, a1, a2, a3, a4, ax, b1, ex, f, fac, fs, fss, p0, &
        p1, p10, p11, p2, p3, p4, p5, p6, p7, p8, p9, s2, s3, s4, thrd, thrd4, &
        vx
!..
!.. Intrinsic Functions ..
      Intrinsic :: exp, log, sqrt
!..
!.. Save statement ..
      Save :: a1, a2, a3, a4, ax, a, b1, thrd, thrd4
!..
!.. Data statements ..
!-----------------------------------------------------------------
      Data a1, a2, a3, a4/0.19645E0_dp, 0.27430E0_dp, 0.15084E0_dp, 100.E0_dp/
      Data ax, a, b1/ -0.7385588E0_dp, 7.7956E0_dp, 0.004E0_dp/
      Data thrd, thrd4/0.333333333333E0_dp, 1.33333333333E0_dp/
!..
!-----------------------------------------------------------------
      fac = ax*d**thrd
      s2 = s*s
      s3 = s2*s
      s4 = s2*s2
      p0 = 1.E0_dp/sqrt(1.E0_dp+a*a*s2)
      p1 = log(a*s+1.E0_dp/p0)
      p2 = exp(-a4*s2)
      p3 = 1.E0_dp/(1.E0_dp+a1*s*p1+b1*s4)
      p4 = 1.E0_dp + a1*s*p1 + (a2-a3*p2)*s2
      f = p3*p4
      ex = fac*f
!  local exchange exl
      exl = fac
      exg = ex - exl

!  energy done. now the potential:
      p5 = b1*s2 - (a2-a3*p2)
      p6 = a1*s*(p1+a*s*p0)
      p7 = 2.E0_dp*(a2-a3*p2) + 2.E0_dp*a3*a4*s2*p2 - 4.E0_dp*b1*s2*f
      fs = p3*(p3*p5*p6+p7)
      p8 = 2.E0_dp*s*(b1-a3*a4*p2)
      p9 = a1*p1 + a*a1*s*p0*(3.E0_dp-a*a*s2*p0*p0)
      p10 = 4.E0_dp*a3*a4*s*p2*(2.E0_dp-a4*s2) - 8.E0_dp*b1*s*f - &
        4.E0_dp*b1*s3*fs
      p11 = -p3*p3*(a1*p1+a*a1*s*p0+4.E0_dp*b1*s3)
      fss = p3*p3*(p5*p9+p6*p8) + 2.E0_dp*p3*p5*p6*p11 + p3*p10 + p7*p11
      vx = fac*(thrd4*f-(u-thrd4*s3)*fss-v*fs)
!  local exchange vxl:
      vxl = fac*thrd4
      vxg = vx - vxl

! in ry and energy density.
      exl = exl*2.E0_dp*d
      exg = exg*2.E0_dp*d
      vxl = vxl*2.E0_dp
      vxg = vxg*2.E0_dp
      Return
    End Subroutine
