    Subroutine cylm02(lmax, cosx, fai, lpot2p, lmmaxd, thet, ylm, dylmt1, &
      dylmt2, dylmf1, dylmf2, dylmtf)
      Use mod_datatypes, Only: dp
!.....------------------------------------------------------------------
!     preparation of cylm0(=ylm(ip,i)), cylmt1(=dylm/dtheta),
!     cylmt2(=d2ylm/dt2),
!     cylmf1, cylmf2 are for fai.
!     cylmtf=d2ylm/dfdt
!     i=1,2,....,(lmax+1)**2
!.....------------------------------------------------------------------
      Implicit None

!.. Parameters ..
      Integer :: ijd
      Parameter (ijd=434)
!..
!.. Scalar Arguments ..
      Integer :: lmax, lmmaxd, lpot2p
!..
!.. Array Arguments ..
      Real (Kind=dp) :: cosx(ijd), dylmf1(ijd, lmmaxd), dylmf2(ijd, lmmaxd), &
        dylmt1(ijd, lmmaxd), dylmt2(ijd, lmmaxd), dylmtf(ijd, lmmaxd), &
        fai(ijd), thet(ijd), ylm(ijd, lmmaxd)
!..
!.. Local Scalars ..
      Complex (Kind=dp) :: ci, em1f, em2f, ep1f, ep2f
      Real (Kind=dp) :: aaa, ccc, di, fi, one, pi, sss
      Integer :: i, ip, l, llmax, lm, lm1, lm1m, lm2, lmm, lmm1, lmm1m, lmm2, &
        m, mm
!..
!.. Local Arrays ..
      Complex (Kind=dp) :: cylm0(lmmaxd), cylmf1(lmmaxd), cylmf2(lmmaxd), &
        cylmt1(lmmaxd), cylmt2(lmmaxd), cylmtf(lmmaxd)
      Real (Kind=dp) :: bb1(lmmaxd), yl(lpot2p)
!..
!.. External Subroutines ..
      External :: spher, trarea
!..
!.. Intrinsic Functions ..
      Intrinsic :: acos, atan, cmplx, conjg, cos, real, sin, sqrt
!..

      ci = cmplx(0.E0_dp, 1.E0_dp, kind=dp)
      one = 1.E0_dp
      pi = 4.E0_dp*atan(one)
      llmax = (lmax+1)**2

      Do ip = 1, ijd

        thet(ip) = acos(cosx(ip))
        fi = fai(ip)
        di = 2*fai(ip)
        ep1f = cmplx(cos(fi), sin(fi), kind=dp)
        em1f = conjg(ep1f)
        ep2f = cmplx(cos(di), sin(di), kind=dp)
        em2f = conjg(ep2f)

        Do l = 0, lmax

          Call spher(yl, l, cosx(ip))
          Do m = -l, l
            mm = l + m + 1
            i = (l+1)**2 - l + m
            aaa = m*fai(ip)
            ccc = cos(aaa)
            sss = sin(aaa)
            cylm0(i) = yl(mm)*cmplx(ccc, sss, kind=dp)
          End Do

          Do m = -l, l
            i = (l+1)**2 - l + m
            cylmt1(i) = 0.E0_dp
            cylmt2(i) = 0.E0_dp
            cylmtf(i) = 0.E0_dp
          End Do

          Do m = -l, l
            i = (l+1)**2 - l + m

            lmm1m = l - m - 1
            lmm = l - m
            lmm1 = l - m + 1
            lmm2 = l - m + 2
            lm1m = l + m - 1
            lm = l + m
            lm1 = l + m + 1
            lm2 = l + m + 2

            cylmt2(i) = cylmt2(i) - (lmm*lm1+lmm1*lm)/4.E0_dp*cylm0(i)

            If (m+2<=l) cylmt2(i) = cylmt2(i) + sqrt(real(lmm1m*lmm*lm1*lm2, &
              kind=dp))/4*cylm0(i+2)*em2f

            If (m+1<=l) cylmt1(i) = cylmt1(i) + sqrt(real(lmm*lm1,kind=dp))/2* &
              cylm0(i+1)*em1f

            If (m-1>=-l) cylmt1(i) = cylmt1(i) - sqrt(real(lm*lmm1,kind=dp))/2 &
              *cylm0(i-1)*ep1f

            If (m-2>=-l) cylmt2(i) = cylmt2(i) + sqrt(real(lmm1*lmm2*lm1m*lm, &
              kind=dp))/4*cylm0(i-2)*ep2f

          End Do

          Do m = -l, l
            i = (l+1)**2 - l + m
            cylmf1(i) = ci*m*cylm0(i)
            cylmf2(i) = -m*m*cylm0(i)
            cylmtf(i) = ci*m*cylmt1(i)
          End Do

        End Do

!        calculate real spherical harmonics differenciated


!        write(6,9005) (cylm0(i),i=1,5)
!9005 format(1x,' cylm0',4f10.5)
        Call trarea(cylm0, bb1, lmax)

        Do m = 1, llmax
          ylm(ip, m) = bb1(m)
        End Do

!        write(6,9006) (ylm(ip,i),i=1,5)
!9006 format(1x,' ylm',10f10.5)


        Call trarea(cylmt1, bb1, lmax)
        Do m = 1, llmax
          dylmt1(ip, m) = bb1(m)
        End Do

        Call trarea(cylmt2, bb1, lmax)
        Do m = 1, llmax
          dylmt2(ip, m) = bb1(m)
        End Do

        Call trarea(cylmf1, bb1, lmax)
        Do m = 1, llmax
          dylmf1(ip, m) = bb1(m)
        End Do

        Call trarea(cylmf2, bb1, lmax)
        Do m = 1, llmax
          dylmf2(ip, m) = bb1(m)
        End Do

        Call trarea(cylmtf, bb1, lmax)
        Do m = 1, llmax
          dylmtf(ip, m) = bb1(m)
        End Do

      End Do
      Return
    End Subroutine
