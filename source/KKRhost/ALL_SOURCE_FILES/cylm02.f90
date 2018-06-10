subroutine cylm02(lmax, cosx, fai, lpot2p, lmmaxd, thet, ylm, dylmt1, dylmt2, &
  dylmf1, dylmf2, dylmtf)
!.....------------------------------------------------------------------
!     preparation of cylm0(=ylm(ip,i)), cylmt1(=dylm/dtheta),
!     cylmt2(=d2ylm/dt2),
!     cylmf1, cylmf2 are for fai.
!     cylmtf=d2ylm/dfdt
!     i=1,2,....,(lmax+1)**2
!.....------------------------------------------------------------------
  implicit none

!.. Parameters ..
  integer :: ijd
  parameter (ijd=434)
!..
!.. Scalar Arguments ..
  integer :: lmax, lmmaxd, lpot2p
!..
!.. Array Arguments ..
  double precision :: cosx(ijd), dylmf1(ijd, lmmaxd), dylmf2(ijd, lmmaxd), &
    dylmt1(ijd, lmmaxd), dylmt2(ijd, lmmaxd), dylmtf(ijd, lmmaxd), fai(ijd), &
    thet(ijd), ylm(ijd, lmmaxd)
!..
!.. Local Scalars ..
  double complex :: ci, em1f, em2f, ep1f, ep2f
  double precision :: aaa, ccc, di, fi, one, pi, sss
  integer :: i, ip, l, llmax, lm, lm1, lm1m, lm2, lmm, lmm1, lmm1m, lmm2, m, &
    mm
!..
!.. Local Arrays ..
  double complex :: cylm0(lmmaxd), cylmf1(lmmaxd), cylmf2(lmmaxd), &
    cylmt1(lmmaxd), cylmt2(lmmaxd), cylmtf(lmmaxd)
  double precision :: bb1(lmmaxd), yl(lpot2p)
!..
!.. External Subroutines ..
  external :: spher, trarea
!..
!.. Intrinsic Functions ..
  intrinsic :: acos, atan, cmplx, conjg, cos, dble, sin, sqrt
!..

  ci = cmplx(0.d0, 1.d0)
  one = 1.d0
  pi = 4.d0*atan(one)
  llmax = (lmax+1)**2

  do ip = 1, ijd

    thet(ip) = acos(cosx(ip))
    fi = fai(ip)
    di = 2*fai(ip)
    ep1f = cmplx(cos(fi), sin(fi))
    em1f = conjg(ep1f)
    ep2f = cmplx(cos(di), sin(di))
    em2f = conjg(ep2f)

    do l = 0, lmax

      call spher(yl, l, cosx(ip))
      do m = -l, l
        mm = l + m + 1
        i = (l+1)**2 - l + m
        aaa = m*fai(ip)
        ccc = cos(aaa)
        sss = sin(aaa)
        cylm0(i) = yl(mm)*cmplx(ccc, sss)
      end do

      do m = -l, l
        i = (l+1)**2 - l + m
        cylmt1(i) = 0.d0
        cylmt2(i) = 0.d0
        cylmtf(i) = 0.d0
      end do

      do m = -l, l
        i = (l+1)**2 - l + m

        lmm1m = l - m - 1
        lmm = l - m
        lmm1 = l - m + 1
        lmm2 = l - m + 2
        lm1m = l + m - 1
        lm = l + m
        lm1 = l + m + 1
        lm2 = l + m + 2

        cylmt2(i) = cylmt2(i) - (lmm*lm1+lmm1*lm)/4.d0*cylm0(i)

        if (m+2<=l) cylmt2(i) = cylmt2(i) + sqrt(dble(lmm1m*lmm*lm1*lm2))/4* &
          cylm0(i+2)*em2f

        if (m+1<=l) cylmt1(i) = cylmt1(i) + sqrt(dble(lmm*lm1))/2*cylm0(i+1)* &
          em1f

        if (m-1>=-l) cylmt1(i) = cylmt1(i) - sqrt(dble(lm*lmm1))/2*cylm0(i-1)* &
          ep1f

        if (m-2>=-l) cylmt2(i) = cylmt2(i) + sqrt(dble(lmm1*lmm2*lm1m*lm))/4* &
          cylm0(i-2)*ep2f

      end do

      do m = -l, l
        i = (l+1)**2 - l + m
        cylmf1(i) = ci*m*cylm0(i)
        cylmf2(i) = -m*m*cylm0(i)
        cylmtf(i) = ci*m*cylmt1(i)
      end do

    end do

!        calculate real spherical harmonics differenciated


!        write(6,9005) (cylm0(i),i=1,5)
!9005 format(1x,' cylm0',4f10.5)
    call trarea(cylm0, bb1, lmax)

    do m = 1, llmax
      ylm(ip, m) = bb1(m)
    end do

!        write(6,9006) (ylm(ip,i),i=1,5)
!9006 format(1x,' ylm',10f10.5)


    call trarea(cylmt1, bb1, lmax)
    do m = 1, llmax
      dylmt1(ip, m) = bb1(m)
    end do

    call trarea(cylmt2, bb1, lmax)
    do m = 1, llmax
      dylmt2(ip, m) = bb1(m)
    end do

    call trarea(cylmf1, bb1, lmax)
    do m = 1, llmax
      dylmf1(ip, m) = bb1(m)
    end do

    call trarea(cylmf2, bb1, lmax)
    do m = 1, llmax
      dylmf2(ip, m) = bb1(m)
    end do

    call trarea(cylmtf, bb1, lmax)
    do m = 1, llmax
      dylmtf(ip, m) = bb1(m)
    end do

  end do
  return
end subroutine
