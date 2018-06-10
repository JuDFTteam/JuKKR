subroutine wfint(qns, cder, dder, qzekdr, pzekdr, vnspll, nsra, irmind, irmd, &
  lmmaxd, irmin, irmax) ! Added IRMIN,IRMAX 1.7.2014
!-----------------------------------------------------------------------
!      determines the integrands CDER, DDER or ADER, BDER in the
!        integral equations for the non-spherical wavefunctions from
!        the non-spherical contributions of the potential vinsPLL.

!      R. Zeller      Aug. 1994
!-----------------------------------------------------------------------
  implicit none
!.. Scalar Arguments ..
  integer :: irmd, irmind, lmmaxd, nsra, irmin, irmax
!..
!.. Array Arguments ..
  double complex :: cder(lmmaxd, lmmaxd, irmind:irmd), &
    dder(lmmaxd, lmmaxd, irmind:irmd), pzekdr(lmmaxd, irmind:irmd, 2), &
    qns(lmmaxd, lmmaxd, irmind:irmd, 2), qzekdr(lmmaxd, irmind:irmd, 2)
  double precision :: vnspll(lmmaxd, lmmaxd, irmind:irmd)
!..
!.. Local Scalars ..
  integer :: ir, lm1, lm2
!..
!.. External Subroutines ..
  external :: dgemm
!..
!.. Local Arrays ..
  double precision :: qnsi(lmmaxd, lmmaxd), qnsr(lmmaxd, lmmaxd), &
    vtqnsi(lmmaxd, lmmaxd), vtqnsr(lmmaxd, lmmaxd)
!..
!.. Intrinsic Functions ..
  intrinsic :: dcmplx, dimag, dble
!     ..

  do ir = irmin, irmax
    do lm2 = 1, lmmaxd
      do lm1 = 1, lmmaxd
        qnsr(lm1, lm2) = dble(qns(lm1,lm2,ir,1))
        qnsi(lm1, lm2) = dimag(qns(lm1,lm2,ir,1))
      end do
    end do
    call dgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, 1.d0, vnspll(1,1,ir), lmmaxd, &
      qnsr, lmmaxd, 0.d0, vtqnsr, lmmaxd)
    call dgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, 1.d0, vnspll(1,1,ir), lmmaxd, &
      qnsi, lmmaxd, 0.d0, vtqnsi, lmmaxd)
    do lm1 = 1, lmmaxd
      do lm2 = 1, lmmaxd
        cder(lm1, lm2, ir) = qzekdr(lm1, ir, 1)*dcmplx(vtqnsr(lm1,lm2), vtqnsi &
          (lm1,lm2))
        dder(lm1, lm2, ir) = pzekdr(lm1, ir, 1)*dcmplx(vtqnsr(lm1,lm2), vtqnsi &
          (lm1,lm2))
      end do
    end do
    if (nsra==2) then
      do lm2 = 1, lmmaxd
        do lm1 = 1, lmmaxd
          qnsr(lm1, lm2) = dble(qns(lm1,lm2,ir,2))
          qnsi(lm1, lm2) = dimag(qns(lm1,lm2,ir,2))
        end do
      end do
      call dgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, 1.d0, vnspll(1,1,ir), &
        lmmaxd, qnsr, lmmaxd, 0.d0, vtqnsr, lmmaxd)
      call dgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, 1.d0, vnspll(1,1,ir), &
        lmmaxd, qnsi, lmmaxd, 0.d0, vtqnsi, lmmaxd)
      do lm2 = 1, lmmaxd
        do lm1 = 1, lmmaxd
          cder(lm1, lm2, ir) = cder(lm1, lm2, ir) + qzekdr(lm1, ir, 2)*dcmplx( &
            vtqnsr(lm1,lm2), vtqnsi(lm1,lm2))
          dder(lm1, lm2, ir) = dder(lm1, lm2, ir) + pzekdr(lm1, ir, 2)*dcmplx( &
            vtqnsr(lm1,lm2), vtqnsi(lm1,lm2))
        end do
      end do
    end if

  end do
end subroutine
