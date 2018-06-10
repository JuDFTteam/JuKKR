subroutine wftsca(drdi, efac, pz, qz, fz, sz, nsra, pzlm, qzlm, pzekdr, &
  qzekdr, ek, loflm, irmind, irmd, irmin, irmax, & ! Added IRMIN,IRMAX 1.7.2014  &
  lmaxd, lmmaxd)
!-----------------------------------------------------------------------
!                 R. Zeller      Oct. 1993
!-----------------------------------------------------------------------
  implicit none
!.. Parameters ..
  double complex :: cone
  parameter (cone=(1.d0,0.d0))
!..
!.. Scalar Arguments ..
  double complex :: ek
  integer :: irmd, irmind, lmaxd, lmmaxd, nsra, irmin, irmax
!..
!.. Array Arguments ..
  double complex :: efac(lmmaxd), fz(irmd, 0:lmaxd), pz(irmd, 0:lmaxd), &
    pzekdr(lmmaxd, irmind:irmd, 2), pzlm(lmmaxd, irmind:irmd, 2), &
    qz(irmd, 0:lmaxd), qzekdr(lmmaxd, irmind:irmd, 2), &
    qzlm(lmmaxd, irmind:irmd, 2), sz(irmd, 0:lmaxd)
  double precision :: drdi(*)
  integer :: loflm(*)
!..
!.. Local Scalars ..
  double complex :: efac1, v1
  integer :: ir, j, l, l1, lm, lm1, m
!..
!.. Intrinsic Functions ..
  intrinsic :: dble
!..


!---> set up array efac : efac(lm) = sqrt(e)**l/(2l - 1)!!

  efac(1) = cone
  v1 = cone
  do l = 1, lmaxd
    v1 = v1*ek/dble(2*l-1)
    do m = -l, l
      lm = l*(l+1) + m + 1
      efac(lm) = v1
    end do
  end do


!---> get wfts of same magnitude by scaling with efac

  do lm1 = 1, lmmaxd
    l1 = loflm(lm1)
    efac1 = efac(lm1)
    do ir = irmin, irmax
      pzlm(lm1, ir, 1) = pz(ir, l1)/efac1
      qzlm(lm1, ir, 1) = qz(ir, l1)*efac1
    end do
    if (nsra==2) then
      do ir = irmin, irmax
        pzlm(lm1, ir, nsra) = fz(ir, l1)/efac1
        qzlm(lm1, ir, nsra) = sz(ir, l1)*efac1
      end do
    end if

    do j = 1, nsra
      do ir = irmin, irmax
        pzekdr(lm1, ir, j) = pzlm(lm1, ir, j)*ek*drdi(ir)
        qzekdr(lm1, ir, j) = qzlm(lm1, ir, j)*ek*drdi(ir)
      end do
    end do
  end do


end subroutine
