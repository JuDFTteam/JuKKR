module mod_wfint0

contains

subroutine wfint0(cder, dder, qzlm, qzekdr, pzekdr, vnspll, nsra, irmind, &
  irmd, lmmaxd, irmin, irmax)      ! Added IRMIN,IRMAX 1.7.2014
  use :: mod_datatypes, only: dp
  ! -----------------------------------------------------------------------
  ! determines the integrands CDER, DDER or ADER, BDER in the
  ! integral equations for the non-spherical wavefunctions from
  ! the non-spherical contributions of the potential vinsPLL.
  ! (This subroutine is used in zeroth order Born approximation,
  ! otherwise subroutine WFINT must be used)
  ! R. Zeller      Aug. 1994
  ! -----------------------------------------------------------------------
  implicit none
  ! .. Scalar Arguments ..
  integer :: irmd, irmind, lmmaxd, nsra, irmin, irmax
  ! ..
  ! .. Array Arguments ..
  complex (kind=dp) :: cder(lmmaxd, lmmaxd, irmind:irmd), &
    dder(lmmaxd, lmmaxd, irmind:irmd), pzekdr(lmmaxd, irmind:irmd, 2), &
    qzekdr(lmmaxd, irmind:irmd, 2), qzlm(lmmaxd, irmind:irmd, 2)
  real (kind=dp) :: vnspll(lmmaxd, lmmaxd, irmind:irmd)
  ! ..
  ! .. Local Scalars ..
  complex (kind=dp) :: v1
  integer :: ir, lm1, lm2

  do ir = irmin, irmax
    do lm2 = 1, lmmaxd
      do lm1 = 1, lmmaxd
        v1 = vnspll(lm1, lm2, ir)*qzlm(lm2, ir, 1)
        cder(lm1, lm2, ir) = qzekdr(lm1, ir, 1)*v1
        dder(lm1, lm2, ir) = pzekdr(lm1, ir, 1)*v1
      end do
    end do
    if (nsra==2) then
      do lm2 = 1, lmmaxd
        do lm1 = 1, lmmaxd
          v1 = vnspll(lm1, lm2, ir)*qzlm(lm2, ir, 2)
          cder(lm1, lm2, ir) = cder(lm1, lm2, ir) + qzekdr(lm1, ir, 2)*v1
          dder(lm1, lm2, ir) = dder(lm1, lm2, ir) + pzekdr(lm1, ir, 2)*v1
        end do
      end do
    end if

  end do
end subroutine wfint0

end module mod_wfint0
