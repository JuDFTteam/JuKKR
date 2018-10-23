!------------------------------------------------------------------------------------
!> Summary: Determines the integrands `CDER`, `DDER` or `ADER`, `BDER`
!> Author: R. Zeller
!> Determines the integrands `CDER`, `DDER` or `ADER`, `BDER` in the integral equations 
!> for the non-spherical wavefunctions from the non-spherical contributions of the 
!> potential `vinsPLL`.
!------------------------------------------------------------------------------------
module mod_wfint
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Determines the integrands `CDER`, `DDER` or `ADER`, `BDER`
  !> Author: R. Zeller
  !> Category: numerical-tools, single-site, KKRhost 
  !> Deprecated: False 
  !> Determines the integrands `CDER`, `DDER` or `ADER`, `BDER` in the integral equations 
  !> for the non-spherical wavefunctions from the non-spherical contributions of the 
  !> potential `vinsPLL`. (This subroutine is used in zeroth order Born approximation,
  !> otherwise subroutine WFINT must be used)
  !-------------------------------------------------------------------------------
  subroutine wfint0(cder,dder,qzlm,qzekdr,pzekdr,vnspll,nsra,irmind,irmd,lmmaxd,    &
    irmin,irmax) ! Added IRMIN,IRMAX 1.7.2014

    implicit none
    ! .. Input variables
    integer, intent(in) :: irmd !! Maximum number of radial points
    integer, intent(in) :: nsra
    integer, intent(in) :: irmin !! Max R for spherical treatment
    integer, intent(in) :: irmax
    integer, intent(in) :: irmind !! irmd - irnsd
    integer, intent(in) :: lmmaxd !! (KREL+KORBIT+1)*(LMAX+1)**2
    real (kind=dp), dimension(lmmaxd, lmmaxd, irmind:irmd), intent(in) :: vnspll
    complex (kind=dp), dimension(lmmaxd, irmind:irmd, 2), intent(in) :: pzekdr
    complex (kind=dp), dimension(lmmaxd, irmind:irmd, 2), intent(in) :: qzekdr
    complex (kind=dp), dimension(lmmaxd, irmind:irmd, 2), intent(in) :: qzlm
    ! .. Output variables
    complex (kind=dp), dimension(lmmaxd, lmmaxd, irmind:irmd), intent(out) :: cder
    complex (kind=dp), dimension(lmmaxd, lmmaxd, irmind:irmd), intent(out) :: dder
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

    !-------------------------------------------------------------------------------
  !> Summary: Determines the integrands `CDER`, `DDER` or `ADER`, `BDER`
  !> Author: R. Zeller
  !> Category: numerical-tools, single-site, KKRhost 
  !> Deprecated: False 
  !> determines the integrands `CDER`, `DDER` or `ADER`, `BDER` in the integral equations 
  !> for the non-spherical wavefunctions from the non-spherical contributions of the 
  !> potential `vinsPLL`. 
  !-------------------------------------------------------------------------------
  subroutine wfint(qns,cder,dder,qzekdr,pzekdr,vnspll,nsra,irmind,irmd,lmmaxd,irmin,&
    irmax) ! Added IRMIN,IRMAX 1.7.2014

    implicit none

    ! .. Input variables
    integer, intent(in) :: irmd !! Maximum number of radial points
    integer, intent(in) :: nsra
    integer, intent(in) :: irmin !! Max R for spherical treatment
    integer, intent(in) :: irmax
    integer, intent(in) :: irmind !! irmd - irnsd
    integer, intent(in) :: lmmaxd !! (KREL+KORBIT+1)*(LMAX+1)**2
    real (kind=dp), dimension(lmmaxd, lmmaxd, irmind:irmd), intent(in) :: vnspll
    complex (kind=dp), dimension(lmmaxd, irmind:irmd, 2), intent(in) :: pzekdr
    complex (kind=dp), dimension(lmmaxd, irmind:irmd, 2), intent(in) :: qzekdr
    complex (kind=dp), dimension(lmmaxd, lmmaxd, irmind:irmd, 2), intent(in) :: qns
    ! .. Output variables
    complex (kind=dp), dimension(lmmaxd, lmmaxd, irmind:irmd), intent(out) :: cder
    complex (kind=dp), dimension(lmmaxd, lmmaxd, irmind:irmd), intent(out) :: dder
    ! .. Local Scalars ..
    integer :: ir, lm1, lm2
    ! .. Local Arrays ..
    real (kind=dp), dimension(lmmaxd, lmmaxd) :: qnsi, qnsr, vtqnsi, vtqnsr

    do ir = irmin, irmax
      do lm2 = 1, lmmaxd
        do lm1 = 1, lmmaxd
          qnsr(lm1, lm2) = real(qns(lm1,lm2,ir,1))
          qnsi(lm1, lm2) = aimag(qns(lm1,lm2,ir,1))
        end do
      end do
      call dgemm('N','N',lmmaxd,lmmaxd,lmmaxd,1.e0_dp,vnspll(1,1,ir),lmmaxd,qnsr,   &
        lmmaxd,0.e0_dp,vtqnsr,lmmaxd)
      call dgemm('N','N',lmmaxd,lmmaxd,lmmaxd,1.e0_dp,vnspll(1,1,ir),lmmaxd,qnsi,   &
        lmmaxd,0.e0_dp,vtqnsi,lmmaxd)
      do lm1 = 1, lmmaxd
        do lm2 = 1, lmmaxd
          cder(lm1, lm2, ir) = qzekdr(lm1, ir, 1)*cmplx(vtqnsr(lm1,lm2), vtqnsi(lm1,lm2), kind=dp)
          dder(lm1, lm2, ir) = pzekdr(lm1, ir, 1)*cmplx(vtqnsr(lm1,lm2), vtqnsi(lm1,lm2), kind=dp)
        end do
      end do
      if (nsra==2) then
        do lm2 = 1, lmmaxd
          do lm1 = 1, lmmaxd
            qnsr(lm1, lm2) = real(qns(lm1,lm2,ir,2))
            qnsi(lm1, lm2) = aimag(qns(lm1,lm2,ir,2))
          end do
        end do
        call dgemm('N','N',lmmaxd,lmmaxd,lmmaxd,1.e0_dp,vnspll(1,1,ir),lmmaxd,qnsr, &
          lmmaxd,0.e0_dp,vtqnsr,lmmaxd)
        call dgemm('N','N',lmmaxd,lmmaxd,lmmaxd,1.e0_dp,vnspll(1,1,ir),lmmaxd,qnsi, &
          lmmaxd,0.e0_dp,vtqnsi,lmmaxd)
        do lm2 = 1, lmmaxd
          do lm1 = 1, lmmaxd
            cder(lm1, lm2, ir) = cder(lm1, lm2, ir) + qzekdr(lm1, ir, 2)*cmplx(vtqnsr(lm1,lm2), vtqnsi(lm1,lm2), kind=dp)
            dder(lm1, lm2, ir) = dder(lm1, lm2, ir) + pzekdr(lm1, ir, 2)*cmplx(vtqnsr(lm1,lm2), vtqnsi(lm1,lm2), kind=dp)
          end do
        end do
      end if

    end do
  end subroutine wfint

end module mod_wfint
