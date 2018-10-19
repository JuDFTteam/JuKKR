!------------------------------------------------------------------------------------
!> Summary: Get wavefunctions of same magnitude by scaling with `efac` 
!> Author: R. Zeller
!> Get wavefunctions of same magnitude by scaling with `efac`
!------------------------------------------------------------------------------------
module mod_wftsca
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Get wavefunctions of same magnitude by scaling with `efac`
  !> Author: R. Zeller
  !> Category: numerical-tools, single-site, KKRhost
  !> Deprecated: False 
  !> Get wavefunctions of same magnitude by scaling with `efac`
  !-------------------------------------------------------------------------------
  subroutine wftsca(drdi,efac,pz,qz,fz,sz,nsra,pzlm,qzlm,pzekdr,qzekdr,ek,loflm,    &
    irmind,irmd,irmin,irmax,lmaxd,lmmaxd)

    use :: constants, only: cone

    implicit none
    ! .. Input variables
    integer, intent(in) :: irmd !! Maximum number of radial points
    integer, intent(in) :: nsra
    integer, intent(in) :: irmin !! Max R for spherical treatment
    integer, intent(in) :: irmax
    integer, intent(in) :: lmaxd  !! Maximum l component in wave function expansion
    integer, intent(in) :: irmind !! irmd - irnsd
    integer, intent(in) :: lmmaxd !! (KREL+KORBIT+1)*(LMAX+1)**2
    complex (kind=dp), intent(in) :: ek
    integer, dimension(*), intent(in) :: loflm !! l of lm=(l,m) (GAUNT)
    real (kind=dp), dimension(*), intent(in) :: drdi  !! Derivative dr/di
    complex (kind=dp), dimension(irmd, 0:lmaxd), intent(in) :: fz
    complex (kind=dp), dimension(irmd, 0:lmaxd), intent(in) :: pz
    complex (kind=dp), dimension(irmd, 0:lmaxd), intent(in) :: sz
    complex (kind=dp), dimension(irmd, 0:lmaxd), intent(in) :: qz
    ! .. Output variables
    complex (kind=dp), dimension(lmmaxd), intent(out) :: efac !! efac(lm) = sqrt(e)**l/(2l - 1)
    complex (kind=dp), dimension(lmmaxd, irmind:irmd, 2), intent(out) :: pzlm
    complex (kind=dp), dimension(lmmaxd, irmind:irmd, 2), intent(out) :: qzlm
    complex (kind=dp), dimension(lmmaxd, irmind:irmd, 2), intent(out) :: qzekdr
    complex (kind=dp), dimension(lmmaxd, irmind:irmd, 2), intent(out) :: pzekdr
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: efac1, v1
    integer :: ir, j, l, l1, lm, lm1, m
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: real

    ! ---> set up array efac : efac(lm) = sqrt(e)**l/(2l - 1)!!

    efac(1) = cone
    v1 = cone
    do l = 1, lmaxd
      v1 = v1*ek/real(2*l-1, kind=dp)
      do m = -l, l
        lm = l*(l+1) + m + 1
        efac(lm) = v1
      end do
    end do

    ! ---> get wfts of same magnitude by scaling with efac

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

  end subroutine wftsca

end module mod_wftsca
