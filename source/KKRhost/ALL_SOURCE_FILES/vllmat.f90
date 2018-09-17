module mod_vllmat
  use :: mod_datatypes, only: dp
  private :: dp

contains

  ! -------------------------------------------------------------------------------
  ! SUBROUTINE: VLLMAT
  ! > @brief
  ! -------------------------------------------------------------------------------
  subroutine vllmat(irmin, nrmaxd, irc, lmmax, lmmaxso, vnspll0, vins, lmpot, cleb, icleb, iend, nspin, zat, rnew, use_sratrick, ncleb)

    implicit none

    ! .. Input variables
    integer, intent (in) :: irc    ! < r point for potential cutting
    integer, intent (in) :: iend   ! < maximal number of non-vanishing Gaunt coefficients
    integer, intent (in) :: ncleb  ! < Number of Gaunt coefficients
    integer, intent (in) :: irmin  ! < max r for spherical treatment, afterwards non-spherical contribution
    integer, intent (in) :: lmmax  ! < (LMAX+1)^2
    integer, intent (in) :: nspin  ! < spin-degree of freedom
    integer, intent (in) :: lmpot  ! < (LPOT+1)**2
    integer, intent (in) :: nrmaxd ! < NTOTD*(NCHEBD+1), maximal number of radial points in Chebychev mesh
    integer, intent (in) :: lmmaxso ! < 2*(LMAX+1)^2, for SOC L=(l,m,s) instead of L=(l,m)
    integer, intent (in) :: use_sratrick ! < switch to use SRA trick (see routine rllsll) or not
    real (kind=dp), intent (in) :: zat ! < atomic charge
    integer, dimension (ncleb, 4), intent (in) :: icleb ! < index array for Gaunt coefficients
    real (kind=dp), dimension (ncleb), intent (in) :: cleb ! < GAUNT coefficients
    real (kind=dp), dimension (irmin:irc, lmpot, nspin), intent (in) :: vins ! < Non-spherical part of the potential
    real (kind=dp), dimension (irmin:nrmaxd), intent (in) :: rnew ! < radial mesh points of Chebychev mesh
    complex (kind=dp), dimension (lmmaxso, lmmaxso, irmin:irc), intent (out) :: vnspll0 ! < output potential in Chebychev mesh and (l,m,s)-space

    ! .. Local variables
    integer :: isp
    integer :: i, ir, j, lm1, lm2, lm3
    real (kind=dp), dimension (lmmax, lmmax, irmin:irc, nspin) :: vnspll

    do isp = 1, nspin
      do lm1 = 1, lmmax
        do lm2 = 1, lm1
          do ir = irmin, irc
            vnspll(lm1, lm2, ir, isp) = 0.0e0_dp
          end do                   ! IR
        end do                     ! LM2
      end do                       ! LM11

      do j = 1, iend
        lm1 = icleb(j, 1)
        lm2 = icleb(j, 2)
        lm3 = icleb(j, 3)
        do i = irmin, irc
          vnspll(lm1, lm2, i, isp) = vnspll(lm1, lm2, i, isp) + cleb(j)*vins(i, lm3, isp)
        end do                     ! I
      end do                       ! J
      ! -------------------------------------------------------------------------
      ! Use symmetry of the gaunt coef.
      ! -------------------------------------------------------------------------
      do lm1 = 1, lmmax
        do lm2 = 1, lm1 - 1
          do i = irmin, irc
            vnspll(lm2, lm1, i, isp) = vnspll(lm1, lm2, i, isp)
          end do                   ! I
        end do                     ! LM2
      end do                       ! LM1

      if (use_sratrick==0) then
        do lm1 = 1, lmmax
          do i = irmin, irc
            vnspll(lm1, lm1, i, isp) = vnspll(lm1, lm1, i, isp) + vins(i, 1, isp) - 2e0_dp*zat/rnew(i)
          end do
        end do
      end if

    end do                         ! NSPIN

    ! Set vnspll as twice as large
    vnspll0(1:lmmax, 1:lmmax, irmin:irc) = cmplx(vnspll(1:lmmax,1:lmmax,irmin:irc,1), 0e0_dp, kind=dp)

    if (nspin==2) then             ! hack to make routine work for Bxc-field
      vnspll0(lmmax+1:lmmaxso, lmmax+1:lmmaxso, irmin:irc) = cmplx(vnspll(1:lmmax,1:lmmax,irmin:irc,nspin), 0e0_dp, kind=dp)
    end if

  end subroutine vllmat

end module mod_vllmat
