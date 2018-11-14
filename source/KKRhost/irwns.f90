!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_irwns
  
  private
  public :: irwns

contains

  !-------------------------------------------------------------------------------
  !> Summary: Determines the irregular non spherical wavefunctions in the n-th.
  !> Born approximation (n given by input parameter icst).
  !> Author: B. Drittler, R. Zeller
  !> Date: Mar. 1989
  !> Category: KKRhost, single-site
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Using the wave functions pz and qz ( regular and irregular
  !> solution ) of the spherically averaged potential, the ir-
  !> regular wavefunction qns is determined by
  !>
  !> qns(ir,lm1,lm2) = cr(ir,lm1,lm2)*pz(ir,l1) + dr(ir,lm1,lm2)*qz(ir,l1)
  !>
  !> the matrices cr and dr are determined by integral equations
  !> containing qns and only the non spherical contributions of
  !> the potential, stored in vinspll. These integral equations
  !> are solved iteratively with Born approximation up to given n.
  !>
  !> The original way of writing the cr and dr matrices in the equa-
  !> tion above caused numerical troubles. Therefore here are used
  !> rescaled cr and dr matrices (compare subroutine wftsca):
  !>
  !> ~
  !> cr(ir,lm1,lm2) = sqrt(e)**(l1+l2) * cr(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)
  !>
  !> ~
  !> dr(ir,lm1,lm2) = sqrt(e)**(l2-l1) * dr(ir,lm1,lm2)*((2*l1-1)!!/(2*l2-1)!!)
  !>
  !> Attention :  the sign of the dr matrix is changed to reduce the
  !> ===========  number of floating point operations
  !>
  !> Modified for the use of shape functions (see notes by B. Drittler)
  !> B. Drittler   Mar.  1989
  !> 
  !> modified by R. Zeller      Aug. 1994
  !-------------------------------------------------------------------------------
  subroutine irwns(cr, dr, efac, qns, vnspll, icst, ipan, ircut, nsra, pzlm, qzlm, pzekdr, qzekdr, cder, cmat, dder, dmat, irmind, irmd, irmin, irmax, ipand, lmmaxd) ! Added IRMIN,IRMAX 1.7.2014

    use :: mod_datatypes, only: dp
    use :: mod_wfint, only: wfint, wfint0
    use :: mod_csinwd, only: csinwd
    use :: mod_constants, only: cone
    implicit none
    ! ..
    ! .. Scalar Arguments ..
    integer :: icst, ipan, ipand, irmd, irmind, lmmaxd, nsra, irmin, irmax
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: cder(lmmaxd, lmmaxd, irmind:irmd), cmat(lmmaxd, lmmaxd, irmind:irmd), cr(lmmaxd, lmmaxd), dder(lmmaxd, lmmaxd, irmind:irmd), &
      dmat(lmmaxd, lmmaxd, irmind:irmd), dr(lmmaxd, lmmaxd), efac(lmmaxd), pzekdr(lmmaxd, irmind:irmd, 2), pzlm(lmmaxd, irmind:irmd, 2), qns(lmmaxd, lmmaxd, irmind:irmd, 2), &
      qzekdr(lmmaxd, irmind:irmd, 2), qzlm(lmmaxd, irmind:irmd, 2)
    real (kind=dp) :: vnspll(lmmaxd, lmmaxd, irmind:irmd)
    integer :: ircut(0:ipand)
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: efac2
    integer :: i, ir, irc1, j, lm1, lm2


    irc1 = ircut(ipan)
    do i = 0, icst
      ! ---> set up integrands for i-th born approximation
      if (i==0) then
        call wfint0(cder, dder, qzlm, qzekdr, pzekdr, vnspll, nsra, irmind, irmd, lmmaxd, irmin, irmax) ! Added IRMIN,IRMAX 1.7.2014
      else
        call wfint(qns, cder, dder, qzekdr, pzekdr, vnspll, nsra, irmind, irmd, lmmaxd, irmin, irmax) ! Added IRMIN,IRMAX 1.7.2014
      end if
      ! ---> call integration subroutines
      call csinwd(cder, cmat, lmmaxd**2, irmind, irmd, irmin, ipan, ircut) ! Added IRMIN 1.7.2014
      call csinwd(dder, dmat, lmmaxd**2, irmind, irmd, irmin, ipan, ircut) ! Added IRMIN 1.7.2014
      do ir = irmin, irc1
        do lm2 = 1, lmmaxd
          dmat(lm2, lm2, ir) = dmat(lm2, lm2, ir) - cone
        end do
      end do
      ! ---> calculate non sph. wft. in i-th born approximation
      do j = 1, nsra
        do ir = irmin, irc1
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              qns(lm1, lm2, ir, j) = cmat(lm1, lm2, ir)*pzlm(lm1, ir, j) - dmat(lm1, lm2, ir)*qzlm(lm1, ir, j)
            end do
          end do
        end do
      end do
    end do
    do lm2 = 1, lmmaxd
      ! ---> store c - and d - matrix
      do lm1 = 1, lmmaxd
        cr(lm1, lm2) = cmat(lm1, lm2, irmin)
        dr(lm1, lm2) = -dmat(lm1, lm2, irmin)
      end do
    end do
    ! ---> rescale with efac
    do j = 1, nsra
      do lm2 = 1, lmmaxd
        efac2 = 1.e0_dp/efac(lm2)
        do ir = irmin, irc1
          do lm1 = 1, lmmaxd
            qns(lm1, lm2, ir, j) = qns(lm1, lm2, ir, j)*efac2
          end do
        end do
      end do
    end do
  end subroutine irwns

end module mod_irwns
