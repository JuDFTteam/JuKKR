module mod_calcrotmat
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine calcrotmat(nk, irel, alfdeg, betdeg, gamdeg, rot, fact, nkmmax)
    ! ********************************************************************
    ! *                                                                  *
    ! *   SETS UP THE ROTATION-MATRICES FOR THE EULER ANGLES             *
    ! *           ( ALFDEG, BETDEG, GAMDEG )                             *
    ! *                                                                  *
    ! *   SEE:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
    ! *            EQS. (4.8), (4.12) AND (4.13)                         *
    ! *                                                                  *
    ! *   for IREL=0,1   NK == NL           non-relativistic (l,m_l)     *
    ! *       IREL=3     NK == odd          relativistic (kappa,mue)     *
    ! *                                                                  *
    ! *   12/11/96  HE  deal with beta = 0                               *
    ! ********************************************************************

    use :: mod_errortrap
    implicit none

    complex (kind=dp) :: ci, c0
    parameter (ci=(0.0e0_dp,1.0e0_dp), c0=(0.0e0_dp,0.0e0_dp))
    real (kind=dp) :: pi
    parameter (pi=3.141592653589793238462643e0_dp)

    real (kind=dp) :: num, msb05, msb05sq, msb05pw, j, m1, m2, rfac, dom, x, cb05, cb05sq, alfdeg, betdeg, gamdeg, sum, cb05pw
    real (kind=dp) :: fact(0:100)

    integer :: s, slow, shigh, off, nkmmax, irel, nk, i1, i2, k, l, im1, im2, nmue
    complex (kind=dp) :: emim2a, emim1g, rot(nkmmax, nkmmax)

    ! INLINE FUNCTION    FACTORIAL FOR REAL ARGUMENT
    rfac(x) = fact(nint(x))

    if (irel==2) call errortrap('calcrotmat', 12, 1)
    if (irel==3 .and. mod(nk,2)==0) call errortrap('CALCROTMAT', 13, 1)

    do i2 = 1, nkmmax
      do i1 = 1, nkmmax
        rot(i1, i2) = c0
      end do
    end do

    cb05 = cos(betdeg*0.5e0_dp*pi/180.0e0_dp)
    cb05sq = cb05*cb05
    msb05 = -sin(betdeg*0.5e0_dp*pi/180.0e0_dp)
    msb05sq = msb05*msb05

    off = 0
    do k = 1, nk
      if (irel<2) then
        l = k - 1
        j = l
      else
        l = k/2
        if (l*2==k) then
          j = l - 0.5e0_dp
        else
          j = l + 0.5e0_dp
        end if
      end if

      nmue = nint(2*j+1)

      do im2 = 1, nmue
        m2 = -j + (im2-1.0e0_dp)
        emim2a = exp(-ci*m2*alfdeg*pi/180.0e0_dp)

        do im1 = 1, nmue
          m1 = -j + (im1-1.0e0_dp)
          emim1g = exp(-ci*m1*gamdeg*pi/180.0e0_dp)

          if (abs(betdeg)<1e-8_dp) then
            if (im1==im2) then
              sum = 1.0e0_dp
            else
              sum = 0.0e0_dp
            end if
          else
            slow = max(0, nint(m1-m2))
            shigh = min(nint(j-m2), nint(j+m1))
            cb05pw = cb05**nint(2*j+m1-m2-2*slow+2)
            msb05pw = msb05**nint(m2-m1+2*slow-2)
            dom = (-1.0e0_dp)**(slow-1)*sqrt(rfac(j+m1)*rfac(j-m1)*rfac(j+m2)*rfac(j-m2))
            sum = 0.0e0_dp

            do s = slow, shigh
              dom = -dom
              num = fact(s)*rfac(j-m2-s)*rfac(j+m1-s)*rfac(m2-m1+s)
              cb05pw = cb05pw/cb05sq
              msb05pw = msb05pw*msb05sq
              sum = sum + (dom/num)*cb05pw*msb05pw
            end do
          end if

          rot(off+im2, off+im1) = emim1g*sum*emim2a
        end do

      end do

      off = off + nmue
    end do

    return
  end subroutine calcrotmat

end module mod_calcrotmat
