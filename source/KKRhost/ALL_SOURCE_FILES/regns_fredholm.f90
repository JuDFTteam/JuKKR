subroutine regns(ar, br, efac, pns, vnspll, icst, ipan, ircut, pzlm, qzlm, &
  pzekdr, qzekdr, ek, ader, amat, bder, bmat, nsra, irmind, irmd, ipand, &
  lmmaxd)
  use :: mod_datatypes, only: dp
  ! -----------------------------------------------------------------------
  ! determines the regular non spherical wavefunctions , the
  ! alpha matrix and the t - matrix in the n-th. born appro-
  ! ximation ( n given by input parameter icst )


  ! using the wave functions pz and qz ( regular and irregular
  ! solution ) of the spherically averaged potential , the
  ! regular wavefunction pns is determined by

  ! pns(ir,lm1,lm2) = ar(ir,lm1,lm2)*pz(ir,l1)
  ! + br(ir,lm1,lm2)*qz(ir,l1)

  ! the matrices ar and br are determined by integral equations
  ! containing pns and only the non spherical contributions of
  ! the potential , stored in vinspll . these integral equations
  ! are  solved iteratively with born approximation up to given n.

  ! the original way of writing the cr and dr matrices in the equa-
  ! tions above caused numerical troubles . therefore here are used
  ! rescaled ar and br matrices :

  ! ~
  ! ar(ir,lm1,lm2) = sqrt(e)**(l1-l2)
  ! * ar(ir,lm1,lm2)*((2*l2-1)!!/(2*l1-1)!!)

  ! ~
  ! br(ir,lm1,lm2) = sqrt(e)**(-l1-l2)
  ! * br(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)

  ! for lloyd's formular is only the determinant of the alpha -
  ! matrix is needed which is identical with the determinant
  ! of the rescaled ar - matrix at the innerst point .

  ! the non spherical t - matrix is the br matrix at r(irc)

  ! modified for the use of shape functions

  ! (see notes by b.drittler)

  ! b.drittler   mar.  1989
  ! -----------------------------------------------------------------------
  ! modified by R. Zeller      Aug. 1994
  ! -----------------------------------------------------------------------
  ! .. Scalar Arguments ..
  complex (kind=dp) :: ek
  integer :: icst, ipan, ipand, irmd, irmind, lmmaxd, nsra
  ! ..
  ! .. Array Arguments ..
  complex (kind=dp) :: ader(lmmaxd, lmmaxd, irmind:irmd), &
    amat(lmmaxd, lmmaxd, irmind:irmd), ar(lmmaxd, lmmaxd), &
    bder(lmmaxd, lmmaxd, irmind:irmd), bmat(lmmaxd, lmmaxd, irmind:irmd), &
    br(lmmaxd, lmmaxd), efac(*), pns(lmmaxd, lmmaxd, irmind:irmd, 2), &
    pzekdr(lmmaxd, irmind:irmd, 2), pzlm(lmmaxd, irmind:irmd, 2), &
    qzekdr(lmmaxd, irmind:irmd, 2), qzlm(lmmaxd, irmind:irmd, 2)
  real (kind=dp) :: vnspll(lmmaxd, lmmaxd, irmind:irmd)
  integer :: ircut(0:ipand)
  ! ..
  ! .. Local Scalars ..
  complex (kind=dp) :: efac1, efac2
  integer :: i, ir, irc1, j, lm1, lm2
  ! ..
  ! .. External Subroutines ..
  external :: csinwd, csout, wfint, wfint0
  ! ..
  ! .. Parameters ..
  complex (kind=dp) :: cone
  parameter (cone=(1.0e0_dp,0.0e0_dp))
  ! ..
  irc1 = ircut(ipan)
  do i = 0, icst
    ! ---> set up integrands for i-th born approximation
    if (i==0) then
      call wfint0(ader, bder, pzlm, qzekdr, pzekdr, vnspll, nsra, irmind, &
        irmd, lmmaxd)
    else
      call wfint(pns, ader, bder, qzekdr, pzekdr, vnspll, nsra, irmind, irmd, &
        lmmaxd)
    end if
    ! ---> call integration subroutines
    call csinwd(ader, amat, lmmaxd**2, irmind, irmd, ipan, ircut)
    call csout(bder, bmat, lmmaxd**2, irmind, irmd, ipan, ircut)
    do ir = irmind, irc1
      do lm2 = 1, lmmaxd
        amat(lm2, lm2, ir) = cone + amat(lm2, lm2, ir)
      end do
    end do
    ! ---> calculate non sph. wft. in i-th born approximation
    do j = 1, nsra
      do ir = irmind, irc1
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            pns(lm1, lm2, ir, j) = (amat(lm1,lm2,ir)*pzlm(lm1,ir,j)+bmat(lm1, &
              lm2,ir)*qzlm(lm1,ir,j))
          end do
        end do
      end do
    end do
  end do
  do lm2 = 1, lmmaxd
    efac2 = efac(lm2)
    ! ---> store alpha and t - matrix
    do lm1 = 1, lmmaxd
      efac1 = efac(lm1)
      ar(lm1, lm2) = amat(lm1, lm2, irmind)
      ! ---> t-matrix
      br(lm1, lm2) = bmat(lm1, lm2, irc1)*efac1*efac2/ek
    end do
  end do
  ! ---> rescale with efac
  do j = 1, nsra
    do lm2 = 1, lmmaxd
      efac2 = efac(lm2)
      do ir = irmind, irc1
        do lm1 = 1, lmmaxd
          pns(lm1, lm2, ir, j) = pns(lm1, lm2, ir, j)*efac2
        end do
      end do
    end do
  end do
end subroutine regns
