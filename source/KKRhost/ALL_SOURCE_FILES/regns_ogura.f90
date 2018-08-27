module mod_regns_ogura

contains

subroutine regns(ar, br, efac, pns, vnspll, icst, ipan, ircut, pzlm, qzlm, &
  pzekdr, qzekdr, ek, ader, amat, bder, bmat, nsra, irmind, irmd, ipand, &
  lmmaxd)
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
  ! added Volterra equation by M. Ogura      Jan. 2006
  ! FRED: true -> use fredholm equation
  ! false -> volterra equation
  ! -----------------------------------------------------------------------
  use :: mod_types, only: t_inc
  use :: mod_datatypes, only: dp
  use mod_regns, only: zgeinv1
  use mod_csout
  use mod_wfint0
  use mod_wfint
  implicit none
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
  real (kind=dp) :: err
  integer :: i, ir, irc1, j, lm1, lm2, lm3
  ! ..
  ! .. Local Arrays ..
  complex (kind=dp) :: pns0(lmmaxd, lmmaxd, irmind:irmd, 2), &
    pns1(lmmaxd, lmmaxd, irmind:irmd)
  integer :: ipiv(lmmaxd)
  ! ..
  ! .. Parameters ..
  complex (kind=dp) :: cone
  parameter (cone=(1.0d0,0.0d0))
  ! ..
  logical :: fred
  data fred/.false./
  ! ..
  ! write(*,*)ek
  irc1 = ircut(ipan)
  ! DO 1 J = 1,NSRA
  ! DO 2 IR = IRMIND,IRC1
  ! DO 3 LM1 = 1,LMMAXD
  ! DO 4 LM2 = 1,LMMAXD
  ! IF(LM1.EQ.LM2)THEN
  ! PNS0(LM1,LM2,IR,J) =  AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)
  ! ELSE
  ! PNS0(LM1,LM2,IR,J) = (0D0,0D0)
  ! ENDIF
  ! 4          CONTINUE
  ! 3        CONTINUE
  ! 2      CONTINUE
  ! 1    CONTINUE
  if (fred) then
    do i = 0, icst
      ! ---> set up integrands for i-th born approximation
      if (i==0) then
        call wfint0(ader, bder, pzlm, qzekdr, pzekdr, vnspll, nsra, irmind, &
          irmd, lmmaxd)
      else
        call wfint(pns, ader, bder, qzekdr, pzekdr, vnspll, nsra, irmind, &
          irmd, lmmaxd)
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
              pns(lm1, lm2, ir, j) = (amat(lm1,lm2,ir)*pzlm(lm1,ir,j)+bmat(lm1 &
                ,lm2,ir)*qzlm(lm1,ir,j))
            end do
          end do
        end do
      end do
      ! -----------------------------------------------------------------------
      ! check convergence
      do j = 1, nsra
        do ir = irmind, irc1
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              pns0(lm1, lm2, ir, j) = pns0(lm1, lm2, ir, j) - &
                pns(lm1, lm2, ir, j)
            end do
          end do
        end do
      end do
      err = 0d0
      do j = 1, nsra
        call csout(pns0(1,1,irmind,j), pns1, lmmaxd**2, irmind, irmd, ipan, &
          ircut)
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            err = max(err, abs(pns1(lm1,lm2,irc1)))
          end do
        end do
      end do
      if (t_inc%i_write>0) write (1337, *) 'Born_Fred', i, err
      ! IF(I.EQ.ICST.AND.ERR.GT.1D-3)WRITE(*,*)'NOT CONVERGENT',ERR
      do j = 1, nsra
        do ir = irmind, irc1
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              pns0(lm1, lm2, ir, j) = pns(lm1, lm2, ir, j)
            end do
          end do
        end do
      end do
      ! -----------------------------------------------------------------------
    end do
  else
    ! -----------------------------------------------------------------------
    ! Volterra equation
    do i = 0, icst
      ! ---> set up integrands for i-th born approximation
      if (i==0) then
        call wfint0(ader, bder, pzlm, qzekdr, pzekdr, vnspll, nsra, irmind, &
          irmd, lmmaxd)
      else
        call wfint(pns, ader, bder, qzekdr, pzekdr, vnspll, nsra, irmind, &
          irmd, lmmaxd)
      end if
      ! ---> call integration subroutines
      call csout(ader, amat, lmmaxd**2, irmind, irmd, ipan, ircut)
      call csout(bder, bmat, lmmaxd**2, irmind, irmd, ipan, ircut)
      do ir = irmind, irc1
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            if (lm1==lm2) then
              amat(lm1, lm2, ir) = cone - amat(lm1, lm2, ir)
            else
              amat(lm1, lm2, ir) = -amat(lm1, lm2, ir)
            end if
          end do
        end do
      end do
      ! ---> calculate non sph. wft. in i-th born approximation
      do j = 1, nsra
        do ir = irmind, irc1
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              pns(lm1, lm2, ir, j) = amat(lm1, lm2, ir)*pzlm(lm1, ir, j) + &
                bmat(lm1, lm2, ir)*qzlm(lm1, ir, j)
            end do
          end do
        end do
      end do
      ! -----------------------------------------------------------------------
      ! check convergence
      do j = 1, nsra
        do ir = irmind, irc1
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              pns0(lm1, lm2, ir, j) = pns0(lm1, lm2, ir, j) - &
                pns(lm1, lm2, ir, j)
            end do
          end do
        end do
      end do
      err = 0d0
      do j = 1, nsra
        call csout(pns0(1,1,irmind,j), pns1, lmmaxd**2, irmind, irmd, ipan, &
          ircut)
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            err = max(err, abs(pns1(lm1,lm2,irc1)))
          end do
        end do
      end do
      if (t_inc%i_write>0) write (1337, *) 'Born', i, err
      ! IF(I.EQ.ICST.AND.ERR.GT.1D-3)WRITE(*,*)'NOT CONVERGENT',ERR
      do j = 1, nsra
        do ir = irmind, irc1
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              pns0(lm1, lm2, ir, j) = pns(lm1, lm2, ir, j)
            end do
          end do
        end do
      end do
      ! -----------------------------------------------------------------------
    end do
    call zgeinv1(amat(1,1,irc1), ar, br, ipiv, lmmaxd)
    do lm1 = 1, lmmaxd
      do lm2 = 1, lmmaxd
        do ir = irmind, irc1
          ader(lm1, lm2, ir) = (0d0, 0d0)
          bder(lm1, lm2, ir) = (0d0, 0d0)
        end do
        do lm3 = 1, lmmaxd
          do ir = irmind, irc1
            ader(lm1, lm2, ir) = ader(lm1, lm2, ir) + &
              amat(lm1, lm3, ir)*ar(lm3, lm2)
            bder(lm1, lm2, ir) = bder(lm1, lm2, ir) + &
              bmat(lm1, lm3, ir)*ar(lm3, lm2)
          end do
        end do
      end do
    end do
    do lm1 = 1, lmmaxd
      do lm2 = 1, lmmaxd
        do ir = irmind, irc1
          amat(lm1, lm2, ir) = ader(lm1, lm2, ir)
          bmat(lm1, lm2, ir) = bder(lm1, lm2, ir)
        end do
      end do
    end do
    do j = 1, nsra
      do ir = irmind, irc1
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            pns(lm1, lm2, ir, j) = amat(lm1, lm2, ir)*pzlm(lm1, ir, j) + &
              bmat(lm1, lm2, ir)*qzlm(lm1, ir, j)
          end do
        end do
      end do
    end do
    ! Volterra equation
    ! -----------------------------------------------------------------------
  end if
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
  ! if(dreal(ek).gt.3.96d-2.and.dreal(ek).lt.3.97d-2)then
  ! if(dimag(ek).gt.0.5018d0.and.dreal(ek).lt.0.5019d0)then
  ! write(*,*)ek
  ! do 250 lm1=5,7,2
  ! write(*,*)'l=',lm1
  ! do 250 ir=irmind,irc1
  ! 250  write(*,*)ir,dreal(pns(lm1,lm1,ir,1)),dimag(pns(lm1,lm1,ir,1))
  ! endif
  ! endif
end subroutine regns
! ************************************************************************

end module mod_regns_ogura
