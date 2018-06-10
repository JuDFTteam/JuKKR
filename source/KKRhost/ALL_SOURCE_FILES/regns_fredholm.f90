SUBROUTINE regns(ar,br,efac,pns,vnspll,icst,ipan,ircut,pzlm,  &
        qzlm,pzekdr,qzekdr,ek,ader,amat,bder,bmat,nsra,  &
        irmind,irmd,ipand,lmmaxd)
!-----------------------------------------------------------------------
!     determines the regular non spherical wavefunctions , the
!       alpha matrix and the t - matrix in the n-th. born appro-
!       ximation ( n given by input parameter icst )


!     using the wave functions pz and qz ( regular and irregular
!       solution ) of the spherically averaged potential , the
!       regular wavefunction pns is determined by

!           pns(ir,lm1,lm2) = ar(ir,lm1,lm2)*pz(ir,l1)
!                                   + br(ir,lm1,lm2)*qz(ir,l1)

!      the matrices ar and br are determined by integral equations
!        containing pns and only the non spherical contributions of
!        the potential , stored in vinspll . these integral equations
!        are  solved iteratively with born approximation up to given n.

!     the original way of writing the cr and dr matrices in the equa-
!        tions above caused numerical troubles . therefore here are used
!        rescaled ar and br matrices :

!              ~
!              ar(ir,lm1,lm2) = sqrt(e)**(l1-l2)
!                             * ar(ir,lm1,lm2)*((2*l2-1)!!/(2*l1-1)!!)

!              ~
!              br(ir,lm1,lm2) = sqrt(e)**(-l1-l2)
!                             * br(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)

!     for lloyd's formular is only the determinant of the alpha -
!        matrix is needed which is identical with the determinant
!        of the rescaled ar - matrix at the innerst point .

!     the non spherical t - matrix is the br matrix at r(irc)

!     modified for the use of shape functions

!                              (see notes by b.drittler)

!                                b.drittler   mar.  1989
!-----------------------------------------------------------------------
!     modified by R. Zeller      Aug. 1994
!-----------------------------------------------------------------------
!.. Scalar Arguments ..
DOUBLE COMPLEX EK
INTEGER ICST,IPAN,IPAND,IRMD,IRMIND,LMMAXD,NSRA
!..
!.. Array Arguments ..
DOUBLE COMPLEX ADER(LMMAXD,LMMAXD,IRMIND:IRMD), &
               AMAT(LMMAXD,LMMAXD,IRMIND:IRMD),AR(LMMAXD,LMMAXD), &
               BDER(LMMAXD,LMMAXD,IRMIND:IRMD), &
               BMAT(LMMAXD,LMMAXD,IRMIND:IRMD),BR(LMMAXD,LMMAXD), &
               EFAC(*),PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2), &
               PZEKDR(LMMAXD,IRMIND:IRMD,2), &
               PZLM(LMMAXD,IRMIND:IRMD,2), &
               QZEKDR(LMMAXD,IRMIND:IRMD,2), &
               QZLM(LMMAXD,IRMIND:IRMD,2)
DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
INTEGER IRCUT(0:IPAND)
!..
!.. Local Scalars ..
DOUBLE COMPLEX EFAC1,EFAC2
INTEGER I,IR,IRC1,J,LM1,LM2
!..
!.. External Subroutines ..
EXTERNAL CSINWD,CSOUT,WFINT,WFINT0
!..
!.. Parameters ..
DOUBLE COMPLEX CONE
PARAMETER (CONE= (1.0D0,0.0D0))
!..
irc1 = ircut(ipan)
DO  i = 0,icst
!---> set up integrands for i-th born approximation
  IF (i == 0) THEN
    CALL wfint0(ader,bder,pzlm,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
  ELSE
    CALL wfint(pns,ader,bder,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
  END IF
!---> call integration subroutines
  CALL csinwd(ader,amat,lmmaxd**2,irmind,irmd,ipan,ircut)
  CALL csout(bder,bmat,lmmaxd**2,irmind,irmd,ipan,ircut)
  DO  ir = irmind,irc1
    DO  lm2 = 1,lmmaxd
      amat(lm2,lm2,ir) = cone + amat(lm2,lm2,ir)
    END DO
  END DO
!---> calculate non sph. wft. in i-th born approximation
  DO  j = 1,nsra
    DO  ir = irmind,irc1
      DO  lm1 = 1,lmmaxd
        DO  lm2 = 1,lmmaxd
          pns(lm1,lm2,ir,j) = (amat(lm1,lm2,ir)*pzlm(lm1,ir,j)+  &
              bmat(lm1,lm2,ir)*qzlm(lm1,ir,j))
        END DO
      END DO
    END DO
  END DO
END DO
DO  lm2 = 1,lmmaxd
  efac2 = efac(lm2)
!---> store alpha and t - matrix
  DO  lm1 = 1,lmmaxd
    efac1 = efac(lm1)
    ar(lm1,lm2) = amat(lm1,lm2,irmind)
!---> t-matrix
    br(lm1,lm2) = bmat(lm1,lm2,irc1)*efac1*efac2/ek
  END DO
END DO
!---> rescale with efac
DO  j = 1,nsra
  DO  lm2 = 1,lmmaxd
    efac2 = efac(lm2)
    DO  ir = irmind,irc1
      DO  lm1 = 1,lmmaxd
        pns(lm1,lm2,ir,j) = pns(lm1,lm2,ir,j)*efac2
      END DO
    END DO
  END DO
END DO
END SUBROUTINE regns
