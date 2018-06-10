SUBROUTINE irwns(cr,dr,efac,qns,vnspll,icst,ipan,ircut,nsra,  &
        pzlm,qzlm,pzekdr,qzekdr,cder,cmat,dder,dmat,  &
        irmind,irmd,irmin,irmax,ipand,lmmaxd)       ! Added IRMIN,IRMAX 1.7.2014
!-----------------------------------------------------------------------
!     determines the irregular non spherical wavefunctions in the n-th.
!       born approximation ( n given by input parameter icst ) .


!     using the wave functions pz and qz ( regular and irregular
!       solution ) of the spherically averaged potential , the ir-
!       regular wavefunction qns is determined by

!           qns(ir,lm1,lm2) = cr(ir,lm1,lm2)*pz(ir,l1)

!                                   + dr(ir,lm1,lm2)*qz(ir,l1)

!      the matrices cr and dr are determined by integral equations
!        containing qns and only the non spherical contributions of
!        the potential , stored in vinspll . these integral equations
!        are solved iteratively with born approximation up to given n.

!     the original way of writing the cr and dr matrices in the equa-
!        tion above caused numerical troubles . therefore here are used
!        rescaled cr and dr matrices (compare subroutine wftsca):

!              ~
!              cr(ir,lm1,lm2) = sqrt(e)**(l1+l2)
!                             * cr(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)

!              ~
!              dr(ir,lm1,lm2) = sqrt(e)**(l2-l1)
!                             * dr(ir,lm1,lm2)*((2*l1-1)!!/(2*l2-1)!!)

!     attention :  the sign of the dr matrix is changed to reduce the
!     ===========  number of floating point operations

!     modified for the use of shape functions

!                              (see notes by b.drittler)

!                                b.drittler   mar.  1989
!-----------------------------------------------------------------------
!     modified by R. Zeller      Aug. 1994
!-----------------------------------------------------------------------
implicit none
!.. Parameters ..
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.D0,0.D0))
!..
!.. Scalar Arguments ..
      INTEGER ICST,IPAN,IPAND,IRMD,IRMIND,LMMAXD,NSRA,IRMIN,IRMAX
!..
!.. Array Arguments ..
DOUBLE COMPLEX CDER(LMMAXD,LMMAXD,IRMIND:IRMD), &
               CMAT(LMMAXD,LMMAXD,IRMIND:IRMD),CR(LMMAXD,LMMAXD), &
               DDER(LMMAXD,LMMAXD,IRMIND:IRMD), &
               DMAT(LMMAXD,LMMAXD,IRMIND:IRMD),DR(LMMAXD,LMMAXD), &
               EFAC(LMMAXD),PZEKDR(LMMAXD,IRMIND:IRMD,2), &
               PZLM(LMMAXD,IRMIND:IRMD,2), &
               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2), &
               QZEKDR(LMMAXD,IRMIND:IRMD,2), &
               QZLM(LMMAXD,IRMIND:IRMD,2)
DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
INTEGER IRCUT(0:IPAND)
!..
!.. Local Scalars ..
      DOUBLE COMPLEX EFAC2
      INTEGER I,IR,IRC1,J,LM1,LM2
!..
!.. External Subroutines ..
      EXTERNAL CSINWD,WFINT,WFINT0
!..
irc1 = ircut(ipan)
DO  i = 0,icst
!---> set up integrands for i-th born approximation
  IF (i == 0) THEN
    CALL wfint0(cder,dder,qzlm,qzekdr,pzekdr,vnspll,nsra,irmind,  &
        irmd,lmmaxd,irmin,irmax)                         ! Added IRMIN,IRMAX 1.7.2014
  ELSE
    CALL wfint(qns,cder,dder,qzekdr,pzekdr,vnspll,nsra,irmind,  &
        irmd,lmmaxd,irmin,irmax)                          ! Added IRMIN,IRMAX 1.7.2014
  END IF
!---> call integration subroutines
  CALL csinwd(cder,cmat,lmmaxd**2,irmind,irmd,irmin,ipan,ircut)     ! Added IRMIN 1.7.2014
  CALL csinwd(dder,dmat,lmmaxd**2,irmind,irmd,irmin,ipan,ircut)     ! Added IRMIN 1.7.2014
  DO  ir = irmin,irc1
    DO  lm2 = 1,lmmaxd
      dmat(lm2,lm2,ir) = dmat(lm2,lm2,ir) - cone
    END DO
  END DO
!---> calculate non sph. wft. in i-th born approximation
  DO  j = 1,nsra
    DO  ir = irmin,irc1
      DO  lm1 = 1,lmmaxd
        DO  lm2 = 1,lmmaxd
          qns(lm1,lm2,ir,j) = cmat(lm1,lm2,ir)*pzlm(lm1,ir,j) -  &
              dmat(lm1,lm2,ir)*qzlm(lm1,ir,j)
        END DO
      END DO
    END DO
  END DO
END DO
DO  lm2 = 1,lmmaxd
!---> store c - and d - matrix
  DO  lm1 = 1,lmmaxd
    cr(lm1,lm2) = cmat(lm1,lm2,irmin)
    dr(lm1,lm2) = -dmat(lm1,lm2,irmin)
  END DO
END DO
!---> rescale with efac
DO  j = 1,nsra
  DO  lm2 = 1,lmmaxd
    efac2 = 1.d0/efac(lm2)
    DO  ir = irmin,irc1
      DO  lm1 = 1,lmmaxd
        qns(lm1,lm2,ir,j) = qns(lm1,lm2,ir,j)*efac2
      END DO
    END DO
  END DO
END DO
END SUBROUTINE irwns
