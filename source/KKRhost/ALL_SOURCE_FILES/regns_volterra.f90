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
!     added Volterra equation by M. Ogura      Jan. 2006
!     FRED: true -> use fredholm equation
!           false -> volterra equation
!-----------------------------------------------------------------------
use mod_types, only: t_inc
implicit none
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
DOUBLE PRECISION ERR
INTEGER I,IR,IRC1,J,LM1,LM2,LM3
!..
!.. Local Arrays ..
DOUBLE COMPLEX PNS0(LMMAXD,LMMAXD,IRMIND:IRMD,2), &
               PNS1(LMMAXD,LMMAXD,IRMIND:IRMD)
INTEGER IPIV(LMMAXD)
!..
!.. External Subroutines ..
EXTERNAL CSINWD,CSOUT,WFINT,WFINT0,ZGEINV1
!..
!.. Parameters ..
DOUBLE COMPLEX CONE
PARAMETER (CONE= (1.0D0,0.0D0))
!..
LOGICAL FRED
DATA FRED/.false./
!..
!      write(*,*)ek
irc1 = ircut(ipan)
!      DO 1 J = 1,NSRA
!        DO 2 IR = IRMIND,IRC1
!          DO 3 LM1 = 1,LMMAXD
!            DO 4 LM2 = 1,LMMAXD
!              IF(LM1.EQ.LM2)THEN
!              PNS0(LM1,LM2,IR,J) =  AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)
!              ELSE
!              PNS0(LM1,LM2,IR,J) = (0D0,0D0)
!              ENDIF
! 4          CONTINUE
! 3        CONTINUE
! 2      CONTINUE
! 1    CONTINUE
IF(fred)THEN
  DO  i = 0,icst
!---> set up integrands for i-th born approximation
    IF (i == 0) THEN
      CALL wfint0(ader,bder,pzlm,qzekdr,pzekdr,vnspll,nsra,irmind,  &
          irmd,lmmaxd)
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
!-----------------------------------------------------------------------
! check convergence
    DO  j = 1,nsra
      DO  ir = irmind,irc1
        DO  lm1 = 1,lmmaxd
          DO  lm2 = 1,lmmaxd
            pns0(lm1,lm2,ir,j) = pns0(lm1,lm2,ir,j)-pns(lm1,lm2,ir,j)
          END DO
        END DO
      END DO
    END DO
    ERR=0D0
    DO  j=1,nsra
      CALL csout(pns0(1,1,irmind,j),pns1,lmmaxd**2,irmind,irmd,ipan, ircut)
      DO  lm1=1,lmmaxd
        DO  lm2=1,lmmaxd
          ERR=MAX(ERR,ABS(pns1(lm1,lm2,irc1)))
        END DO
      END DO
    END DO
    IF(t_inc%i_write>0) WRITE(1337,*) 'Born_Fred',i,ERR
!      IF(I.EQ.ICST.AND.ERR.GT.1D-3)WRITE(*,*)'NOT CONVERGENT',ERR
    DO  j = 1,nsra
      DO  ir = irmind,irc1
        DO  lm1 = 1,lmmaxd
          DO  lm2 = 1,lmmaxd
            pns0(lm1,lm2,ir,j) = pns(lm1,lm2,ir,j)
          END DO
        END DO
      END DO
    END DO
!-----------------------------------------------------------------------
  END DO
ELSE
!-----------------------------------------------------------------------
! Volterra equation
  DO  i = 0,icst
!---> set up integrands for i-th born approximation
    IF (i == 0) THEN
      CALL wfint0(ader,bder,pzlm,qzekdr,pzekdr,vnspll,nsra,irmind,  &
          irmd,lmmaxd)
    ELSE
      CALL wfint(pns,ader,bder,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
    END IF
!---> call integration subroutines
    CALL csout(ader,amat,lmmaxd**2,irmind,irmd,ipan,ircut)
    CALL csout(bder,bmat,lmmaxd**2,irmind,irmd,ipan,ircut)
    DO  ir = irmind,irc1
      DO  lm1 = 1,lmmaxd
        DO  lm2 = 1,lmmaxd
          IF(lm1 == lm2)THEN
            amat(lm1,lm2,ir) = cone - amat(lm1,lm2,ir)
          ELSE
            amat(lm1,lm2,ir) = - amat(lm1,lm2,ir)
          END IF
        END DO
      END DO
    END DO
!---> calculate non sph. wft. in i-th born approximation
    DO  j = 1,nsra
      DO  ir = irmind,irc1
        DO  lm1 = 1,lmmaxd
          DO  lm2 = 1,lmmaxd
            pns(lm1,lm2,ir,j) =  amat(lm1,lm2,ir)*pzlm(lm1,ir,j)  &
                +bmat(lm1,lm2,ir)*qzlm(lm1,ir,j)
          END DO
        END DO
      END DO
    END DO
!-----------------------------------------------------------------------
! check convergence
    DO  j = 1,nsra
      DO  ir = irmind,irc1
        DO  lm1 = 1,lmmaxd
          DO  lm2 = 1,lmmaxd
            pns0(lm1,lm2,ir,j) = pns0(lm1,lm2,ir,j)-pns(lm1,lm2,ir,j)
          END DO
        END DO
      END DO
    END DO
    ERR=0D0
    DO  j=1,nsra
      CALL csout(pns0(1,1,irmind,j),pns1,lmmaxd**2,irmind,irmd,ipan, ircut)
      DO  lm1=1,lmmaxd
        DO  lm2=1,lmmaxd
          ERR=MAX(ERR,ABS(pns1(lm1,lm2,irc1)))
        END DO
      END DO
    END DO
    IF(t_inc%i_write>0) WRITE(1337,*) 'Born',i,ERR
!      IF(I.EQ.ICST.AND.ERR.GT.1D-3)WRITE(*,*)'NOT CONVERGENT',ERR
    DO  j = 1,nsra
      DO  ir = irmind,irc1
        DO  lm1 = 1,lmmaxd
          DO  lm2 = 1,lmmaxd
            pns0(lm1,lm2,ir,j) = pns(lm1,lm2,ir,j)
          END DO
        END DO
      END DO
    END DO
!-----------------------------------------------------------------------
  END DO
  CALL zgeinv1(amat(1,1,irc1),ar,br,ipiv,lmmaxd)
  DO  lm1=1,lmmaxd
    DO  lm2=1,lmmaxd
      DO  ir=irmind,irc1
        ader(lm1,lm2,ir)=(0D0,0D0)
        bder(lm1,lm2,ir)=(0D0,0D0)
      END DO
      DO  lm3=1,lmmaxd
        DO  ir=irmind,irc1
          ader(lm1,lm2,ir)=ader(lm1,lm2,ir)+amat(lm1,lm3,ir)*ar(lm3,lm2)
          bder(lm1,lm2,ir)=bder(lm1,lm2,ir)+bmat(lm1,lm3,ir)*ar(lm3,lm2)
        END DO
      END DO
    END DO
  END DO
  DO  lm1=1,lmmaxd
    DO  lm2=1,lmmaxd
      DO  ir=irmind,irc1
        amat(lm1,lm2,ir)=ader(lm1,lm2,ir)
        bmat(lm1,lm2,ir)=bder(lm1,lm2,ir)
      END DO
    END DO
  END DO
  DO  j = 1,nsra
    DO  ir = irmind,irc1
      DO  lm1 = 1,lmmaxd
        DO  lm2 = 1,lmmaxd
          pns(lm1,lm2,ir,j) =  amat(lm1,lm2,ir)*pzlm(lm1,ir,j)  &
              +bmat(lm1,lm2,ir)*qzlm(lm1,ir,j)
        END DO
      END DO
    END DO
  END DO
! Volterra equation
!-----------------------------------------------------------------------
END IF
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
!      if(dreal(ek).gt.3.96d-2.and.dreal(ek).lt.3.97d-2)then
!      if(dimag(ek).gt.0.5018d0.and.dreal(ek).lt.0.5019d0)then
!      write(*,*)ek
!      do 250 lm1=5,7,2
!      write(*,*)'l=',lm1
!      do 250 ir=irmind,irc1
! 250  write(*,*)ir,dreal(pns(lm1,lm1,ir,1)),dimag(pns(lm1,lm1,ir,1))
!      endif
!      endif
END SUBROUTINE regns
! ************************************************************************

SUBROUTINE zgeinv1(a,u,aux,ipiv,dim)
! ************************************************************************
!   - inverts a general double complex matrix A,
!   - the result is return in U,
!   - input matrix A is returned unchanged,
!   - AUX is a auxiliary matrix,
!   - A,U and AUX are of dimension (DIM,DIM),
! ------------------------------------------------------------------------

REAL, INTENT(IN OUT)                     :: a
DOUBLE COMPLEX a(dim,*, INTENT(OUT)      :: u(dim,*)
DOUBLE COMPLEX a(dim,*, INTENT(IN OUT)   :: aux(dim,*)
INTEGER, INTENT(IN OUT)                  :: ipiv(*)
INTEGER, INTENT(IN)                      :: dim

DOUBLE COMPLEX a(dim,*)

!     .. PARAMETER

DOUBLE COMPLEX cone
REAL, PARAMETER :: cone=(1.d0,0.d0)

INTEGER :: lm1,info
EXTERNAL zcopy,zgetrs,zgetrf
! ------------------------------------------------------------------------
CALL cinit(dim*dim,u)
DO  lm1=1,dim
  u(lm1,lm1) = cone
END DO
END DO

CALL zcopy(dim*dim,a,1,aux,1)
CALL zgetrf(dim,dim,aux,dim,ipiv,info)
CALL zgetrs('N',dim,dim,aux,dim,ipiv,u,dim,info)

RETURN
END SUBROUTINE zgeinv1
