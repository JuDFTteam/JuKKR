c 05.10.10 ***************************************************************
      SUBROUTINE NORMCOEFF_SO_SPINFLUX(IRCUT,
     +                   LMMAX,PNS,
     +                   KSRA,DRDI)
c ************************************************************************
c     Calculates the KKR matrix elements for the spin flux operator, i.e.,
c     
c           INT dr [R^{mu}_{Ls}]^dagger Q_s R^{mu}_{L's'}.
c
c     Details are in http://arxiv.org/pdf/1602.03417v1.pdf
c     
c     This subroutine was adapted from NORMCOEFF_SO.
c
c                                            Guillaume Géranton, 2016
c-----------------------------------------------------------------------
C     .. Parameters ..
      IMPLICIT NONE

      include 'inc.p'
      INTEGER          LMMAXD
      PARAMETER        (LMMAXD= (LMAXD+1)**2)
      INTEGER          LMPOTD
      PARAMETER        (LMPOTD= (LPOTD+1)**2)
      INTEGER          IRMIND,IRLMD
      PARAMETER        (IRMIND=IRMD-IRNSD,IRLMD= (IRNSD+1)*LMMAXD)
      INTEGER, PARAMETER :: NSPD=NSPIND
C     ..
C     .. Scalar Arguments ..
      INTEGER          LMMAX,KSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   PNS(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2,NATYPD)
      DOUBLE PRECISION DRDI(IRMD,NATYPD)
      INTEGER          IRCUT(0:IPAND,NATYPD)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   CZERO
      INTEGER          LM1,LM2,LM1P,
     +                 IR,I1,I1SP1,I1SP2,
     +                 LMSP1,LMSP2,ISIGMA,I2SP1,I2SP2,INSRA,NSRA
      INTEGER          LMMAXSO
C     ..
C     .. External Subroutines ..
      EXTERNAL         ZGEMM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        DATAN,DIMAG,DSQRT
C     ..
C     .. Save statement ..
      SAVE             CZERO
C     ..
C     ..Local Arrays..
      DOUBLE COMPLEX, ALLOCATABLE  ::   SPINFLUX(:,:,:,:),
     +                                  RLL_12(:),
     +                                  DENS(:,:,:,:,:,:,:),
     +                                  RLL(:,:,:,:,:,:,:)
      DOUBLE COMPLEX :: DELTA1, DELTA2
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
C     ..
c
      LMMAXSO=2*LMMAXD

      WRITE(6,*) "KSRA",KSRA
      IF (KSRA.GE.1) THEN    ! previously this was .GT. which is wrong for kvrel=1
         NSRA = 2
      ELSE
         NSRA = 1
      END IF

      WRITE(6,*) "NSRA",NSRA
c      WRITE(60,*)"ALAT=",ALAT
      WRITE(6,*) "LMMAX",LMMAX
      WRITE(6,*) "LMMAXD",LMMAXD
      WRITE(6,*) "LMMAXSO",LMMAXSO


      ALLOCATE(RLL(IRMD,LMMAX,LMMAX,2,2,2,NATYPD))
      ALLOCATE(RLL_12(LMMAX))
      ALLOCATE(DENS(LMMAXD,LMMAXD,2,2,2,2,NATYPD))
      ALLOCATE(SPINFLUX(LMMAXSO,LMMAXSO,NATYPD,3))

      RLL =CZERO
      DENS=CZERO

c    rewrite the wavefunctions in RLL arrays of 1,2*LMMAXD

      DO I1=1,NATYPD

        DO INSRA=1,NSRA
          DO IR=1,IRMD

            DO I1SP1=1,2
              DO I1SP2=1,2
                DO LM1 =1,LMMAXD
                  LMSP1=(I1SP1-1)*LMMAXD+LM1
                  DO LM2 =1,LMMAXD
                    LMSP2=(I1SP2-1)*LMMAXD+LM2
                    RLL(IR,LM2,LM1,I1SP2,I1SP1,INSRA,I1)=
     +                        PNS(LMSP2,LMSP1,IR,INSRA,I1)
                  END DO      !LM1=1,LMMAXD
                END DO      !LM1=1,LMMAXD
              END DO      !ISP1=1,2
            END DO      !ISP1=1,2

          END DO      !IR
        END DO      !INSRA


c set up the array R*_L1L2 R_L3L4 
        IF (I1.EQ.4) write(55,99)DRDI(:,I1)
        IF (I1.EQ.4) write(56,99)REAL(RLL(:,4,6,1,
     +                          1,1,I1))
        IF (I1.EQ.4) write(57,99)AIMAG(RLL(:,4,6,1,
     +                          1,1,I1))
        IF (I1.EQ.4) write(58,99)REAL(RLL(:,4,6,1,
     +                          1,2,I1))
        IF (I1.EQ.4) write(59,99)AIMAG(RLL(:,4,6,1,
     +                          1,2,I1))
                          
 99     FORMAT(1d20.10)

        DO I1SP1=1,2
          DO I1SP2=1,2
            DO I2SP1=1,2
              DO I2SP2=1,2

c                WRITE(6,*) "I1SP1,I1SP2,I2SP1,I2SP2"
c                WRITE(6,*) I1SP1,I1SP2,I2SP1,I2SP2

                DO LM1 = 1,LMMAX
                  DO LM2 = 1,LMMAX

                    DO INSRA=1,NSRA

                      RLL_12=CZERO

                      DO LM1P = 1,LMMAX

                        DELTA1=(RLL(IRCUT(1,I1),LM1P,LM2,I2SP1,
     +                          I2SP2,INSRA,I1)-RLL(IRCUT(1,I1)-1,
     +                          LM1P,LM2,I2SP1,I2SP2,INSRA,I1))
     +                          /DRDI(IRCUT(1,I1),I1)
                        DELTA2=(RLL(IRCUT(1,I1),LM1P,LM1,I1SP1,
     +                          I1SP2,INSRA,I1)-RLL(IRCUT(1,I1)-1,
     +                          LM1P,LM1,I1SP1,I1SP2,INSRA,I1))
     +                          /DRDI(IRCUT(1,I1),I1)

                        RLL_12(LM1P)=
     +                  DCONJG(RLL(IRCUT(1,I1)-1,LM1P,LM1,I1SP1,I1SP2,
     +                  INSRA,I1))*DELTA1-
     +                  RLL(IRCUT(1,I1)-1,LM1P,LM2,I2SP1,I2SP2,INSRA,I1)
     +                  *DCONJG(DELTA2)

                      END DO!LM1P

!                      DO LM1P = 1,LMMAX
!
!                        DELTA1=(RLL(361,LM1P,LM2,I2SP1,
!     +                          I2SP2,INSRA,I1)-RLL(360,
!     +                          LM1P,LM2,I2SP1,I2SP2,INSRA,I1))
!     +                          /(DRDI(361,I1))
!                        DELTA2=(RLL(361,LM1P,LM1,I1SP1,
!     +                          I1SP2,INSRA,I1)-RLL(360,
!     +                          LM1P,LM1,I1SP1,I1SP2,INSRA,I1))
!     +                          /(DRDI(361,I1))
!
!                        RLL_12(LM1P)=
!     +                  DCONJG(RLL(360,LM1P,LM1,I1SP1,I1SP2,
!     +                  INSRA,I1))*DELTA1-
!     +                  RLL(360,LM1P,LM2,I2SP1,I2SP2,INSRA,I1)*
!     +                  DCONJG(DELTA2)
!
!                      END DO!LM1P

                      DO LM1P = 1,LMMAX
                        DENS(LM1,LM2,I1SP1,I1SP2,I2SP1,I2SP2,I1) =
     +                  DENS(LM1,LM2,I1SP1,I1SP2,I2SP1,I2SP2,I1)
     +                  + RLL_12(LM1P)
                      END DO!LM1P

                    END DO   !NSRA

                  END DO         !LM2
                END DO           !LM1

              END DO            !I2SP2
            END DO              !I2SP1
          END DO         !I1SP2
        END DO           !I1SP1 

      END DO             !I1


      SPINFLUX=CZERO

      DO ISIGMA=1,3  !ISIGMA == 1 --> Q_x
                     !ISIGMA == 2 --> Q_y
                     !ISIGMA == 3 --> Q_z

        WRITE(6,*) "ISIGMA",ISIGMA
        DO I1=1,NATYPD

          IF (ISIGMA==1) THEN  !Q_x 

            DO I1SP1=1,2
              DO I1SP2=1,2
                DO LM1=1,LMMAX
                  DO LM2=1,LMMAX
            SPINFLUX((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +               -(0d0,1d0)*(DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)+
     +               DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1))/2
                  END DO !LM2
                END DO !LM1
              END DO !I1SP2
            END DO !I1SP1

          ELSE IF (ISIGMA==2) THEN !Q_y

            DO I1SP1=1,2
              DO I1SP2=1,2
                DO LM1=1,LMMAX
                  DO LM2=1,LMMAX
            SPINFLUX((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +      -(0d0,1d0)*(-1)*(0d0,1d0)*(DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)-
     +      DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1))/2
                  END DO !LM2
                END DO !LM1
              END DO !I1SP2
            END DO !I1SP1

          ELSE IF (ISIGMA==3) THEN !Q_z

            DO I1SP1=1,2
              DO I1SP2=1,2
                DO LM1=1,LMMAX
                  DO LM2=1,LMMAX
            SPINFLUX((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +               (0d0,1d0)*(DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)-
     +               DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1))/2
                  END DO !LM2
                END DO !LM1
              END DO !I1SP2
            END DO !I1SP1

          END IF

        END DO  !I1
      END DO    !ISIGMA

      open(unit=12, file='TBkkr_spinflux.txt', form='formatted',
     +     action='write')
      DO ISIGMA=1,3
        DO I1=1,NATYPD
          DO LM2=1,LMMAXSO
            DO LM1=1,LMMAXSO
c                IF(LM1.lt.LM2) THEN
c minus sign to get the spin flux into the sphere :
                   WRITE(12,'(2ES25.16)') -SPINFLUX(LM1,LM2,I1,ISIGMA)
c                   WRITE(12,'(2ES25.16)') SPINFLUX(LM2,LM1,I1,ISIGMA)
c                END IF
            END DO
          END DO
        END DO
      END DO
      close(12)

      DEALLOCATE(RLL)
      DEALLOCATE(DENS)
      DEALLOCATE(RLL_12)
      DEALLOCATE(SPINFLUX)

      END




