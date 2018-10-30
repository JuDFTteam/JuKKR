c 12.05.2016 *************************************************************
      SUBROUTINE NORMCOEFF_SO_TORQ(IRMINSO,IRCUT,LMAX,
     +                   LMMAX,PNS,THETAS,NTCELL,
     +                   IFUNM,IPAN,LMSP,KSRA,CLEB,ICLEB,IEND,DRDI,
     +                   IRWS,VISP,NSPIN,VINS,IRMIN)
c ************************************************************************
c     Calculates the KKR matrix elements for the torque operator, i.e.,
c     
c           INT dr [R^{mu}_{Ls}]^dagger T^{mu}(r) R^{mu}_{L's'}.
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
C     ..
C     .. Scalar Arguments ..
      INTEGER          IEND,LMAX,LMMAX,KSRA,IRWS(*),NSPIN,IRMIN(*)
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   PNS(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2,NATYPD)
      DOUBLE PRECISION CLEB(*),
     +                 THETAS(IRID,NFUND,*),
     +                 DRDI(IRMD,NATYPD),
     +                 VISP(IRMD,*),              ! spherical part of the potential
     +                 VINS(IRMIND:IRMD,LMPOTD,*), ! non-spher. part of the potential
     +                 THETA,PHI,THETA_TMP,PHI_TMP,
     +                 SQA(3)
      INTEGER          ICLEB(NCLEB,4),IFUNM(NATYPD,LMPOTD),
     +                 LMSP(NATYPD,*),IRMINSO,IRCUT(0:IPAND,NATYPD),
     +                 IPAN(NATYPD),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   CONE,CZERO,NORM
      INTEGER          LM1,LM2,LM1P,LM2P,
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
      SAVE             CZERO,CONE
C     ..
C     ..Local Arrays..
      DOUBLE COMPLEX, ALLOCATABLE  ::   TORQ(:,:,:,:),
     +                                  RLL_12(:,:,:),
     +                                  DENS(:,:,:,:,:,:,:),
     +                                  RLL(:,:,:,:,:,:,:)
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA CONE/ (1.0D0,0.0D0)/
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

      WRITE(6,*) "LMMAX",LMMAX
      WRITE(6,*) "LMMAXD",LMMAXD
      WRITE(6,*) "LMMAXSO",LMMAXSO


      ALLOCATE(RLL(IRMD,LMMAX,LMMAX,2,2,2,NATYPD))

      ALLOCATE(RLL_12(IRMD,LMMAX,LMMAX))
      ALLOCATE(DENS(LMMAXD,LMMAXD,2,2,2,2,NATYPD))
      ALLOCATE(TORQ(LMMAXSO,LMMAXSO,NATYPD,3))

      RLL =CZERO
      DENS=CZERO

c    rewrite the wavefunctions in RLL arrays of 1,2*LMMAXD

      DO I1=1,NATYPD

        DO INSRA=1,NSRA
          DO IR=IRMINSO,IRMD

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

        DO I1SP1=1,2
          DO I1SP2=1,2
            DO I2SP1=1,2
              DO I2SP2=1,2

                DO LM1 = 1,LMMAX
                  DO LM2 = 1,LMMAX

                    DO INSRA=1,NSRA

                      RLL_12=CZERO

                      DO LM1P = 1,LMMAX
                        DO LM2P = 1,LMMAX

                          DO IR=1,IRMD
                              RLL_12(IR,LM1P,LM2P)=
     +                    DCONJG(RLL(IR,LM1P,LM1,I1SP1,I1SP2,INSRA,I1))*
     +                           RLL(IR,LM2P,LM2,I2SP1,I2SP2,INSRA,I1)
                          END DO       !IR

                        END DO         !LM2P
                      END DO           !LM1P

                       CALL CALC_TORQ_LL_SS(LMAX,LMMAX,RLL_12,
     +                     IRCUT(0:IPAND,I1),IPAN(I1),NTCELL(I1),THETAS,
     +                     CLEB,ICLEB,IEND,IFUNM,LMSP,IRMINSO,IRWS(I1),
     +                    DRDI(:,I1),NORM,VISP,NSPIN,I1,VINS,IRMIN(I1))

                       DENS(LM1,LM2,I1SP1,I1SP2,I2SP1,I2SP2,I1) =
     +                 DENS(LM1,LM2,I1SP1,I1SP2,I2SP1,I2SP2,I1) + NORM
                    END DO   !NSRA

                  END DO         !LM2
                END DO           !LM1

              END DO            !I2SP2
            END DO              !I2SP1
          END DO         !I1SP2
        END DO           !I1SP1 

      END DO             !I1


c reads spin quantization axis from file
      OPEN(UNIT=11,FILE='nonco_angle.dat',FORM='FORMATTED')

      DO I1=1,NATYPD 
        READ(11,*) THETA_TMP,PHI_TMP
        IF ( I1 > 1) THEN
          IF ((THETA_TMP .ne. THETA ) .or. (PHI_TMP .ne. PHI))
     +      stop "It seems you want to compute the torque 
     +      in a non-colinear system. This is not implemented yet."
        END IF        
        THETA = THETA_TMP
        PHI = PHI_TMP
      END DO
      CLOSE(11)

      SQA(1) = sin(THETA)*cos(PHI)
      SQA(2) = sin(THETA)*sin(PHI)
      SQA(3) = cos(THETA)      

      TORQ=CZERO

      DO ISIGMA=1,3  !ISIGMA == 1 --> T_x
                     !ISIGMA == 2 --> T_y
                     !ISIGMA == 3 --> T_z

        WRITE(6,*) "ISIGMA",ISIGMA
        DO I1=1,NATYPD

c       TEMPORARY IMPLEMENTATION OF THE TORQUE OPERATOR FOR B // z
          IF (ISIGMA==1) THEN  !T_x 

            DO I1SP1=1,2
              DO I1SP2=1,2
                DO LM1=1,LMMAX
                  DO LM2=1,LMMAX
               TORQ((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +               (0d0,1d0)*(DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1)-
     +                          DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1))*SQA(3)
     +                       -(-DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)+
     +                          DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1))*SQA(2)
                  END DO !LM2
                END DO !LM1
              END DO !I1SP2
            END DO !I1SP1

          ELSE IF (ISIGMA==2) THEN !T_y

            DO I1SP1=1,2
              DO I1SP2=1,2
                DO LM1=1,LMMAX
                  DO LM2=1,LMMAX
               TORQ((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +                        (-DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)+
     +                          DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1))*SQA(1)
     +                        -(DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)+
     +                          DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1))*SQA(3)
                  END DO !LM2
                END DO !LM1
              END DO !I1SP2
            END DO !I1SP1

          ELSE IF (ISIGMA==3) THEN !T_z

            DO I1SP1=1,2
              DO I1SP2=1,2
                DO LM1=1,LMMAX
                  DO LM2=1,LMMAX
            TORQ((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +                         (DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)+
     +                          DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1))*SQA(2)
     +             -(0d0,1d0)*(-DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)+
     +                          DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1))*SQA(1)
                  END DO !LM2
                END DO !LM1
              END DO !I1SP2
            END DO !I1SP1

          END IF

        END DO  !I1
      END DO    !ISIGMA

        open(unit=12, file='TBkkr_torq.txt',
     +     form='formatted',action='write')
        DO ISIGMA=1,3
          DO I1=1,NATYPD
            DO LM2=1,LMMAXSO
              DO LM1=1,LMMAXSO
                WRITE(12,'(2ES25.16)') TORQ(LM1,LM2,I1,ISIGMA)
              END DO
            END DO
          END DO
        END DO
        close(12)

      DEALLOCATE(RLL)
      DEALLOCATE(DENS)
      DEALLOCATE(RLL_12)
      DEALLOCATE(TORQ)

      END

