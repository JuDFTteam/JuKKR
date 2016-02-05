      SUBROUTINE NORMCOEFF_SO(IRMINSO,IRCUT,LMAX,
     +                   LMMAX,PNS,THETAS,NTCELL,
     +                   IFUNM,IPAN,LMSP,KSRA,CLEB,ICLEB,IEND,DRDI,
     +                   IRWS,INTERF,ISP,IMPLAYER,NSPOH)
c ************************************************************************
c
c     calculates the norm of the wavefunctions with full potential and
c     spin orbit coupling. this is needed for the normalization of the
c     coefficients c_Lks .
c
c     attention : the gaunt coeffients are stored in index array
c                   (see subroutine gaunt)
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      IMPLICIT NONE

      include 'inc.p'
      INTEGER          LMMAXD
      PARAMETER        (LMMAXD= (LMAXD+1)**2)
      INTEGER          LMPOTD
      PARAMETER        (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER      ::   IEND,LMAX,LMMAX,KSRA,IRWS(*),ISP,LMMAXSO
      LOGICAL      ::   INTERF
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   PNS(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2,NATYPD)   ! non-sph. eigen states of single pot 
      DOUBLE PRECISION CLEB(*),THETAS(IRID,NFUND,*),
     +                 DRDI(IRMD,NATYPD)                            ! derivative dr/di
      INTEGER          ICLEB(NCLEB,4),IFUNM(NATYPD,LMPOTD),
     +                 LMSP(NATYPD,*),IRMINSO,IRCUT(0:IPAND,NATYPD),
     +                 IPAN(NATYPD),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   CONE,CZERO,PROD,NORM,SZHILF,NORM1
      DOUBLE PRECISION PI
      INTEGER          I,IFUN,IR,J,LM1,LM2,LM3,LM1P,LM2P,ICELL,ENT,
     +                 I1,LM3P,I1SP1,I1SP2,I1SP,
     +                 LMSP1,LMSP2,I2SP1,I2SP2,INSRA,NSRA,
     +                 NBET,NA,IMPLAYER,NSPOH,
     +                 LMDOS,ISIGMA
      INTEGER STATUS_READ
      LOGICAL          OPT,TEST
C     .. External Subroutines ..
      EXTERNAL         ZGEMM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        DATAN,DIMAG,DSQRT
C     ..Local Arrays..
      DOUBLE COMPLEX, ALLOCATABLE  :: DENS(:,:,:,:,:,:,:)
      DOUBLE COMPLEX, ALLOCATABLE  ::   RLL_12(:,:,:),
     +                                  RLL(:,:,:,:,:,:,:),
     +                                  RHOD(:,:,:,:)
c      DOUBLE COMPLEX RLL(IRMD,LMMAX,LMMAX,NSPD,NSPD,NSPD,NATYPD),
c     +                RLL_12(IRMD,LMMAX,LMMAX)
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA CONE/ (0.0D0,1.0D0)/
C     ..
c
      PI=4.d0*DATAN(1.d0)
      LMMAXSO=2*LMMAXD
      WRITE(6,*) 'in norm coefficient'
      WRITE(6,*) "KSRA",KSRA
      IF (KSRA.GE.1) THEN    ! previously this was .GT. which is wrong for kvrel=1
         NSRA = 2
      ELSE
         NSRA = 1
      END IF
c      NSRA=1
      WRITE(6,*) "NSRA",NSRA

      ALLOCATE(RLL(IRMD,LMMAX,LMMAX,NSPOH,NSPOH,NSPOH,NATYPD))
      ALLOCATE(RLL_12(IRMD,LMMAX,LMMAX))
      ALLOCATE(DENS(LMMAXD,LMMAXD,NSPD,NSPD,NSPD,NSPD,NATYPD))
c      ALLOCATE(DENS1(LMMAXD,LMMAXD,LMMAXD,2,8,NATYPD))

      RLL=CZERO
      DENS=CZERO
c      DENS1=CZERO

c    rewrite the wavefunctions in RLL arrays of 1,2*LMMAXD
      DO I1=1,NATYPD
       WRITE(6,*) 'ATOM',I1

        DO INSRA=1,NSRA
          DO IR=IRMINSO,IRMD

            DO I1SP1=1,NSPOH
              DO I1SP2=1,NSPOH
                DO LM1 =1,LMMAX
                  LMSP1=(I1SP1-1)*LMMAXD+LM1
                  DO LM2 =1,LMMAX
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

        DO I1SP1=1,NSPOH
          DO I1SP2=1,NSPOH
            DO I2SP1=1,NSPOH
              DO I2SP2=1,NSPOH

                DO LM1 = 1,LMMAX
                  DO LM2 = 1,LMMAX
         
                    DO INSRA=1,NSRA
 
                    RLL_12=CZERO
         
                      DO LM1P = 1,LMMAX
                        DO LM2P = 1,LMMAX
           
                          DO IR=1,IRMD
           
                            RLL_12(IR,LM1P,LM2P)=
     +                    DCONJG(RLL(IR,LM1P,LM1,I1SP1,I1SP2,INSRA,I1))*
     +                            RLL(IR,LM2P,LM2,I2SP1,I2SP2,INSRA,I1)
                          END DO       !IR
           
                        END DO         !LM2P
                      END DO           !LM1P


                      CALL CALC_RHO_LL_SS(LMAX,LMMAX,RLL_12,
     +                     IRCUT(0:IPAND,I1),IPAN(I1),NTCELL(I1),THETAS,
     +                     CLEB,ICLEB,IEND,IFUNM,LMSP,IRMINSO,IRWS(I1),
     +                    DRDI(:,I1),NORM)
                      DENS(LM1,LM2,I1SP1,I1SP2,I2SP1,I2SP2,I1)=
     +                DENS(LM1,LM2,I1SP1,I1SP2,I2SP1,I2SP2,I1) + NORM
                    END DO   !NSRA

                  END DO         !LM2
                END DO           !LM1

              END DO            !I2SP2
            END DO              !I2SP1

          END DO         !I1SP2
        END DO           !I1SP1 

c      IF (TEST('LMDOS   ').AND.I1.EQ.IMPLAYER) THEN
      IF (TEST('LMDOS   ')) THEN
        DO I1SP1=1,2
           I1SP=0
          DO I1SP2=1,NSPOH
            DO I2SP1=1,NSPOH
               DO I2SP2=1,NSPOH
                 I1SP=I1SP+1

                DO LM1 = 1,LMMAX
                  DO LM2 = 1,LMMAX
         
                    DO INSRA=1,NSRA
 
                    RLL_12=CZERO
         
                      DO LM1P = 1,LMMAX
                        DO LM2P = 1,LMMAX
           
                          DO IR=1,IRMD
           
c                            RLL_12(IR,LM1P,LM2P)=
c     +                    DCONJG(RLL(IR,LM1P,LM1,I1SP1,I1SP2,INSRA,I1))*
c     +                            RLL(IR,LM2P,LM2,I2SP1,I2SP2,INSRA,I1)
                          END DO       !IR
           
                        END DO         !LM2P
                      END DO           !LM1P

                      DO LMDOS=1,LMMAX
c                      CALL CALC_RHO_LL_SS_LMDOS(LMAX,LMMAX,RLL_12,
c     +                     IRCUT(0:IPAND,I1),IPAN(I1),NTCELL(I1),THETAS,
c     +                     CLEB,ICLEB,IEND,IFUNM,LMSP,IRMINSO,IRWS(I1),
c     +                    DRDI(:,I1),NORM1,LMDOS)
          
c                      DENS1(LMDOS,LM1,LM2,I1SP1,I1SP,I1)=
c     +                DENS1(LMDOS,LM1,LM2,I1SP1,I1SP,I1) + NORM1
                      ENDDO
                    END DO   !NSRA

                  END DO         !LM2
                END DO           !LM1
              END DO            !I2SP2
            END DO              !I2SP1

          END DO         !I1SP2
        END DO           !I1SP1 

      ENDIF


      END DO             !I1
      DEALLOCATE(RLL)
      DEALLOCATE(RLL_12) 

c calculate rho
      ALLOCATE(RHOD(LMMAXSO,LMMAXSO,NATYPD,4))
        IF (NSPOH.NE.1) THEN
          DO ISIGMA=1,4
           DO I1=1,NATYPD
            IF (ISIGMA.EQ.1) THEN
             DO I1SP1=1,NSPOD
              DO I1SP2=1,NSPOD
               DO LM1=1,LMMAX
                DO LM2=1,LMMAX
          RHOD((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +           DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)+ 
     +           DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1) 
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ELSEIF (ISIGMA.EQ.2) THEN
             DO I1SP1=1,NSPOD
              DO I1SP2=1,NSPOD
               DO LM1=1,LMMAX
                DO LM2=1,LMMAX
          RHOD((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +             DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)+ 
     +             DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ELSEIF (ISIGMA.EQ.3) THEN
             DO I1SP1=1,NSPOD
              DO I1SP2=1,NSPOD
               DO LM1=1,LMMAX
                DO LM2=1,LMMAX
          RHOD((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
c     +          (0D0,1D0)*(DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)- 
c     +                     DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1)) 
     +          -(0D0,1D0)*(DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)-
     +                     DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1)) 
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ELSEIF (ISIGMA.EQ.4) THEN
              DO I1SP1=1,NSPOD
               DO I1SP2=1,NSPOD
                DO LM1=1,LMMAX
                 DO LM2=1,LMMAX
          RHOD((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
c     +            DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)- 
c     +            DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1) 
     +           -DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)+
     +            DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1) 
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ENDIF
            ENDDO
           ENDDO

c write to the file
        OPEN(UNIT=12,FILE='TBkkr_rhod.txt',FORM='formatted',
     +       ACTION='write')
          DO ISIGMA=1,4
           DO I1=1,NATYPD
            DO LM2=1,LMMAXSO
             DO LM1=1,LMMAXSO
              WRITE(12,'(2ES25.16)') RHOD(LM1,LM2,I1,ISIGMA)
             ENDDO
            ENDDO
           ENDDO
          ENDDO        
        CLOSE(12)
       ENDIF
      DEALLOCATE(DENS)
      DEALLOCATE(RHOD)
      END




