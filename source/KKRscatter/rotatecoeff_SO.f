      SUBROUTINE ROTATECOEFF_SO(LMMAX,INTERF,IMPLAYER,ILAYERS)
c ************************************************************************
c
c    rotate the degenerated coefficients along the certain quantization axis
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
      INTEGER      ::   LMMAX,LMMAXSO
      LOGICAL      ::   INTERF
C     ..
C     .. Array Arguments ..
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   CONE,CZERO,PROD,NORM,SZHILF,NORM1
      DOUBLE PRECISION PI,NUM_TAN,
     +                 THETA,PHI,SDE0,SDE1,SDE2,
     +                 SDE4(NATYPD),SDE5(NATYPD),
     +                 MAXSdn,x0(3),y0(3),z0(3),NUM_TAN1(2),
     +                 ALPHA1(4),BETA1(2),Sdn(4),STOT(2),
     +                 RA,DPHI,DTHETA,SUMBSUR,SUMBBULK
      INTEGER          I,IFUN,IR,J,LM1,LM2,LM3,LM1P,LM2P,ICELL,ENT,ENT0,
     +                 NKSUM,NK,NKDUM,ENTDUM,I1,LM3P,I1SP1,I1SP2,I1SP,
     +                 LMSP1,LMSP2,ISIGMA,I2SP1,I2SP2,INSRA,NSRA,
     +                 NBET,NA,NKF,KSUM,JTAKE,IMPLAYER,ILAYERS,
     +                 IPHI,NPHI,NTHETA,ITHETA,NCONDI,NKDUM1,NKDUM2
      LOGICAL          OPT,TEST
C     .. External Subroutines ..
      EXTERNAL         ZGEMM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        DATAN,DIMAG,DSQRT
C     ..Local Arrays..
      INTEGER, ALLOCATABLE         ::   ENTG(:)
      DOUBLE COMPLEX, ALLOCATABLE  ::   RHOD(:,:,:,:),
     +                                  CKHELP(:),
     +                                  CKGES(:,:,:,:),
     +                                  CKGES_INI(:,:,:,:),
     +                                  DENS_GESAMT(:,:,:),
     +                                  DENS_GESAMT_I1(:,:,:,:),
     +                                  S_12(:,:),S_12n(:,:),
     +                                  DENS_GESAMTn(:,:,:),
     +                                  DENS_GESAMT_I1n(:,:,:)
      DOUBLE PRECISION,ALLOCATABLE ::   KPOINT(:,:),BZKP(:,:),
     +                                  DENS_GESAMT_INI(:,:),
     +                                  ALPHA(:),BETA(:),V_FERMI(:,:),
     +                                  WEIGHTK(:),FERMIVL(:),
     +                                  A_SQUARE(:),B_SQUARE(:),
     +                                  A_SQUARE_I1(:,:),
     +                                  B_SQUARE_I1(:,:)
      DOUBLE COMPLEX               ::   BETAHILF
      CHARACTER*9                  ::   FILENAME
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA CONE/ (0.0D0,1.0D0)/
C     ..
c
      PI=4.d0*DATAN(1.d0)
      LMMAXSO=NSPD*LMMAXD

      WRITE(6,*) 'in rotation coeff'
c  read in the coefficients

c      OPEN(UNIT=19, FILE="coefficients_0",FORM="FORMATTED")
      OPEN(UNIT=19, FILE="coefficients_0",FORM="UNFORMATTED")
      OPEN(UNIT=21, FILE="S_12",FORM="FORMATTED")
      OPEN(UNIT=22, FILE="S_xyz",FORM="FORMATTED")
c      OPEN(UNIT=23, FILE="rhod",FORM="FORMATTED")
      OPEN(UNIT=23, FILE="rhod",FORM="UNFORMATTED")
      OPEN(UNIT=27, FILE="kpoints_FS",form='FORMATTED')
      OPEN(UNIT=28, FILE="vectorFERMIVL",form='FORMATTED')

      ALLOCATE(RHOD(LMMAXSO,LMMAXSO,NATYPD,4))
      RHOD=CZERO
      DO ISIGMA=1,4
       DO I1=1,NATYPD
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
c          READ(23,'(2e17.9)') RHOD(LM1,LM2,I1,ISIGMA)
          READ(23) RHOD(LM1,LM2,I1,ISIGMA)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

c      READ(19,'(I9)') NKSUM
      READ(19) NKSUM
      WRITE(6,*) 'NKSUM',NKSUM 


      ALLOCATE(KPOINT(3,NKSUM))
      ALLOCATE(BZKP(3,NKSUM))
      ALLOCATE(CKGES_INI(2*LMMAXD,NATYPD,4,NKSUM))
      ALLOCATE(ENTG(NKSUM))
      ALLOCATE(S_12(NKSUM,3))
      ALLOCATE(DENS_GESAMT_INI(NKSUM,3))
      ALLOCATE(WEIGHTK(NKSUM))
      ALLOCATE(FERMIVL(NKSUM))
      ALLOCATE(V_FERMI(3,NKSUM))

      KPOINT=0d0
      CKGES_INI=CZERO
      DENS_GESAMT_INI=0d0
      ENTG=0
            

      DO NK=1,NKSUM
        
c        READ(19,'((I9),(3E17.9))') NKDUM,(KPOINT(J,NK), J=1,3)
c        READ(19,'(I5)') ENTG(NK)
        READ(19) NKDUM,(KPOINT(J,NK), J=1,3)
        READ(19) ENTG(NK)

        DO ENT=1,ENTG(NK)
c          READ(19,"(I5)") ENTDUM
          READ(19) ENTDUM
          DO I1=1,NATYPD
            DO LM1 =1,LMMAXSO 
c             READ(19,"(4(E20.12,2X))") CKGES_INI(LM1,I1,ENT,NK)
             READ(19) CKGES_INI(LM1,I1,ENT,NK)

            END DO      !LM1=1,2*LMMAXD
          END DO        !I1=1,NATYP
        END DO          !ENT=1,ENTG
      
        READ(21,'((I5),(7e17.9))') NKDUM1,(S_12(NK,ISIGMA),ISIGMA=1,3)
        READ(22,'((I9),(7e17.9))') NKDUM2,(BZKP(J,NK),J=1,2),
     +                        (DENS_GESAMT_INI(NK,ISIGMA),ISIGMA=1,3)
        READ(27,"((4E17.9))") (BZKP(J,NK),J=1,3),WEIGHTK(NK)
        READ(28,"((3E17.9))") (V_FERMI(J,NK),J=1,3)
        FERMIVL(NK)=SQRT(V_FERMI(1,NK)**2+V_FERMI(2,NK)**2+
     +                   V_FERMI(3,NK)**2)

      END DO            !NK=1,NKSUM

      CLOSE(19)
      CLOSE(21)
      CLOSE(22)
      CLOSE(27)
      CLOSE(28)

      ALLOCATE(CKGES(2*LMMAXD,NATYPD,4,NKSUM))
      ALLOCATE(CKHELP(2*LMMAXD))
      ALLOCATE(DENS_GESAMT(4,NKSUM,4))
      ALLOCATE(DENS_GESAMT_I1(NATYPD,4,NKSUM,4))
      ALLOCATE(S_12n(NKSUM,3))
      ALLOCATE(DENS_GESAMTn(2,NKSUM,3))
      ALLOCATE(DENS_GESAMT_I1n(NATYPD,2,NKSUM))
      ALLOCATE(ALPHA(NKSUM))
      ALLOCATE(BETA(NKSUM))
      ALLOCATE(A_SQUARE(NKSUM))
      ALLOCATE(B_SQUARE(NKSUM))
      ALLOCATE(A_SQUARE_I1(NATYPD,NKSUM))
      ALLOCATE(B_SQUARE_I1(NATYPD,NKSUM))
      CKHELP=CZERO
      CKGES=CZERO
      DENS_GESAMTn=CZERO
      DENS_GESAMT_I1n=CZERO
      S_12n=CZERO


      RA=1d0
      NTHETA=10
      DTHETA=(PI/2d0)/(NTHETA-1)
c      NTHETA=180
c      DTHETA=(2d0*PI)/(NTHETA)
      NPHI=10
      DPHI=(PI/2d0)/(NPHI-1)
c      NPHI=180
c      DPHI=(2d0*PI)/(NPHI)
c      DO ITHETA=1,NTHETA
c      DO ITHETA=1,1
       DO ITHETA=NTHETA,NTHETA
       THETA=(ITHETA-1)*DTHETA
c       THETA=0d0
       z0(1)=RA*SIN(THETA)
       z0(2)=0d0
       z0(3)=-RA*COS(THETA)
c        DO IPHI=1,NPHI
c        DO IPHI=NPHI,NPHI
        DO IPHI=1,1
         IF (THETA.EQ.NTHETA.AND.IPHI.GT.1) GO TO 100
         PHI=(IPHI-1)*DPHI
c         PHI=PI/4d0
         x0(1)=RA*COS(THETA)*COS(PHI)
         y0(1)=RA*COS(THETA)*SIN(PHI)
         x0(2)=RA*SIN(PHI)
         y0(2)=-RA*COS(PHI)
         x0(3)=RA*SIN(THETA)*COS(PHI)
         y0(3)=RA*SIN(THETA)*SIN(PHI)
       WRITE(*,*) 'IPHI',IPHI,'ITHETA',ITHETA
       WRITE(*,"(2E17.9)") PHI/PI*180,THETA/PI*180
c       x0(1)=1d0
c       y0(1)=0d0
c       z0(1)=0d0 
c       x0(2)=0d0
c       y0(2)=1d0
c       z0(2)=0d0 
c       x0(3)=0d0
c       y0(3)=0d0
c       z0(3)=1d0 
c      DO NCONDI=1,2
c      DO NCONDI=1,1  ! S_x=S_y=0
      DO NCONDI=2,2   ! S_z=max
      DO NK =1,NKSUM

        IF (ENTG(NK) .NE. 1) THEN

          BETA(NK)=0d0
          ALPHA(NK)=0d0

c          BETAHILF= -(1d0)*( DENS_GESAMT(1,NK,4)*CONJG(S_12(NK,3,1))-
c     +                       DENS_GESAMT(1,NK,3)*CONJG(S_12(NK,4,1)))/
c     +                     ( DENS_GESAMT(1,NK,4)*S_12(NK,3,1)-
c     +                       DENS_GESAMT(1,NK,3)*S_12(NK,4,1))
c          BETA(NK)=REAL((0d0,-1d0)/2d0*LOG(BETAHILF))
c          NUM_TAN=REAL(EXP((0d0,1d0)*BETA(NK))*S_12(NK,4,1))
c          ALPHA(NK)= DATAN2(REAL(DENS_GESAMT(1,NK,4)),-NUM_TAN)
c alpha,beta for quantized axis along z direction (S_x and S_y =0)  
c          BETAHILF= -(1d0)*( DENS_GESAMT(1,NK,2)*CONJG(S_12(NK,3,1))-
c     +                       DENS_GESAMT(1,NK,3)*CONJG(S_12(NK,2,1)))/
c     +                     ( DENS_GESAMT(1,NK,2)*S_12(NK,3,1)-
c     +                       DENS_GESAMT(1,NK,3)*S_12(NK,2,1))
c          BETA(NK)=REAL((0d0,-1d0)/2d0*LOG(BETAHILF))
c          NUM_TAN=REAL(EXP((0d0,1d0)*BETA(NK))*S_12(NK,2,1))
c          ALPHA(NK)= DATAN2(REAL(DENS_GESAMT(1,NK,2)),-NUM_TAN)

c alpha,beta for quatizied axis along z direction (S_z max) 
c          BETA(NK)= -DATAN2(AIMAG(S_12(NK,4,1)),REAL(S_12(NK,4,1)))
c          NUM_TAN=REAL(EXP((0d0,1d0)*BETA(NK))*S_12(NK,4,1))
c          ALPHA(NK)= DATAN2(NUM_TAN,REAL(DENS_GESAMT(1,NK,4)))

          DO J=1,3
          S_12n(NK,J)=x0(J)*S_12(NK,1)+
     +              y0(J)*S_12(NK,2)+
     +              z0(J)*S_12(NK,3)
          DENS_GESAMTn(1,NK,J)=x0(J)*DENS_GESAMT_INI(NK,1)+
     +                     y0(J)*DENS_GESAMT_INI(NK,2)+
     +                     z0(J)*DENS_GESAMT_INI(NK,3)
          ENDDO
c alpha,beta for quatizied axis along n direction (S_x=S_y=0) 
          IF (NCONDI.EQ.1) THEN
          BETAHILF= -(1d0)*( DENS_GESAMTn(1,NK,3)*CONJG(S_12n(NK,2))-
     +                       DENS_GESAMTn(1,NK,2)*CONJG(S_12n(NK,3)))/
     +                     ( DENS_GESAMTn(1,NK,3)*S_12n(NK,2)-
     +                       DENS_GESAMTn(1,NK,2)*S_12n(NK,3))
          BETA(NK)=REAL((0d0,-1d0)/2d0*LOG(BETAHILF))
          NUM_TAN=REAL(EXP((0d0,1d0)*BETA(NK))*S_12n(NK,3))
          ALPHA(NK)= DATAN2(REAL(DENS_GESAMTn(1,NK,3)),-NUM_TAN)

c alpha,beta for quatizied axis along n direction (S_n max) 
          ELSE IF (NCONDI.EQ.2) THEN

          BETA1=0d0
          ALPHA1=0d0
          NUM_TAN1=0d0
          Sdn=0d0

          BETA1(1)= -DATAN2(AIMAG(S_12n(NK,1)),REAL(S_12n(NK,1)))
          BETA1(2)= BETA1(1)+PI
          NUM_TAN1(1)=REAL(EXP((0d0,1d0)*BETA1(1))*S_12n(NK,1))
          NUM_TAN1(2)=REAL(EXP((0d0,1d0)*BETA1(2))*S_12n(NK,1))
          ALPHA1(1)= DATAN2(NUM_TAN1(1),REAL(DENS_GESAMTn(1,NK,1)))
          ALPHA1(2)= ALPHA1(1)+PI
          ALPHA1(3)= DATAN2(NUM_TAN1(2),REAL(DENS_GESAMTn(1,NK,1)))
          ALPHA1(4)= ALPHA1(3)+PI
          Sdn(1)=COS(ALPHA1(1))*REAL(DENS_GESAMTn(1,NK,1))+
     +              SIN(ALPHA1(1))*NUM_TAN1(1)          
          Sdn(2)=COS(ALPHA1(2))*REAL(DENS_GESAMTn(1,NK,1))+
     +              SIN(ALPHA1(2))*NUM_TAN1(1)          
          Sdn(3)=COS(ALPHA1(3))*REAL(DENS_GESAMTn(1,NK,1))+
     +              SIN(ALPHA1(3))*NUM_TAN1(2)          
          Sdn(4)=COS(ALPHA1(4))*REAL(DENS_GESAMTn(1,NK,1))+
     +              SIN(ALPHA1(4))*NUM_TAN1(2)          
          MAXSdn=0d0
          DO J=1,4
           IF (ABS(Sdn(J)).GT.MAXSdn) THEN
             MAXSdn=ABS(Sdn(J))
             JTAKE=J
           ENDIF
          ENDDO
          IF (JTAKE.GE.3) THEN
          ALPHA(NK)=ALPHA1(JTAKE)
          BETA(NK)=BETA1(2)
          ELSE
          ALPHA(NK)=ALPHA1(JTAKE)
          BETA(NK)=BETA1(1)          
          ENDIF
         ENDIF ! NCONDI
          DO I1=1,NATYPD
            DO LM1=1,LMMAXSO
              CKGES(LM1,I1,1,NK)=
     +           COS(ALPHA(NK)/2d0)*CKGES_INI(LM1,I1,1,NK)+
     +           SIN(ALPHA(NK)/2d0)*EXP((0d0,1d0)*BETA(NK))*
     +                                       CKGES_INI(LM1,I1,2,NK)
              CKGES(LM1,I1,2,NK)=
     +          -SIN(ALPHA(NK)/2d0)*CKGES_INI(LM1,I1,1,NK)+
     +           COS(ALPHA(NK)/2d0)*EXP((0d0,1d0)*BETA(NK))*
     +                                        CKGES_INI(LM1,I1,2,NK)
             IF (ENTG(NK).GT.2) THEN
              CKGES(LM1,I1,3,NK)=
     +           COS(ALPHA(NK)/2d0)*CKGES_INI(LM1,I1,3,NK)+
     +           SIN(ALPHA(NK)/2d0)*EXP((0d0,1d0)*BETA(NK))*
     +                                       CKGES_INI(LM1,I1,4,NK)
              CKGES(LM1,I1,4,NK)=
     +          -SIN(ALPHA(NK)/2d0)*CKGES_INI(LM1,I1,3,NK)+
     +           COS(ALPHA(NK)/2d0)*EXP((0d0,1d0)*BETA(NK))*
     +                                        CKGES_INI(LM1,I1,4,NK)
             ENDIF  
            END DO
          END DO        !I1=1,NATYP
        END IF ! ENTG

      END DO    !NK


      DENS_GESAMT=CZERO
      DENS_GESAMT_I1=CZERO


      DO NK=1,NKSUM                    

        DO ENT=1,ENTG(NK)

          DO I1=1,NATYPD

            CKHELP=CZERO
            CALL ZGEMM('N','N',LMMAXSO,1,LMMAXSO,1d0,
     +                                           RHOD(:,:,I1,1),
     +          LMMAXSO,CKGES(:,I1,ENT,NK),LMMAXSO,0d0,CKHELP,LMMAXSO)

            DENS_GESAMT_I1(I1,ENT,NK,1)=
     +                        DOT_PRODUCT(CKGES(:,I1,ENT,NK),CKHELP)
            DENS_GESAMT(ENT,NK,1)=DENS_GESAMT(ENT,NK,1)+
     +                                 DENS_GESAMT_I1(I1,ENT,NK,1)

          END DO  !I1

          DO I1=1,NATYPD
            DO LM1=1,LMMAXSO
              CKGES(LM1,I1,ENT,NK)=
     +           CKGES(LM1,I1,ENT,NK)/CDSQRT((DENS_GESAMT(ENT,NK,1)))
            END DO           !LM1
 
          END DO             !I1
        END DO               !ENT
      END DO                 !NK


c      IF (NCONDI.EQ.1.AND.z0(1).EQ.1d0.AND.IPHI.EQ.1) THEN     
      IF (NCONDI.EQ.2.AND.z0(1).EQ.1d0.AND.IPHI.EQ.1) THEN     
c      OPEN(UNIT=171,FILE="Density_k_after_rot",FORM='FORMATTED')
c      OPEN(UNIT=172,FILE="Density_k_layers",FORM='FORMATTED')
      OPEN(UNIT=152,FILE="S_xyz_alpha",FORM='FORMATTED')
      OPEN(UNIT=160,FILE="Parameter_a_b",FORM='FORMATTED')
c      OPEN(UNIT=181,FILE="coefficients_0_rot",FORM='FORMATTED')
      OPEN(UNIT=181,FILE="coefficients_0_rot",FORM='UNFORMATTED')



      DENS_GESAMT=CZERO
      DENS_GESAMT_I1=CZERO
      DENS_GESAMTn=CZERO
      DENS_GESAMT_I1n=CZERO

c      WRITE(181,"(I9)") NKSUM 
      WRITE(181) NKSUM 
c      WRITE(171,"(I9)") NKSUM 

      DO NK=1,NKSUM                    

c        WRITE(181,"((I9),(3(E17.9)))") NK,(KPOINT(J,NK), J=1,3)
c        WRITE(181,"(I5)") ENTG(NK)
        WRITE(181) NK,(KPOINT(J,NK), J=1,3)
        WRITE(181) ENTG(NK)

        DO ENT=1,ENTG(NK)


          DO ISIGMA=1,4
            DO I1=1,NATYPD

              CKHELP=CZERO
              CALL ZGEMM('N','N',LMMAXSO,1,LMMAXSO,1d0,
     +                                           RHOD(:,:,I1,ISIGMA),
     +          LMMAXSO,CKGES(:,I1,ENT,NK),LMMAXSO,0d0,CKHELP,LMMAXSO)

              DENS_GESAMT_I1(I1,ENT,NK,ISIGMA)=
     +                        DOT_PRODUCT(CKGES(:,I1,ENT,NK),CKHELP)
              DENS_GESAMT(ENT,NK,ISIGMA)=DENS_GESAMT(ENT,NK,ISIGMA)+
     +                                 DENS_GESAMT_I1(I1,ENT,NK,ISIGMA)


            END DO  !I1

         END DO  !ISIGMA

c          WRITE(171,"((2I9),(9(e17.9)))") ENT,NK,
c     +             (KPOINT(J,NK), J=1,3), 
c     +         DENS_GESAMT(ENT,NK,1),CDSQRT(DENS_GESAMT(ENT,NK,1))
        END DO               !ENT
c          DO I1=1,NATYPD
c          WRITE(172,"(1I9,4e17.9)") NK,(KPOINT(J,NK),J=1,2),
c     +         DENS_GESAMT_I1(I1,1,NK,1)
c          ENDDO
        IF (ENTG(NK) .GE. 2 ) THEN 
         DO ENT=1,ENTG(NK)
          DO J=1,3
           DENS_GESAMTn(ENT,NK,J)=x0(J)*DENS_GESAMT(ENT,NK,2)+
     +                       y0(J)*DENS_GESAMT(ENT,NK,3)+
     +                       z0(J)*DENS_GESAMT(ENT,NK,4)
          ENDDO 
          DO I1=1,NATYPD
           DENS_GESAMT_I1n(I1,ENT,NK)=x0(1)*DENS_GESAMT_I1(I1,ENT,NK,2)+
     +                                y0(1)*DENS_GESAMT_I1(I1,ENT,NK,3)+
     +                                z0(1)*DENS_GESAMT_I1(I1,ENT,NK,4)
          ENDDO
         ENDDO
          IF (REAL(DENS_GESAMTn(1,NK,1)) .GT. 0d0 ) THEN
            DO ENT=1,ENTG(NK)
c              WRITE(181,"(I5)") ENT
              WRITE(181) ENT
              DO I1=1,NATYPD
                DO LM1=1,LMMAXSO
c                  WRITE(181,"(4(E20.12,2x))") CKGES(LM1,I1,ENT,NK)
                  WRITE(181) CKGES(LM1,I1,ENT,NK)
                END DO           !LM1
              END DO             !I1
            END DO
          ELSE
            ENT0=0
            DO ENT=ENTG(NK),1,-1
              ENT0=ENT0+1
c              WRITE(181,"(I5)") ENT0
              WRITE(181) ENT0
              DO I1=1,NATYPD
                DO LM1=1,LMMAXSO
c                  WRITE(181,"(4(E20.12,2x))") CKGES(LM1,I1,ENT,NK)
                  WRITE(181) CKGES(LM1,I1,ENT,NK)

                END DO           !LM1
              END DO             !I1
            END DO
            SZHILF=CZERO
            SZHILF=DENS_GESAMTn(2,NK,1)
            DENS_GESAMTn(2,NK,1)=DENS_GESAMTn(1,NK,1)
            DENS_GESAMTn(1,NK,1)=SZHILF
            DO J=2,4
             SZHILF=CZERO
             SZHILF=DENS_GESAMT(2,NK,J)
             DENS_GESAMT(2,NK,J)=DENS_GESAMT(1,NK,J)
             DENS_GESAMT(1,NK,J)=SZHILF
            ENDDO
           DO I1=1,NATYPD
            SZHILF=CZERO
            SZHILF=DENS_GESAMT_I1n(I1,2,NK)
            DENS_GESAMT_I1n(I1,2,NK)=DENS_GESAMT_I1n(I1,1,NK)
            DENS_GESAMT_I1n(I1,1,NK)=SZHILF
           ENDDO
          ENDIF
        ELSE
          DO ENT=1,ENTG(NK)
c            WRITE(181,"(I5)") ENT
            WRITE(181) ENT
            DO I1=1,NATYPD
              DO LM1=1,LMMAXSO
c                WRITE(181,"(4(E20.12,2x))") CKGES(LM1,I1,ENT,NK)
                WRITE(181) CKGES(LM1,I1,ENT,NK)
              END DO           !LM1
            END DO             !I1
          END DO
        END IF

        WRITE(152,"((I9),(7(e17.9)))") NK,(KPOINT(J,NK),J=1,2),
     +               REAL(DENS_GESAMT(1,NK,2)),
     +               REAL(DENS_GESAMT(1,NK,3)),
     +               REAL(DENS_GESAMT(1,NK,4))
       A_SQUARE(NK)=0.5d0*(REAL(DENS_GESAMT(1,NK,1))+
     +             ABS(REAL(DENS_GESAMTn(1,NK,1))))
       B_SQUARE(NK)=0.5d0*(REAL(DENS_GESAMT(1,NK,1))-
     +             ABS(REAL(DENS_GESAMTn(1,NK,1))))
        WRITE(160,"((I9),(7(e17.9)))") NK,(KPOINT(J,NK),J=1,2),
     +       A_SQUARE(NK),B_SQUARE(NK)
       DO I1=1,NATYPD
       A_SQUARE_I1(I1,NK)=0.5d0*(REAL(DENS_GESAMT_I1(I1,1,NK,1))+
     +             ABS(REAL(DENS_GESAMT_I1n(I1,1,NK))))
       B_SQUARE_I1(I1,NK)=0.5d0*(REAL(DENS_GESAMT_I1(I1,1,NK,1))-
     +             ABS(REAL(DENS_GESAMT_I1n(I1,1,NK))))
       ENDDO
      END DO              !NK

      OPEN(UNIT=29,file='SUM_ALL_FS',form='formatted')

      SDE0=0D0
      SDE1=0D0
      SDE2=0D0
      SDE4=0D0
      SDE5=0D0
      DO NKF=1,NKSUM
        SDE0=SDE0+WEIGHTK(NKF)/FERMIVL(NKF)
        SDE1=SDE1+A_SQUARE(NKF)*WEIGHTK(NKF)/FERMIVL(NKF)
        SDE2=SDE2+B_SQUARE(NKF)*WEIGHTK(NKF)/FERMIVL(NKF)
       DO I1=1,NATYPD
        SDE4(I1)=SDE4(I1)+A_SQUARE_I1(I1,NKF)*WEIGHTK(NKF)/FERMIVL(NKF)
        SDE5(I1)=SDE5(I1)+B_SQUARE_I1(I1,NKF)*WEIGHTK(NKF)/FERMIVL(NKF)
       ENDDO
      ENDDO
      WRITE(29,*) 'A_SQUARE',SDE1/SDE0,1d0-SDE1/SDE0
      WRITE(29,*) 'B_SQUARE',SDE2/SDE0,1d0-SDE2/SDE0
      WRITE(29,*) 'B_SQUARE LAYER'
      DO I1=1,NATYPD
      WRITE(29,*) I1,SDE5(I1)/SDE0,SDE4(I1)/SDE0
      ENDDO
      CLOSE(29)

c      CLOSE(171)
c      CLOSE(172)
      CLOSE(152)
      CLOSE(160)
      CLOSE(181)

      ENDIF      ! writing at certain direction

c anisotropy of b square 
      IF (TEST('RASHBA1  ')) THEN
        WRITE(300,"(4E17.9)") THETA/PI*180,PHI/PI*180
        WRITE(301,"(4E17.9)") THETA/PI*180,PHI/PI*180
        WRITE(400,"(4E17.9)") THETA/PI*180,PHI/PI*180
        WRITE(401,"(4E17.9)") THETA/PI*180,PHI/PI*180
        WRITE(500,"(4E17.9)") THETA/PI*180,PHI/PI*180
      ENDIF

      IF (NCONDI.EQ.2) THEN
      IF (TEST('COEFF   ')) THEN
      WRITE(FILENAME,"((A5),(I3.3))") 'coeff',(ITHETA-1)*10+IPHI
      OPEN(318,FILE=FILENAME,FORM='UNFORMATTED')
c      WRITE(200+(ITHETA-1)*10+IPHI,"(I9)") NKSUM 
      WRITE(318) NKSUM 
      ENDIF
      ENDIF

        DENS_GESAMT=CZERO
        DENS_GESAMT_I1=CZERO
        DENS_GESAMTn=CZERO
        DENS_GESAMT_I1n=CZERO

       DO NK=1,NKSUM

        DO ENT=1,ENTG(NK)

          DO ISIGMA=1,4
            DO I1=1,NATYPD

              CKHELP=CZERO
              CALL ZGEMM('N','N',LMMAXSO,1,LMMAXSO,1d0,
     +                                           RHOD(:,:,I1,ISIGMA),
     +          LMMAXSO,CKGES(:,I1,ENT,NK),LMMAXSO,0d0,CKHELP,LMMAXSO)

              DENS_GESAMT_I1(I1,ENT,NK,ISIGMA)=
     +                        DOT_PRODUCT(CKGES(:,I1,ENT,NK),CKHELP)
              DENS_GESAMT(ENT,NK,ISIGMA)=DENS_GESAMT(ENT,NK,ISIGMA)+
     +                                 DENS_GESAMT_I1(I1,ENT,NK,ISIGMA)

            END DO  !I1

         END DO  !ISIGMA

        END DO               !ENT

        IF (ENTG(NK) .GE. 2 ) THEN 
         DO ENT=1,ENTG(NK)
          DO J=1,3
          DENS_GESAMTn(ENT,NK,J)=x0(J)*DENS_GESAMT(ENT,NK,2)+
     +                       y0(J)*DENS_GESAMT(ENT,NK,3)+
     +                       z0(J)*DENS_GESAMT(ENT,NK,4)
          ENDDO
          DO I1=1,NATYPD
           DENS_GESAMT_I1n(I1,ENT,NK)=x0(1)*DENS_GESAMT_I1(I1,ENT,NK,2)+
     +                              y0(1)*DENS_GESAMT_I1(I1,ENT,NK,3)+
     +                              z0(1)*DENS_GESAMT_I1(I1,ENT,NK,4)
          ENDDO
        ENDDO
          
        
c write down coefficients after rotation
      IF (NCONDI.EQ.2) THEN
      IF (TEST('COEFF   ')) THEN
c        WRITE(200+(ITHETA-1)*10+IPHI,"((I9),(3(E17.9)))") NK,
c     +               (KPOINT(J,NK), J=1,3)
c        WRITE(200+(ITHETA-1)*10+IPHI,"(I5)") ENTG(NK)
        WRITE(318) NK,
     +               (KPOINT(J,NK), J=1,3)
        WRITE(318) ENTG(NK)
          IF (REAL(DENS_GESAMTn(1,NK,1)) .GT. 0d0 ) THEN
            DO ENT=1,ENTG(NK)
c              WRITE(200+(ITHETA-1)*10+IPHI,"(I5)") ENT
              WRITE(318) ENT
              DO I1=1,NATYPD
                DO LM1=1,LMMAXSO
c                  WRITE(200+(ITHETA-1)*10+IPHI,"(4(E20.12,2x))") 
c     +                   CKGES(LM1,I1,ENT,NK)
                  WRITE(318) 
     +                   CKGES(LM1,I1,ENT,NK)
                END DO           !LM1
              END DO             !I1
            END DO
          ELSE
            ENT0=0
            DO ENT=ENTG(NK),1,-1
              ENT0=ENT0+1
c              WRITE(200+(ITHETA-1)*10+IPHI,"(I5)") MOD(ENT,2) +1
              WRITE(318) ENT0
              DO I1=1,NATYPD
                DO LM1=1,LMMAXSO
c                  WRITE(200+(ITHETA-1)*10+IPHI,"(4(E20.12,2x))") 
c     +                    CKGES(LM1,I1,ENT,NK)
                  WRITE(318) 
     +                    CKGES(LM1,I1,ENT,NK)
                END DO           !LM1
              END DO             !I1
            END DO
         ENDIF
       ENDIF
       ENDIF
        IF (TEST('RASHBA1  ')) THEN
          WRITE(400,"(4e17.9)") (KPOINT(J,NK), J=1,2),
     +          REAL(DENS_GESAMT_I1n(IMPLAYER,1,NK)),
     +          REAL(DENS_GESAMT_I1n(IMPLAYER,2,NK))
          WRITE(401,"(4e17.9)") (KPOINT(J,NK), J=1,2),
     +          REAL(DENS_GESAMT_I1n(IMPLAYER+ILAYERS-1,1,NK)),
     +          REAL(DENS_GESAMT_I1n(IMPLAYER+ILAYERS-1,2,NK))
        ENDIF
          IF (REAL(DENS_GESAMTn(1,NK,1)) .LT. 0d0 ) THEN
            SZHILF=CZERO
            SZHILF=DENS_GESAMTn(2,NK,1)
            DENS_GESAMTn(2,NK,1)=DENS_GESAMTn(1,NK,1)
            DENS_GESAMTn(1,NK,1)=SZHILF
          DO I1=1,NATYPD
            SZHILF=CZERO
            SZHILF=DENS_GESAMT_I1n(I1,2,NK)
            DENS_GESAMT_I1n(I1,2,NK)=DENS_GESAMT_I1n(I1,1,NK)
            DENS_GESAMT_I1n(I1,1,NK)=SZHILF
          ENDDO
         ENDIF
        ENDIF

        IF (TEST('RASHBA1  ')) THEN
          WRITE(300,"(4e17.9)") (KPOINT(J,NK), J=1,2),
     +          REAL(DENS_GESAMT_I1n(IMPLAYER,1,NK)),
     +          REAL(DENS_GESAMT_I1n(IMPLAYER,2,NK))
          WRITE(301,"(4e17.9)") (KPOINT(J,NK), J=1,2),
     +          REAL(DENS_GESAMT_I1n(IMPLAYER+ILAYERS-1,1,NK)),
     +          REAL(DENS_GESAMT_I1n(IMPLAYER+ILAYERS-1,2,NK))
          WRITE(500,"(4e17.9)") (KPOINT(J,NK), J=1,2),
     +          REAL(DENS_GESAMTn(1,NK,1)),
     +          REAL(DENS_GESAMTn(2,NK,1))
        ENDIF

! write out the bsquare and asquare

       A_SQUARE(NK)=0.5d0*(REAL(DENS_GESAMT(1,NK,1))+
     +             REAL(DENS_GESAMTn(1,NK,1)))
       B_SQUARE(NK)=0.5d0*(REAL(DENS_GESAMT(1,NK,1))-
     +             REAL(DENS_GESAMTn(1,NK,1)))
       DO I1=1,NATYPD
       A_SQUARE_I1(I1,NK)=0.5d0*(REAL(DENS_GESAMT_I1(I1,1,NK,1))+
     +             REAL(DENS_GESAMT_I1n(I1,1,NK)))
       B_SQUARE_I1(I1,NK)=0.5d0*(REAL(DENS_GESAMT_I1(I1,1,NK,1))-
     +             REAL(DENS_GESAMT_I1n(I1,1,NK)))
       ENDDO
       
       END DO              !NK

      IF (NCONDI.EQ.1) THEN
      OPEN(UNIT=30,file='anibsqutot_zeroz',form='formatted')
c      OPEN(UNIT=31,file='anibsqusurf_zeroz',form='formatted')
c      OPEN(UNIT=32,file='anibsqusurf_topbottom_zeroz',form='formatted')
c      OPEN(UNIT=33,file='anibsqusurflayer_zeroz',form='formatted')
c      OPEN(UNIT=34,file='anibsqubulklayer_zeroz',form='formatted')
      ELSE
      OPEN(UNIT=40,file='anibsqutot_maxz',form='formatted')
c      OPEN(UNIT=41,file='anibsqusurf_maxz',form='formatted')
c      OPEN(UNIT=42,file='anibsqusurf_topbottom_maxz',form='formatted')
c      OPEN(UNIT=43,file='anibsqusurflayer_maxz',form='formatted')
c      OPEN(UNIT=44,file='anibsqubulklayer_maxz',form='formatted')
      ENDIF
    
      SDE0=0D0
      SDE1=0D0
      SDE2=0D0
      SDE4=0D0
      SDE5=0D0
      DO NKF=1,NKSUM
        SDE0=SDE0+WEIGHTK(NKF)/FERMIVL(NKF)
        SDE1=SDE1+A_SQUARE(NKF)*WEIGHTK(NKF)/FERMIVL(NKF)
        SDE2=SDE2+B_SQUARE(NKF)*WEIGHTK(NKF)/FERMIVL(NKF)
       DO I1=1,NATYPD
        SDE4(I1)=SDE4(I1)+A_SQUARE_I1(I1,NKF)*WEIGHTK(NKF)/FERMIVL(NKF)
        SDE5(I1)=SDE5(I1)+B_SQUARE_I1(I1,NKF)*WEIGHTK(NKF)/FERMIVL(NKF)
       ENDDO
      ENDDO
      WRITE(30+(NCONDI-1)*10,"(4E17.9)") x0(1),y0(1),z0(1),SDE2/SDE0
      IF (IPHI.EQ.NPHI) WRITE(30+(NCONDI-1)*10,*) ''
c      WRITE(31+(NCONDI-1)*10,"(4E17.9)") x0(1),y0(1),z0(1),
c     +                   (SDE5(IMPLAYER)+SDE5(IMPLAYER+ILAYERS-1))/SDE0
c      IF (IPHI.EQ.NPHI) WRITE(31+(NCONDI-1)*10,*) ''
c      WRITE(32+(NCONDI-1)*10,"(5E17.9)")x0(1),y0(1),z0(1),
c     +                SDE5(IMPLAYER)/SDE0,SDE5(IMPLAYER+ILAYERS-1)/SDE0
      SUMBSUR=0d0
      SUMBBULK=0d0
      DO I1=1,NATYPD
       IF (I1.LE.IMPLAYER.OR.I1.GE.(IMPLAYER+ILAYERS-1)) THEN
       SUMBSUR=SUMBSUR+SDE5(I1)/SDE0
       ELSE
       SUMBBULK=SUMBBULK+SDE5(I1)/SDE0
       ENDIF
      ENDDO
c      WRITE(33+(NCONDI-1)*10,"(4E17.9)") x0(1),y0(1),z0(1),
c     +                   SUMBSUR
c      WRITE(34+(NCONDI-1)*10,"(4E17.9)") x0(1),y0(1),z0(1),
c     +                   SUMBBULK
      ENDDO ! NCONDI
      ENDDO ! THETA
      ENDDO ! PHI
  100 CONTINUE
      CLOSE(30)
      CLOSE(31)
      CLOSE(32)
      CLOSE(33)
      CLOSE(34)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      CLOSE(43)
      CLOSE(44)
      CLOSE(318)
  

      DEALLOCATE(KPOINT)
      DEALLOCATE(CKGES)
      DEALLOCATE(ENTG)

      DEALLOCATE(CKGES_INI)
      DEALLOCATE(CKHELP)
      DEALLOCATE(DENS_GESAMT)
      DEALLOCATE(DENS_GESAMT_INI)
      DEALLOCATE(RHOD)
      DEALLOCATE(DENS_GESAMT_I1)
      DEALLOCATE(S_12)

      DEALLOCATE(WEIGHTK)
      DEALLOCATE(FERMIVL)
      DEALLOCATE(S_12n)
      DEALLOCATE(DENS_GESAMTn)
      DEALLOCATE(DENS_GESAMT_I1n)
      DEALLOCATE(ALPHA)
      DEALLOCATE(BETA)
      DEALLOCATE(A_SQUARE)
      DEALLOCATE(B_SQUARE)
      DEALLOCATE(A_SQUARE_I1)
      DEALLOCATE(B_SQUARE_I1)


      STOP "after norm"

      END




