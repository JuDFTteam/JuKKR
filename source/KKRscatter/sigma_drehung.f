c 18.10.95 ***************************************************************
      SUBROUTINE SIGMA(IRMINSO,IRCUT,LMAX,
     +                   LMMAX,PNS,THETAS,NTCELL,
     +                   IFUNM,IPAN,LMSP,NSRA,CLEB,ICLEB,
     +                   IEND,DRDI,RM,IRWS,INTERF,ISP,IROTMAX)
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
      INTEGER          LMMAXSO
      PARAMETER        (LMMAXSO= 2*(LMAXD+1)**2)
      INTEGER          LMPOTD
      PARAMETER        (LMPOTD= (LPOTD+1)**2)
      INTEGER          IRMIND,IRLMD
      PARAMETER        (IRMIND=IRMD-IRNSD,IRLMD= (IRNSD+1)*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      INTEGER          IEND,LMAX,LMMAX,NSRA,NSPO,IRWS(*),ISP
      LOGICAL,INTENT(IN)      ::   INTERF
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   PNS(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2,NATYPD)   ! non-sph. eigen states of single pot 
      DOUBLE PRECISION CLEB(*),
     +                 THETAS(IRID,NFUND,*),
     +                 DRDI(IRMD,NATYPD),                            ! derivative dr/di
     +                 RM(IRMD,NATYPD)
      INTEGER          ICLEB(NCLEB,4),IFUNM(NATYPD,LMPOTD),
     +                 LMSP(NATYPD,*),IRMINSO,IRCUT(0:IPAND,NATYPD),
     +                 IPAN(NATYPD),NTCELL(*),IROTMAX
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   CLT,CONE,CZERO,LAMBDA,LINEAR,FAC1,FAC2,DIFF,
     +                 DENS_KR,DENS_SURF
      DOUBLE PRECISION C0LL,THETA,PHI
      INTEGER          I,IFUN,IR,J,LM1,LM2,LM3,LM1P,LM2P,ICELL,ENT,
     +                 NKSUM,NK,NKDUM,ENTDUM,I1,LM3P,ISP1,ISP2,
     +                 LMSP1,LMSP2,ISIGMA
      INTEGER          IRCUTM(0:IPAND),IPAN1,NK_IRREDUC
      LOGICAL          OPT
      CHARACTER (LEN=16)      ::   FILEN
      CHARACTER (LEN=2)       ::   NUM_LAY
      CHARACTER (LEN=1)       ::   ISPCHAR
      CHARACTER (LEN=18)      ::   FILENAME
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
      INTEGER, ALLOCATABLE         ::   ENTG(:)
      DOUBLE COMPLEX, ALLOCATABLE  ::   CKGES(:,:,:,:),ROTSPIN(:,:),
     +                                  RL1P12P2(:,:,:),CKHELP(:),
     +                                  DENS(:,:,:),RLL(:,:,:,:,:,:),
     +                                  DENS_GESAMT(:,:,:),
     +                                  DENS_GESAMT_I1(:,:,:,:),
     +                                  STOT(:,:),STOT_I1(:,:,:)!,DENSGES(:,:,:)
      DOUBLE PRECISION,ALLOCATABLE ::   KPOINT(:,:),DIREC_VEL_FS(:,:)
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA CONE/ (1.0D0,0.0D0)/
C     ..
c
      C0LL = 1.0d0/DSQRT(16.0D0*DATAN(1.0D0))
     

c  read in the normalized coefficients...

      WRITE(6,*) "in Sigma_natypd"
      OPEN(UNIT=149,FILE="coefficients_0",
     +         FORM="FORMATTED",STATUS="OLD", ACTION="READ")

      READ(149,"(I9)") NKSUM
      WRITE(6,*) "NKSUM:",NKSUM 

      IF (INTERF .EQ. .TRUE.) THEN
        NK_IRREDUC=NKSUM/IROTMAX    !  number of k points in the irreducible
                             !  part of the BZ
      ELSE
        NK_IRREDUC=NKSUM/48    
      END IF

      ALLOCATE(KPOINT(3,NKSUM))
      ALLOCATE(CKGES(2*LMMAXD,NATYPD,2,NKSUM))
      ALLOCATE(CKHELP(2*LMMAXD))
      ALLOCATE(ENTG(NKSUM))

      ALLOCATE(RLL(IRMD,LMMAX,LMMAX,2,2,NATYPD))
      ALLOCATE(RL1P12P2(IRMD,LMMAX,LMMAX))
      
      ALLOCATE(DENS(2*LMMAXD,2*LMMAXD,NATYPD))
      ALLOCATE(DENS_GESAMT(2,NKSUM,4))
      ALLOCATE(DENS_GESAMT_I1(NATYPD,2,NKSUM,4))

      ALLOCATE(STOT(2,NKSUM))
      ALLOCATE(ROTSPIN(LMMAXSO,LMMAXSO))

      WRITE(6,*) "k values eigenvec direct"
      WRITE(6,"(I5)") NKSUM 

      DENS_GESAMT=0d0
      DENS_GESAMT_I1=0d0

      DO NK=1,NKSUM
                      
        READ(149,"((I9),(3(E17.9)))") NKDUM,(KPOINT(J,NK), J=1,3)
        READ(149,"(I9)") ENTG(NK)

        DO ENT=1,ENTG(NK)
          READ(149,"(I5)") ENTDUM

          DO I1=1,NATYPD

            DO LM1 =1,LMMAXSO
              READ(149,"(2(E20.12,2X))") CKGES(LM1,I1,ENT,NK)
            END DO      !LM1=1,LMMAXD

          END DO        !I1=1,NATYP
        END DO          !ENT=1,ENTG
      END DO            !NK=1,NKSUM

      CLOSE(149)

      DO I1=1,NATYPD

        DO ISP1=1,2
          DO ISP2=1,2
            DO LM1 =1,LMMAXD
              LMSP1=(ISP1-1)*LMMAXD+LM1
              DO LM2 =1,LMMAXD
                LMSP2=(ISP2-1)*LMMAXD+LM2
                DO IR=IRMINSO,IRMD
                  RLL(IR,LM2,LM1,ISP2,ISP1,I1)=
     +                    PNS(LMSP2,LMSP1,IR,1,I1)
                END DO
              END DO      !LM1=1,LMMAXD
            END DO      !LM1=1,LMMAXD
          END DO      !ISP1=1,2
        END DO      !ISP1=1,2

      END DO        !I1=1,NATYP




c loop over the three expectation values M_x,y,z (given by the Pauli matrices)
c if sigma is 4, then simply the density is calculated 
c                              (-> test, whether the normalization is correct)

      DO ISIGMA=1,4
        WRITE(6,*) "ISIGMA",ISIGMA

c set up the array R*_L1L2 R_L3L4 

        DENS=0d0
       
        DO I1=1,NATYPD
       
          DO ISP2=1,2
            DO ISP1=1,2
       
              IF (ISIGMA==1) THEN
c             eigenvalue of the Pauli-Matrix sigma_x
                DO LM2 = 1,LMMAX
                  LMSP2=(ISP2-1)*LMMAXD+LM2

                  DO LM1 = 1,LMMAX
                    LMSP1=(ISP1-1)*LMMAXD+LM1
               
                    RL1P12P2=0d0
               
                    DO LM2P = 1,LMMAX
                      DO LM1P = 1,LMMAX
               
                        DO IR=1,IRMD
               
                          RL1P12P2(IR,LM1P,LM2P)=
     +                     CONJG(RLL(IR,LM1P,LM1,1,ISP1,I1))*
     +                              RLL(IR,LM2P,LM2,2,ISP2,I1)
     +                   + CONJG(RLL(IR,LM1P,LM1,2,ISP1,I1))*
     +                              RLL(IR,LM2P,LM2,1,ISP2,I1)
               
                        END DO       !IR
               
                      END DO         !LM1P
                    END DO           !LM2P
         
c calculate the 
                    CALL CALC_RHO_LL_SS(LMAX,LMMAX,RL1P12P2,
     +                       IRCUT(0:IPAND,I1),IPAN(I1),NTCELL(I1),
     +                       THETAS,CLEB,ICLEB,IEND,IFUNM,LMSP,IRMINSO,
     +                       IRWS(I1),DRDI(:,I1),DENS(LMSP1,LMSP2,I1))
         
                  END DO            !LM1 
                END DO              !LM2 
              ELSE IF (ISIGMA==2) THEN
c             eigenvalue of the Pauli-Matrix sigma_y
                DO LM2 = 1,LMMAX
                  LMSP2=(ISP2-1)*LMMAXD+LM2

                  DO LM1 = 1,LMMAX
                    LMSP1=(ISP1-1)*LMMAXD+LM1
               
c                WRITE(6,*) "LM1,LM2",LM1,LM2
             
                    RL1P12P2=0d0
             
                    DO LM2P = 1,LMMAX
                      DO LM1P = 1,LMMAX
             
                        DO IR=1,IRMD

                          RL1P12P2(IR,LM1P,LM2P)=(0d0,-1d0)*
     +                    (CONJG(RLL(IR,LM1P,LM1,1,ISP1,I1))*
     +                              RLL(IR,LM2P,LM2,2,ISP2,I1)
     +                   - CONJG(RLL(IR,LM1P,LM1,2,ISP1,I1))*
     +                              RLL(IR,LM2P,LM2,1,ISP2,I1))

                        END DO       !IR
             
                      END DO         !LM2P
                    END DO           !LM1P

c calculate the 
                    CALL CALC_RHO_LL_SS(LMAX,LMMAX,RL1P12P2,
     +                     IRCUT(0:IPAND,I1),IPAN(I1),NTCELL(I1),THETAS,
     +                     CLEB,ICLEB,IEND,IFUNM,LMSP,IRMINSO,IRWS(I1),
     +                     DRDI(:,I1),DENS(LMSP1,LMSP2,I1))

                  END DO            !LM2 
                END DO              !LM1 

              ELSE IF (ISIGMA==3) THEN
c             eigenvalue of the Pauli-Matrix sigma_z
                DO LM2 = 1,LMMAX
                  LMSP2=(ISP2-1)*LMMAXD+LM2

                  DO LM1 = 1,LMMAX
                    LMSP1=(ISP1-1)*LMMAXD+LM1
             
c              WRITE(6,*) "LM1,LM2",LM1,LM2
             
                    RL1P12P2=0d0
             
                    DO LM2P = 1,LMMAX
                      DO LM1P = 1,LMMAX
             
                        DO IR=1,IRMD
             
                          RL1P12P2(IR,LM1P,LM2P)=              
     +                      CONJG(RLL(IR,LM1P,LM1,1,ISP1,I1))*
     +                            RLL(IR,LM2P,LM2,1,ISP2,I1)
     +                    - CONJG(RLL(IR,LM1P,LM1,2,ISP1,I1))*
     +                            RLL(IR,LM2P,LM2,2,ISP2,I1)
             
                        END DO       !IR
             
                      END DO         !LM2P
                    END DO           !LM1P

c calculate the 
                    CALL CALC_RHO_LL_SS(LMAX,LMMAX,RL1P12P2,
     +                   IRCUT(0:IPAND,I1),IPAN(I1),NTCELL(I1),THETAS,
     +                   CLEB,ICLEB,IEND,IFUNM,LMSP,IRMINSO,IRWS(I1),
     +                   DRDI(:,I1),DENS(LMSP1,LMSP2,I1))

                  END DO            !LM2 
                END DO              !LM1 

              ELSE IF (ISIGMA==4) THEN
c             test of the right normalization...
                DO LM2 = 1,LMMAX
                  LMSP2=(ISP2-1)*LMMAXD+LM2

                  DO LM1 = 1,LMMAX
                    LMSP1=(ISP1-1)*LMMAXD+LM1
              
c                WRITE(6,*) "LM1,LM2",LM1,LM2
              
                    RL1P12P2=0d0
              
                    DO LM2P = 1,LMMAX
                      DO LM1P = 1,LMMAX
              
                        DO IR=1,IRMD
              
                          RL1P12P2(IR,LM1P,LM2P)=
     +                     CONJG(RLL(IR,LM1P,LM1,1,ISP1,I1))*
     +                              RLL(IR,LM2P,LM2,1,ISP2,I1)
     +                   + CONJG(RLL(IR,LM1P,LM1,2,ISP1,I1))*
     +                              RLL(IR,LM2P,LM2,2,ISP2,I1)
              
                        END DO       !IR
              
                      END DO         !LM2P
                    END DO           !LM1P

c calculate the 
                    CALL CALC_RHO_LL_SS(LMAX,LMMAX,RL1P12P2,
     +                     IRCUT(0:IPAND,I1),IPAN(I1),NTCELL(I1),THETAS,
     +                     CLEB,ICLEB,IEND,IFUNM,LMSP,IRMINSO,IRWS(I1),
     +                     DRDI(:,I1),DENS(LMSP1,LMSP2,I1))

                  END DO            !LM1 
                END DO              !LM2 
              END IF  ! SIGMA 1,2,3,4

            END DO         !ISP2
          END DO           !ISP1 

        END DO             !I1
        

        IF (ISIGMA == 1 ) THEN
          OPEN(UNIT=49,FILE="S_x",FORM='FORMATTED')
c          OPEN(UNIT=50,FILE="S_x_atom",FORM='FORMATTED')
        ELSE IF (ISIGMA == 2 ) THEN
          OPEN(UNIT=49,FILE="S_y",FORM='FORMATTED')
c          OPEN(UNIT=50,FILE="S_y_atom",FORM='FORMATTED')
        ELSE IF (ISIGMA == 3 ) THEN
          OPEN(UNIT=49,FILE="S_z",FORM='FORMATTED')
c          OPEN(UNIT=50,FILE="S_z_atom",FORM='FORMATTED')
        ELSE IF (ISIGMA == 4 ) THEN
          OPEN(UNIT=49,FILE="Dens_k_check",FORM='FORMATTED')
c          OPEN(UNIT=50,FILE="Dens_k_check_atom",FORM='FORMATTED')
        END IF
      
        WRITE(49,"(I5)") NKSUM 
c        WRITE(50,"(I5)") NKSUM 

c        DENS_GESAMT=0d0
c        DENS_GESAMT_I1=0d0

        DO NK=1,NKSUM                    
   
          DO ENT=1,ENTG(NK)
           
            DO I1=1,NATYPD
           
              CKHELP=0d0
              CALL ZGEMM('N','N',LMMAXSO,1,LMMAXSO,1d0,DENS(:,:,I1),
     +            LMMAXSO,CKGES(:,I1,ENT,NK),LMMAXSO,0d0,CKHELP,LMMAXSO)
            
            
              DENS_GESAMT_I1(I1,ENT,NK,ISIGMA)=
     +                            DOT_PRODUCT(CKGES(:,I1,ENT,NK),CKHELP)
              DENS_GESAMT(ENT,NK,ISIGMA)=
     +                    DENS_GESAMT(ENT,NK,ISIGMA)+
     +                                 DENS_GESAMT_I1(I1,ENT,NK,ISIGMA)

            END DO

          END DO

c         DO I1=1,NATYPD
c         
c           WRITE(50,"((2I9),(7(e17.9)))") NK,I1,
c    +           (KPOINT(J,NK), J=1,3), 
c    +       (DENS_GESAMT_I1(I1,ENT,NK),ENT=1,ENTG(NK))
c         
c         END DO             !I1

          WRITE(49,"((I9),(7(e17.9)))") NK,(KPOINT(J,NK),J=1,3),
     +                 (DENS_GESAMT(ENT,NK,ISIGMA),ENT=1,ENTG(NK))

        END DO                 !NK

c        CLOSE(50)
        CLOSE(49)

      END DO                   !ISIGMA=1,4


      STOT=0d0

      DO NK =1,NKSUM
        DO ENT=1,ENTG(NK)
          DO ISIGMA=1,3
            STOT(ENT,NK) =STOT(ENT,NK)+
     +        (DENS_GESAMT(ENT,NK,ISIGMA)**2) 
          END DO
          STOT(ENT,NK)=SQRT((STOT(ENT,NK)))
        END DO

c   rotate such that S_z is maximum... 
c   COS(THETA)= S_z / |S_tot|

c        THETA=ACOS(REAL(DENS_GESAMT(1,NK,3)/ABS(STOT(1,NK))))
        THETA=REAL(SQRT(DENS_GESAMT(1,NK,1)**2+DENS_GESAMT(1,NK,2)**2))
        THETA=ATAN2(THETA,REAL(DENS_GESAMT(1,NK,3)))
        PHI=ATAN2(REAL(DENS_GESAMT(1,NK,2)),REAL(DENS_GESAMT(1,NK,1)))
        WRITE(226,"((I5),(4e17.9))") NK,THETA,PHI

c   construct rotation matrix in spin space
        ROTSPIN=0d0
        DO LM1=1,LMMAX
c         ROTSPIN(LM1,LM1)=COS(THETA/2d0)        
c         ROTSPIN(LMMAX+LM1,LMMAX+LM1)=COS(THETA/2d0)        
c         ROTSPIN(LM1,LMMAX+LM1)=SIN(THETA/2d0)        
c         ROTSPIN(LMMAX+LM1,LM1)=-SIN(THETA/2d0)        
c         ROTSPIN(LM1,LM1)=COS(THETA/2d0)        
c         ROTSPIN(LMMAX+LM1,LMMAX+LM1)=COS(THETA/2d0)        
c         ROTSPIN(LM1,LMMAX+LM1)=-SIN(THETA/2d0)        
c         ROTSPIN(LMMAX+LM1,LM1)=SIN(THETA/2d0)        
          ROTSPIN(LM1,LM1)=COS(THETA/2d0)*EXP((0d0,0.5d0)*PHI)        
          ROTSPIN(LMMAX+LM1,LMMAX+LM1)=
     +                     COS(THETA/2d0)*EXP(-(0d0,0.5d0)*PHI)
c         ROTSPIN(LM1+LMMAX,LM1)=
c    +                    -SIN(THETA/2d0)*EXP(-(0d0,0.5d0)*PHI)        
c         ROTSPIN(LM1,LMMAX+LM1)=
c    +                     SIN(THETA/2d0)*EXP((0d0,0.5d0)*PHI)        
          ROTSPIN(LM1,LMMAX+LM1)=
     +                     SIN(THETA/2d0)*EXP(-(0d0,0.5d0)*PHI)        
          ROTSPIN(LMMAX+LM1,LM1)=
     +                    -SIN(THETA/2d0)*EXP((0d0,0.5d0)*PHI)        
c         ROTSPIN(LM1+LMMAX,LM1)=
c    +                    -SIN(THETA/2d0)*EXP(-(0d0,0.5d0)*PHI)        
c         ROTSPIN(LM1,LMMAX+LM1)=
c    +                     SIN(THETA/2d0)*EXP((0d0,0.5d0)*PHI)        
        END DO
        
        IF (NK==3) THEN
          DO LM2=1,LMMAXSO
            DO LM1=1,LMMAXSO
              WRITE(227,"((2I5),(4e17.9))") LM2,LM1,
     +           ROTSPIN(LM1,LM2),ROTSPIN(LM2,LM1) 
            END DO
          END DO
        END IF
c       ALLOCATE(CKGES(2*LMMAXD,NATYPD,2,NKSUM))

        DO I1=1,NATYPD
          DO ENT=1,ENTG(NK)
            CKHELP=0d0
            CALL ZGEMM('N','N',LMMAXSO,1,LMMAXSO,1d0,ROTSPIN,
     +            LMMAXSO,CKGES(:,I1,ENT,NK),LMMAXSO,0d0,CKHELP,LMMAXSO)

            DO LM1=1,LMMAXSO
              CKGES(LM1,I1,ENT,NK)=CKHELP(LM1)
            END DO

          END DO 
        END DO 
            
      END DO

c     write out rotated coefficients      
      OPEN(UNIT=49,FILE="coefficients_0_rot",FORM='FORMATTED')
      WRITE(49,"(I9)") NKSUM 

      DO NK=1,NKSUM                    

        WRITE(49,"((I9),(3(E17.9)))") NK,(KPOINT(J,NK), J=1,3)
        WRITE(49,"(I5)") ENTG(NK)

        DO ENT=1,ENTG(NK)
          WRITE(49,"(I5)") ENT
          DO I1=1,NATYPD
            DO LM1=1,LMMAXSO
              WRITE(49,"(4(E20.12,2x))") CKGES(LM1,I1,ENT,NK)
            END DO           !LM1
          END DO             !I1
        END DO               !ENT

      END DO                 !NK
      CLOSE(49)

      OPEN(UNIT=49,FILE="S_tot",FORM='FORMATTED')
      DO NK =1,NKSUM
        WRITE(49,"((I9),(7e17.9))") NK,(KPOINT(J,NK), J=1,3), 
     +                   (STOT(ENT,NK),ENT=1,ENTG(NK))
      END DO
      CLOSE(49)

      IF (NATYPD .NE. 1 ) THEN
        ALLOCATE(STOT_I1(2,NKSUM,NATYPD))
        STOT_I1=0d0

        DO I1=1,NATYPD
          DO NK =1,NKSUM
            DO ENT=1,ENTG(NK)
              DO ISIGMA=1,3
                STOT_I1(ENT,NK,I1) =STOT_I1(ENT,NK,I1)+
     +          (DENS_GESAMT_I1(I1,ENT,NK,ISIGMA)**2) 
              END DO
              STOT_I1(ENT,NK,I1)=SQRT((STOT_I1(ENT,NK,I1)))
            END DO
          END DO
        END DO

        OPEN(UNIT=49,FILE="S_tot_atom",FORM='FORMATTED')
        DO NK =1,NKSUM
          DO ENT=1,ENTG(NK)
            WRITE(49,"((I9),(7e17.9))") NK,(KPOINT(J,NK), J=1,3),ENT, 
     +                 (STOT_I1(ENT,NK,I1),I1=1,NATYPD)
          END DO
        END DO
        CLOSE(49)

        DEALLOCATE(STOT_I1)

      END IF 

      DEALLOCATE(RLL)
      DEALLOCATE(DENS)
      DEALLOCATE(DENS_GESAMT_I1)

      DEALLOCATE(ROTSPIN)

      DEALLOCATE(STOT)
      DEALLOCATE(CKGES)
      DEALLOCATE(ENTG)
      DEALLOCATE(DENS_GESAMT)
      DEALLOCATE(CKHELP)
      DEALLOCATE(RL1P12P2)

c     ALLOCATE(DIREC_VEL_FS(3,NK_IRREDUC))

c     CALL CREATE_GRID_FS_3D(KPOINT(:,1:NK_IRREDUC),NK_IRREDUC,
c    +                                       KPOINT,DIREC_VEL_FS)

      DEALLOCATE(KPOINT)
c      DEALLOCATE(DIREC_VEL_FS)

      END




