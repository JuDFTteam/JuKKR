C ************************************************************************
      SUBROUTINE CHECK_REAL_RL(YR,PNS,I1,ALPHA,DAR,NSPD,NSP,LOFLM,LM2D,
     +                         LMAX,LMMAX,LMPOT,IRMIN,IRMAX,IJEND,RR)

      IMPLICIT NONE

C ************************************************************************
      INTEGER,INTENT(IN)            ::   LMMAX,IJEND,IRMIN,IRMAX,LMPOT,
     +                                   NSPD,I1,LMAX,NSP,
     +                                   LM2D,LOFLM(LM2D)
C     .. Array Arguments ..
      DOUBLE PRECISION,INTENT(IN)   ::   YR(IJEND,LMPOT),RR(IRMAX,*)
      DOUBLE COMPLEX,INTENT(IN)     ::   PNS(NSPD*LMMAX,NSPD*LMMAX,
     +                                                  IRMIN:IRMAX,2),
     +                                   DAR(NSPD*LMMAX,NSPD*LMMAX),
     +                                   ALPHA(0:LMAX,NSP)


C     .. Local Scalars ..
      INTEGER                       ::   LM1,LM2,NSRAN,IJ,IWF,ISP
      CHARACTER (LEN=4)             ::   FILEN
      CHARACTER (LEN=2)             ::   NUM_LAY
      CHARACTER (LEN=6)             ::   FILENAME

C     .. Local Arrays ..
      DOUBLE COMPLEX,ALLOCATABLE    ::   ARINV(:,:),
     +                                   URL(:,:,:),RRL(:,:,:),
     +                                   ALPHA_SP(:,:),AR(:,:)

C     .. External Subroutines ..
      EXTERNAL                  ::   ZGEMM,INVERS_MATRIX
C     ..
c
c multiply AR with the ALPHA matrix in order to obtain the full AR
c   ( the AR calculated in TMATRX01_SO is only delta AR ) 

      ALLOCATE(ALPHA_SP(NSPD*LMMAX,NSPD*LMMAX))
      ALLOCATE(AR(NSPD*LMMAX,NSPD*LMMAX))
      ALPHA_SP=0d0
      AR=0d0

      DO ISP=1,NSP
        DO LM1=1,LMMAX
          ALPHA_SP((ISP-1)*LMMAX+LM1,(ISP-1)*LMMAX+LM1)=
     +                                 ALPHA(LOFLM(LM1),ISP)
        END DO
      END DO

      IF (NSP==1) THEN
        DO LM1=1,LMMAX
          ALPHA_SP(LMMAX+LM1,LMMAX+LM1)=ALPHA(LOFLM(LM1),1)
        END DO
      END IF

      DO LM1=1,NSPD*LMMAX
        DO LM2=1,NSPD*LMMAX
          AR(LM2,LM1)=ALPHA_SP(LM2,LM2)*DAR(LM2,LM1)
          WRITE(123,"((2I5),(2e17.9))") LM2,LM1,AR(LM2,LM1)
        END DO
      END DO

      DEALLOCATE(ALPHA_SP)

c  check whether the AR transformation transforms the complex radial
c  wavefunctions to real wavefunctions
c        
c   U_L (r)=sum _L' AR^-1 _LL' R_L' (r)

c first calculate R_L'(vec r)=sum_L'' R_L''L'(r) Y_L''(vec r)

        FILEN="PNS_"
        WRITE(NUM_LAY,'(i2)') I1
        NUM_LAY=ADJUSTL(NUM_LAY)
        FILENAME=FILEN//NUM_LAY
        WRITE(6,*) "FILENAME:",FILENAME
    
        OPEN(UNIT=13,FILE=FILENAME,FORM='FORMATTED')

        DO IWF=IRMIN,IRMAX        !radial mesh 
          DO LM1 =1,NSPD*LMMAX
            WRITE(13,"((I5),(4e17.9))") LM1,RR(IWF,1),PNS(LM1,1,IWF,1)
          END DO
        END DO

        CLOSE(13)

        ALLOCATE(ARINV(NSPD*LMMAX,NSPD*LMMAX))
        ALLOCATE(URL(NSPD*LMMAX,IRMIN:IRMAX,IJEND))
        ALLOCATE(RRL(NSPD*LMMAX,IRMIN:IRMAX,IJEND))

        URL=0d0
        RRL=0d0
        ARINV=0d0

        DO IJ = 1,IJEND              !angular mesh
          DO IWF=IRMIN,IRMAX         !radial mesh 
            DO LM2 =1,NSPD*LMMAX
              DO LM1=1,NSPD*LMMAX
                RRL(LM2,IWF,IJ)=RRL(LM2,IWF,IJ)+
     +                PNS(LM1,LM2,IWF,1)*YR(IJ,LM1)
              END DO
            END DO
          END DO
        END DO

        FILEN="WF_C"
        WRITE(NUM_LAY,'(i2)') I1
        NUM_LAY=ADJUSTL(NUM_LAY)
        FILENAME=FILEN//NUM_LAY
        WRITE(6,*) "FILENAME:",FILENAME
    
        OPEN(UNIT=13,FILE=FILENAME,FORM='FORMATTED')

        DO IJ = 1,IJEND              !angular mesh
          DO IWF=IRMIN,IRMAX         !radial mesh 
            DO LM1 =1,NSPD*LMMAX
              WRITE(13,"((2I5),(4e17.9))") IJ,LM1,RR(IWF,1),
     +                                       RRL(LM1,IWF,IJ)
            END DO
          END DO
        END DO

        CLOSE(13)

        CALL INVERT_MATRIX(AR,ARINV,NSPD*LMMAX)

c       OPEN(UNIT=13,FILE="ARINV_CHECK",FORM="FORMATTED")

c       DO LM2=1,NSPD*LMMAX             !radial mesh 
c         DO LM1 =1,NSPD*LMMAX
c           WRITE(13,"((2I5),(4e17.9))") LM1,LM2,ARINV(LM1,LM2)     
c         END DO
c       END DO

c       CLOSE(13)


        FILEN="WF_R"
        WRITE(NUM_LAY,'(i2)') I1
        NUM_LAY=ADJUSTL(NUM_LAY)
        FILENAME=FILEN//NUM_LAY
        WRITE(6,*) "FILENAME:",FILENAME
        OPEN(UNIT=13,FILE=FILENAME,FORM='FORMATTED')

        DO IJ = 1,IJEND              !angular mesh
          DO IWF=IRMIN+1,IRMAX         !radial mesh 

            CALL ZGEMM('T','N',NSPD*LMMAX,1,NSPD*LMMAX,1d0,ARINV,
     +             NSPD*LMMAX,RRL(:,IWF,IJ),NSPD*LMMAX,
     +                          0d0,URL(:,IWF,IJ),NSPD*LMMAX)

            DO LM1 =1,NSPD*LMMAX
              WRITE(13,"((2I5),(4e17.9))") IJ,LM1,RR(IWF,1),
     +                                      URL(LM1,IWF,IJ)
            END DO
          END DO

        END DO

        CLOSE(13)

        DEALLOCATE(ARINV)
        DEALLOCATE(URL)
        DEALLOCATE(RRL)
        DEALLOCATE(AR)


        END SUBROUTINE CHECK_REAL_RL
