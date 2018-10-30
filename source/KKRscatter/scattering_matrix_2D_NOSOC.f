      SUBROUTINE  SCATTERING_MATRIX_2D_NOSOC(SCAT_MAT,NKPOID,NAEZ,LMAX,
     +   LMMAX,DTMATLM,DELTAMAT,GLL,COEFF_0,KPOINT,NCL,IMPLAYER,
     +   INDEX_LAYER,ISPIN,NSPO,NSPOH,ENTG,VOLBZ,WEIGHT,TAU_K_0,
     +   RN,ATOMIMP)


      IMPLICIT NONE

      INTEGER,INTENT(IN)   ::  NKPOID,LMAX,LMMAX,NSPO,NSPOH,NAEZ,ISPIN,
     +                         NCL,IMPLAYER,INDEX_LAYER(NCL),
     +                         ENTG(NKPOID)


      COMPLEX,INTENT(IN)   ::  COEFF_0(NSPOH*LMMAX,NAEZ,NSPOH,NKPOID),
     +                         GLL(NSPOH*LMMAX*NCL,NSPO*LMMAX*NCL),
     +                         DTMATLM(NSPOH*LMMAX*NCL,NSPOH*LMMAX*NCL),
     +                         DELTAMAT(NSPOH*LMMAX*NCL,NSPOH*LMMAX*NCL)


      DOUBLE PRECISION,INTENT(IN)   ::  VOLBZ,
     +                                  KPOINT(3,NKPOID),WEIGHT(NKPOID)

      DOUBLE PRECISION,INTENT(OUT) ::SCAT_MAT(NKPOID,NSPOH,NKPOID,NSPO),
     +                                  TAU_K_0(NKPOID,NSPOH,NSPO)
      DOUBLE PRECISION                  TAU_0(NSPOH,NSPO),
     +                                  RN(3,NCL),
     +                                  LEFT(2,NKPOID),RIGHT(2,NKPOID)
      INTEGER                      ::   ATOMIMP(NCL)

c   local variables

      INTEGER                       ::  N_CL,J,KA,KSL,ENT,NK,
     +                                  LM1,LM2,NF,LMCLSO,LMMAXSO,
     +                                  ISH,ISIMP,I

      COMPLEX, ALLOCATABLE          ::  HILF_2(:,:),HILF_CK_1(:),
     +                                  HILF(:,:),HILF_CK_2(:),
     +                                  C_kLn(:,:,:),
     +                                  SCATTER_COMPLEX_1(:,:,:,:),
     +                                  SCATTER_COMPLEX_2(:,:,:,:)
      COMPLEX                           FACTOR
      DOUBLE PRECISION              ::  COSTH,PI,RSLAB(3,NAEZ),
     +                                  CONTR_POINT,DIFF(3,NCL),
     +                                  BZKP(3,NKPOID),FERMI_VL(NKPOID),
     +                                  FACT,TAU,T1,CTR1(NSPOH,NSPO),
     +                                  V_FERMI(NKPOID,3),CONTR_POINT1,
     +                                  BZKP1(3,NKPOID),STOT(NKPOID,3)  
      DOUBLE PRECISION, ALLOCATABLE ::  LAMDA0(:,:),LAMDA1(:,:),
     +                                  LAMDA2(:,:)
      DOUBLE PRECISION              ::  DELTAPL,DELTAMN,SIGMA(3,3),
     +                                  SIGMAMN(3,3),SIGMASPL(3,3),
     +                                  SIGMASMN(3,3),HALLANGLE1,FACT2,
     +                                  HALLANGLE2,DELTA,MAXDT,FACT1,
     +                                  INTEGRAL
      INTEGER                       ::  ITR,NITR,ITRMAX,NKMAX
      
      PI = 4.0D0*DATAN(1.0D0)
      FACT= 206.7084995*2*PI
c       OPEN(unit=77,file="FERMI_VELOCITY",status="old",
c     +              action="read")
       OPEN(unit=78,file="vectorFERMIVL",status="old",
     +              action="read")
       DO NK=1,NKPOID 
c        READ(77,"(4(E17.9))")
c     +    (BZKP(J,NK),J=1,3),FERMI_VL(NK) 
        READ(78,"(3(E17.9))")
     +    (V_FERMI(NK,J),J=1,3) 
       FERMI_VL(NK)=sqrt(V_FERMI(NK,1)**2+V_FERMI(NK,2)**2+
     +                   V_FERMI(NK,3)**2)
       END DO 
c       CLOSE(77)
       CLOSE(78)

c      LMMAXSO=LMMAX*NSPOH
c      LMCLSO=LMMAX*NCL*NSPOH
      LMMAXSO=LMMAX*NSPO
      LMCLSO=LMMAX*NCL*NSPO
      SCAT_MAT=0d0
      
      ALLOCATE(SCATTER_COMPLEX_1(NKPOID,NSPOH,NKPOID,NSPO))
      ALLOCATE(SCATTER_COMPLEX_2(NKPOID,NSPOH,NKPOID,NSPO))
      ALLOCATE(HILF(LMCLSO,LMCLSO))
      ALLOCATE(HILF_2(LMCLSO,LMCLSO))
      ALLOCATE(HILF_CK_1(LMCLSO))
      ALLOCATE(HILF_CK_2(LMCLSO))


      write (6,*) "calculate the coefficients c_kL_n"
c     for a cluster, the coefficients c_kL_n_0 have to be multiplied by a phase shift,

      ALLOCATE(C_kLn(LMCLSO,NSPO,NKPOID))

      C_kLn=0d0
      DO N_CL=1,NCL
        write(6,*) N_CL,(RN(J,N_CL),J=1,3),ATOMIMP(N_CL)

        IF (NSPO == 2 .AND. NSPOH == 1) THEN
          DO NK=1,NKPOID          
                  FACTOR=exp((0d0,1d0)*2d0*PI*
     +             dot_product(KPOINT(:,NK),RN(:,N_CL)))
                  DO LM1=1,LMMAX
                    C_kLn(LM1+LMMAXSO*(N_CL-1),1,NK)=
     +               FACTOR*COEFF_0(LM1,ATOMIMP(N_CL),1,NK)
                    C_kLn(LM1+LMMAX+LMMAXSO*(N_CL-1),1,NK)=
     +               FACTOR*COEFF_0(LM1,ATOMIMP(N_CL),1,NK)
                  END DO 
          END DO
        ELSE   ! CASE SOC in host and imp or 
               ! No SOC in host and NO SOC in the imp
          DO NK=1,NKPOID          
              DO ENT=1,ENTG(NK)
                  FACTOR=exp((0d0,1d0)*2d0*PI*
     +             dot_product(KPOINT(:,NK),RN(:,N_CL)))
                  DO LM1=1,LMMAXSO
                    C_kLn(LM1+LMMAXSO*(N_CL-1),ENT,NK)=
     +               FACTOR*COEFF_0(LM1,ATOMIMP(N_CL),ENT,NK)
                  END DO 
              END DO
          END DO

        END IF

      END DO  ! 1,NCL

      HILF=0.d0

c     G_LL'* delta t_l'
      write(6,*) "G_LL'* delta t_l'"

      CALL ZGEMM('N','N',LMCLSO,LMCLSO,LMCLSO,1d0,GLL,LMCLSO,DTMATLM,
     +                         LMCLSO,0d0,HILF,LMCLSO)

c     1 + G_LL'* delta t_l'
      write(6,*) "1+ G_LL'* delta t_L'L''"

      do LM1=1,LMCLSO
        HILF(LM1,LM1)=1d0+HILF(LM1,LM1)
      end do

      HILF_2=0.d0

      write(6,*) "Delta_LL' (1+G_L'L''* delta t_L''L''')"

      CALL ZGEMM('N','N',LMCLSO,LMCLSO,LMCLSO,1d0,DELTAMAT,LMCLSO,HILF,
     +                         LMCLSO,0d0,HILF_2,LMCLSO)

      OPEN(unit=78,file="Tau_k_integrated",FORM="FORMATTED")

      DO ISH=1,NSPOH

         DO NK=1,NKPOID ! k points loop
          
          HILF_CK_1=0.d0
          HILF_CK_2=0.d0

c   c_kL* ^0 * delta_LL'* (1 + G t)
          CALL ZGEMM('N','N',1,LMCLSO,LMCLSO,1d0,
     +       conjg(C_kLn(:,ISH,NK)),1,
     +            HILF_2,LMCLSO,0d0,HILF_CK_1,1)

c  delta_LL'* (1 + G t)  c_kL ^0
c          CALL ZGEMM('N','N',LMCLSO,1,LMCLSO,1d0,
c     +         HILF_2,LMCLSO,
c     +         C_kLn(:,ISIMP,NK),LMCLSO,0d0,HILF_CK_2,LMCLSO)

       DO ISIMP=1,NSPO
          TAU_K_0(NK,ISH,ISIMP)=0d0
           DO NF=1,NKPOID
              CONTR_POINT=0d0

c   c_kL* ^0 * delta_LL'* (1 + G t) c_k'L'' ^0
                SCATTER_COMPLEX_1(NF,ISH,NK,ISIMP)=
     +          DOT_PRODUCT(CONJG(HILF_CK_1),C_KLN(:,ISIMP,NF))

c   c_kL* ^0 * delta_LL'* (1 + G t) c_k'L'' ^0
c                SCATTER_COMPLEX_2(NF,ISH,NK,ISIMP)=
c     +          DOT_PRODUCT(C_KLN(:,ISH,NF),HILF_CK_2)

                SCAT_MAT(NF,ISH,NK,ISIMP)=
     +                REAL(SCATTER_COMPLEX_1(NF,ISH,NK,ISIMP)
     +                     *CONJG(SCATTER_COMPLEX_1(NF,ISH,NK,ISIMP)))
c     +                REAL(SCATTER_COMPLEX_2*CONJG(SCATTER_COMPLEX_2))
                

c integrate |T_kk'|**2             
             CONTR_POINT=WEIGHT(NF)*SCAT_MAT(NF,ISH,NK,ISIMP)
     +                   /FERMI_VL(NF)/VOLBZ
             TAU_K_0(NK,ISH,ISIMP)=TAU_K_0(NK,ISH,ISIMP)+CONTR_POINT
          END DO             !NF
        WRITE(78,"(2I5,4E17.9)") ISH,ISIMP,
     +    (KPOINT(J,NK),J=1,3),TAU_K_0(NK,ISH,ISIMP) 
        END DO               !ISIMP
      
          END DO             !NK

      END DO                 !ISH
      CLOSE(78)
c optical-theorem ImT_kk^sigmasigma=-pi*sum_sigma'(tau_k_0^sigmasigma')
      OPEN(unit=88,file="Optical-theorem",FORM="FORMATTED")
      DO NK=1,NKPOID
      LEFT(1,NK)=DIMAG(SCATTER_COMPLEX_1(NK,1,NK,1))
      RIGHT(1,NK)=-PI*(TAU_K_0(NK,1,1))
      WRITE(88,'((I5),(3E17.9))') NK,LEFT(1,NK),RIGHT(1,NK),
     +                          LEFT(1,NK)-RIGHT(1,NK)
      ENDDO
      CLOSE(88)
      OPEN(unit=60,file='scatt_matrix',form='unformatted')
      DO NK=1,NKPOID
       DO NF=1,NKPOID
        WRITE(60) SCAT_MAT(NF,1,NK,1)
       ENDDO
      ENDDO
      CLOSE(60)
      STOP
c calculate conductivity and Hall angle by Boltzman equation
c first calculate mean free path by self-consistent equation
      OPEN(unit=90,file="meanfreepath",FORM="FORMATTED")
      OPEN(unit=91,file="conductivity",FORM="FORMATTED")
      ALLOCATE(LAMDA0(NKPOID,3))
      ALLOCATE(LAMDA1(NKPOID,3))
      ALLOCATE(LAMDA2(NKPOID,3))
      DO NK=1,NKPOID
       DO NF=NK+1,NKPOID
c         SCAT_MAT(NF,1,NK,1)=SCAT_MAT(NK,1,NF,1)
       ENDDO
      ENDDO
      NITR=50
      FACT1=2D0*PI
      FACT2=1D0
      DO J=1,1
        DO NK=1,NKPOID
         LAMDA0(NK,J)=0d0
          DO NF=1,NKPOID
           CONTR_POINT=0d0
           CONTR_POINT=FACT1*WEIGHT(NF)*SCAT_MAT(NF,1,NK,1)
     +                   /FERMI_VL(NF)/VOLBZ
c           CONTR_POINT=FACT1*WEIGHT(NF)*SCAT_MAT(NF,1,NK,1)
c     +                   /ABS(V_FERMI(NK,J))/VOLBZ
          LAMDA0(NK,J)=LAMDA0(NK,J)+CONTR_POINT
          ENDDO
          LAMDA0(NK,J)=LAMDA0(NK,J)*FACT2
          LAMDA1(NK,J)=(1D0/LAMDA0(NK,J))*V_FERMI(NK,J)
c          LAMDA1(NK,J)=1D0/LAMDA0(NK,J)
c         LAMDA1(NK,J)=0d0
        ENDDO
        DO ITR=1,NITR
        DO NK=1,NKPOID
          LAMDA2(NK,J)=0d0
          INTEGRAL=0d0
          DO NF=1,NKPOID
           CONTR_POINT=0d0
           CONTR_POINT=FACT1*WEIGHT(NF)*SCAT_MAT(NF,1,NK,1)*LAMDA1(NF,J)
     +                   /FERMI_VL(NF)/VOLBZ
c           CONTR_POINT=FACT1*WEIGHT(NF)*SCAT_MAT(NF,1,NK,1)*LAMDA1(NF,J)
c     +                   *V_FERMI(NK,J)*V_FERMI(NF,J)/ABS(V_FERMI(NK,J))
c     +                   /FERMI_VL(NF)/VOLBZ
c           CONTR_POINT=FACT1*WEIGHT(NF)*SCAT_MAT(NF,1,NK,1)*LAMDA1(NF,J)
c     +      *(V_FERMI(NK,1)*V_FERMI(NF,1)+V_FERMI(NK,2)*V_FERMI(NF,2)*
c     +                  V_FERMI(NK,3)*V_FERMI(NF,3))/FERMI_VL(NK)
c     +                   /FERMI_VL(NF)/VOLBZ
c           CONTR_POINT=FACT1*WEIGHT(NF)*SCAT_MAT(NF,1,NK,1)*LAMDA1(NF,J)
c     +                   /ABS(V_FERMI(NF,J))/VOLBZ
           INTEGRAL=INTEGRAL+CONTR_POINT
          ENDDO ! NF
           LAMDA2(NK,J)=(INTEGRAL+V_FERMI(NK,J))*(1D0/LAMDA0(NK,J))
c           LAMDA2(NK,J)=(LAMDA2(NK,J)+FERMI_VL(NK))*(1D0/LAMDA0(NK,J))*
c     +                    (1D0/FERMI_VL(NK))
c           LAMDA2(NK,J)=(LAMDA2(NK,J)*(V_FERMI(NK,J)**2/FERMI_VL(NK))+
c     _                         V_FERMI(NK,J))*(1D0/LAMDA0(NK,J))
        ENDDO ! NK
        MAXDT=1D-07
        DO NK=1,NKPOID
          DELTA=0d0
           DELTA=LAMDA2(NK,J)-LAMDA1(NK,J)
           IF (ABS(DELTA).GT.MAXDT) THEN
           MAXDT=ABS(DELTA)
           NKMAX=NK
           ENDIF
        IF (ITR.EQ.NITR) WRITE(55,*) NK,DELTA,LAMDA2(NK,J)
          LAMDA1(NK,J)=LAMDA2(NK,J)
        ENDDO
          WRITE(6,*) ITR,MAXDT,NKMAX
          IF (MAXDT.LE.1D-06) THEN
          ITRMAX=ITR
          GO TO 100
          ENDIF
        IF (ITR.EQ.NITR) THEN 
        WRITE(6,*) 'not converged',MAXDT
        STOP
        ENDIF
        ENDDO ! NITR
100     CONTINUE
      ENDDO ! dimension J
        DO NK=1,NKPOID
        WRITE(90,'((1I5),(3e17.9))') NK,(LAMDA1(NK,J),J=1,3)
        ENDDO
      DO I=1,3
c       DO J=1,3
         SIGMA(I,1)=0d0
          DO NK=1,NKPOID
           CONTR_POINT=0d0
c           CONTR_POINT=FACT1*WEIGHT(NK)*LAMDA1(NK,I)*V_FERMI(NK,J)
c     +                   /FERMI_VL(NK)/VOLBZ
           CONTR_POINT=FACT1*WEIGHT(NK)*LAMDA1(NK,I)
     +                   *FERMI_VL(NK)/VOLBZ
           SIGMA(I,1)=SIGMA(I,1)+CONTR_POINT
          ENDDO
        WRITE(91,*) I,1,SIGMA(I,1)
c       ENDDO
      ENDDO
            
      CLOSE(90)
      STOP
      OPEN(unit=81,file="Tau_k_averaged",FORM="FORMATTED")
      OPEN(unit=82,file="Relaxationtime",FORM="FORMATTED")
      DO ISH=1,NSPOH
       DO ISIMP=1,NSPO
        TAU_0(ISH,ISIMP)=0d0
        COSTH=0d0
        CTR1=0d0
        DO NK=1,NKPOID
         CONTR_POINT=0d0
         CONTR_POINT=WEIGHT(NK)*TAU_K_0(NK,ISH,ISIMP)/FERMI_VL(NK)/VOLBZ
         COSTH=COSTH+WEIGHT(NK)/FERMI_VL(NK)/VOLBZ
         CTR1(ISH,ISIMP)=CTR1(ISH,ISIMP)+WEIGHT(NK)*
     +              (1d0/TAU_K_0(NK,ISH,ISIMP))/FERMI_VL(NK)/VOLBZ
         TAU_0(ISH,ISIMP)=TAU_0(ISH,ISIMP)+CONTR_POINT
        ENDDO
         TAU_0(ISH,ISIMP)=TAU_0(ISH,ISIMP)/COSTH
        WRITE(81,'((2I5),(2E17.9))') ISH,ISIMP,TAU_0(ISH,ISIMP),
     +                               TAU_0(ISH,ISIMP)*FACT
       ENDDO
      ENDDO
        WRITE(6,*) 'DOS (normalization factor)',COSTH
        TAU=1d0/(TAU_0(1,1))
        WRITE(82,*) 'spin lifetimes in Ryd and ps unit'
        WRITE(82,'(2E17.9)') 2d0*TAU,2d0*TAU/FACT
         
      CLOSE(81)
      CLOSE(82)

        DEALLOCATE(LAMDA0)
        DEALLOCATE(LAMDA1)
        DEALLOCATE(LAMDA2)
        



      DEALLOCATE(SCATTER_COMPLEX_1)
      DEALLOCATE(SCATTER_COMPLEX_2)
      DEALLOCATE(HILF)
      DEALLOCATE(HILF_2)
      DEALLOCATE(HILF_CK_1)
      DEALLOCATE(HILF_CK_2)
      DEALLOCATE(C_kLn)


      WRITE(6,*) "end of scattering matrix 2D"

      END SUBROUTINE SCATTERING_MATRIX_2D_NOSOC

