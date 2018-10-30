      SUBROUTINE  SCATTERING_MATRIX_2D(SCAT_MAT,NKPOID,NAEZ,LMAX,LMMAX,
     +   DTMATLM,DELTAMAT,GLL,COEFF_0,KPOINT,NCL,IMPLAYER,
     +   INDEX_LAYER,ISPIN,NSPO,NSPOH,ENTG,VOLBZ,WEIGHT,TAU_K_0,
     +   RN,ATOMIMP,DOSW,FERMI_VL)


      IMPLICIT NONE

      INTEGER,INTENT(IN)   ::  NKPOID,LMAX,LMMAX,NSPO,NSPOH,NAEZ,ISPIN,
     +                         NCL,IMPLAYER,INDEX_LAYER(NCL),
     +                         ENTG(NKPOID)


      COMPLEX,INTENT(IN)   ::  COEFF_0(NSPOH*LMMAX,NAEZ,4,NKPOID),
     +                         GLL(NSPOH*LMMAX*NCL,NSPO*LMMAX*NCL),
     +                         DTMATLM(NSPOH*LMMAX*NCL,NSPOH*LMMAX*NCL),
     +                         DELTAMAT(NSPOH*LMMAX*NCL,NSPOH*LMMAX*NCL)


      DOUBLE PRECISION,INTENT(IN)   ::  VOLBZ,
     +                                  KPOINT(3,NKPOID),WEIGHT(NKPOID),
     +                                  FERMI_VL(NKPOID)

c      DOUBLE PRECISION,INTENT(OUT) ::SCAT_MAT(NKPOID,NSPOH,NKPOID,NSPO),
c     +                                  TAU_K_0(NKPOID,NSPOH,NSPO)
      DOUBLE PRECISION,INTENT(OUT) ::SCAT_MAT(NKPOID,2,NKPOID,2),
     +                                  TAU_K_0(NKPOID,2,2)
      DOUBLE PRECISION                  TAU_0(2,2),
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
c     +                                  SCATTER_COMPLEX_1(:,:,:,:),
     +                                  SCATTER_COMPLEX_2(:,:,:,:)
      COMPLEX                           FACTOR
      DOUBLE PRECISION              ::  DOSW,PI,RSLAB(3,NAEZ),
     +                                  CONTR_POINT,DIFF(3,NCL),
     +                                  FACT,TAU,T1,CTR1(NSPOH,NSPO),
     +                                  CONTR_POINT1
      DOUBLE PRECISION, ALLOCATABLE ::  V_FERMI(:,:),BZKP(:,:),Sz(:,:),
     +                                  LAMDAPL0(:),LAMDAMN0(:),
     +                                  LAMDAPL1(:),LAMDAMN1(:),
     +                                  LAMDAPL2(:,:),LAMDAMN2(:,:),
     +                                  DELTAPL(:),DELTAMN(:)
      DOUBLE PRECISION              ::  SIGMA(3,3),SIGMAS(3,3),
     +                                  ANGLE(2),MAXLD,CTR
      INTEGER                       ::  JDIM,NKMAX,NK1
      INTEGER                       ::  ITR,ITRMAX,NITR
      PARAMETER (NITR=500)
      LOGICAL                       ::  OPT
      EXTERNAL                      ::  OPT
      
      PI = 4.0D0*DATAN(1.0D0)
      FACT= 206.7084995*2*PI

c      LMMAXSO=LMMAX*NSPOH
c      LMCLSO=LMMAX*NCL*NSPOH
      LMMAXSO=LMMAX*NSPO
      LMCLSO=LMMAX*NCL*NSPO
      SCAT_MAT=0d0
      
c      ALLOCATE(SCATTER_COMPLEX_1(NKPOID,NSPOH,NKPOID,NSPO))
c      ALLOCATE(SCATTER_COMPLEX_2(NKPOID,NSPOH,NKPOID,NSPO))
c      ALLOCATE(SCATTER_COMPLEX_1(NKPOID,2,NKPOID,2))
      ALLOCATE(SCATTER_COMPLEX_2(NKPOID,2,NKPOID,2))
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
                    C_kLn(LM1+LMMAX+LMMAXSO*(N_CL-1),2,NK)=
     +               FACTOR*COEFF_0(LM1,ATOMIMP(N_CL),1,NK)
                  END DO 
          END DO
        ELSE   ! CASE SOC in host and imp or 
               ! No SOC in host and NO SOC in the imp
          DO NK=1,NKPOID          
c              DO ENT=1,ENTG(NK)
              DO ENT=1,2
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
      OPEN(unit=84,file="Tau_k_upup",FORM="FORMATTED")
      OPEN(unit=86,file="Tau_k_updown",FORM="FORMATTED")
      OPEN(unit=85,file="Tau_k_downup",FORM="FORMATTED")
      OPEN(unit=87,file="Tau_k_downdown",FORM="FORMATTED")
 
      DO ISIMP=1,NSPO

         DO NK=1,NKPOID ! k points loop

          
          HILF_CK_1=0.d0
          HILF_CK_2=0.d0

c   c_kL* ^0 * delta_LL'* (1 + G t)
c          CALL ZGEMM('N','N',1,LMCLSO,LMCLSO,1d0,
c     +       conjg(C_kLn(:,ISH,NK)),1,
c     +            HILF_2,LMCLSO,0d0,HILF_CK_1,1)

c  delta_LL'* (1 + G t)  c_kL ^0
          CALL ZGEMM('N','N',LMCLSO,1,LMCLSO,1d0,
     +         HILF_2,LMCLSO,
     +         C_kLn(:,ISIMP,NK),LMCLSO,0d0,HILF_CK_2,LMCLSO)

       DO ISH=1,NSPOH
          TAU_K_0(NK,ISIMP,ISH)=0d0
           DO NF=1,NKPOID
              CONTR_POINT=0d0

c   c_kL* ^0 * delta_LL'* (1 + G t) c_k'L'' ^0
c                SCATTER_COMPLEX_1(NK,ISH,NF,ISIMP)=
c     +          DOT_PRODUCT(CONJG(HILF_CK_1),C_KLN(:,ISIMP,NF))

c   c_kL* ^0 * delta_LL'* (1 + G t) c_k'L'' ^0
                SCATTER_COMPLEX_2(NK,ISIMP,NF,ISH)=
     +          DOT_PRODUCT(C_KLN(:,ISH,NF),HILF_CK_2)

                SCAT_MAT(NK,ISIMP,NF,ISH)=
c     +                REAL(SCATTER_COMPLEX_1(NK,ISH,NF,ISIMP)
c     +                     *CONJG(SCATTER_COMPLEX_1(NK,ISH,NF,ISIMP)))
     +                REAL(SCATTER_COMPLEX_2(NK,ISIMP,NF,ISH)
     +                     *CONJG(SCATTER_COMPLEX_2(NK,ISIMP,NF,ISH)))
                

c integrate |T_kk'|**2             
             CONTR_POINT=WEIGHT(NF)*SCAT_MAT(NK,ISIMP,NF,ISH)
     +                   /FERMI_VL(NF)/VOLBZ
             TAU_K_0(NK,ISIMP,ISH)=TAU_K_0(NK,ISIMP,ISH)+CONTR_POINT
          END DO             !NF
        WRITE(78,"(2I5,4E17.9)") ISIMP,ISH,
     +    (KPOINT(J,NK),J=1,3),TAU_K_0(NK,ISIMP,ISH) 
        END DO               !ISH
        WRITE(83+ISIMP,"(4E17.9)") 
     +    (KPOINT(J,NK),J=1,3),TAU_K_0(NK,ISIMP,1) 
        WRITE(85+ISIMP,"(4E17.9)") 
     +    (KPOINT(J,NK),J=1,3),TAU_K_0(NK,ISIMP,2) 
         
      
          END DO             !NK

      END DO                 !ISIMP
      CLOSE(78)
      CLOSE(84)
      CLOSE(85)
      CLOSE(86)
      CLOSE(87)
c calculate spin Hall conductivity
      IF (OPT('HALLCOND')) THEN
      ALLOCATE(BZKP(3,NKPOID))
      ALLOCATE(V_FERMI(3,NKPOID))
      ALLOCATE(Sz(3,NKPOID))
      ALLOCATE(LAMDAPL0(NKPOID))
      ALLOCATE(LAMDAMN0(NKPOID))
      ALLOCATE(LAMDAPL1(NKPOID))
      ALLOCATE(LAMDAMN1(NKPOID))
      ALLOCATE(LAMDAPL2(3,NKPOID))
      ALLOCATE(LAMDAMN2(3,NKPOID))
      ALLOCATE(DELTAPL(NKPOID))
      ALLOCATE(DELTAMN(NKPOID))
      OPEN(UNIT=70,FILE='vectorFERMIVL',form='formatted')
      OPEN(UNIT=71,FILE='S_xyz_alpha',form='formatted')
      OPEN(UNIT=72,FILE='conductivity',form='formatted')
      OPEN(UNIT=73,FILE='meanpath_x',form='formatted')
      OPEN(UNIT=74,FILE='meanpath_y',form='formatted')
      DO NK=1,NKPOID
       READ(70,'(3e17.9)') (V_FERMI(J,NK),J=1,3)
       READ(71,'((I9),(5e17.9))') NK1,(BZKP(J,NK),J=1,2),
     +                          (Sz(J,NK),J=1,3)
      ENDDO
      CLOSE(70)
      CLOSE(71)
c write down scattering matrix
c      OPEN(unit=61,file="scatt_matrix_upup",FORM="UNFORMATTED")
c      OPEN(unit=63,file="scatt_matrix_updn",FORM="UNFORMATTED")
c      OPEN(unit=62,file="scatt_matrix_dnup",FORM="UNFORMATTED")
c      OPEN(unit=64,file="scatt_matrix_dndn",FORM="UNFORMATTED")
c       DO ISIMP=1,2
c        DO NK=1,NKPOID
c         DO ISH=1,2
c            DO NF=1,NKPOID
c             IF (ISH.EQ.1) THEN
c             WRITE(60+ISIMP) 
c     +                 SCAT_MAT(NK,ISIMP,NF,ISH)
c             ELSE
c             WRITE(62+ISIMP) 
c     +                 SCAT_MAT(NK,ISIMP,NF,ISH)
c             ENDIF
c            ENDDO
c          ENDDO
c        ENDDO
c       ENDDO
c       CLOSE(61)
c       CLOSE(62)
c       CLOSE(63)
c       CLOSE(64)
       DO JDIM=1,2
        DO NK=1,NKPOID
         LAMDAPL0(NK)=0d0
         LAMDAMN0(NK)=0d0
         LAMDAPL1(NK)=0d0
         LAMDAMN1(NK)=0d0
         LAMDAPL0(NK)=TAU_K_0(NK,1,1)+TAU_K_0(NK,1,2)       
         LAMDAPL1(NK)=(1D0/LAMDAPL0(NK))*V_FERMI(JDIM,NK)
         LAMDAMN0(NK)=TAU_K_0(NK,2,2)+TAU_K_0(NK,2,1)       
         LAMDAMN1(NK)=(1D0/LAMDAMN0(NK))*V_FERMI(JDIM,NK)
        ENDDO
        DO ITR=1,NITR
         DO NK=1,NKPOID
          LAMDAPL2(JDIM,NK)=0d0
          LAMDAMN2(JDIM,NK)=0d0
          DO NF=1,NKPOID
           CTR=0d0
           CTR=WEIGHT(NF)*(SCAT_MAT(NF,1,NK,1)*LAMDAPL1(NF)+
     +                     SCAT_MAT(NF,2,NK,1)*LAMDAMN1(NF))
     +                     /FERMI_VL(NF)/VOLBZ
           LAMDAPL2(JDIM,NK)=LAMDAPL2(JDIM,NK)+CTR
           CTR=0d0
           CTR=WEIGHT(NF)*(SCAT_MAT(NF,2,NK,2)*LAMDAMN1(NF)+
     +                     SCAT_MAT(NF,1,NK,2)*LAMDAPL1(NF))
     +                     /FERMI_VL(NF)/VOLBZ
           LAMDAMN2(JDIM,NK)=LAMDAMN2(JDIM,NK)+CTR
          ENDDO
          LAMDAPL2(JDIM,NK)=(LAMDAPL2(JDIM,NK)+
     +                       V_FERMI(JDIM,NK))*(1D0/LAMDAPL0(NK))
          LAMDAMN2(JDIM,NK)=(LAMDAMN2(JDIM,NK)+
     +                       V_FERMI(JDIM,NK))*(1D0/LAMDAMN0(NK))
         ENDDO ! NK
         MAXLD=1D-07
         DO NK=1,NKPOID
          DELTAPL(NK)=LAMDAPL2(JDIM,NK)-LAMDAPL1(NK)
          DELTAMN(NK)=LAMDAMN2(JDIM,NK)-LAMDAMN1(NK)
           IF (ABS(DELTAPL(NK)).GT.MAXLD) THEN
            MAXLD=ABS(DELTAPL(NK))
            NKMAX=NK
           ENDIF
           IF (ABS(DELTAMN(NK)).GT.MAXLD) THEN
            MAXLD=ABS(DELTAMN(NK))
            NKMAX=NK
           ENDIF
          LAMDAPL1(NK)=LAMDAPL2(JDIM,NK)
          LAMDAMN1(NK)=LAMDAMN2(JDIM,NK)
         ENDDO ! NK
         WRITE(6,*) ITR,MAXLD,NKMAX
         IF (MAXLD.LE.1D-07) THEN
         ITRMAX=ITR
         GO TO 100
         ENDIF
        ENDDO ! ITR
100     CONTINUE
       DO NK=1,NKPOID
       WRITE(72+JDIM,'((I5),(2e17.9))') NK,
     +                        LAMDAPL2(JDIM,NK),LAMDAMN2(JDIM,NK)
       ENDDO
       ENDDO ! dimension JDIM
      WRITE(6,*) 'conductivity'
c calculate conductivity
      WRITE(72,*) 'conductivity'
      DO J=1,2
       DO I=1,2
        SIGMA(I,J)=0d0
         DO NK=1,NKPOID
          CTR=0d0
          CTR=WEIGHT(NK)*V_FERMI(I,NK)*
     +            (LAMDAPL2(J,NK)+LAMDAMN2(J,NK))
     +              /FERMI_VL(NK)/VOLBZ
          SIGMA(I,J)=SIGMA(I,J)+CTR
         ENDDO
        WRITE(72,*) I,J,SIGMA(I,J)
       ENDDO
       ENDDO
c calculate spin polarized conductivity
      WRITE(72,*) 'spin polarized conductivity'
      DO J=1,2
       DO I=1,2
        SIGMAS(I,J)=0d0
         DO NK=1,NKPOID
          CTR=0d0
          CTR=WEIGHT(NK)*V_FERMI(I,NK)*
     +             (LAMDAPL2(J,NK)-LAMDAMN2(J,NK))
     +              *Sz(3,NK)/FERMI_VL(NK)/VOLBZ
          SIGMAS(I,J)=SIGMAS(I,J)+CTR
         ENDDO
        WRITE(72,*) I,J,SIGMAS(I,J)
       ENDDO
       ENDDO
       WRITE(72,*) 'spin Hall angle'
       ANGLE(1)=SIGMAS(2,1)/SIGMA(1,1)
       ANGLE(2)=SIGMAS(1,2)/SIGMA(2,2)
       WRITE(72,*) ANGLE(1),ANGLE(2)
      CLOSE(72)
      CLOSE(73)
      CLOSE(74)
      DEALLOCATE(V_FERMI)
      DEALLOCATE(Sz)
      DEALLOCATE(BZKP)
      DEALLOCATE(LAMDAPL0)
      DEALLOCATE(LAMDAMN0)
      DEALLOCATE(LAMDAPL1)
      DEALLOCATE(LAMDAMN1)
      DEALLOCATE(LAMDAPL2)
      DEALLOCATE(LAMDAMN2)
      DEALLOCATE(DELTAPL)
      DEALLOCATE(DELTAMN)
      ENDIF
c       STOP
c optical-theorem ImT_kk^sigmasigma=-pi*sum_sigma'(tau_k_0^sigmasigma')
      OPEN(unit=88,file="Optical-theorem_up",FORM="FORMATTED")
      OPEN(unit=89,file="Optical-theorem_down",FORM="FORMATTED")
      DO NK=1,NKPOID
c      LEFT(1,NK)=DIMAG(SCATTER_COMPLEX_1(NK,1,NK,1))
c      LEFT(2,NK)=DIMAG(SCATTER_COMPLEX_1(NK,2,NK,2))
      LEFT(1,NK)=DIMAG(SCATTER_COMPLEX_2(NK,1,NK,1))
      LEFT(2,NK)=DIMAG(SCATTER_COMPLEX_2(NK,2,NK,2))
      RIGHT(1,NK)=-PI*(TAU_K_0(NK,1,1)+TAU_K_0(NK,1,2))
      RIGHT(2,NK)=-PI*(TAU_K_0(NK,2,2)+TAU_K_0(NK,2,1))
c      RIGHT(1,NK)=-PI*(TAU_K_0(NK,1,1)+TAU_K_0(NK,2,1))
c      RIGHT(2,NK)=-PI*(TAU_K_0(NK,2,2)+TAU_K_0(NK,1,2))
      WRITE(88,'((I5),(3E17.9))') NK,LEFT(1,NK),RIGHT(1,NK),
     +                          LEFT(1,NK)-RIGHT(1,NK)
      WRITE(89,'((I5),(3E17.9))') NK,LEFT(2,NK),RIGHT(2,NK),
     +                          LEFT(2,NK)-RIGHT(2,NK)
      ENDDO
      CLOSE(88)
      CLOSE(89)
      OPEN(unit=79,file="Tau_k_inverse",FORM="FORMATTED")
      OPEN(unit=80,file="T_k_inverse",FORM="FORMATTED")
      DO NK=1,NKPOID
        WRITE(79,"(4E17.9)") (KPOINT(J,NK),J=1,3),
     +     TAU_K_0(NK,1,1)+TAU_K_0(NK,2,2) 
        WRITE(80,"(4E17.9)") (KPOINT(J,NK),J=1,3),
     +     TAU_K_0(NK,1,2)+TAU_K_0(NK,2,1) 
      END DO
      CLOSE(79)
      CLOSE(80)
      OPEN(unit=81,file="Tau_k_averaged",FORM="FORMATTED")
      OPEN(unit=82,file="Relaxationtime",FORM="FORMATTED")
      DO ISIMP=1,2
       DO ISH=1,2
        TAU_0(ISIMP,ISH)=0d0
c        CTR1=0d0
        DO NK=1,NKPOID
         CONTR_POINT=0d0
         CONTR_POINT=WEIGHT(NK)*TAU_K_0(NK,ISIMP,ISH)/FERMI_VL(NK)/VOLBZ
c         CTR1(ISIMP,ISH)=CTR1(ISIMP,ISH)+WEIGHT(NK)*
c     +              (1d0/TAU_K_0(NK,ISIMP,ISH))/FERMI_VL(NK)/VOLBZ
         TAU_0(ISIMP,ISH)=TAU_0(ISIMP,ISH)+CONTR_POINT
        ENDDO
         TAU_0(ISIMP,ISH)=TAU_0(ISIMP,ISH)/DOSW
        WRITE(81,'((2I5),(2E17.9))') ISIMP,ISH,TAU_0(ISIMP,ISH),
     +                               TAU_0(ISIMP,ISH)*FACT
c        WRITE(55,'((2I5),(2E17.9))') ISH,ISIMP,1d0/CTR1(ISH,ISIMP),
c     +                               (1d0/CTR1(ISH,ISIMP))*FACT
       ENDDO
      ENDDO
        TAU=1d0/(TAU_0(1,1)+TAU_0(2,2))
        T1=1d0/(TAU_0(1,2)+TAU_0(2,1))
        WRITE(82,*) 'spin conserving lifetimes in Ryd and ps unit'
        WRITE(82,'(2E17.9)') 2d0*TAU,2d0*TAU/FACT
        WRITE(82,*) 'spin flip lifetimes in Ryd and ps unit'
        WRITE(82,'(2E17.9)') 2d0*T1,2d0*T1/FACT
        WRITE(82,*) 'spin relaxation times in Ryd and ps unit'
        WRITE(82,'(2E17.9)') T1,T1/FACT
        WRITE(82,*) 'momentum relaxation times in Ryd and ps unit'
        WRITE(82,'(2E17.9)') 2d0*1d0/(1d0/T1+1d0/TAU),
     +                       2d0*1d0/(1d0/T1+1d0/TAU)/FACT
        WRITE(82,*) 'spin flip rate/momentum relaxation rate'
        WRITE(82,'(1E17.9)') (1d0/(T1))/(1d0/(T1)+1d0/(TAU))
         
      CLOSE(81)
      CLOSE(82)

        



c      DEALLOCATE(SCATTER_COMPLEX_1)
      DEALLOCATE(SCATTER_COMPLEX_2)
      DEALLOCATE(HILF)
      DEALLOCATE(HILF_2)
      DEALLOCATE(HILF_CK_1)
      DEALLOCATE(HILF_CK_2)
      DEALLOCATE(C_kLn)


      WRITE(6,*) "end of scattering matrix 2D"

      END SUBROUTINE SCATTERING_MATRIX_2D

