      SUBROUTINE  SCATTERING_MATRIX_2D_ANI(SCAT_MAT,NKPOID,NAEZ,LMAX,
     +   LMMAX,DTMATLM,DELTAMAT,GLL,COEFF_0,KPOINT,NCL,IMPLAYER,
     +   INDEX_LAYER,ISPIN,NSPO,NSPOH,ENTG,VOLBZ,WEIGHT,TAU_K_0,
     +   RN,ATOMIMP,ITHETA,IPHI,NTHETA,NPHI,DOSW,FERMI_VL)


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
      DOUBLE PRECISION                  TAU_0(NSPOH,NSPO),
     +                                  RN(3,NCL),
     +                                  LEFT(2,NKPOID),RIGHT(2,NKPOID)
      INTEGER                      ::   ATOMIMP(NCL)

c   local variables

      INTEGER                       ::  N_CL,J,KA,KSL,ENT,NK,
     +                                  LM1,LM2,NF,LMCLSO,LMMAXSO,
     +                                  ISH,ISIMP

      COMPLEX, ALLOCATABLE          ::  HILF_2(:,:),HILF_CK_1(:),
     +                                  HILF(:,:),HILF_CK_2(:),
     +                                  C_kLn(:,:,:),
c     +                                  SCATTER_COMPLEX_1(:,:,:,:),
     +                                  SCATTER_COMPLEX_2(:,:,:,:)
      COMPLEX                           FACTOR
      DOUBLE PRECISION              ::  DOSW,PI,RSLAB(3,NAEZ),
     +                                  CONTR_POINT,DIFF(3,NCL),
     +                                  FACT,TAU,T1
      INTEGER                       ::  ITHETA,IPHI,NTHETA,NPHI
      DOUBLE PRECISION              ::  x0,y0,z0,DTHETA,DPHI,PHI,THETA
    
      PI = 4.0D0*DATAN(1.0D0)
      FACT= 206.7084995*2D0*PI
      DTHETA=(PI/2d0)/(NTHETA-1)
      DPHI=(PI/2d0)/(NPHI-1)
      THETA=(ITHETA-1)*DTHETA
      PHI=(IPHI-1)*DPHI
      x0=COS(THETA)*COS(PHI)
      y0=COS(THETA)*SIN(PHI)
      z0=SIN(THETA)

c      LMMAXSO=LMMAX*NSPOH
c      LMCLSO=LMMAX*NCL*NSPOH
      LMMAXSO=LMMAX*NSPO
      LMCLSO=LMMAX*NCL*NSPO
      SCAT_MAT=0d0
      
c      ALLOCATE(SCATTER_COMPLEX_1(NKPOID,NSPOH,NKPOID,NSPO))
c      ALLOCATE(SCATTER_COMPLEX_2(NKPOID,NSPOH,NKPOID,NSPO))
c      ALLOCATE(SCATTER_COMPLEX_1(NKPOID,NSPOH,NKPOID,NSPO))
      ALLOCATE(SCATTER_COMPLEX_2(NKPOID,2,NKPOID,2))
      ALLOCATE(HILF(LMCLSO,LMCLSO))
      ALLOCATE(HILF_2(LMCLSO,LMCLSO))
      ALLOCATE(HILF_CK_1(LMCLSO))
      ALLOCATE(HILF_CK_2(LMCLSO))


c     for a cluster, the coefficients c_kL_n_0 have to be multiplied by a phase shift,

      ALLOCATE(C_kLn(LMCLSO,NSPO,NKPOID))

      C_kLn=0d0
      DO N_CL=1,NCL
c        write(6,*) N_CL,(RN(J,N_CL),J=1,3),ATOMIMP(N_CL)

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
c      write(6,*) "G_LL'* delta t_l'"

      CALL ZGEMM('N','N',LMCLSO,LMCLSO,LMCLSO,1d0,GLL,LMCLSO,DTMATLM,
     +                         LMCLSO,0d0,HILF,LMCLSO)

c     1 + G_LL'* delta t_l'
c      write(6,*) "1+ G_LL'* delta t_L'L''"

      do LM1=1,LMCLSO
        HILF(LM1,LM1)=1d0+HILF(LM1,LM1)
      end do

      HILF_2=0.d0

c      write(6,*) "Delta_LL' (1+G_L'L''* delta t_L''L''')"

      CALL ZGEMM('N','N',LMCLSO,LMCLSO,LMCLSO,1d0,DELTAMAT,LMCLSO,HILF,
     +                         LMCLSO,0d0,HILF_2,LMCLSO)

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
c                SCATTER_COMPLEX_1(NF,ISH,NK,ISIMP)=
c     +          DOT_PRODUCT(CONJG(HILF_CK_1),C_KLN(:,ISIMP,NF))

c   c_kL* ^0 * delta_LL'* (1 + G t) c_k'L'' ^0
                SCATTER_COMPLEX_2(NK,ISIMP,NF,ISH)=
     +          DOT_PRODUCT(C_KLN(:,ISH,NF),HILF_CK_2)

                SCAT_MAT(NK,ISIMP,NF,ISH)=
c     +                REAL(SCATTER_COMPLEX_1(NF,ISH,NK,ISIMP)
c     +                     *CONJG(SCATTER_COMPLEX_1(NF,ISH,NK,ISIMP)))
     +                REAL(SCATTER_COMPLEX_2(NK,ISIMP,NF,ISH)
     +                     *CONJG(SCATTER_COMPLEX_2(NK,ISIMP,NF,ISH)))
                

c integrate |T_kk'|**2             
             CONTR_POINT=WEIGHT(NF)*SCAT_MAT(NK,ISIMP,NF,ISH)
     +                   /FERMI_VL(NF)/VOLBZ
             TAU_K_0(NK,ISIMP,ISH)=TAU_K_0(NK,ISIMP,ISH)+CONTR_POINT
          END DO             !NF
         END DO               !ISIMP
         
      
          END DO             !NK

      END DO                 !ISH

      DO ISIMP=1,NSPO
       DO ISH=1,NSPOH
        TAU_0(ISIMP,ISH)=0d0
        DO NK=1,NKPOID
         CONTR_POINT=0d0
         CONTR_POINT=WEIGHT(NK)*TAU_K_0(NK,ISIMP,ISH)/FERMI_VL(NK)/VOLBZ
         TAU_0(ISIMP,ISH)=TAU_0(ISIMP,ISH)+CONTR_POINT
        ENDDO
         TAU_0(ISIMP,ISH)=TAU_0(ISIMP,ISH)/DOSW
       ENDDO
      ENDDO
        TAU=1d0/(TAU_0(1,1)+TAU_0(2,2))
        T1=1d0/(TAU_0(1,2)+TAU_0(2,1))
c        spin conserving lifetime in Ryd and ps unit
        WRITE(81,'(5E17.9)') x0,y0,z0,2D0*TAU,2D0*TAU/FACT
        IF (IPHI.EQ.NPHI) WRITE(81,*) ''
c        spin relaxation time in Ryd and ps unit
        WRITE(82,'(5E17.9)') x0,y0,z0,T1,T1/FACT
        IF (IPHI.EQ.NPHI) WRITE(82,*) ''
c        momentum relaxation time in Ryd and ps unit
        WRITE(83,'(5E17.9)') x0,y0,z0,2D0*1d0/(1d0/T1+1d0/TAU),
     +                                2D0*1d0/(1d0/T1+1D0/TAU)/FACT
        IF (IPHI.EQ.NPHI) WRITE(83,*) ''
         

c      DEALLOCATE(SCATTER_COMPLEX_1)
      DEALLOCATE(SCATTER_COMPLEX_2)
      DEALLOCATE(HILF)
      DEALLOCATE(HILF_2)
      DEALLOCATE(HILF_CK_1)
      DEALLOCATE(HILF_CK_2)
      DEALLOCATE(C_kLn)


      WRITE(6,*) "end of scattering matrix 2D"

      END SUBROUTINE SCATTERING_MATRIX_2D_ANI

