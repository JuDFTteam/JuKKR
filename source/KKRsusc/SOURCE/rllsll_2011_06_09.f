      MODULE MOD_RLLSLL
        CONTAINS
      SUBROUTINE RLLSLL(RINTV,RNEW,VLL,EIN,
     +                 C1,SRLLP,RLL,ULL,SLL,TLLP,
     +                 N,MMAX,LMAX,LMAXP1,LMMAX,LMMAX2,NMMAX,
     +                 use_fullgmat,nvec)
C radial wave functions by the integral equation method of
C Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
C Problem: the chosen signs lead to the matrix -A instead of A
C          therefore the alpha matrix is given (A-iB)^(-1) instead of
C          by -(A+iB)^(-1)
      USE MOD_TIMING
      USE MOD_BESHANK
      USE MOD_CHEBINT
      USE MOD_CONFIG, only: config_testflag
      USE MOD_RLLSLLTOOLS
      USE mod_physic_params,only: cvlight
      IMPLICIT NONE
C---------------------------------------------------------------------
C     .. Parameters ..
      INTEGER N,MMAX,LMAX,LMAXP1,LMMAX,LMMAX2
!       PARAMETER(N=16,MMAX=64,LMAX=4,LMAXP1=LMAX+1,LMMAX=(LMAX+1)**2)
C Explanation of parameters:
C LMAX = highest angular momentum
C N = order of the Chebyshev expansion
C MMAX = number of subintervals
C---------------------------------------------------------------------
      INTEGER NMMAX
      INTEGER IVEC, NVEC
!       PARAMETER (NMMAX= (N+1)*MMAX)
      DOUBLE COMPLEX CI,CONE,CZERO
      PARAMETER (CI= (0.0D0,1.0D0),CONE=(1.0D0,0.0D0))
      PARAMETER (CZERO=(0.0D0,0.0D0))
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX ERYD,EIN,EK,EK2
      DOUBLE COMPLEX AU,BU,AS,BS,F1U,F2U,F1S,F2S
      DOUBLE PRECISION PLLM
      INTEGER INFO,J,K,K1,M,MN,NM,NPLM
      INTEGER L1,LM1,LM2,LM3
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX HLK(0:(LMAX+1)*NVEC-1,0:N),
     +               JLK(0:(LMAX+1)*NVEC-1,0:N),
     +               NLK(0:(LMAX+1)*NVEC-1,0:N)
      DOUBLE COMPLEX 
     +               ULL(LMMAX2,LMMAX,NMMAX),SLL(LMMAX2,LMMAX,NMMAX),
     +               RLL(LMMAX2,LMMAX,NMMAX), TLLP(LMMAX,LMMAX),
     +               SRLLP(LMMAX,LMMAX),
     +               VLL(LMMAX*nvec,LMMAX*nvec,NMMAX)

      DOUBLE COMPLEX,allocatable :: 
!      +               SRLLP(:,:),
     +               WORK(:,:),
     +               ALLP(:,:,:),BLLP(:,:,:),
     +               CLLP(:,:,:),DLLP(:,:,:),
     +               SLV(:,:,:,:),SRV(:,:,:,:),
     +               MRNVY(:,:,:),MRNVZ(:,:,:),
     +               MRJVY(:,:,:),MRJVZ(:,:,:),
     +               MIHVY(:,:,:),MIHVZ(:,:,:),
     +               MIJVY(:,:,:),MIJVZ(:,:,:),
     +               YILL(:,:,:),ZILL(:,:,:),
     +               YRLL(:,:,:),ZRLL(:,:,:),
!      +               ULL(:,:,:),SLL(:,:,:),
!      +               RLL(:,:,:), !HLL(LMMAX,LMMAX,NMMAX),
     +               VNL(:,:,:),VHL(:,:,:),
     +               VJL(:,:,:)
!      +               VLL(:,:,:)
      DOUBLE COMPLEX,allocatable :: YIF(:,:,:,:),
     +               YRF(:,:,:,:),
     +               ZIF(:,:,:,:),
     +               ZRF(:,:,:,:)
      DOUBLE COMPLEX ZSLC1SUM(0:N)
      DOUBLE PRECISION C1(0:N,0:N),RINTV(0:MMAX)
      DOUBLE PRECISION CSLC1(0:N,0:N),CSRC1(0:N,0:N),
     +                 TAU(0:N,0:MMAX),
     +                 SLC1SUM(0:N),RNEW(NMMAX)
      INTEGER LOFLM(LMMAX2),IPIV(0:N,LMMAX2)
      LOGICAL TEST

      INTEGER use_fullgmat,ISPINfullgmat

!       DOUBLE PRECISION,parameter :: CVLIGHT = 274.0720442D0
C     ..
C     .. External Subroutines ..
      EXTERNAL ZGETRF,ZGETRS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,COS,DIMAG,EXP,MAX,MIN,SIN,SQRT
C     ..
         call timing_start('rllsll')
!           call timing_start('prestuff2')

      ALLOCATE( 
     +          WORK(LMMAX,LMMAX),
     +          ALLP(LMMAX,LMMAX,0:MMAX),BLLP(LMMAX,LMMAX,0:MMAX),
     +          CLLP(LMMAX,LMMAX,0:MMAX),DLLP(LMMAX,LMMAX,0:MMAX),
     +          SLV(0:N,LMMAX2,0:N,LMMAX2),SRV(0:N,LMMAX2,0:N,LMMAX2),
     +          MRNVY(LMMAX,LMMAX,MMAX),MRNVZ(LMMAX,LMMAX,MMAX),
     +          MRJVY(LMMAX,LMMAX,MMAX),MRJVZ(LMMAX,LMMAX,MMAX),
     +          MIHVY(LMMAX,LMMAX,MMAX),MIHVZ(LMMAX,LMMAX,MMAX),
     +          MIJVY(LMMAX,LMMAX,MMAX),MIJVZ(LMMAX,LMMAX,MMAX),
     +          YILL(0:N,LMMAX2,LMMAX),ZILL(0:N,LMMAX2,LMMAX),
     +          YRLL(0:N,LMMAX2,LMMAX),ZRLL(0:N,LMMAX2,LMMAX),
!      +               ULL(LMMAX,LMMAX,NMMAX),SLL(LMMAX,LMMAX,NMMAX),
!      +               RLL(LMMAX,LMMAX,NMMAX), !HLL(LMMAX,LMMAX,NMMAX),
     +          VNL(LMMAX,LMMAX2,0:N),VHL(LMMAX,LMMAX2,0:N),
     +          VJL(LMMAX,LMMAX2,0:N))


      ALLOCATE(
     +               YIF(LMMAX2,LMMAX,0:N,MMAX),
     +               YRF(LMMAX2,LMMAX,0:N,MMAX),
     +               ZIF(LMMAX2,LMMAX,0:N,MMAX),
     +               ZRF(LMMAX2,LMMAX,0:N,MMAX) )

      ERYD = EIN

!       LM1 = 1
!       DO L1 = 0,LMAX
!         DO M = -L1,L1
!           LOFLM(LM1) = L1
!           LM1 = LM1 + 1
!         END DO   
!       END DO  
!          write(*,*) nvec,lmmax2,lmmax
!           stop
C---------------------------------------------------------------------
C Bauer added for Spinorbit Coupling
C LOFLM is an (LMAX+1)**2 array. LOFLM = (/ LOFLM, LOFLM /) 
C---------------------------------------------------------------------

! !       VLL=czero
      LM1 = 1
      DO IVEC=1,NVEC
        DO ISPINfullgmat=0,use_fullgmat
          DO L1 = 0,LMAX
            DO M = -L1,L1
              LOFLM(LM1) = L1+(IVEC-1)*(LMAX+1)
!               print *, lm1,loflm(lm1)
              LM1 = LM1 + 1
            END DO   
          END DO  
        END DO!ISPINORBIT=0,use_fullgmat
      END DO !NVEC



C---------------------------------------------------------------------
      DO M = 1,MMAX
        DO K = 0,N
          MN = M*N + M - K
          TAU(K,M) = RNEW(MN)
        END DO
      END DO

      IF (nvec==1) THEN 
        EK = SQRT(ERYD)
        EK2 = SQRT(ERYD)
      ELSEIF (nvec==2) THEN
        EK = SQRT(ERYD+(ERYD/CVLIGHT)**2)
        EK2 = SQRT(ERYD+(ERYD/CVLIGHT)**2) *(1.0D0+ERYD/CVLIGHT**2)
      ELSE
        stop'[rllsll] wrong value for nvec'
      END IF

C---------------------------------------------------------------------
      CALL CHEBINT(CSLC1,CSRC1,SLC1SUM,C1,N)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call timing_start('local')

C loop over subintervals
      DO M = 1,MMAX
C initialization
!       call timing_start('local1')
        DO LM2 = 1,LMMAX
          DO LM1 = 1,LMMAX2
            DO K = 0,N
              YRLL(K,LM1,LM2) = CZERO
              ZRLL(K,LM1,LM2) = CZERO
              YILL(K,LM1,LM2) = CZERO
              ZILL(K,LM1,LM2) = CZERO
            END DO
          END DO
        END DO
        DO LM2 = 1,LMMAX2
          DO LM1 = 1,LMMAX2
            DO K = 0,N
              DO K1 = 0,N
                SLV(K1,LM1,K,LM2) = 0.0D0
                SRV(K1,LM1,K,LM2) = 0.0D0
              END DO
            END DO
          END DO
        END DO
        DO LM1 = 1,LMMAX2
          DO K = 0,N
            SLV(K,LM1,K,LM1) = 1.0D0
            SRV(K,LM1,K,LM1) = 1.0D0
          END DO
        END DO
        VNL=czero
        VHL=czero
        VJL=czero
C---------------------------------------------------------------------
C 1. prepare VJL, VNL, VHL, which appear in the integrands
C TAU(K,M) is used instead of TAU(K,M)**2, which directly gives
C RLL(r) and SLL(r) multiplied with r
C
C 2. prepare the source terms YR, ZR, YI, ZI
C because of the conventions used by
C Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
C a factor sqrt(E) is included in the source terms
C this factor is removed by the definition of ZSLC1SUM given below
C
        DO K = 0,N
          MN = M*N + M - K
          CALL BESHANK(HLK(0,K),JLK(0,K),NLK(0,K),EK*TAU(K,M),LMAX)

          IF (NVEC==2) THEN
            CALL BESHANK_SMALLCOMP(HLK(0,K),JLK(0,K),NLK(0,K),
     +                        EK*TAU(K,M),TAU(K,M),ERYD,LMAX)
          END IF

!           write(6546,'(5000E)') TAU(K,M),EK*TAU(K,M),hlk(:,k) 
!           write(6547,'(5000E)') TAU(K,M),EK*TAU(K,M),jlk(:,k) 
!           write(6548,'(5000E)') TAU(K,M),EK*TAU(K,M),nlk(:,k) 


!           IF (NVEC/=2) THEN
            DO L1 = 0,NVEC*(LMAX+1)-1
              HLK(L1,K) = -CI*HLK(L1,K)
            END DO
!           END IF

          DO IVEC=1,NVEC
            DO LM2 = 1,LMMAX2
              DO LM1 = 1,LMMAX
                L1 = LOFLM( LM1+LMMAX*(IVEC-1) )
!                  print *, VJL(LM1,LM2,K)
                VJL(LM1,LM2,K) = VJL(LM1,LM2,K) +
     +             EK2*TAU(K,M)*JLK(L1,K)*VLL(LM1+LMMAX*(IVEC-1),LM2,MN)
                VNL(LM1,LM2,K) = VNL(LM1,LM2,K) +
     +            -EK2*TAU(K,M)*NLK(L1,K)*VLL(LM1+LMMAX*(IVEC-1),LM2,MN)
                VHL(LM1,LM2,K) = VHL(LM1,LM2,K) +
     +            -EK2*TAU(K,M)*HLK(L1,K)*VLL(LM1+LMMAX*(IVEC-1),LM2,MN)

              END DO
            END DO
          END DO !NVEC
!           DO IVEC=1,NVEC
!             DO LM2 = 1,LMMAX
!               DO LM1 = 1,LMMAX
!                 L1 = LOFLM( LM1+LMMAX*(IVEC-1) )
!                 VJL(LM1,LM2+LMMAX*(IVEC-1),K) =
!      +             EK*TAU(K,M)*JLK(L1,K)*VLL(LM1,LM2,MN)
!                 VNL(LM1,LM2+LMMAX*(IVEC-1),K) = 
!      +            -EK*TAU(K,M)*NLK(L1,K)*VLL(LM1,LM2,MN)
!                 VHL(LM1,LM2+LMMAX*(IVEC-1),K) = 
!      +            -EK*TAU(K,M)*HLK(L1,K)*VLL(LM1,LM2,MN)
!               END DO
!             END DO
!           END DO !IVEC=1,NVEC


          DO IVEC=1,NVEC
            DO LM1 = 1,LMMAX
              L1 = LOFLM( LM1+LMMAX*(IVEC-1) )
              YRLL(K,LM1+LMMAX*(IVEC-1),LM1) =  EK*TAU(K,M)*JLK(L1,K) 
              ZRLL(K,LM1+LMMAX*(IVEC-1),LM1) = -EK*TAU(K,M)*NLK(L1,K) 
              YILL(K,LM1+LMMAX*(IVEC-1),LM1) = -EK*TAU(K,M)*HLK(L1,K)
              ZILL(K,LM1+LMMAX*(IVEC-1),LM1) =  EK*TAU(K,M)*JLK(L1,K)
            END DO
          END DO !IVEC=1,NVEC
        END DO
C-------------------------------------------------------
C determine the matrices in equations (4.5a) and (4.5b)

        DO J = 0,N
          DO K = 0,N
            DO LM2 = 1,LMMAX2
              DO IVEC=1,NVEC
                DO LM3 = 1,LMMAX
                  LM1=LM3+(IVEC-1)*LMMAX
                  L1 = LOFLM(LM1)
                  SLV(K,LM1,J,LM2) = 
     +                   (-TAU(K,M)*JLK(L1,K)*CSLC1(K,J)*VNL(LM3,LM2,J)
     +                    -TAU(K,M)*NLK(L1,K)*CSLC1(K,J)*VJL(LM3,LM2,J))
     +                   *(RINTV(M)-RINTV(M-1))/ 2.D0
                  SRV(K,LM1,J,LM2) =
     +                    (TAU(K,M)*JLK(L1,K)*CSRC1(K,J)*VNL(LM3,LM2,J)
     +                    +TAU(K,M)*NLK(L1,K)*CSRC1(K,J)*VJL(LM3,LM2,J))
     +                   *(RINTV(M)-RINTV(M-1))/ 2.D0
                END DO
              END DO
            END DO
          END DO
        END DO
        DO LM1 = 1,LMMAX2
          DO J = 0,N
            SLV(J,LM1,J,LM1) = SLV(J,LM1,J,LM1) + 1.D0
            SRV(J,LM1,J,LM1) = SRV(J,LM1,J,LM1) + 1.D0
          END DO
        END DO
!         call timing_stop('local1')
!         call timing_start('local2')

C-------------------------------------------------------
C determine the local solutions
C solve the equations SLV*YRLL=S and SLV*ZRLL=C 
C                 and SRV*YILL=C and SRV*ZILL=S
!         DO LM1 = 1,LMMAX
!         DO J = 0,N
!                write(453,'(5000F)') SLV(:,:,J,LM1)
!         END DO
!         END DO

!         DO LM1 = 1,LMMAX2
!           DO J = 0,N
!         write(3883,'(5000E)') slv(:,:,j,lm1)
!         end do
!         end do
!         stop
        NPLM = (N+1)*LMMAX2
        IF (.true.) THEN
          CALL ZGETRF(NPLM,NPLM,SLV,NPLM,IPIV,INFO)
          CALL ZGETRS('N',NPLM,LMMAX,SLV,NPLM,IPIV,YRLL,NPLM,INFO)
          CALL ZGETRS('N',NPLM,LMMAX,SLV,NPLM,IPIV,ZRLL,NPLM,INFO)
          CALL ZGETRF(NPLM,NPLM,SRV,NPLM,IPIV,INFO)
          CALL ZGETRS('N',NPLM,LMMAX,SRV,NPLM,IPIV,YILL,NPLM,INFO)
          CALL ZGETRS('N',NPLM,LMMAX,SRV,NPLM,IPIV,ZILL,NPLM,INFO)
        ELSE
          call rllslltools(NPLM,n+1,lmmax,SLV,YRLL,ZRLL)
          call rllslltools(NPLM,n+1,lmmax,SRV,YILL,ZILL)
        END IF

!         DO LM1 = 1,LMMAX
!         DO J = 0,N
!                write(453,'(5000F)') YRLL(:,:,LM1)
!         END DO
!         END DO
!           stop




C-------------------------------------------------------
        DO K = 0,N
          DO LM2 = 1,LMMAX
            DO LM1 = 1,LMMAX2
              YRF(LM1,LM2,K,M) = YRLL(K,LM1,LM2)
              ZRF(LM1,LM2,K,M) = ZRLL(K,LM1,LM2)
              YIF(LM1,LM2,K,M) = YILL(K,LM1,LM2)
              ZIF(LM1,LM2,K,M) = ZILL(K,LM1,LM2)
            END DO
          END DO
        END DO
C-------------------------------------------------------
C determine the left hand sides of equations (3.5), (3,6)
C (3.8), and (3.9) with the equations on page 143.
!       call timing_stop('local2')
!       call timing_start('local3')

      DO K = 0,N
        ZSLC1SUM(K) = SLC1SUM(K) * (RINTV(M)-RINTV(M-1))/ (2.D0*EK)
      END DO
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(0),VNL(1,1,0),
     +              LMMAX,YRF(1,1,0,M),LMMAX2,CZERO,MRNVY(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(0),VJL(1,1,0),
     +              LMMAX,YRF(1,1,0,M),LMMAX2,CZERO,MRJVY(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(0),VNL(1,1,0),
     +              LMMAX,ZRF(1,1,0,M),LMMAX2,CZERO,MRNVZ(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(0),VJL(1,1,0),
     +              LMMAX,ZRF(1,1,0,M),LMMAX2,CZERO,MRJVZ(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(0),VHL(1,1,0),
     +              LMMAX,YIF(1,1,0,M),LMMAX2,CZERO,MIHVY(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(0),VJL(1,1,0),
     +              LMMAX,YIF(1,1,0,M),LMMAX2,CZERO,MIJVY(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(0),VHL(1,1,0),
     +              LMMAX,ZIF(1,1,0,M),LMMAX2,CZERO,MIHVZ(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(0),VJL(1,1,0),
     +              LMMAX,ZIF(1,1,0,M),LMMAX2,CZERO,MIJVZ(1,1,M),LMMAX)
        DO K = 1,N
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(K),VNL(1,1,K),
     +              LMMAX,YRF(1,1,K,M),LMMAX2,CONE,MRNVY(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(K),VJL(1,1,K),
     +              LMMAX,YRF(1,1,K,M),LMMAX2,CONE,MRJVY(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(K),VNL(1,1,K),
     +              LMMAX,ZRF(1,1,K,M),LMMAX2,CONE,MRNVZ(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(K),VJL(1,1,K),
     +              LMMAX,ZRF(1,1,K,M),LMMAX2,CONE,MRJVZ(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(K),VHL(1,1,K),
     +              LMMAX,YIF(1,1,K,M),LMMAX2,CONE,MIHVY(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(K),VJL(1,1,K),
     +              LMMAX,YIF(1,1,K,M),LMMAX2,CONE,MIJVY(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(K),VHL(1,1,K),
     +              LMMAX,ZIF(1,1,K,M),LMMAX2,CONE,MIHVZ(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX2,ZSLC1SUM(K),VJL(1,1,K),
     +              LMMAX,ZIF(1,1,K,M),LMMAX2,CONE,MIJVZ(1,1,M),LMMAX)
        END DO
!         call timing_stop('local3')

      END DO
C end of loop over the subintervals

          call timing_stop('local')
           call timing_start('afterlocal')

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C calculate A(M), B(M), C(M), D(M) for m from 1 to MMAX
C starting from A(0) = 1, B(0) = 0, C(MMAX) = 0 and D(MMAX) = 1
      DO LM2 = 1,LMMAX
        DO LM1 = 1,LMMAX
          BLLP(LM1,LM2,0) = 0.D0
          ALLP(LM1,LM2,0) = 0.D0
        END DO
      END DO
      DO LM1 = 1,LMMAX
        ALLP(LM1,LM1,0) = 1.D0
      END DO
      DO M = 1,MMAX
        CALL ZCOPY(LMMAX*LMMAX,ALLP(1,1,M-1),1,ALLP(1,1,M),1)
        CALL ZCOPY(LMMAX*LMMAX,BLLP(1,1,M-1),1,BLLP(1,1,M),1)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,MRNVY(1,1,M),
     +              LMMAX,ALLP(1,1,M-1),LMMAX,CONE,ALLP(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,MRNVZ(1,1,M),
     +              LMMAX,BLLP(1,1,M-1),LMMAX,CONE,ALLP(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,-CONE,MRJVY(1,1,M),
     +              LMMAX,ALLP(1,1,M-1),LMMAX,CONE,BLLP(1,1,M),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,-CONE,MRJVZ(1,1,M),
     +              LMMAX,BLLP(1,1,M-1),LMMAX,CONE,BLLP(1,1,M),LMMAX)
      END DO
C
      DO LM2 = 1,LMMAX
        DO LM1 = 1,LMMAX
          DLLP(LM1,LM2,MMAX) = 0.D0
          CLLP(LM1,LM2,MMAX) = 0.D0
        END DO
      END DO
      DO LM1 = 1,LMMAX
        DLLP(LM1,LM1,MMAX) = 1.D0
      END DO
      DO M = MMAX,1,-1
        CALL ZCOPY(LMMAX*LMMAX,CLLP(1,1,M),1,CLLP(1,1,M-1),1)
        CALL ZCOPY(LMMAX*LMMAX,DLLP(1,1,M),1,DLLP(1,1,M-1),1)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,-CONE,MIHVZ(1,1,M),
     +              LMMAX,CLLP(1,1,M),LMMAX,CONE,CLLP(1,1,M-1),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,-CONE,MIHVY(1,1,M),
     +              LMMAX,DLLP(1,1,M),LMMAX,CONE,CLLP(1,1,M-1),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,MIJVZ(1,1,M),
     +              LMMAX,CLLP(1,1,M),LMMAX,CONE,DLLP(1,1,M-1),LMMAX)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,MIJVY(1,1,M),
     +              LMMAX,DLLP(1,1,M),LMMAX,CONE,DLLP(1,1,M-1),LMMAX)
      END DO
C---------------------------------------------------------------------
C determine the regular solution ULL and the irregular solution SLL
      DO M = 1,MMAX
        DO K = 0,N
          MN = M*N + M - K
        CALL ZGEMM('N','N',LMMAX2,LMMAX,LMMAX,CONE,YRF(1,1,K,M),
     +              LMMAX2,ALLP(1,1,M-1),LMMAX,CZERO,ULL(1,1,MN),LMMAX2)
        CALL ZGEMM('N','N',LMMAX2,LMMAX,LMMAX,CONE,ZRF(1,1,K,M),
     +              LMMAX2,BLLP(1,1,M-1),LMMAX,CONE,ULL(1,1,MN),LMMAX2)
        CALL ZGEMM('N','N',LMMAX2,LMMAX,LMMAX,CONE,ZIF(1,1,K,M),
     +              LMMAX2,CLLP(1,1,M),LMMAX,CZERO,SLL(1,1,MN),LMMAX2)
        CALL ZGEMM('N','N',LMMAX2,LMMAX,LMMAX,CONE,YIF(1,1,K,M),
     +              LMMAX2,DLLP(1,1,M),LMMAX,CONE,SLL(1,1,MN),LMMAX2)
        END DO
      END DO

          call timing_stop('afterlocal')
!             call timing_start('endstuff')

C-----------------------------------------------------------------
C replace regular wave function in the first subinterval by a
C linear function times r**(l+1)
C replace irregular wave function in the first subinterval by a
C linear function divided by r**l
      IF(config_testflag('wforigin')) THEN
      DO LM2 =1,LMMAX
      DO LM1 =1,LMMAX
      PLLM = 0.5D0*(LOFLM(LM1)+LOFLM(LM2))
      PLLM = LOFLM(LM1)
      F1U = ULL(LM1,LM2,2*N+2-N)/TAU(N,2)**(PLLM+1.D0)
      F2U = ULL(LM1,LM2,2*N+2-N/2)/TAU(N/2,2)**(PLLM+1.D0)
      AU = (F1U*TAU(N/2,2)-F2U*TAU(N,2))/(TAU(N/2,2)-TAU(N,2))
      BU = (F1U-F2U)/(TAU(N,2)-TAU(N/2,2))
      F1S = SLL(LM1,LM2,2*N+2-N)*TAU(N,2)**PLLM
      F2S = SLL(LM1,LM2,2*N+2-N/2)*TAU(N/2,2)**PLLM
      AS = (F1S*TAU(N/2,2)-F2S*TAU(N,2))/(TAU(N/2,2)-TAU(N,2))
      BS = (F1S-F2S)/(TAU(N,2)-TAU(N/2,2))
      DO K = 0,N
        MN = N + 1 - K
        ULL(LM1,LM2,MN) = (AU+BU*TAU(K,1))*TAU(K,1)**(PLLM+1.D0)
        SLL(LM1,LM2,MN) = (AS+BS*TAU(K,1))/TAU(K,1)**PLLM
      END DO
      END DO
      END DO
      END IF
C-----------------------------------------------------------------
C transform from Volterra solution to Fredholm solution
C calculate alpha and t matrices
      CALL ZAXPY(LMMAX*LMMAX,-CI,BLLP(1,1,MMAX),1,ALLP(1,1,MMAX),1)       ! calculate the transformation matrix alpha
                                                                          ! assuming A is calculated with a neuman function
                                                                          ! n=h+ij (?)
!David
!       DO NM = 1,NMMAX
!       CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,SLL(1,1,NM),
!      +            LMMAX,ALLP(1,1,MMAX),LMMAX,CZERO,HLL(1,1,NM),LMMAX)
!       END DO
! end David

      CALL ZGETRF(LMMAX,LMMAX,ALLP(1,1,MMAX),LMMAX,IPIV,INFO)              !invert alpha
      CALL ZGETRI(LMMAX,ALLP(1,1,MMAX),LMMAX,IPIV,WORK,LMMAX*LMMAX,INFO)   !invert alpha -> transformation matrix RLL=alpha^-1*RLL
      CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,-CONE/EK,BLLP(1,1,MMAX),        ! calc t-matrix TLL = BLL*alpha^-1
     +            LMMAX,ALLP(1,1,MMAX),LMMAX,CZERO,TLLP,LMMAX)
      DO LM2 = 1,LMMAX
      DO LM1 = 1,LMMAX
      SRLLP(LM1,LM2) = 2.D0*EK*TLLP(LM1,LM2)
      END DO
      SRLLP(LM2,LM2) = SRLLP(LM2,LM2) + CI
      END DO
      CALL ZGETRF(LMMAX,LMMAX,SRLLP,LMMAX,IPIV,INFO)
      CALL ZGETRI(LMMAX,SRLLP,LMMAX,IPIV,WORK,LMMAX*LMMAX,INFO)
C
      DO NM = 1,NMMAX
      CALL ZGEMM('N','N',LMMAX2,LMMAX,LMMAX,CONE,ULL(1,1,NM),
     +            LMMAX2,ALLP(1,1,MMAX),LMMAX,CZERO,RLL(1,1,NM),LMMAX2)
      END DO
!       call timing_stop('endstuff')
      call timing_stop('rllsll')


      RETURN
      END SUBROUTINE
      END MODULE MOD_RLLSLL
