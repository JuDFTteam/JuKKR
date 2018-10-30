c 28.9.99 ***************************************************************
      SUBROUTINE BANDSTRSO(EZ,LMMAXSO,NSPO,TMATLL,
     +                  IE,ALAT,
     +                  NAEZ,CLS,EQINV,NACLS,RR,EZOA,ATOM,KAOEZ,
     +                  GINP,RCLS,MY_RANK)
      implicit none
c ************************************************************************
c It calculates the band structure along symmetry directions.
c (optionally calculates also complex band structure)
c use OPT('COMPLEX')
c ------------------------------------------------------------------------
c     .. parameters ..
      include 'inc.p'
      include 'inc.cls'
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAX+1)**2)
      INTEGER ALM
      PARAMETER (ALM = NAEZD*LMAXSQ)
      INTEGER ALMSO
      PARAMETER (ALMSO = NAEZD*NSPD*LMAXSQ)
      INTEGER NKPOID
      PARAMETER(NKPOID=1000)
      INTEGER NAUX
      PARAMETER(NAUX = 2*ALMSO**2+5*ALMSO)
      DOUBLE PRECISION ZERO
      DOUBLE COMPLEX CI,CZERO,CONE,CONEM
      PARAMETER (ZERO  = 0.0D0)
      PARAMETER (CI=(0.D0,1.D0),CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0),
     +           CONEM=(-1d0,0d0))
      DOUBLE PRECISION FILTER_EIGVAL
      PARAMETER (FILTER_EIGVAL=1d0)
      
c     ..
c     .. scalar arguments ..
      DOUBLE COMPLEX EZ(IEMXD),FAC,NORM,NORMT
      DOUBLE PRECISION ALAT,RFCTOR,NORMABS,PHI,MINVALUE,DIFF
      INTEGER IC,IE,NAEZ
      INTEGER LMMAXSO,NSPO
c     ..
c     .. array arguments ..
      DOUBLE COMPLEX
     +     GINP(LMAXSQ*NACLSD,LMAXSQ,NCLSD),
     +     TMATLL(LMMAXSO,LMMAXSO,NAEZD)
      DOUBLE PRECISION
     +     RR(3,0:NRD),
     +     RCLS(3,NACLSD,NCLSD),
     +     DAUX(2*ALMSO)
      INTEGER
     +     ATOM(NACLSD,NAEZD),
     +     CLS(NATYPD),
     +     EQINV(NAEZD),
     +     EZOA(NACLSD,NAEZD),
     +     KAOEZ(NAEZD+NEMBD),
     +     NACLS(NCLSD)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,INFO,J,IL1,IL2,K,IS1,IS2,IL1T,IL2T,
     +        LM1,LM2,KSUM,I1,M,JN,IM,LM3,LM4,LMTAKE
C     .. LOCAL ARRAYS ..
      DOUBLE PRECISION XINTERSECT,REINTERSECT
      DOUBLE COMPLEX
     +     AUX(NAUX)
      DOUBLE PRECISION BZKP(6,NKPOID),DK(6,NKPOID),KPINT(6),
     +                 KP(6,NKPOID)
      INTEGER ENT(ALMSO),ENTW(ALMSO,ALMSO) 
      INTEGER MY_RANK
      LOGICAL TEST,OPT
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL CINIT,OPT,
     +         TEST,
     +         ZGETRF,ZGETRS
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DATAN,EXP
C
      DOUBLE COMPLEX, ALLOCATABLE :: GLLKE(:,:),GLLKE1(:,:),
     +                        GLLKEDK(:,:),DGLLKEDK(:,:),
     +                        HILF_1(:,:),HILF_2(:,:),
     +                        TMAT_1(:,:),TMAT_2(:,:),
     +                        RV(:,:),LV(:,:),W(:),GLLKES(:,:)
      DOUBLE COMPLEX       :: WAPR(ALMSO),WK(ALMSO)
c      
c     ..
c ------------------------------------------------------------------------
c

c
c     reading from the inputcard the points defining the symmetry
c     directions for a band structure calculation.

      CALL BANDINPUT(IE,KSUM,BZKP,NKPOID) 

      ALLOCATE(GLLKE1(NAEZ*LMAXSQ,LMAXSQ))
      ALLOCATE(GLLKE(ALMSO,ALMSO))
      ALLOCATE(GLLKES(ALMSO,ALMSO))
      ALLOCATE(GLLKEDK(ALMSO,ALMSO))
      ALLOCATE(DGLLKEDK(ALMSO,ALMSO))
      ALLOCATE(HILF_1(ALMSO,ALMSO))
      ALLOCATE(HILF_2(ALMSO,ALMSO))
      ALLOCATE(TMAT_1(ALMSO,ALMSO))
      ALLOCATE(TMAT_2(ALMSO,ALMSO))
      ALLOCATE(LV(ALMSO,ALMSO))
      ALLOCATE(RV(ALMSO,ALMSO))
      ALLOCATE(W(ALMSO))
             
      RFCTOR = ALAT/(8D0*DATAN(1D0))

c      IF (IE==1) THEN
c        WRITE(6,*) "DK"
c        WRITE(6,"(6e17.9)") (DK(J,K),J=1,6)
c      END IF

      DO 300 K = 1,KSUM-1         ! K-POINT-LOOP

      DO J=1,6
        DK(J,K)=BZKP(J,K+1)-BZKP(J,K)
      ENDDO

c        IF (IE.EQ.1) THEN
           IF (K.EQ.1) THEN 
c           WRITE(6,*) 
c     &'      Real       and         imaginary parts of k-point '
           ENDIF
c           WRITE (6,9090)k,(BZKP(J,K),J=1,3), (BZKP(J,K),J=4,6)
           WRITE (6,9090)MY_RANK+1,IE,k,(BZKP(J,K),J=1,3)
c        ENDIF
 9090   FORMAT('proc',I5,'  ener',I5,'  kpts',I5,3F10.5,3X,3F10.5)
c     
c --->  fourier transformation
c
        GLLKE=CZERO
        GLLKEDK=CZERO
        DGLLKEDK=CZERO
       
        DO I=1,NAEZ
         IC=CLS(KAOEZ(I))
         FAC=CONE
           IF (I.NE.EQINV(I)) FAC=CONEM
           CALL CINIT(NAEZD*LMAXSQ*LMAXSQ,GLLKE1)
           CALL DLKE1(GLLKE1,ALAT,NACLS(IC),RR,EZOA(1,I),ATOM(1,I),
     +           BZKP(1,K),IE,IC,FAC,GINP(1,1,IC),RCLS(1,1,IC))
           DO M=1,LMAXSQ
            IM=(I-1)*LMAXSQ+M
             DO JN=1,LMAXSQ*NAEZ
               GLLKE(JN,IM) = GLLKE(JN,IM)+GLLKE1(JN,M)
             ENDDO
           ENDDO
         ENDDO ! atom loop
c  for spin-orbit coupling, G_LL'(k) 
        IF (NSPD == 2) THEN
          DO LM1=1,ALM
            DO LM2=1,ALM
              GLLKE(ALM+LM2,ALM+LM1) = GLLKE(LM2,LM1) 
            END DO
          END DO
        END IF
        DO LM1=1,ALMSO
          DO LM2=1,ALMSO
            GLLKES(LM1,LM2) = GLLKE(LM1,LM2) 
          END DO
        END DO

c       the matrix M=[1 - G^r*delta t] is constructed
        IS1=1
        IS2=1
        TMAT_1=CZERO
        HILF_1=CZERO
        DO I1=1,NAEZ
          DO IS1=1,NSPD
            DO LM1=1,LMAXSQ
              DO IS2=1,NSPD
                DO LM2=1,LMAXSQ
                  IL1=(IS1-1)*ALM+LMAXSQ*(I1-1)+LM1
                  IL2=(IS2-1)*ALM+LMAXSQ*(I1-1)+LM2
                  IL1T=(IS1-1)*LMAXSQ+LM1
                  IL2T=(IS2-1)*LMAXSQ+LM2
                  TMAT_1(IL2,IL1)=TMATLL(IL2T,IL1T,I1)/RFCTOR
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
       
        CALL ZGEMM('N','N',ALMSO,ALMSO,ALMSO,CONE,GLLKES,ALMSO,
     +                           TMAT_1,ALMSO,CZERO,HILF_1,ALMSO)
        DO LM1=1,ALMSO
          DO LM2=1,ALMSO
            IF (LM1.EQ.LM2) THEN
              HILF_1(LM1,LM2)=CONE-HILF_1(LM1,LM2)
            ELSE
              HILF_1(LM1,LM2)=-HILF_1(LM1,LM2)
            ENDIF
          ENDDO
        ENDDO
c       eigenvalues determination
c       use lapack routine to calculate the eigenvalues of the
c       constructred matrix
c LV and RV are left or right eigen values; A*RV=lamda*RV;LV**H*A=lamda*LV**H 
        LV=CZERO
        RV=CZERO
        W=CZERO

        CALL ZGEEV('V','V',ALMSO,HILF_1,ALMSO,W,LV,ALMSO,RV,ALMSO,
     +              AUX,NAUX,DAUX,INFO)

        ENT(:)=1
        ENTW(:,:)=1
 
        DO LM2=1,ALMSO
          ENTW(1,LM2)=LM2
          DO LM1=LM2+1,ALMSO
            IF (ABS(W(LM1)-W(LM2)) .LT. 1.D-4 ) THEN 
              ENT(LM2)=ENT(LM2)+1
              ENT(LM1)=ENT(LM1)+1
              ENTW(ENT(LM2),LM2)=LM1
              ENTW(ENT(LM1),LM1)=LM2
            END IF
          END DO
        END DO

c  ---> check whether LV and RV are normed correctly
        
          DO LM2=1,ALMSO
            DO LM3=1,ALMSO
              NORM=0d0
              CALL ZGEMM('C','N',1,1,ALMSO,(1d0,0d0),LV(:,LM3),ALMSO,
     +                     RV(:,LM2),ALMSO,(0d0,0d0),NORM,1)
              IF (LM2 .EQ. LM3 ) THEN
                PHI=atan(aimag(NORM)/real(NORM))
                NORMT=NORM*exp(-(0,1.d0)*PHI)
                NORMABS=SQRT(REAL(NORMT)**2+AIMAG(NORMT)**2)
                IF(REAL(NORMT)<0.D0) THEN
                  DO LM4=1,ALMSO
                    RV(LM4,LM2)=-EXP(-(0,1.D0)*PHI)/NORMABS*RV(LM4,LM2) 
                  END DO
                ELSE
                  DO LM4=1,ALMSO
                    RV(LM4,LM2)=EXP(-(0,1.D0)*PHI)/NORMABS*RV(LM4,LM2) 
                  END DO
                END IF
              END IF
            END DO  ! LM3=1,ALMSO
          END DO    ! LM2=1,ALMSO

        IF (K.GT.1) THEN
c find the nearest band
        DO LM1=1,ALMSO
          IF (ABS(WAPR(LM1)).LT.1d6)  THEN
           MINVALUE=1d6
           LMTAKE=0
           DO LM2=1,ALMSO
            DIFF=ABS(W(LM2)-WAPR(LM1))
            IF (DIFF.LT.MINVALUE) THEN
             MINVALUE=DIFF
             LMTAKE=LM2
            ENDIF
           ENDDO
c find the connection if change size
       IF (LMTAKE.NE.0) THEN
       IF (AIMAG(W(LMTAKE))*AIMAG(WK(LM1)).LT.0d0) THEN
         XINTERSECT=-AIMAG(WK(LM1))/(AIMAG(W(LMTAKE))-AIMAG(WK(LM1)))
         REINTERSECT= (1d0-XINTERSECT)*REAL(WK(LM1))+
     &                XINTERSECT*REAL(W(LMTAKE))
          DO J=1,3
           KPINT(J)=(1d0-XINTERSECT)*BZKP(J,K-1)+XINTERSECT*BZKP(J,K)
          ENDDO
         IF (ABS(REINTERSECT).LE.1d-02) THEN
          WRITE(61,"(7e17.9)") (KPINT(J),J=1,3),EZ(IE),REINTERSECT      
         ELSEIF (ABS(REINTERSECT).GT.1D-02.AND.
     +           ABS(REINTERSECT).LE.1D0) THEN
          WRITE(62,"(7e17.9)") (KPINT(J),J=1,3),EZ(IE),REINTERSECT      
         ENDIF
       ENDIF
       ENDIF
       ENDIF
       ENDDO
       ENDIF

c save eigenvector
        DO LM1=1,ALMSO
          WK(LM1)=W(LM1)
        ENDDO

c --->  calculate the gradient of the Green function at k times dk
c       DGLLKE for extrapolating the eigenvalue at k+dk
         DO J=1,3
          KP(J,K) = BZKP(J,K)+DK(J,K)
         ENDDO
         DO I=1,NAEZ
          IC=CLS(KAOEZ(I))
          FAC=CONE
           IF (I.NE.EQINV(I)) FAC = CONEM      
           CALL CINIT(NAEZD*LMAXSQ*LMAXSQ,GLLKE1)
           CALL DLKE1(GLLKE1,ALAT,NACLS(IC),RR,EZOA(1,I),ATOM(1,I),
     +                KP(1,K),IE,IC,FAC,GINP(1,1,IC),RCLS(1,1,IC))   
             DO M=1,LMAXSQ
               IM=(I-1)*LMAXSQ+M
                DO JN=1,LMAXSQ*NAEZ
                 GLLKEDK(JN,IM) = GLLKEDK(JN,IM)+ GLLKE1(JN,M)
                ENDDO
             ENDDO
         ENDDO ! atom loop
          DO LM2=1,ALM
            DO LM1=1,ALM
              DGLLKEDK(LM1,LM2) = GLLKEDK(LM1,LM2)-GLLKE(LM1,LM2)
            END DO
          END DO
       
c  for spin-orbit coupling, dG_LL'(k) 
        IF (NSPD == 2) THEN
          DO LM1=1,ALM
            DO LM2=1,ALM
              DGLLKEDK(ALM+LM2,ALM+LM1) = DGLLKEDK(LM2,LM1) 
            END DO
          END DO
        END IF


c       the matrix dM=[ - dG*delta t] is constructed
        IS1=1
        IS2=1
        TMAT_2=CZERO
        HILF_2=CZERO
        DO I1=1,NAEZ
          DO IS1=1,NSPD
            DO LM1=1,LMAXSQ
              DO IS2=1,NSPD
                DO LM2=1,LMAXSQ
                  IL1=(IS1-1)*ALM+LMAXSQ*(I1-1)+LM1
                  IL2=(IS2-1)*ALM+LMAXSQ*(I1-1)+LM2
                  IL1T=(IS1-1)*LMAXSQ+LM1
                  IL2T=(IS2-1)*LMAXSQ+LM2
                  TMAT_2(IL2,IL1)=TMATLL(IL2T,IL1T,I1)/RFCTOR
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

c calculate delta(lamda)/dk

        CALL ZGEMM('N','N',ALMSO,ALMSO,ALMSO,CONE,DGLLKEDK,ALMSO,
     +                           TMAT_2,ALMSO,CZERO,HILF_2,ALMSO)
        DO LM1=1,ALMSO
          DO LM2=1,ALMSO
              HILF_2(LM1,LM2)=-HILF_2(LM1,LM2)
          ENDDO
        ENDDO

c pertubation method for dM --> dlamda/dk
c         write(55,*) IE,EZ(IE),K
         CALL CALC_WAPR_BSTR(W,ALMSO,LV,RV,HILF_2,
     +       ENT,ENTW,MAXVAL(ENT),EZ(IE),WAPR,FILTER_EIGVAL)
 300  END DO                    ! K = 1,KSUM

      DEALLOCATE(GLLKE1)
      DEALLOCATE(GLLKE)
      DEALLOCATE(GLLKES)
      DEALLOCATE(GLLKEDK)
      DEALLOCATE(DGLLKEDK)
      DEALLOCATE(HILF_1)
      DEALLOCATE(HILF_2)
      DEALLOCATE(TMAT_1)
      DEALLOCATE(TMAT_2)
      DEALLOCATE(RV)
      DEALLOCATE(LV)
      DEALLOCATE(W)


c      write(6,*) "end of KKRMATEIGEN_SO"
      RETURN

 9000 FORMAT(20I4)
 9010 FORMAT(3f12.5)
 9011 FORMAT('K = ',i4,3f12.5)
 9030 FORMAT(A1,4x,':',I6,2d12.3,2f10.4)
 9040 FORMAT(' NUMTAU =',I6)
 9050 FORMAT(' NUM_LR =',I6)

      END















