      SUBROUTINE CALC_WAPR_SOF(W,ALMSO,LV,RV,HILF_2,ENT,ENTW,
     +           ENTMAX,PMUP,PMDOWN,RHOD,NATYPD,LMMAXSO,LMAXSQ,DENS1,
     +           CKGES)

      implicit none

      DOUBLE COMPLEX        :: CZERO,CONE
      PARAMETER (CZERO=(0d0,0d0),CONE=(1d0,0d0))
      INTEGER               :: ALMSO,ENT(ALMSO),ENTW(ALMSO,ALMSO),
     +                         ENTMAX
      DOUBLE COMPLEX        :: W(ALMSO),LV(ALMSO,ALMSO),
     +                         RV(ALMSO,ALMSO),HILF_2(ALMSO,ALMSO)      

      DOUBLE COMPLEX        :: WAPR(ALMSO)
      DOUBLE COMPLEX        :: RHOD(LMMAXSO,LMMAXSO,NATYPD,4),DENS(2),
     +                         DENS1(2,4),TMP
      DOUBLE COMPLEX        :: CKGES(LMMAXSO,NATYPD,2),CKHELP(LMMAXSO)
      INTEGER I1,LMMAXSO,LMAXSQ,IENT,NATYPD,ISIGMA

c ---> local variables

      INTEGER               :: LM1,LM2,LM3,IENT1,IENT2,
     +                         LMENT1,LMENT2,J,I,
     +                         INFO,LM1TAKE,LM2TAKE,LM3TAKE
      DOUBLE PRECISION      :: MINLM1,MINLM2,DELTAK


      DOUBLE COMPLEX, allocatable  :: 
     +                         APPR(:),
     +                         W_CONSTR(:,:),
     +                         ALPHA(:),BETA(:),
     +                         LVENT(:,:),RVENT(:,:),
     +                         SOVERLAP(:,:)
      DOUBLE COMPLEX        :: TR_APPR,TR_SOVLP
      DOUBLE COMPLEX        :: PMUP,PMDOWN,PROD
      DOUBLE COMPLEX        :: NORM,NORMT,NORMABS,SOVER,C1,
     +                         RV2(ALMSO,2),LV2(ALMSO,2)
      DOUBLE COMPLEX        :: M(2,2),C(2),W_CONSTR1(2,2)
      DOUBLE PRECISION      :: PHI
      DOUBLE COMPLEX, ALLOCATABLE :: AUX(:)
      DOUBLE PRECISION, ALLOCATABLE :: DAUX(:)
      INTEGER NAUX
      NAUX=2*ENTMAX**2+5*ENTMAX
      ALLOCATE(AUX(NAUX))
      ALLOCATE(DAUX(8*ENTMAX))

      ALLOCATE(APPR(ALMSO))
     

      MINLM1=1D6
      DO LM1=1,ALMSO
       IF (ABS(W(LM1)).LE.MINLM1) THEN
        MINLM1=ABS(W(LM1))
        LM1TAKE=LM1
       ENDIF
      ENDDO
      LM2TAKE=ENTW(1,LM1TAKE)
      LM3TAKE=ENTW(2,LM1TAKE)
      
         IF (ABS(AIMAG(W(LM1TAKE))). LE .1d-02 
     +        .AND. ABS(REAL(W(LM1TAKE))). LE. 1D-02) THEN
        IF (ENT(LM1TAKE) .EQ. 1) THEN 

          APPR=CZERO

          CALL ZGEMM('N','N', ALMSO,1,ALMSO,CONE,HILF_2,ALMSO,
     +                 RV(:,LM1TAKE),ALMSO,CZERO,APPR,ALMSO)

          TR_APPR=CZERO
          CALL ZGEMM('C','N', 1,1,ALMSO,CONE,LV(:,LM1TAKE),ALMSO,
     +                 APPR,ALMSO,CZERO,TR_APPR,1)

          TR_SOVLP=CONE
          CALL ZGEMM('C','N', 1,1,ALMSO,CONE,LV(:,LM1TAKE),ALMSO,
     +                 RV(:,LM1TAKE),ALMSO,CZERO,TR_SOVLP,1)

          WAPR(LM1TAKE)=W(LM1TAKE)+TR_APPR/TR_SOVLP

c        ELSE IF (ENT(LM1TAKE) > 1 .AND. 
c     +           ENTW(1,LM1TAKE) .LT. ENTW(2,LM1TAKE) ) THEN 
        ELSE IF (ENT(LM1TAKE) > 1) THEN
          ALLOCATE(W_CONSTR(ENT(LM1TAKE),ENT(LM1TAKE)))
          ALLOCATE(SOVERLAP(ENT(LM1TAKE),ENT(LM1TAKE)))
          ALLOCATE(LVENT(ENT(LM1TAKE),ENT(LM1TAKE)))
          ALLOCATE(RVENT(ENT(LM1TAKE),ENT(LM1TAKE)))
          ALLOCATE(ALPHA(ENT(LM1TAKE)))
          ALLOCATE(BETA(ENT(LM1TAKE)))
           
          DO IENT1=1,ENT(LM1TAKE) 
           LMENT1=ENTW(IENT1,LM1TAKE)
           DO IENT2=1,ENT(LM1TAKE)
            LMENT2=ENTW(IENT2,LM1TAKE)
            APPR=CZERO
            CALL ZGEMM('N','N', ALMSO,1,ALMSO,CONE,HILF_2,ALMSO,
     +                 RV(:,LMENT2),ALMSO,CZERO,APPR,ALMSO)

            W_CONSTR(IENT1,IENT2)=CZERO
            CALL ZGEMM('C','N', 1,1,ALMSO,CONE,LV(:,LMENT1),ALMSO,
     +                 APPR,ALMSO,CZERO,W_CONSTR(IENT1,IENT2),1)

            SOVERLAP(IENT1,IENT2)=CZERO
            CALL ZGEMM('C','N',1,1,ALMSO,CONE,LV(:,LMENT1),ALMSO,
     +                 RV(:,LMENT2),ALMSO,CZERO,SOVERLAP(IENT1,IENT2),1)
            END DO
           END DO

           LVENT=CZERO
           RVENT=CZERO
           ALPHA=CZERO
           BETA=CONE
         CALL ZGGEV('V','V',ENT(LM1TAKE),W_CONSTR,ENT(LM1TAKE),SOVERLAP,
     +         ENT(LM1TAKE),ALPHA,BETA,LVENT,ENT(LM1TAKE),RVENT,
     +        ENT(LM1TAKE),AUX,NAUX,DAUX,INFO)
          PMUP= ALPHA(1)/BETA(1)
          PMDOWN= ALPHA(2)/BETA(2)

          DO IENT1=1,2
           DO IENT2=1,2
c            WRITE(125,'((2I5),(4e17.9))') IENT1,IENT2,
c     +                    RVENT(IENT1,IENT2),LVENT(IENT1,IENT2)
           ENDDO
          ENDDO

c test coefficient
c          TMP=DCONJG(RVENT(2,2))/
c     +        SQRT(DBLE(RVENT(2,2))**2+DIMAG(RVENT(2,2)**2))
c          RVENT(2,2)=RVENT(2,2)*TMP
c          WRITE(6,*) RVENT(2,2)

          DO LM1=1,ALMSO
          RV2(LM1,1)=RVENT(1,1)*RV(LM1,LM2TAKE)+
     +               RVENT(2,1)*RV(LM1,LM3TAKE)
          LV2(LM1,1)=LVENT(1,1)*LV(LM1,LM2TAKE)+
     +               LVENT(2,1)*LV(LM1,LM3TAKE)
          RV2(LM1,2)=RVENT(1,2)*RV(LM1,LM2TAKE)+
     +               RVENT(2,2)*RV(LM1,LM3TAKE)
          LV2(LM1,2)=LVENT(1,2)*LV(LM1,LM2TAKE)+
     +               LVENT(2,2)*LV(LM1,LM3TAKE)
          ENDDO
          DO IENT=1,2
           DO LM1=1,ALMSO
c            RV2(LM1,IENT)=RV2(LM1,IENT)
c     +                    /SQRT(DOT_PRODUCT(LV2(:,IENT),RV2(:,IENT)))
c            LV2(LM1,IENT)=LV2(LM1,IENT)
c     +                    /SQRT(DOT_PRODUCT(LV2(:,IENT),RV2(:,IENT)))
           ENDDO
          ENDDO

          CKGES=CZERO
          DENS1=CZERO
          DENS=CZERO
          DO IENT=1,2
           LM2=0
            DO I1=1,NATYPD
             DO LM1=1,LMAXSQ
              LM2=LM2+1
              CKGES(LM1,I1,IENT)=RV2(LM2,IENT)
             ENDDO
            ENDDO
            DO I1=1,NATYPD
             DO LM1=1,LMAXSQ
              LM2=LM2+1
              CKGES(LM1+LMAXSQ,I1,IENT)=RV2(LM2,IENT)
             ENDDO
            ENDDO
           DO I1=1,NATYPD
            CKHELP=CZERO
             CALL ZGEMM('N','N',LMMAXSO,1,LMMAXSO,CONE,RHOD(:,:,I1,1),
     +           LMMAXSO,CKGES(:,I1,IENT),LMMAXSO,CZERO,CKHELP,LMMAXSO) 
            DENS(IENT)=DENS(IENT)+
     +             DOT_PRODUCT(CKGES(:,I1,IENT),CKHELP)
           ENDDO
            DO I1=1,NATYPD
             DO LM1=1,LMMAXSO
              CKGES(LM1,I1,IENT)=CKGES(LM1,I1,IENT)/
     +                    CDSQRT(DENS(IENT))
             ENDDO
            ENDDO
          DO ISIGMA=1,4
           DO I1=1,NATYPD
            CKHELP=CZERO
          CALL ZGEMM('N','N',LMMAXSO,1,LMMAXSO,CONE,RHOD(:,:,I1,ISIGMA),
     +           LMMAXSO,CKGES(:,I1,IENT),LMMAXSO,CZERO,CKHELP,LMMAXSO) 
            DENS1(IENT,ISIGMA)=DENS1(IENT,ISIGMA)+
     +           DOT_PRODUCT(CKGES(:,I1,IENT),CKHELP)
           ENDDO
          ENDDO
c           WRITE(145,'((I5),(2e17.9))') IENT,REAL(DENS1(IENT,1))
c           WRITE(144+IENT,'((6e17.9))') 
c     +              (REAL(DENS1(IENT,J)),J=2,4)
          ENDDO
          DEALLOCATE(W_CONSTR)
          DEALLOCATE(SOVERLAP)
          DEALLOCATE(LVENT)
          DEALLOCATE(RVENT)
          DEALLOCATE(ALPHA)
          DEALLOCATE(BETA)
        END IF     ! ENT .EQ. 1

      ENDIF
      DEALLOCATE(APPR)
      DEALLOCATE(AUX,DAUX)

      END SUBROUTINE CALC_WAPR_SOF
