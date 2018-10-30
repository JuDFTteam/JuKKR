      SUBROUTINE CALC_WAPR_FV(W,ALMSO,LV,RV,DGLLKE,
     +           ENT,ENTW,ENTMAX,EZ,DELTALAMDAK)

      implicit none

      DOUBLE COMPLEX        :: CZERO,CONE
      PARAMETER (CZERO=(0d0,0d0),CONE=(1d0,0d0))
      
      INTEGER, INTENT(IN)   :: ALMSO,ENT(ALMSO),ENTW(ALMSO,ALMSO),
     +                         ENTMAX
      COMPLEX,INTENT(IN)    :: W(ALMSO),LV(ALMSO,ALMSO),EZ,
     +                         RV(ALMSO,ALMSO),DGLLKE(ALMSO,ALMSO)

      DOUBLE COMPLEX        :: WAPR(ALMSO)

c ---> local variables

      INTEGER               :: LM1,LM2,LM3,IENT1,IENT2,
     +                         LMENT1,LMENT2,J,
     +                         NAUX,INFO,LM1TAKE,LM2TAKE
      DOUBLE PRECISION      :: MINLM1,MINLM2

      DOUBLE PRECISION,allocatable :: DAUX(:)

      COMPLEX, allocatable  :: APPR(:),
     +                         W_CONSTR(:,:),WENT(:),
     +                         LVENT(:,:),RVENT(:,:),
     +                         AUX(:)

      COMPLEX               :: DEW,TR_APPR,TR_APPR2
      DOUBLE COMPLEX        :: DELTALAMDAK

      NAUX=2*ENTMAX**2+5*ENTMAX
      ALLOCATE(AUX(NAUX))
      ALLOCATE(DAUX(2*ENTMAX))

      ALLOCATE(APPR(ALMSO))

      MINLM1=1D6
      DO LM1=1,ALMSO
       IF (ABS(W(LM1)).LE.MINLM1) THEN
        MINLM1=ABS(W(LM1))
        LM1TAKE=LM1
       ENDIF
      ENDDO
      
         IF (ABS(AIMAG(W(LM1TAKE))). LE .1d-02 
     +        .AND. ABS(REAL(W(LM1TAKE))). LE. 1D-02) THEN

        IF (ENT(LM1TAKE) .EQ. 1) THEN 

          APPR=CZERO
          CALL ZGEMM('N','N', ALMSO,1,ALMSO,CONE,DGLLKE,ALMSO,
     +                 RV(:,LM1TAKE),ALMSO,CZERO,APPR,ALMSO)

          TR_APPR=CZERO
          CALL ZGEMM('C','N', 1, 1,ALMSO,CONE,LV(:,LM1TAKE),ALMSO,
     +                 APPR,ALMSO,CZERO,TR_APPR,1)

          WAPR(LM1TAKE)=W(LM1TAKE)+TR_APPR

        ELSE IF (ENT(LM1TAKE) > 1) THEN
c     +           ENTW(1,LM1TAKE) .LT. ENTW(2,LM1TAKE) ) THEN 

          ALLOCATE(W_CONSTR(ENT(LM1TAKE),ENT(LM1TAKE)))
          ALLOCATE(WENT(ENT(LM1TAKE)))
          ALLOCATE(LVENT(ENT(LM1TAKE),ENT(LM1TAKE)))
          ALLOCATE(RVENT(ENT(LM1TAKE),ENT(LM1TAKE)))

          DO IENT1=1,ENT(LM1TAKE) 

            LMENT1=ENTW(IENT1,LM1TAKE)

            DO IENT2=1,ENT(LM1TAKE) 

              LMENT2=ENTW(IENT2,LM1TAKE)

              APPR=CZERO

              CALL ZGEMM('N','N', ALMSO,1,ALMSO,CONE,DGLLKE,ALMSO,
     +                      RV(:,LMENT2),ALMSO,CZERO,APPR,ALMSO)

              W_CONSTR(IENT1,IENT2)=CZERO
              CALL ZGEMM('C','N', 1, 1,ALMSO,CONE,LV(:,LMENT1),ALMSO,
     +                     APPR,ALMSO,CZERO,W_CONSTR(IENT1,IENT2),1)

            END DO
          END DO
                       
          LVENT=CZERO
          RVENT=CZERO
          WENT=CZERO

          CALL ZGEEV('V','V',ENT(LM1TAKE),W_CONSTR,ENT(LM1TAKE),WENT,
     +        LVENT,ENT(LM1TAKE),RVENT,ENT(LM1TAKE),AUX,NAUX,DAUX,INFO)

          DO IENT1=1,ENT(LM1TAKE) 
            WAPR(ENTW(IENT1,LM1TAKE))=W(ENTW(IENT1,LM1TAKE))+WENT(IENT1)
          END DO

          DEALLOCATE(W_CONSTR)
          DEALLOCATE(WENT)
          DEALLOCATE(LVENT)
          DEALLOCATE(RVENT)

        END IF     ! ENT .EQ. 1

c ---> calculate delta(lamda)/dk
c         IF (ABS(AIMAG(W(LM1))). LE .1d-04 
c     +        .AND. ABS(REAL(W(LM1))). LE. 1D-04) THEN
          DELTALAMDAK= WAPR(LM1TAKE)-W(LM1TAKE)
c         ENDIF
      ENDIF
      DEALLOCATE(APPR)

      DEALLOCATE(AUX)
      DEALLOCATE(DAUX)

      END SUBROUTINE CALC_WAPR_FV
