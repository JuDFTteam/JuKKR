      SUBROUTINE CALC_WAPR_BSTR(W,ALMSO,LV,RV,DGLLKE,
     +                        ENT,ENTW,ENTMAX,ENER,WAPR,FILTER_EIGVAL)

      implicit none

      DOUBLE COMPLEX        :: CZERO,CONE
      PARAMETER (CZERO=(0d0,0d0),CONE=(1d0,0d0))

      INTEGER           :: ALMSO,ENT(ALMSO),ENTW(ALMSO,ALMSO),
     +                         ENTMAX
      DOUBLE COMPLEX    :: W(ALMSO),LV(ALMSO,ALMSO),ENER,
     +                         RV(ALMSO,ALMSO),DGLLKE(ALMSO,ALMSO)

      DOUBLE COMPLEX    :: WAPR(ALMSO)

c ---> local variables

      INTEGER           :: LM1,LM2,LM3,LM4,IENT1,IENT2,
     +                     LMENT1,LMENT2,J,
     +                     NAUX,INFO
      DOUBLE PRECISION  :: FILTER_EIGVAL

      DOUBLE PRECISION,allocatable :: DAUX(:)

      DOUBLE COMPLEX, allocatable  :: APPR(:),AUX(:)
      DOUBLE COMPLEX, allocatable ::  W_CONSTR(:,:),WENT(:),
     +                                LVENT(:,:),RVENT(:,:)
c     +                                SOVERLAP(:,:),
c     +                                ALPHA(:),BETA(:)

      DOUBLE COMPLEX               :: TR_APPR,TR_SOVLP
c      DOUBLE COMPLEX               :: W_CONSTR(2,2),WENT(2),
c     +                                LVENT(2,2),RVENT(2,2),
c     +                                SOVERLAP(2,2),ALPHA(2),BETA(2)

      NAUX=2*ENTMAX**2+5*ENTMAX
      ALLOCATE(AUX(NAUX))
      ALLOCATE(DAUX(2*ENTMAX))
      

      ALLOCATE(APPR(ALMSO))
      
      WAPR=CZERO

      DO LM1=1,ALMSO
 
       IF (ABS(W(LM1)). LE .FILTER_EIGVAL) THEN

        IF (ENT(LM1) .EQ. 1) THEN 

          APPR=CZERO

          CALL ZGEMM('N','N', ALMSO,1,ALMSO,CONE,DGLLKE,ALMSO,
     +                 RV(:,LM1),ALMSO,CZERO,APPR,ALMSO)

          TR_APPR=CZERO
          CALL ZGEMM('C','N',1,1,ALMSO,CONE,LV(:,LM1),ALMSO,
     +                 APPR,ALMSO,CZERO,TR_APPR,1)
c          TR_SOVLP=CONE
c          CALL ZGEMM('C','N',1,1,ALMSO,CONE,LV(:,LM1),ALMSO,
c     +             RV(:,LM1),ALMSO,CZERO,TR_SOVLP,1)

          WAPR(LM1)=W(LM1)+TR_APPR

        ELSE IF (ENT(LM1).GT.1.AND.
     +                ENTW(1,LM1).LT.ENTW(2,LM1)) THEN

          ALLOCATE(W_CONSTR(ENT(LM1),ENT(LM1)))
c          ALLOCATE(SOVERLAP(ENT(LM1),ENT(LM1)))
          ALLOCATE(WENT(ENT(LM1)))
          ALLOCATE(LVENT(ENT(LM1),ENT(LM1)))
          ALLOCATE(RVENT(ENT(LM1),ENT(LM1)))
c          ALLOCATE(ALPHA(ENT(LM1)))
c          ALLOCATE(BETA(ENT(LM1)))

          DO IENT1=1,ENT(LM1) 

            LMENT1=ENTW(IENT1,LM1)

            DO IENT2=1,ENT(LM1) 

              LMENT2=ENTW(IENT2,LM1)

              APPR=CZERO

              CALL ZGEMM('N','N', ALMSO,1,ALMSO,CONE,DGLLKE,ALMSO,
     +                 RV(:,LMENT2),ALMSO,CZERO,APPR,ALMSO)

              W_CONSTR(IENT1,IENT2)=CZERO
              CALL ZGEMM('C','N',1,1,ALMSO,CONE,LV(:,LMENT1),ALMSO,
     +                 APPR,ALMSO,CZERO,W_CONSTR(IENT1,IENT2),1)

c              SOVERLAP(IENT1,IENT2)=CZERO
c              CALL ZGEMM('C','N',1,1,ALMSO,CONE,LV(:,LMENT1),ALMSO,
c     +             RV(:,LMENT2),ALMSO,CZERO,SOVERLAP(IENT1,IENT2),1)
          

            END DO
          END DO
          
                      
          LVENT=CZERO
          RVENT=CZERO
          WENT=CZERO
c          ALPHA=CZERO
c          BETA=CONE
          
         CALL ZGEEV('V','V',ENT(LM1),W_CONSTR,ENT(LM1),WENT,LVENT,
     +              ENT(LM1),RVENT,ENT(LM1),AUX,NAUX,DAUX,INFO)
                
c         CALL ZGGEV('N','N',ENT(LM1),W_CONSTR,ENT(LM1),SOVERLAP,
c     +         ENT(LM1),ALPHA,BETA,LVENT,ENT(LM1),RVENT,
c     +        ENT(LM1),AUX,NAUX,DAUX,INFO)
c         IF (INFO.NE.0) STOP "Problems with ZGGEV"
c          CALL ZEIGENOVERLAP(W_CONSTR,SOVERLAP,ALPHA,LVENT,RVENT,
c     +                      ENT(LM1))
          DO IENT1=1,ENT(LM1) 
c            WENT(IENT1)=ALPHA(IENT1)/BETA(IENT1)
c            WENT(IENT1)=ALPHA(IENT1)
            WAPR(ENTW(IENT1,LM1))=W(ENTW(IENT1,LM1))+WENT(IENT1)
          END DO

          DEALLOCATE(W_CONSTR)
c          DEALLOCATE(SOVERLAP)
          DEALLOCATE(WENT)
          DEALLOCATE(LVENT)
          DEALLOCATE(RVENT)
c          DEALLOCATE(ALPHA)
c          DEALLOCATE(BETA)

        END IF     ! ENT .EQ. 1

c ---> calculate delta(lamda)/dk
       ELSE
       WAPR(LM1)=(1d6,1d6)
       ENDIF
      ENDDO
      DEALLOCATE(APPR)

      DEALLOCATE(AUX)
      DEALLOCATE(DAUX)

      END SUBROUTINE CALC_WAPR_BSTR
