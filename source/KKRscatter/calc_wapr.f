      SUBROUTINE CALC_WAPR(W,WAPR,ALMSO,LV,RV,WSIGN,I_WSIGN,DGLLKE,
     +            ENT,ENTW,ENTMAX,DK,K0,EZ,KINTER,WREALAPP)

      implicit none

c      INTEGER               :: ENTMAX
c      PARAMETER (ENTMAX=MAXVAL(ENT))

      INTEGER, INTENT(IN)   :: ALMSO,ENT(ALMSO),ENTW(ALMSO,ALMSO),
     +                         ENTMAX
      REAL, INTENT(IN)      :: DK(6),K0(6)
      COMPLEX,INTENT(IN)    :: W(ALMSO),LV(ALMSO,ALMSO),EZ,
     +                         RV(ALMSO,ALMSO),DGLLKE(ALMSO,ALMSO)

      INTEGER,INTENT(OUT)   :: WSIGN,I_WSIGN(ALMSO)
      COMPLEX,INTENT(OUT)   :: WAPR(ALMSO)
      REAL,INTENT(OUT)      :: KINTER(6,ALMSO),WREALAPP(ALMSO)

c ---> local variables

      INTEGER               :: LM1,LM2,LM3,IENT1,IENT2,
     +                         LMENT1,LMENT2,J,
     +                         NAUX,INFO
c      PARAMETER(NAUX = 2*ENTMAX**2+5*ENTMAX)

c      DOUBLE PRECISION      :: DAUX(2*ENTMAX)
      DOUBLE PRECISION,allocatable :: DAUX(:)

      COMPLEX, allocatable  :: PROJ(:,:),APPR(:,:),
     +                         W_CONSTR(:,:),WENT(:),
     +                         LVENT(:,:),RVENT(:,:),
     +                         AUX(:)

      REAL                  :: DK_INTER,DKABS!,KINTER(3)
      COMPLEX               :: DEW,TR_APPR!,AUX(NAUX)

      NAUX=2*ENTMAX**2+5*ENTMAX
      KINTER=0d0
      ALLOCATE(AUX(NAUX))
      ALLOCATE(DAUX(2*ENTMAX))

      ALLOCATE(PROJ(ALMSO,ALMSO))
      ALLOCATE(APPR(ALMSO,ALMSO))

c     WRITE(6,*) "DGLKKE"

c     DO LM2=1,ALMSO
c       DO LM1=1,ALMSO
c         WRITE(6,"((2I5),(2e17.9))") LM2,LM1,DGLLKE(LM2,LM1)
c       END DO
c     END DO

c     WRITE(6,*) "after DGLKKE"

c      WRITE(147,"(6e17.9)") (K0(J),J=1,6)
c      write (6,*) "IN CALC"

      PROJ=0d0
      APPR=0d0
      WSIGN=0
      I_WSIGN(:)=0
      WAPR=0d0

      DKABS=0d0
      DO J=1,3
        DKABS=DKABS+DK(J)**2
      END DO
      DKABS=SQRT(DKABS)

      DO LM1=1,ALMSO

c       write (6,*) "RV,LV"
c        write (6,*) "LM1",LM1

c       DO LM2=1,ALMSO
c         WRITE(6,"((I5),(4e17.9))") LM2,
c    +                 RV(LM2,LM1),LV(LM2,LM1)
c       END DO

c       write (6,*) " "
c       write (6,*) "before first ZGEMM"
c       write (6,*) " "
        PROJ=0d0
        APPR=0d0

        IF (ENT(LM1) .EQ. 1) THEN 

          PROJ=0d0
          APPR=0d0

c          WRITE(6,*) "before ZGEMM"

          CALL ZGEMM('N','C', ALMSO,ALMSO,1,1d0,RV(:,LM1),ALMSO,
     +                 LV(:,LM1),ALMSO,0d0,PROJ,ALMSO)

c          WRITE(6,*) "after ZGEMM"
c         WRITE(6,*) "PROJ"

c         DO IENT1=1,ALMSO
c           DO IENT2=1,ALMSO
c             WRITE(6,"((2I5),(4e17.9))") IENT2,IENT1,
c    +                 PROJ(IENT2,IENT1)!,DGLLKE(LM3,LM2)
c           END DO
c         END DO

c        write (6,*) " "
c          write (6,*) "after proj"
c        write (6,*) " "

c         STOP "N A N " 

          CALL ZGEMM('N','N', ALMSO,ALMSO,ALMSO,1d0,DGLLKE,ALMSO,
     +                 PROJ,ALMSO,0d0,APPR,ALMSO)

          TR_APPR=0d0

          DO LM2=1,ALMSO
            TR_APPR=TR_APPR+APPR(LM2,LM2)
c            WRITE(6,"((I5),(2e17.9))") LM2,TR_APPR
          END DO
                       
c          write (6,*) "after TR_APPR"
          WAPR(LM1)=W(LM1)+TR_APPR
c          write (6,*) "LM1,WAPR(LM1)",LM1,WAPR(LM1)

        ELSE IF (ENT(LM1) > 1 .AND. 
     +           ENTW(1,LM1) .LT. ENTW(2,LM1) ) THEN 

          ALLOCATE(W_CONSTR(ENT(LM1),ENT(LM1)))
          ALLOCATE(WENT(ENT(LM1)))
          ALLOCATE(LVENT(ENT(LM1),ENT(LM1)))
          ALLOCATE(RVENT(ENT(LM1),ENT(LM1)))

          DO IENT1=1,ENT(LM1) 

            LMENT1=ENTW(IENT1,LM1)

            DO IENT2=1,ENT(LM1) 

              LMENT2=ENTW(IENT2,LM1)
              PROJ=0d0
              APPR=0d0

c             DO LM2=1,ALMSO
c               WRITE(147,"((5I5),(4e17.9))") IENT1,LMENT1,IENT2,LMENT2,
c    +                   LM2,RV(LM2,LMENT2),LV(LM2,LMENT1)
c             END DO
c             WRITE(147,*) "" 

              CALL ZGEMM('N','C', ALMSO,ALMSO,1,1d0,RV(:,LMENT2),ALMSO,
     +                 LV(:,LMENT1),ALMSO,0d0,PROJ,ALMSO)

              CALL ZGEMM('N','N', ALMSO,ALMSO,ALMSO,1d0,DGLLKE,ALMSO,
     +                 PROJ,ALMSO,0d0,APPR,ALMSO)

              W_CONSTR(IENT1,IENT2)=0d0

              DO LM2=1,ALMSO
                W_CONSTR(IENT1,IENT2) = W_CONSTR(IENT1,IENT2) + 
     +                                                 APPR(LM2,LM2)
              END DO

            END DO
          END DO
                       
          LVENT=0d0
          RVENT=0d0
          WENT=0d0

c          WRITE(6,*) "W_CONSTR "

c         DO IENT1=1,ENT(LM1)
c           DO IENT2=1,ENT(LM1)
c             WRITE(6,"((2I5),(2e17.9))") IENT2,IENT1,
c    +                           W_CONSTR(IENT1,IENT2)
c           END DO
c         END DO

c          WRITE(6,*) "before ZGEEV calc_wapr"

          CALL ZGEEV('V','V',ENT(LM1),W_CONSTR,ENT(LM1),WENT,LVENT,
     +          ENT(LM1),RVENT,ENT(LM1),AUX,NAUX,DAUX,INFO)

c          WRITE(6,*) "After  ZGEEV calc_wapr"

          DO IENT1=1,ENT(LM1) 
c            WRITE(146, "((2I5),(2e17.9))") IENT1,ENTW(IENT1,LM1),
c     +                                                 WENT(IENT1)
            WAPR(ENTW(IENT1,LM1))=W(ENTW(IENT1,LM1))+WENT(IENT1)
          END DO

          DEALLOCATE(W_CONSTR)
          DEALLOCATE(WENT)
          DEALLOCATE(LVENT)
          DEALLOCATE(RVENT)

        END IF     ! ENT .EQ. 1

c ---> check whether one of the eigenvalues changes sign

        IF( ABS(AIMAG(WAPR(LM1))  + AIMAG(W(LM1))) .NE. 
     +        ABS(AIMAG(WAPR(LM1))) + ABS(AIMAG(W(LM1))) ) THEN 
          WSIGN=WSIGN+1
          I_WSIGN(WSIGN)=LM1

c ---> calculate the intersection with the x-axis
          DEW=(WAPR(LM1)-W(LM1))/DKABS
          DK_INTER=-AIMAG(W(LM1))/AIMAG(DEW)
          WREALAPP(WSIGN)=REAL(W(LM1)) + REAL(DEW)*DK_INTER            
          DO J=1,3
            KINTER(J,WSIGN)=K0(J)+DK_INTER/DKABS*DK(J)
          END DO
c         IF (ABS(WREALAPP) .LT. 1.D-2) THEN 
c           WRITE(39,"(7e17.9)") (KINTER(J),J=1,3),EZ,WREALAPP
c         ELSEIF (ABS(WREALAPP) .LT. 0.5D-1) THEN 
c           WRITE(140,"(7e17.9)") (KINTER(J),J=1,3),EZ,WREALAPP
c         ELSEIF (ABS(WREALAPP) .LT. 1.0D-1) THEN 
c           WRITE(141,"(7e17.9)") (KINTER(J),J=1,3),EZ,WREALAPP
c         END IF
        END IF             

      END DO

      DEALLOCATE(PROJ)
      DEALLOCATE(APPR)

      DEALLOCATE(AUX)
      DEALLOCATE(DAUX)

      END SUBROUTINE CALC_WAPR
