      SUBROUTINE ROT_RLL(ROTSPIN,RLL,RLLROT,LMMAX,IRMD,NATYPD)

      IMPLICIT NONE

      INTEGER,INTENT(IN)       :: LMMAX,IRMD,NATYPD
      COMPLEX,INTENT(IN)       :: RLL(IRMD,LMMAX,LMMAX,2,2,NATYPD),
     +                            ROTSPIN(2*LMMAX,2*LMMAX)
      COMPLEX,INTENT(OUT)      :: RLLROT(IRMD,LMMAX,LMMAX,2,2,NATYPD)
      COMPLEX,ALLOCATABLE      :: RLLHELP(:,:),UDRLL(:,:),UDRLLU(:,:)

      INTEGER                  :: IR,I1,LMMAXSO,LM1,LM2,LM1S,LM2S,
     +                            ISP1,ISP2

      LMMAXSO=2*LMMAX

      ALLOCATE(UDRLL(2*LMMAX,2*LMMAX))
      ALLOCATE(UDRLLU(2*LMMAX,2*LMMAX))
      ALLOCATE(RLLHELP(2*LMMAX,2*LMMAX))

      RLLROT=0d0

      DO I1=1,NATYPD
        DO IR = 1,IRMD

          RLLHELP=0d0
          UDRLLU=0d0
          UDRLL=0d0

          DO ISP2=1,2
            DO LM2=1,LMMAX
              LM2S=(ISP2-1)*LMMAX+LM2
              DO ISP1=1,2
                DO LM1=1,LMMAX
                  LM1S=(ISP1-1)*LMMAX+LM1
                  RLLHELP(LM1S,LM2S)=RLL(IR,LM1,LM2,ISP1,ISP2,I1)
                END DO
              END DO
            END DO
          END DO

c         CALL ZGEMM('C','N',LMMAXSO,LMMAXSO,LMMAXSO,1d0,ROTSPIN,
c    +                    LMMAXSO,RLLHELP,LMMAXSO,0d0,UDRLL,LMMAXSO)

c         CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,1d0,UDRLL,
c    +                    LMMAXSO,ROTSPIN,LMMAXSO,0d0,UDRLLU,LMMAXSO)

          CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,1d0,ROTSPIN,
     +                    LMMAXSO,RLLHELP,LMMAXSO,0d0,UDRLL,LMMAXSO)

          CALL ZGEMM('N','C',LMMAXSO,LMMAXSO,LMMAXSO,1d0,UDRLL,
     +                    LMMAXSO,ROTSPIN,LMMAXSO,0d0,UDRLLU,LMMAXSO)

          DO ISP2=1,2
            DO LM2=1,LMMAX
              LM2S=(ISP2-1)*LMMAX+LM2
              DO ISP1=1,2
                DO LM1=1,LMMAX
                  LM1S=(ISP1-1)*LMMAX+LM1
                  RLLROT(IR,LM1,LM2,ISP1,ISP2,I1)=UDRLLU(LM1S,LM2S)
                END DO
              END DO
            END DO
          END DO

        END DO
      END DO

      DEALLOCATE(RLLHELP)
      DEALLOCATE(UDRLL)
      DEALLOCATE(UDRLLU)

      END SUBROUTINE ROT_RLL
