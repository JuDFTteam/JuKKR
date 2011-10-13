      SUBROUTINE SREDGM(
     >                  NSPIN,IELAST,
     >                  MYRANK,
     >                  SMPIC,SMYRANK,SRANK,
     >                  EMPIC,EMYRANK,ERANK,EPROC,
     >                  GMATN,LLY_GRDT,
     <                  GMATN_ALL,LLY_GRDT_ALL)
C ======================================================================
C
C this routine preforms an MPI_ALLREDUCE command for
C a given set of ranks
C
C                               Alexander Thiess, 7th of December 2009
C =======================================================================
C
      IMPLICIT NONE
C
      include 'mpif.h'
C     .. Parameters ..
      include 'inc.p'
C 
      INTEGER        LMMAXD
      PARAMETER      (LMMAXD= (LMAXD+1)**2)
C
C     .. global scalars ..
      INTEGER        ISPIN,NSPIN
C     .. global arrays ..
      DOUBLE COMPLEX GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND),
     +               GMATN_ALL(LMMAXD,LMMAXD,IEMXD,NSPIND),
     +               LLY_GRDT(IEMXD,NSPIND),
     +               LLY_GRDT_ALL(IEMXD,NSPIND)
C     .. local scalars ..
      INTEGER        IE,IELAST,
     +               LM1,LM2,SEND,RECV
C     .. local arrays ..
      DOUBLE COMPLEX GSEND(LMMAXD*LMMAXD*IEMXD+IEMXD),
     +               GRECV(LMMAXD*LMMAXD*IEMXD+IEMXD),
     +               GESEND(LMMAXD*LMMAXD+1,NSPIND),
     +               GERECV(LMMAXD*LMMAXD+1,NSPIND)
C
C     .. MPI ..
      INTEGER, dimension(MPI_STATUS_SIZE) :: STATUS
      INTEGER        IERR,MAPBLOCK
C     .. N-MPI
      INTEGER        MYRANK
C     .. LS-MPI
      INTEGER        SMYRANK(SMPID,NAEZD*LMPID*EMPID),
     +               SRANK(SMPID,NAEZD*LMPID*EMPID),
     +               SMPI,SMPIC
C     .. E-MPI
      INTEGER        EMYRANK(EMPID,NAEZD*LMPID*SMPID),
     +               ERANK(EMPID,NAEZD*LMPID*SMPID),
     +               EMPI,EMPIC,
     +               EPROC(IEMXD)
C
C     fist copy the complete array to destination array:
C
      CALL ZCOPY(LMMAXD*LMMAXD*IEMXD*NSPIN,GMATN,1,GMATN_ALL,1)
      DO ISPIN=1,NSPIN
        DO IE=1,IELAST
          LLY_GRDT_ALL(IE,ISPIN) = LLY_GRDT(IE,ISPIN)
        ENDDO
      ENDDO

C
C
C=======================================================================
C     now take care that all processors get the results of
C     all energy points
C=======================================================================
C
      IF (EMPID.GT.1) THEN
C
CIE
      DO IE = 1, IELAST
C
C       identify sending proc:
          SEND  =
     +    EMYRANK(MAPBLOCK(EPROC(IE),1,EMPID,1,0,EMPID-1)+1,EMPIC)
C
C
CEMPI     begin loop over E-parallel processors
          DO EMPI = 1, EMPID
CEMPI
            RECV = EMYRANK(MAPBLOCK(EMPI,1,EMPID,1,0,EMPID-1)+1,EMPIC)


C
C         check if communication is necessary
            IF (RECV.EQ.SEND.AND.RECV.EQ.MYRANK) THEN
C
              DO ISPIN=1,NSPIN
                DO LM1=1,LMMAXD
                  DO LM2=1,LMMAXD
                    GMATN_ALL(LM1,LM2,IE,ISPIN)=GMATN(LM1,LM2,IE,ISPIN)
                  ENDDO
                ENDDO
                LLY_GRDT_ALL(IE,ISPIN)=LLY_GRDT(IE,ISPIN) 
              ENDDO
C
            ELSE
C
              IF (MYRANK.EQ.SEND) THEN
C
                DO ISPIN=1,NSPIN
                  CALL ZCOPY(LMMAXD*LMMAXD,GMATN(1,1,IE,ISPIN),1,
     +                       GESEND(1,ISPIN),1)
                  GESEND(LMMAXD*LMMAXD+1,ISPIN)=
     +                       LLY_GRDT(IE,ISPIN)
                ENDDO
C
                CALL MPI_SEND(GESEND,(LMMAXD*LMMAXD+1)*NSPIND,
     +                        MPI_DOUBLE_COMPLEX,
     +                        RECV,96,MPI_COMM_WORLD,IERR)
C
              ENDIF
C
              IF (MYRANK.EQ.RECV) THEN
C
                CALL MPI_RECV(GERECV,(LMMAXD*LMMAXD+1)*NSPIND,
     +                        MPI_DOUBLE_COMPLEX,
     +                        SEND,96,MPI_COMM_WORLD,STATUS,IERR)
C
                DO ISPIN=1,NSPIN
                  CALL ZCOPY(LMMAXD*LMMAXD,GERECV(1,ISPIN),1,
     +                       GMATN_ALL(1,1,IE,ISPIN),1)
                  LLY_GRDT_ALL(IE,ISPIN)=GERECV(LMMAXD*LMMAXD+1,ISPIN)
                ENDDO
C
              ENDIF
C
            ENDIF
C
CEMPI
          ENDDO
CEMPI
C
C
        IF (SMPID.EQ.2) THEN
          DO ISPIN=1,NSPIN
            DO LM1=1,LMMAXD
              DO LM2=1,LMMAXD
                GMATN(LM1,LM2,IE,ISPIN) = GMATN_ALL(LM1,LM2,IE,ISPIN)
              ENDDO
            ENDDO
            LLY_GRDT(IE,ISPIN) = LLY_GRDT_ALL(IE,ISPIN)
          ENDDO
        ENDIF
CIE
      ENDDO
CIE
C
      ENDIF
C
C=======================================================================
C=======================================================================
C
C
C
C
C=======================================================================
C     as a third step exchange information of different spin channels
C=======================================================================
C
      IF (SMPID.EQ.2) THEN
C
      DO SMPI=1, SMPID
C
        SEND = SMYRANK(MAPBLOCK(SMPI,1,SMPID,1,0,SMPID-1)+1,SMPIC)
        RECV = SMYRANK(
     +         MAPBLOCK(1+MOD(SMPI,2),1,SMPID,1,0,SMPID-1)+1,SMPIC)
C
        IF (MYRANK.EQ.SEND) THEN
C
          CALL ZCOPY(LMMAXD*LMMAXD*IEMXD,GMATN_ALL(1,1,1,SMPI),1,
     +               GSEND,1)
          CALL ZCOPY(IEMXD,LLY_GRDT_ALL(1,SMPI),1,
     +               GSEND(LMMAXD*LMMAXD*IEMXD+1),1)
C
          CALL MPI_SEND(GSEND,LMMAXD*LMMAXD*IEMXD+IEMXD,
     +                  MPI_DOUBLE_COMPLEX,
     +                  RECV,97,MPI_COMM_WORLD,IERR)
C
        ENDIF
C
        IF (MYRANK.EQ.RECV) THEN
C
          CALL MPI_RECV(GRECV,LMMAXD*LMMAXD*IEMXD+IEMXD,
     +                  MPI_DOUBLE_COMPLEX,
     +                  SEND,97,MPI_COMM_WORLD,STATUS,IERR)
C
          CALL ZCOPY(LMMAXD*LMMAXD*IEMXD,GRECV,1,
     +               GMATN_ALL(1,1,1,SMPI),1) 
          CALL ZCOPY(IEMXD,GRECV(LMMAXD*LMMAXD*IEMXD+1),1,
     +               LLY_GRDT_ALL(1,SMPI),1) 
C
        ENDIF
C
      ENDDO
C
      ENDIF
C
      RETURN
C
      END
