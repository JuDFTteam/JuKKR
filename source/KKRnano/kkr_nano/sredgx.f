      SUBROUTINE SREDGX(
     >                  ISPIN,NSPIN,
     >                  MYRANK,NROFNODES,
     >                  SMPIB,SMPIC,SMYRANK,SRANK,
     >                  GMATXIJ,
     <                  GXIJ_ALL)
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
      DOUBLE COMPLEX GMATXIJ(LMMAXD,LMMAXD,NXIJD,NSPIND),
     +               GXIJ_ALL(LMMAXD,LMMAXD,NXIJD,NSPIND)
C     .. local scalars ..
      INTEGER        IX,LM1,LM2,SEND,RECV
C     .. local arrays ..
      DOUBLE COMPLEX GSEND(LMMAXD,LMMAXD,NXIJD),
     +               GRECV(LMMAXD,LMMAXD,NXIJD)
C
C     .. MPI ..
      INTEGER, dimension(MPI_STATUS_SIZE) :: STATUS
      INTEGER        IERR,MAPBLOCK
C     .. N-MPI
      INTEGER        MYRANK,NROFNODES
C     .. LS-MPI
      INTEGER        SMYRANK(SMPID,NAEZD*LMPID*EMPID),
     +               SRANK(SMPID,NAEZD*LMPID*EMPID),
     +               SMPI,SMPIB,SMPIC
C
C
      DO ISPIN=1,NSPIN
        DO IX=1,NXIJD
          DO LM1=1,LMMAXD
            DO LM2=1,LMMAXD
              GXIJ_ALL(LM1,LM2,IX,ISPIN)=GMATXIJ(LM1,LM2,IX,ISPIN)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
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
          DO IX=1,NXIJD
            DO LM1=1,LMMAXD
              DO LM2=1,LMMAXD
                GSEND(LM1,LM2,IX)=GMATXIJ(LM1,LM2,IX,SMPI)
              ENDDO
            ENDDO
          ENDDO
C
          CALL MPI_SEND(GSEND,LMMAXD*LMMAXD*NXIJD,
     +                  MPI_DOUBLE_COMPLEX,
     +                  RECV,97,MPI_COMM_WORLD,IERR)
C
        ENDIF
C
        IF (MYRANK.EQ.RECV) THEN
C
          CALL MPI_RECV(GRECV,LMMAXD*LMMAXD*NXIJD,
     +                  MPI_DOUBLE_COMPLEX,
     +                  SEND,97,MPI_COMM_WORLD,STATUS,IERR)
C
          DO IX=1,NXIJD
            DO LM1=1,LMMAXD
              DO LM2=1,LMMAXD
                GXIJ_ALL(LM1,LM2,IX,SMPI)=GRECV(LM1,LM2,IX)
              ENDDO
            ENDDO
          ENDDO
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
