      SUBROUTINE EPRDIST(
     >                   IELAST,KMESH,NOFKS,
     >                   PRSC,SPRS,CNVFAC,
     >                   MYRANK,EMPIC,EMYRANK,
     >                   EPROC,EPROCO)
C ======================================================================
C
C this routine redistributes on demand the initial guess
C information to the individual processors EMPI=1,...,EMPID
C
C                               Alexander Thiess, 14th of December 2009
C =======================================================================
C
      IMPLICIT NONE
C
      include 'mpif.h'
C     .. Parameters ..
      include 'inc.p'
C
      INTEGER          LMMAXD
      PARAMETER        (LMMAXD= (LMAXD+1)**2)
      INTEGER          MAXMSHD
      PARAMETER        (MAXMSHD=8)
C
C     .. global scalars ..
      INTEGER          IELAST
C     .. global arrays ..
      COMPLEX          PRSC(NGUESSD*LMMAXD,EKMD)
      DOUBLE PRECISION CNVFAC(EKMD)
      INTEGER          SPRS(NGUESSD*LMMAXD+1,EKMD+1)
      INTEGER          EPROC(IEMXD),
     +                 EPROCO(IEMXD),
     +                 KMESH(IEMXD),NOFKS(MAXMSHD)
C     .. local scalars ..
      INTEGER          IDIM,IE,JEKM,RECV,SEND
C     .. local arrays ..

C
C     .. MPI ..
      INTEGER, dimension(MPI_STATUS_SIZE) :: STATUS
      INTEGER          IERR,MAPBLOCK
C     .. N-MPI
      INTEGER          MYRANK
C     .. E-MPI
      INTEGER          EMYRANK(EMPID,NAEZD*LMPID*SMPID),
     +                 EMPIC
C
C
      JEKM = 0
C
      DO IE= 1, IELAST
C
        IF (EPROC(IE).NE.EPROCO(IE)) THEN
C
          SEND  =
     +    EMYRANK(MAPBLOCK(EPROCO(IE),1,EMPID,1,0,EMPID-1)+1,EMPIC)
C
          RECV  =
     +    EMYRANK(MAPBLOCK(EPROC(IE),1,EMPID,1,0,EMPID-1)+1,EMPIC)
C
C
C         first exchange PRSC
C
          IDIM = NGUESSD*LMMAXD*NOFKS(KMESH(IE))
C
          IF (MYRANK.EQ.SEND) THEN
          CALL MPI_SEND(PRSC(1,JEKM+1),IDIM,
     +                  MPI_COMPLEX,
     +                  RECV,92,MPI_COMM_WORLD,IERR)
          ELSEIF (MYRANK.EQ.RECV) THEN
          CALL MPI_RECV(PRSC(1,JEKM+1),IDIM,
     +                  MPI_COMPLEX,
     +                  SEND,92,MPI_COMM_WORLD,STATUS,IERR)
          ENDIF
C
C
C         second step: exchange SPRS
C
          IDIM = (NGUESSD*LMMAXD+1)*NOFKS(KMESH(IE))
C
          IF (MYRANK.EQ.SEND) THEN
          CALL MPI_SEND(SPRS(1,JEKM+1),IDIM,
     +                  MPI_COMPLEX,
     +                  RECV,91,MPI_COMM_WORLD,IERR)
          ELSEIF (MYRANK.EQ.RECV) THEN
          CALL MPI_RECV(SPRS(1,JEKM+1),IDIM,
     +                  MPI_COMPLEX,
     +                  SEND,91,MPI_COMM_WORLD,STATUS,IERR)
          ENDIF
C
C
C         third step: exchange CNVFAC
C
          IDIM = NOFKS(KMESH(IE))
C
          IF (MYRANK.EQ.SEND) THEN
          CALL MPI_SEND(CNVFAC(JEKM+1),IDIM,
     +                  MPI_DOUBLE_PRECISION,
     +                  RECV,93,MPI_COMM_WORLD,IERR)
          ELSEIF (MYRANK.EQ.RECV) THEN
          CALL MPI_RECV(CNVFAC(JEKM+1),IDIM,
     +                  MPI_DOUBLE_PRECISION,
     +                  SEND,93,MPI_COMM_WORLD,STATUS,IERR)
          ENDIF
C
        ENDIF
C
        JEKM = JEKM + NOFKS(KMESH(IE))
C
      ENDDO
C
C
      RETURN
C
      END
