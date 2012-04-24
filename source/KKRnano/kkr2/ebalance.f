      SUBROUTINE EBALANCE(
     >                    INFO,ITER,SCFSTEPS,
     >                    IERLAST,NPNT1,
     >                    MYACTVRANK,ACTVCOMM,
     >                    ETIME,EPROC,EPROCO,
C                         new input parameters after inc.p removal
     &                    empid, iemxd)
C ======================================================================
C                       build MPI_Groups 
C ======================================================================
C
C this routine balances the distribution of energy points
C to the set of processors of EMPI=1,...,EMPID
C
C
C                               Alexander Thiess, 7th of December 2009
C =======================================================================
C
      IMPLICIT NONE
C
      include 'mpif.h'

      INTEGER empid
      INTEGER iemxd
C
C     .. global scalars ..
      INTEGER        ITER,SCFSTEPS,IERLAST,NPNT1
      CHARACTER      INFO
C     .. global arrays ..
      REAL           ETIME(IEMXD)
      INTEGER        EPROC(IEMXD),
     +               EPROCO(IEMXD),
     +               EPOINT(IEMXD,EMPID)
C     .. local scalars ..
      REAL           AIMTM,SUMTM
      INTEGER        EMPI,IER,IERI,BLNCD,BLNCD1,IOS
      LOGICAL        LSAME,STOPIT,READIT,TEST
C     .. local arrays ..
      REAL           PROCTM(EMPID),
     +               MTIME(IEMXD)
C
C     .. MPI ACTV ..
      INTEGER        MYACTVRANK,ACTVCOMM,IERR
C
      EXTERNAL       LSAME
C
C
C=======================================================================
C 1st guess >>>
C=======================================================================
C     Trying to figure out what this code means:
C     *) if there is only one energy process -> it calculates everything
C     *) if there are two processes -> last two points calculated by
C        process 1, rest is calculated by process 2
C     *) 3 processes: 1st process calculates only last point
C                     2nd process calculates first 2/3 of all points
C                     3rd process calculates last 1/3 of points except
C                     last
C     *) MORE than 3 processes: USELESS! proc > 3 do nothing
C        except for DOS
      IF (LSAME(INFO,'I')) THEN
C
        IF (EMPID.EQ.1) THEN
          DO IER=1,IERLAST
            EPROC(IER)     = 1
          ENDDO
        ELSEIF (EMPID.EQ.2) THEN
          EPROC(IERLAST)   = 1
          EPROC(IERLAST-1) = 1
          DO IER=1,(IERLAST-2)
            EPROC(IER)     = 2
          ENDDO
        ELSE
          EPROC(IERLAST)   = 1
          IERI = FLOOR(REAL(IERLAST)/(REAL(3))*REAL(2))
          DO IER=1, IERI
            EPROC(IER) = 2
          ENDDO
          DO IER=(IERI+1), (IERLAST-1)
            EPROC(IER) = 3
          ENDDO
        ENDIF
C
C     if DOS shall be calculated make use of special scheme allowing for
C     broader Energy-parallelization
C
        IF (TEST('DOS     ')) THEN 
          DO IER=1,IERLAST
            EPROC(IER)     = MOD(IER,EMPID) + 1
          ENDDO
        ENDIF
        
C
C     check whether information on energy-point load-balancing is 
C     available
C
      INQUIRE(FILE='balance',EXIST=READIT)
C
        IF ((ITER.EQ.1).AND.(READIT)) THEN
C
          OPEN (50,FILE='balance',FORM='unformatted')
C
          READ (50,IOSTAT=IOS) BLNCD,BLNCD1
C
          IF (IOS.EQ.0.AND.BLNCD.EQ.IERLAST.AND.BLNCD1.EQ.EMPID) THEN
            READ (50) EPROC
          ENDIF
C
          CLOSE(50)
C
        ENDIF
C
      ENDIF
C=======================================================================
C >>> 1st guess
C=======================================================================
C
C
C=======================================================================
C use timing of ITER-1 >>>
C=======================================================================
      IF (LSAME(INFO,'R')) THEN
C
C     gather on all processors all timings
C
      CALL MPI_ALLREDUCE(ETIME,MTIME,IEMXD,MPI_REAL,MPI_MAX,
     +                   ACTVCOMM,IERR)
C
C
C----------------------------------------------------------------------
C     analyze balance of last iteration
C----------------------------------------------------------------------
C
      DO EMPI=1,EMPID
        PROCTM(EMPI)        = 0.0D0
      ENDDO
C
      SUMTM                 = 0.0D0
C
      DO IER=1,IERLAST
        SUMTM               = SUMTM + MTIME(IER)
        PROCTM(EPROC(IER))  = PROCTM(EPROC(IER)) + MTIME(IER)
C
        EPROCO(IER)         = EPROC(IER)
      ENDDO
C
      DO EMPI=1,EMPID
        IF (MYACTVRANK.EQ.0)
     +    WRITE(6,*) 'EMPID: group ',EMPI,' took',PROCTM(EMPI),'%:',
     +             REAL(100)*PROCTM(EMPI)/SUMTM
      ENDDO

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
C
C----------------------------------------------------------------------
C     intitialize
C----------------------------------------------------------------------
C
      SUMTM = 0.0D0
C
      DO IER = 1, IERLAST
        SUMTM         = SUMTM + MTIME(IER)
        DO EMPI = 1, EMPID
          EPOINT(IER,EMPI) = 0
        ENDDO
        EPROC(IER)    = 0
      ENDDO
C
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
C
C     goal for total time on each processor 1,...,EMPID
C
      AIMTM = SUMTM/EMPID
C
C     1st processor (EMPI=1) acts always on IER=IERLAST
C     and if more time is availible on energy points starting with NPNT1
C     2nd ... EMPID processors sum up timing from IER=1 up to IERLAST-1
C     taking time limits into account
C
C     loop over processors EMPI=1,EMPID >>>
      DO EMPI=1,EMPID
C
        PROCTM(EMPI)        = 0.0D0
C
        IF (EMPI.EQ.1) THEN
          PROCTM(1)         = MTIME(IERLAST)
          EPOINT(IERLAST,1) = IERLAST
          EPROC(IERLAST)    = 1
          IERI              = NPNT1
        ELSE
          IERI              = 1
        ENDIF
C
        IF (EMPID.EQ.1) THEN
          IERI              = 1
        ENDIF
C
C
C     Quick fix: if NPNT1=0 then IERI=0 -> leads to read out of bounds
C     at line 212   E.R.
      if (IERI .eq. 0) then
        IERI = 1
      end if

C     loop over processors IER=IERI,IERLAST >>>
        DO IER = IERI, IERLAST
          IF (EMPI.LT.EMPID) THEN
            IF ((PROCTM(EMPI).LT.AIMTM).AND.(EPROC(IER).EQ.0)) THEN
              PROCTM(EMPI)    = PROCTM(EMPI) + MTIME(IER)
              EPOINT(IER,EMPI)= IER
              EPROC(IER)      = EMPI
            ENDIF
          ELSEIF (EMPI.EQ.EMPID) THEN
            IF (EPROC(IER).EQ.0) THEN
              PROCTM(EMPI)    = PROCTM(EMPI) + MTIME(IER)
              EPOINT(IER,EMPI)= IER
              EPROC(IER)      = EMPI
            ENDIF
          ENDIF
        ENDDO
C     >>> loop over processors IER=IERI,IERLAST
C
      IF (MYACTVRANK.EQ.0)
     +  WRITE(6,*) 'EMPID: group ',EMPI,'will take',PROCTM(EMPI),'%:',
     +           REAL(100)*PROCTM(EMPI)/SUMTM
C
C
      ENDDO
C     >>> loop over processors EMPI=1,EMPID
C
C
C     write information on load-balancing to formatted file 'balance'
C
      INQUIRE(FILE='STOP',EXIST=STOPIT)
      IF ((ITER.EQ.SCFSTEPS).OR.(STOPIT)) THEN
C
        IF (MYACTVRANK.EQ.0) THEN
C
          OPEN (50,FILE='balance',FORM='unformatted')
          WRITE(50) IERLAST,EMPID
          WRITE(50) EPROC
          CLOSE(50)
C          DO IER = 1,IERLAST
C          WRITE(6,*) IER,MTIME(IER)
C          END DO
C
        ENDIF
C
      ENDIF
C
C
      ENDIF
C=======================================================================
C >>> use timing of ITER-1
C=======================================================================
C
      RETURN
C
      END
