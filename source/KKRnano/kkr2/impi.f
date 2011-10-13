      SUBROUTINE IMPI(
     >                NAEZ,
     <                MYRANK,NROFNODES,
     <                LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE,
     <                LSMPIB,LSMPIC,LSRANK,LSMYRANK,
     <                SMPIB,SMPIC,SRANK,SMYRANK,
     <                EMPIB,EMPIC,ERANK,EMYRANK,
     <                MYACTVRANK,ACTVGROUP,ACTVCOMM,ACTVSIZE)
C ======================================================================
C                       build MPI_Groups 
C ======================================================================
C
C the required information for l-parallelization are
C intitialized.
C
C all depends on the additional switch in inc.p:
C
C MPI_L = 1  .... still one process per atom 
C MPI_L > 1  .... LMPI number of processes running per atom
C
C parallelization strategy:
C only the solution of the Dyson equation, resp. the inversion of
C the Lloyd-matrix are parallelized of LM. All tasks before are
C executed anyway, to save communication. After INVIT and LLYINVIT
C finished, the results for different L are communicated to the 
C root process by MPI_REDUCE(...,MPI_SUM,...) and from here on only
C the rout processes are running.
C
C                               Alexander Thiess, 19th of November 2009
C =======================================================================
C
      IMPLICIT NONE
C
      include 'mpif.h'
      include 'inc.p'
C
C     .. Parameters ..
      INTEGER      LSED
      PARAMETER    (LSED = LMPID*SMPID*EMPID)
      INTEGER      NSED
      PARAMETER    (NSED = NAEZD*SMPID*EMPID)
      INTEGER      NLSD
      PARAMETER    (NLSD = NAEZD*LMPID*SMPID)
      INTEGER      NLED
      PARAMETER    (NLED = NAEZD*LMPID*EMPID)
C
      INTEGER      I1,NAEZ,ISPIN,IRANK
C     .. MPI ..
      INTEGER      WGROUP,IERR
C     .. N-MPI
      INTEGER      MYRANK,NROFNODES
C     .. L-MPI
      INTEGER      MYLRANK(LSED),
     +             LRANKS(NAEZD,LSED),
     +             LCOMM(LSED),
     +             LGROUP(LSED),
     +             LSIZE(LSED),
     +             LMPI,LMPIC
C     .. LS-MPI
      INTEGER      LSMYRANK(LMPID,NSED),
     +             LSRANK(LMPID,NSED),
     +             LSMPI,LSMPIC,LSMPIB
C     .. S-MPI
      INTEGER      SMYRANK(SMPID,NLED),
     +             SRANK(SMPID,NLED),
     +             SMPI,SMPIC,SMPIB
C     .. E-MPI
      INTEGER      EMYRANK(EMPID,NLSD),
     +             ERANK(EMPID,NLSD),
     +             EMPI,EMPIC,EMPIB
C     .. ACTV-MPI
      INTEGER      MYACTVRANK,ACTVRANKS(NAEZD*LSED),
     +             ACTVCOMM,ACTVGROUP,ACTVSIZE
C
      EXTERNAL     MPI_INIT,
     +             MPI_COMM_RANK,
     +             MPI_COMM_SIZE
C
C intialize MPI and standard groups and communicator
C
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NROFNODES,IERR)
C
C define ranks for all non-standard groups of processors
C
      LMPIC = 0
C
      LSMPIC = 0
      LSMPIB = 0
C
      SMPIC = 0
      SMPIB = 0
C
      EMPIC = 0
      EMPIB = 0
C
C
      DO I1 = 1, NAEZ
        DO ISPIN = 1, SMPID
          DO LMPI = 1, LMPID
            DO EMPI = 1, EMPID
C
            IRANK = (ISPIN-1)*(NAEZ*LMPID*EMPID) +
     +              (LMPI-1)*NAEZ*EMPID +
     +              (I1-1)*EMPID +
     +              EMPI - 1
C
            LRANKS(I1,(ISPIN-1)*LMPID*EMPID+(LMPI-1)*EMPID+EMPI) = IRANK
C
            LSMYRANK(LMPI,(ISPIN-1)*NAEZ*EMPID+(I1-1)*EMPID+EMPI)= IRANK
            LSRANK(LMPI,(ISPIN-1)*NAEZ*EMPID+(I1-1)*EMPID+EMPI)  =LMPI-1
C
            SMYRANK(ISPIN,(LMPI-1)*NAEZ*EMPID+(I1-1)*EMPID+EMPI) = IRANK
            SRANK(ISPIN,(LMPI-1)*NAEZ*EMPID+(I1-1)*EMPID+EMPI)  =ISPIN-1
C
            EMYRANK(EMPI,(ISPIN-1)*NAEZ*LMPID+(I1-1)*LMPID+LMPI)= IRANK
            ERANK(EMPI,(ISPIN-1)*NAEZ*LMPID+(I1-1)*LMPID+LMPI)  =LMPI-1
C
            ACTVRANKS(IRANK+1)                = IRANK
C
            IF (MYRANK.EQ.IRANK) THEN
C
              LMPIC  = (ISPIN-1)*LMPID*EMPID+(LMPI-1)*EMPID + EMPI
C
              LSMPIC = (ISPIN-1)*NAEZ*EMPID +(I1-1)  *EMPID + EMPI
              LSMPIB = LMPI
C
              SMPIC  = (LMPI-1)*NAEZ*EMPID  +(I1-1)*EMPID   + EMPI
              SMPIB  = ISPIN
C
              EMPIC  = (ISPIN-1)*NAEZ*LMPID +(I1-1)*LMPID   + LMPI
              EMPIB  = EMPI
C
            ENDIF
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C build groups and communicators
C
C
C ACTIVE GROUP (ACTVGROUP) ...............................................
C
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD,WGROUP,IERR) ! get default
        CALL MPI_GROUP_INCL(WGROUP,NAEZ*LMPID*SMPID*EMPID,
     +                      ACTVRANKS(1),
     +                      ACTVGROUP,IERR) !create a group ACTVGROUP
        CALL MPI_COMM_CREATE(MPI_COMM_WORLD,ACTVGROUP,ACTVCOMM,IERR)

        MYACTVRANK=-2
        CALL MPI_GROUP_RANK(ACTVGROUP,MYACTVRANK,IERR)
        CALL MPI_GROUP_SIZE(ACTVGROUP,ACTVSIZE,IERR)
C
C LGROUP ...............................................................
C
      DO LMPI=1,LMPID*SMPID*EMPID
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD, WGROUP, IERR) ! get default
        CALL MPI_GROUP_INCL(WGROUP, NAEZ, LRANKS(1,LMPI),
     +                      LGROUP(LMPI), IERR) !create groups LGROUP
        CALL MPI_COMM_CREATE(MPI_COMM_WORLD, LGROUP(LMPI),
     +                       LCOMM(LMPI), IERR)!create communicators for
      ENDDO
C
      DO LMPI=1,LMPID*SMPID*EMPID
        MYLRANK(LMPI) = -2
        CALL MPI_GROUP_RANK(LGROUP(LMPI),MYLRANK(LMPI),IERR)
        CALL MPI_GROUP_SIZE(LGROUP(LMPI),LSIZE(LMPI),IERR)
      ENDDO
C
C.......................................................................
C
      IF (MYRANK.EQ.0) THEN
        WRITE(6,'(79(1H=))')
        WRITE(6,*) '  groups of processes created'
        WRITE(6,*) '  NMPI                    = ',NAEZ
        WRITE(6,*) '  LMPI                    = ',LMPID
        WRITE(6,*) '  SMPI                    = ',SMPID
        WRITE(6,*) '  EMPI                    = ',EMPID
        WRITE(6,*) '  NTHRDS                  = ',NTHRDS
        WRITE(6,*) '  MPI-processes           = ',NAEZ*LSED
        WRITE(6,*) '  total no. of processors = ',NAEZ*LSED*NTHRDS
        WRITE(6,'(79(1H=))')
      ENDIF
C
C
      RETURN
C
      END
