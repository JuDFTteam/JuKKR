      SUBROUTINE KKRJIJ(
     >                  BZKP,VOLCUB,KPT,
     >                  NSYMAT,NAEZ,I3,
     >                  NXIJ,IXCP,ZKRXIJ,
     >                  GLLKE1,
     <                  GSXIJ,
     >                  communicator, comm_size,
C                       new input parameters after inc.p removal
     &                  lmmaxd, nxijd)

      IMPLICIT NONE
C ------------------------------------------------------------------------
C a) performs scattering of the off-diagonal GLLKE-elements
C    required to calculate Jij's
C b) multiplication with appropriate exp-factor, k-dependent
C                                                       A. Thiess Sep'09
C    "simplification": Elias Rabel, 2012
C ------------------------------------------------------------------------
      include 'mpif.h'

      integer communicator
      integer comm_size
      INTEGER lmmaxd
      INTEGER nxijd

      INTEGER          NSYMAXD
      PARAMETER        (NSYMAXD=48)

      DOUBLE COMPLEX   CZERO
      PARAMETER        (CZERO=(0.0D0,0.0D0))
      DOUBLE COMPLEX   CIONE
      PARAMETER        (CIONE  = ( 0.0D0,-1.0D0))
c     ..
c     .. GLOBAL SCALAR ARGUMENTS ..
      INTEGER          NAEZ,          ! number of atoms per unit cell
     +                 NSYMAT,        ! number of active symmetries
     +                 KPT,           ! k-point - run. index in kkrmat01
     +                 NXIJ           ! number of atoms in Jij-cluster
c     ..

c     .. ARRAY ARGUMENTS ..
      DOUBLE COMPLEX   GSXIJ(LMMAXD,LMMAXD,NSYMAXD,NXIJD)

      DOUBLE PRECISION BZKP(3,*),   ! TODO: better pass only 1 k-point
     +                 VOLCUB(*),
     +                 ZKRXIJ(48,3,NXIJD) ! position of atoms and sites
                                          ! connected by symmetry
      INTEGER          IXCP(NXIJD)      ! corresp. to Jij no. XIJ on the
                                        ! real space lattice

      DOUBLE COMPLEX   GLLKE1(NAEZ * lmmaxd, lmmaxd)

C     .. LOCAL ARRAYS ..
C     .. Fortran 90 automatic arrays
      DOUBLE COMPLEX   GSEND(lmmaxd, lmmaxd)
      DOUBLE COMPLEX   GRECV(lmmaxd, lmmaxd)
      DOUBLE COMPLEX   GXIJ (lmmaxd, lmmaxd, NXIJD)
      DOUBLE COMPLEX   EKRXIJ(48, NXIJD)

      INTEGER          IXCPS(NXIJD)       ! copydummy for IXCP used MPI

C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC ATAN,EXP


C     .. LOCAL SCALARS ..
      INTEGER          I3,I5,IV,ISYM,LM1,LM2,LM,ILM,XIJ
      DOUBLE COMPLEX   CARG

C     .. MPI ..
      INTEGER, dimension(MPI_STATUS_SIZE) :: STATUS

C     .. N-MPI
      INTEGER          MAPBLOCK,IERR

c ------------------------------------------------------------------------

C ================================================================
C       XCCPL communicate off-diagonal elements
C ================================================================

        DO I5 = 1, NAEZ
C .....
          DO XIJ = 1, NXIJD
C ....
            IXCPS(XIJ) = 0
            IF (XIJ.LE.NXIJ) IXCPS(XIJ) = IXCP(XIJ)
C ....
          ENDDO
C
C         broadcast IXCPS from processor of atom I5 to all
C
          CALL MPI_BCAST(IXCPS,NXIJD,MPI_INTEGER,
     +                   MAPBLOCK(I5,1,NAEZ,1,0,comm_size-1),
     +                   communicator,IERR)
C
          CALL MPI_BARRIER(communicator,IERR)
C
          DO XIJ = 2, NXIJ
C ....
C           if local GLL(k,E) is required just copy
C
            IF (IXCPS(XIJ).EQ.I5) THEN
C ...
              IF (I3.EQ.I5) THEN
C ..
                DO LM = 1, LMMAXD
C .
                  ILM = LMMAXD*(I3-1) + 1
                  CALL ZCOPY(LMMAXD,GLLKE1(ILM,LM),1,
     +                       GXIJ(1,LM,XIJ),1)
C .
                ENDDO
C ..
              ENDIF
C ...
C           required GLL(k,E) send by IXCPS(XIJ).EQ.I3
C           and recieved by I5.EQ.I3
C
            ELSE
C ...
              IF (IXCPS(XIJ).EQ.I3) THEN
C ..
                DO LM = 1, LMMAXD
C .
                  ILM = LMMAXD*(I5-1) + 1
                  CALL ZCOPY(LMMAXD,GLLKE1(ILM,LM),1,
     +                       GSEND(1,LM),1)
C .
                ENDDO
C
           CALL MPI_SEND(GSEND,LMMAXD*LMMAXD,MPI_DOUBLE_COMPLEX,
     +                   MAPBLOCK(I5,1,NAEZ,1,0,comm_size-1),
     +                   99,communicator,IERR)
C ..
              ENDIF


              IF (I5.EQ.I3) THEN
C ..
           CALL MPI_RECV(GRECV,LMMAXD*LMMAXD,MPI_DOUBLE_COMPLEX,
     +          MAPBLOCK(IXCPS(XIJ),1,NAEZ,1,0,comm_size-1),
     +          99,communicator,STATUS,IERR)

                DO LM = 1, LMMAXD
C .
                  CALL ZCOPY(LMMAXD,GRECV(1,LM),1,
     +                       GXIJ(1,LM,XIJ),1)
C .
                ENDDO
C ..
              ENDIF

           CALL MPI_BARRIER(communicator,IERR)
C ...
            ENDIF
C ....
          ENDDO
C .....
        ENDDO

C ================================================================
C       XCCPL calculate exponential factor and multiply with GXIJ
C ================================================================

        DO XIJ = 2, NXIJ
          DO ISYM  = 1,NSYMAT

            CARG = CZERO
            DO IV = 1,3
              CARG =  CARG + ZKRXIJ(ISYM,IV,XIJ)*BZKP(IV,KPT)
            ENDDO

            EKRXIJ(ISYM,XIJ) =
     +      VOLCUB(KPT) * EXP(CARG*(CIONE*8.D0*ATAN(1.D0)))

          ENDDO
        ENDDO

        DO XIJ = 2, NXIJ
          DO ISYM = 1,NSYMAT
            DO LM2=1,LMMAXD
              DO LM1=1,LMMAXD
                GSXIJ(LM1,LM2,ISYM,XIJ) =
     &          GSXIJ(LM1,LM2,ISYM,XIJ) +
     &          EKRXIJ(ISYM,XIJ) * GXIJ(LM1,LM2,XIJ)
              ENDDO
            ENDDO
          ENDDO             ! ISYM = 1,NSYMAT
        ENDDO

C ================================================================
C       XCCPL calculated integrated over KPT GSXIJ for off-diag-
C       onal elements
C       now back to KLOOPZ1 ..
C ================================================================

      RETURN

      END
