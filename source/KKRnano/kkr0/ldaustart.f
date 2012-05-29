C=======================================================================
C initialization of LDA+U calculation consisting of
C
C    probe file wldau.unf if it matches the given set of parameter
C    NLDAU,LLDAU,ULDAU,JLDAU. In case of agreement of all arguments
C    leave file wldau.unf as it is. If disagreeing information exist
C    initialize WLDAU to zero and write WLDAU to wldau.unf.
C
C    important note: all arrays are hard-coded allowing maximal
C    one s-, one p-, one d-, and one f-orbital to be treated in LDA+U
C    fashion
C
C=======================================================================
C
      SUBROUTINE LDAUSTART(
     >                     I,NLDAU,LLDAU,ULDAU,JLDAU,LWLDAU,
C                          new input parameters after inc.p removal
     &                     LMAXD, NSPIND)
C
      IMPLICIT NONE

      INTEGER LMAXD
      INTEGER NSPIND

C global arrays ..
      INTEGER            LLDAU(LMAXD+1), ! quantum no. of affected orbital
     +                   LLDAU0(LMAXD+1)
      DOUBLE PRECISION   ULDAU(LMAXD+1), ! effective repulsion U
     +                   ULDAU0(LMAXD+1),
     +                   JLDAU(LMAXD+1), ! effective exchange J
     +                   JLDAU0(LMAXD+1)
C ..
C global scalars ..
      INTEGER            I,        ! atomic index
     +                   NLDAU,    ! number of orbitals affected by LDAU
     +                   NLDAU0
      LOGICAL            LWLDAU
C ..
C local arrays ..
C     DOUBLE PRECISION   WMLDAU(MMAXD,MMAXD,LMAXD1,NSPIND)
      DOUBLE PRECISION   WMLDAU(2*LMAXD + 1,2*LMAXD + 1,LMAXD+1, NSPIND)
C ..
C local scalars ..
      INTEGER            ILDAU
      LOGICAL            LINIT

      INTEGER             LMAXD1
      INTEGER             MMAXD

      MMAXD  = 2*LMAXD + 1
      LMAXD1 = LMAXD + 1

      LINIT = .TRUE.
C
      IF (LWLDAU) THEN
        READ(65,REC=I) NLDAU0,LLDAU0,ULDAU0,JLDAU0,WMLDAU
        IF (NLDAU.NE.NLDAU0) LINIT = .FALSE.
        IF (LINIT) THEN
          DO ILDAU=1,NLDAU
            IF (LLDAU(ILDAU).NE.LLDAU0(ILDAU)) LINIT = .FALSE.
            IF (ULDAU(ILDAU).NE.ULDAU0(ILDAU)) LINIT = .FALSE.
            IF (JLDAU(ILDAU).NE.JLDAU0(ILDAU)) LINIT = .FALSE.
          ENDDO
        ENDIF
      ELSE
        LINIT = .FALSE.
      ENDIF
C
C
C
      IF (.NOT.LINIT) THEN
C
        CALL RINIT(MMAXD*MMAXD*NSPIND*LMAXD1,WMLDAU)
        WRITE(6,*) 'set WMLDAU to zero and write to file wldau.unf'
C
        WRITE(65,REC=I) NLDAU,LLDAU,ULDAU,JLDAU,WMLDAU
C
      ELSE
        WRITE(6,*) 'use existing WMLDAU from file wldau.unf'
      ENDIF
C
C
      END
C
