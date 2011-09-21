C     Initialises a LDA+U calculation
C     Read file 'ldauinfo' and write 'wldau.unf',
C     if this file does not already exist with the correct set of
C     parameters (given in ldauinfo)
      SUBROUTINE ldauinfo_read(LMAXD, NSPIND, Z, NAEZD)

      IMPLICIT NONE
c ----------------------------------------------------------------------
c     Parameters
c ----------------------------------------------------------------------

C     all are intent(in)
      INTEGER, INTENT(IN) :: LMAXD
      INTEGER, INTENT(IN) :: NSPIND
      INTEGER, INTENT(IN) :: NAEZD
      DOUBLE PRECISION, INTENT(IN) :: Z(NAEZD)

C     local variables
      DOUBLE PRECISION ULDAU(LMAXD+1),JLDAU(LMAXD+1)
      CHARACTER(len=10) TXTLDAU(0:3)
      INTEGER LRECLDAU
      
      INTEGER LMAXD1
      INTEGER MMAXD
      INTEGER NLDAU
      INTEGER LLDAU(LMAXD+1)
      LOGICAL LWLDAU,LLDAUINFO

      INTEGER I,J

      LMAXD1 = LMAXD + 1
      MMAXD  = 2*LMAXD + 1
c
        TXTLDAU(0) = 's'
        TXTLDAU(1) = 'p'
        TXTLDAU(2) = 'd'
        TXTLDAU(3) = 'f'
c
c verify that file 'ldauinfo' exists
c
        INQUIRE(FILE='ldauinfo',EXIST=LLDAUINFO)
        IF (.NOT.LLDAUINFO) THEN
          WRITE(6,*) 'file ldauinfo has to be provided!'
          CALL RCSTOP('LDAUINFO')
        ENDIF
c
c open 'ldauinfo'
c
        OPEN(77,FILE='ldauinfo',FORM='formatted')
        WRITE(6,2100)
        WRITE(6,*) 'LDA:'
c
c-------------------------------------------------------------
c    read L,U, and J for each atom from ldauinfo
c-------------------------------------------------------------
c
        LWLDAU = .false.
        INQUIRE(FILE='wldau.unf',EXIST=LWLDAU)
c
c determine rec-length for file 'wldau.unf', opened here and accessed
c in LDAUSTART
c
        LRECLDAU = 4*(1+LMAXD1)                   ! NLDAU & LLDAU
     &           + 8*2*LMAXD1                     ! ULDAU & JLDAU
     &           + 8*MMAXD*MMAXD*NSPIND*LMAXD1    ! WMLDAU
c
        OPEN (65,ACCESS='direct',RECL=LRECLDAU,FILE='wldau.unf',
     &        FORM='unformatted')
c
c
        DO I=1,NAEZD
          READ (77,FMT=*) NLDAU,(LLDAU(J), J=1,NLDAU),
     &                   (ULDAU(J),J=1,NLDAU),
     &                   (JLDAU(J), J=1,NLDAU)
          IF (NLDAU.GT.4) CALL RCSTOP('NLDAU')
          IF (NLDAU.GT.0) THEN
            WRITE(6,*) 'atom=',I,' (Z=',INT(Z(I)),') with ',NLDAU,
     &                 ' Coulomb rep. coeff.'
            DO J=1,NLDAU
              WRITE(6,*) TXTLDAU(LLDAU(J)),ULDAU(J),JLDAU(J)
            ENDDO
          ENDIF
c
c initialize or read matrix WMLDAU
c
          CALL LDAUSTART(I,NLDAU,LLDAU,ULDAU,JLDAU,LWLDAU)
c
        ENDDO
c
c
        CLOSE(65)
c
c-------------------------------------------------------------

        WRITE(6,2100)
        CLOSE(77)
C
c
c-------------------------------------------------------------
c-------------------------------------------------------------
c end of initialization of LDAU arrays
c-------------------------------------------------------------
c-------------------------------------------------------------
c

 2100 FORMAT(79(1H-))

      END SUBROUTINE


