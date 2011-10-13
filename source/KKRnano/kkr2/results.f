C
C process with MYLRANK(LMPIC).EQ.0 and LMPIC.EQ.1 writes results
C
      SUBROUTINE RESULTS(LRECRES2,IELAST,ITSCF,LMAX,NAEZ,NPOL,NSPIN,
     +     KPRE,KTE,LPOT,E1,E2,TK,EFERMI,ALAT,ITITLE,CHRGNT,ZAT,EZ,WEZ,
     +     LDAU)

      IMPLICIT NONE
      INCLUDE 'inc.p'

C     .. Parameters ..
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NAEZD)
      INTEGER LMAXD1
      PARAMETER (LMAXD1=LMAXD+1)
C     ..
C     .. Local Scalars ..
      INTEGER IELAST,ITSCF,LMAX,NAEZ,NPOL,NSPIN
      INTEGER KPRE,KTE
      INTEGER I1,ISPIN,LPOT
      INTEGER LRECRES1,LRECRES2
C     ..
      DOUBLE PRECISION E1,E2,TK,EFERMI
      DOUBLE PRECISION CHRGNT,TOTSMOM,ALAT,PI
C     ..
      LOGICAL TEST,LDAU
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX EZ(IEMXD),WEZ(IEMXD)
      DOUBLE COMPLEX DEN(0:LMAXD1,IEMXD,NSPIND)
      DOUBLE PRECISION DOSTOT(0:LMAXD1,2)
      DOUBLE PRECISION ECOU(0:LPOTD),EPOTIN,EULDAU,EDCLDAU,
     +                 ESPC(0:3,NSPIND),ESPV(0:LMAXD1,NSPIND),
     +                 EXC(0:LPOTD)
      DOUBLE PRECISION ECORE(20,2)
      DOUBLE PRECISION ZAT(NAEZD)
      DOUBLE PRECISION CHARGE(0:LMAXD1,2)
C     ..
      DOUBLE PRECISION CATOM(NSPIND),QC
      DOUBLE PRECISION VMAD
C     ,,
      INTEGER ITITLE(20,NPOTD)
      INTEGER LCOREMAX
C
      PI = 4.0D0*ATAN(1.0D0)
C
      LRECRES1 = 8*43 + 16*(LMAXD+2)
C
C
      IF (NPOL.EQ.0 .OR. TEST('DOS     ')) THEN 
        LRECRES1 = LRECRES1 + 32*(LMAXD+2)*IEMXD
      END IF
C
C
      OPEN (71,ACCESS='direct',RECL=LRECRES1,FILE='results1',
     +      FORM='unformatted')
C
C
      DO I1 = 1,NAEZ
        IF (NPOL.EQ.0 .OR. TEST('DOS     ')) THEN 
          READ(71,REC=I1) QC,CATOM,CHARGE,ECORE,DEN
        ELSE
          READ(71,REC=I1) QC,CATOM,CHARGE,ECORE
        END IF
        CALL WRMOMS(NAEZ,NSPIN,CHARGE,I1,LMAXD,LMAXD1)
      END DO
C
C
      IF (NPOL.EQ.0 .OR. TEST('DOS     ')) THEN 
        DO I1 = 1,NAEZ
          READ(71,REC=I1) QC,CATOM,CHARGE,ECORE,DEN
          CALL WRLDOS(DEN,EZ,WEZ,
     +                LMAXD1,IEMXD,NPOTD,ITITLE,EFERMI,E1,E2,ALAT,TK,
     +                NSPIN,NAEZ,IELAST,I1,DOSTOT)
        END DO
      END IF
C
C
      TOTSMOM = 0.0D0
      DO I1 = 1,NAEZ
        IF (NPOL.EQ.0 .OR. TEST('DOS     ')) THEN 
          READ(71,REC=I1) QC,CATOM,CHARGE,ECORE,DEN
        ELSE
          READ(71,REC=I1) QC,CATOM,CHARGE,ECORE
        END IF
        DO ISPIN = 1,NSPIN
          IF (ISPIN.NE.1) THEN
            WRITE (6,FMT=9011) CATOM(ISPIN)
          ELSE
            WRITE (6,FMT=9001) I1,CATOM(ISPIN)
          END IF
        END DO
        WRITE (6,FMT=9041) ZAT(I1),QC
        IF (NSPIN.EQ.2) TOTSMOM = TOTSMOM + CATOM(NSPIN)
      END DO
      WRITE(6,'(79(1H+))')
      WRITE (6,FMT=9021) ITSCF,CHRGNT
      IF (NSPIN.EQ.2) WRITE (6,FMT=9031) TOTSMOM
      WRITE(6,'(79(1H+))')
C
C
 9001 FORMAT ('  Atom ',I4,' charge in wigner seitz cell =',f10.6)
 9011 FORMAT (7X,'spin moment in wigner seitz cell =',f10.6)
 9021 FORMAT ('      ITERATION',I4,
     &     ' charge neutrality in unit cell = ',f12.6)
 9031 FORMAT ('                   ',
     &     ' TOTAL mag. moment in unit cell = ',f12.6)
 9041 FORMAT (4X,'nuclear charge  ',F10.6,9X,'core charge =   ',F10.6)
C        WRITE(6,FMT=99001)
C        WRITE(6,FMT=99002)
C99001 FORMAT (79(1H=),/,18X,' MADELUNG POTENTIALS ',
C     &        '(spherically averaged) ')
C99002 FORMAT (/,25X,' ATOM ','  Delta_Q  ','     VMAD',/,25X,30(1H-))
C99003 FORMAT (25X,I4,2X,F10.6,1X,F12.6)
        CLOSE(71)
C
C=======================================================================
C output of information stored in 'results2'
C set KTE=1 in inputcard for output of energy contributions
C=======================================================================
C
      OPEN (72,ACCESS='direct',RECL=LRECRES2,FILE='results2',
     +     FORM='unformatted')
C
      DO I1 = 1,NAEZ
        READ(72,REC=I1) CATOM,VMAD,ECOU,EPOTIN,ESPC,ESPV,EXC,LCOREMAX,
     +                  EULDAU,EDCLDAU
C        WRITE (6,FMT=99003) I1,(CATOM(1)-ZAT(I1)),VMAD
      END DO
C      WRITE(6,'(25X,30(1H-),/)')
C      WRITE(6,'(79(1H=))')
C
      IF (KTE.EQ.1) THEN 
        DO I1 = 1,NAEZ
          READ(72,REC=I1) CATOM,VMAD,ECOU,EPOTIN,ESPC,ESPV,EXC,LCOREMAX,
     +                    EULDAU,EDCLDAU
          CALL ETOTB1(ECOU,E2,EPOTIN,ESPC,ESPV,EXC,
     &                EULDAU,EDCLDAU,LDAU,
     &                KPRE,LMAX,LPOT,
     &                LCOREMAX,NSPIN,I1,NAEZ)
        END DO
      END IF
C
      CLOSE(72)
C
C=======================================================================
C
      RETURN

      END
