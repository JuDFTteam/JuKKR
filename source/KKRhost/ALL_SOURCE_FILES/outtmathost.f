      SUBROUTINE OUTTMATHOST(ALAT,INS,KREL,KMROT,NSPIN,NAEZ,LMMAX,
     &                       BRAVAIS,RBASIS,QMTET,QMPHI,E2IN,TK,
     &                       NPOL,NPNT1,NPNT2,NPNT3)
C **********************************************************************
C *                                                                    *
C *  Writes out the header of the t-matrices decimation file           *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Arguments ..
      INTEGER INS,KREL,KMROT,NSPIN,NAEZ,LMMAX,NPOL,NPNT1,NPNT2,NPNT3
      DOUBLE PRECISION ALAT,E2IN,TK
      DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*),QMTET(*),QMPHI(*)
C     ..
C     .. Locals ..
      INTEGER I,IH
C ----------------------------------------------------------------------
      WRITE(6,'(5X,A,/)')
     &                 '< DECIOPT > : writing header of decimation file'
      OPEN (37,FILE='decifile',STATUS='unknown')
      WRITE (37,FMT=*) 'INVERSE T-MATRIX AND CMOMS'
      WRITE (37,FMT=99001)
      WRITE (37,FMT=99002) ALAT,NSPIN,NAEZ,LMMAX,INS,KREL,KMROT
      WRITE (37,FMT=99003) BRAVAIS
      IF ( KREL.EQ.0 ) THEN
         WRITE (37,FMT=99004)
         DO IH = 1,NAEZ
            WRITE (37,FMT=99006) (RBASIS(I,IH),I=1,3)
         END DO
      ELSE
         WRITE (37,FMT=99005)
         DO IH = 1,NAEZ
            WRITE (37,FMT=99007) (RBASIS(I,IH),I=1,3),
     &                           QMTET(IH),QMPHI(IH)
         END DO
      END IF
      WRITE (37,FMT=99008) E2IN,TK
      WRITE (37,FMT=99009) NPNT1,NPNT2,NPNT3,NPOL
      CLOSE (37)
C ----------------------------------------------------------------------
99001 FORMAT (' Vectors in lattice constant units',/,
     &        '                                 ')
99002 FORMAT ('ALAT=',F9.6,' NSPIN=',I2,'  NAEZ=',I3,' LMMAX=',I3,
     &        ' INS=',I1,' KREL=',I1,' KMROT=',I1)
99003 FORMAT ('BRAVAIS ',/,3F8.4,/,3F8.4,/,3F8.4)
99004 FORMAT ('RBASIS')
99005 FORMAT ('RBASIS',20X,'MAGNETISATION ANGLES THETA/PHI')
99006 FORMAT (3F8.4)
99007 FORMAT (3F8.4,2F9.4)
99008 FORMAT ('EF=',F10.6,' TEMP=',F10.4,' Kelvin')
99009 FORMAT ('N1=',I3,' N2=',I3,' N3=',I3,' NPOL=',I3)
      END
