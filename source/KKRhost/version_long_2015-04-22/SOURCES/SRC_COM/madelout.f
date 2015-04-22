C **********************************************************************
      SUBROUTINE MADEL2OUT(IPRINT,NAEZ,LRECAMAD,LMPOTD,
     &                     NLEFTOFF,NRIGHTOFF,NLEFTALL,NRIGHTALL)
      IMPLICIT NONE
      INTEGER IPRINT,NAEZ,LMPOTD
      INTEGER LRECAMAD,NLEFTOFF,NRIGHTOFF,NLEFTALL,NRIGHTALL

      INTEGER LFMT,IQ1,IQ2,LM1,LM2
      DOUBLE PRECISION SMAT(6,6)
      DOUBLE PRECISION SMAT1(6,200),SMAT2(6,200)
      DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD)
      CHARACTER*80 FMT
      INTEGER IREC
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,99001) 'A(1,1)',MIN(6,NAEZ),'avmad.dat'
      FMT = ' '
      LFMT = 0
      DO IQ1 = 1,MIN(6,NAEZ)
         FMT = FMT(1:LFMT)//'------------'
         LFMT = LFMT + 12
      END DO
      WRITE (6,'(4X,A,/,8X," Inside the slab ",/,4X,A)') 
     &     (FMT(1:LFMT),IQ1=1,2)
C
      OPEN (69,ACCESS='direct',RECL=LRECAMAD,FILE='avmad.unformatted',
     +     FORM='unformatted')
      DO IQ1 = 1,MIN(6,NAEZ)
         DO IQ2 = 1,MIN(6,NAEZ)
            IREC = IQ2 + NAEZ*(IQ1-1) 
            READ(69,REC=IREC) AVMAD
            SMAT(IQ1,IQ2) = AVMAD(1,1)
         END DO
         DO IQ2 = 1,MIN(200,NLEFTALL)
            IREC = IQ2 + NLEFTALL*(IQ1-1) + NLEFTOFF
            READ(69,REC=IREC) AVMAD
            SMAT1(IQ1,IQ2) = AVMAD(1,1)
         END DO
         DO IQ2 = 1,MIN(200,NRIGHTALL)
            IREC = IQ2 + NRIGHTALL*(IQ1-1) + NRIGHTOFF
            READ(69,REC=IREC) AVMAD
            SMAT2(IQ1,IQ2) = AVMAD(1,1)
         END DO
      END DO
      CLOSE(69)
C
      DO IQ1 = 1,MIN(6,NAEZ)
         WRITE (6,99002) (SMAT(IQ1,IQ2),IQ2=1,MIN(6,NAEZ))
      END DO
      WRITE (6,'(4X,A,/,8X," Slab - left host",/,4X,A)') 
     &     (FMT(1:LFMT),IQ1=1,2)
      DO IQ2 = 1,MIN(200,NLEFTALL)
         WRITE (6,99002) (SMAT1(IQ1,IQ2),IQ1=1,MIN(6,NAEZ))
      END DO
      WRITE (6,'(4X,A,/,8X," Slab - right host",/,4X,A)') 
     &     (FMT(1:LFMT),IQ1=1,2)
      DO IQ2 = 1,MIN(200,NRIGHTALL)
         WRITE (6,99002) (SMAT2(IQ1,IQ2),IQ1=1,MIN(6,NAEZ))
      END DO
      WRITE (6,'(4X,A,/)') FMT(1:LFMT)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      IF ( IPRINT.LT.3 ) RETURN
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      OPEN (78,FILE='avmad.dat',FORM='formatted')
      WRITE (78,99003) 'A(IQ1,IQ2,LM1,LM2) '
C     
      OPEN (69,ACCESS='direct',RECL=LRECAMAD,FILE='avmad.unformatted',
     +     FORM='unformatted')
C
      WRITE (78,99006) ' Inside the slab '
      DO IQ1 = 1,NAEZ
         DO IQ2 = 1,NAEZ
            WRITE (78,99004) IQ1,IQ2
C
            IREC = IQ2 + NAEZ*(IQ1-1)
            READ (69,REC=IREC) AVMAD
C
            DO LM1 = 1,LMPOTD
               DO LM2 = 1,LMPOTD
                  IF ( ABS(AVMAD(LM1,LM2)).GT.1D-10 )WRITE (78,99005)
     &                 LM1,LM2,AVMAD(LM1,LM2)
               END DO
            END DO
            WRITE (78,'(33(1H-))')
         END DO
      END DO
      WRITE (78,99006) ' Slab - Left Host '
      DO IQ1 = 1,NAEZ
         DO IQ2 = 1,NLEFTALL
            WRITE (78,99004) IQ1,IQ2
C
            IREC = IQ2 + NLEFTALL*(IQ1-1) + NLEFTOFF
            READ (69,REC=IREC) AVMAD
C
            DO LM1 = 1,LMPOTD
               DO LM2 = 1,LMPOTD
                  IF ( ABS(AVMAD(LM1,LM2)).GT.1D-10 )WRITE (78,99005)
     &                 LM1,LM2,AVMAD(LM1,LM2)
               END DO
            END DO
            WRITE (78,'(33(1H-))')
         END DO
      END DO
      WRITE (78,99006) ' Slab - Right Host '
      DO IQ1 = 1,NAEZ
         DO IQ2 = 1,NRIGHTALL
            WRITE (78,99004) IQ1,IQ2
C
            IREC = IQ2 + NRIGHTALL*(IQ1-1) + NRIGHTOFF
            READ (69,REC=IREC) AVMAD
C
            DO LM1 = 1,LMPOTD
               DO LM2 = 1,LMPOTD
                  IF ( ABS(AVMAD(LM1,LM2)).GT.1D-10 )WRITE (78,99005)
     &                 LM1,LM2,AVMAD(LM1,LM2)
               END DO
            END DO
            WRITE (78,'(33(1H-))')
         END DO
      END DO
      CLOSE (69)
      CLOSE (78)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
99001 FORMAT (5X,'Madelung potential coefficients ',A,' up to NAEZ =',
     &        I2,/,5X,'(full matrix (nonzeros) in file ',A,
     &        ' if IPRINT.GT.2)')
99002 FORMAT (4X,6(D12.4))
99003 FORMAT (' Madelung potential coefficients ',A,/)
99004 FORMAT ('IQ1 =',I3,' IQ2 =',I3,/,17('-'),/,5X,
     &        '  LM1  LM2  A(LM1,LM2)')
99005 FORMAT (5X,2I5,1P,D18.10)
99006 FORMAT (70(1H*),/,10X,A,/,70(1H*))
      END 
C **********************************************************************
C
C **********************************************************************
      SUBROUTINE MADEL3OUT(IPRINT,NAEZ,LRECABMAD,SMAT1,SMAT2,LMPOTD)
      IMPLICIT NONE
      INTEGER IPRINT,NAEZ,LMPOTD
      INTEGER LRECABMAD
      INTEGER LFMT,IQ1,IQ2,LM1,LM2
      DOUBLE PRECISION SMAT1(6,6),SMAT2(6,6)
      DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD),BVMAD(LMPOTD)
      CHARACTER*80 FMT
      INTEGER IREC
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,99001) 'A(1,1)',MIN(6,NAEZ),'avmad.dat'
      FMT = ' '
      LFMT = 0
      DO IQ1 = 1,MIN(6,NAEZ)
         FMT = FMT(1:LFMT)//'------------'
         LFMT = LFMT + 12
      END DO
      WRITE (6,'(4X,A)') FMT(1:LFMT)
      DO IQ1 = 1,MIN(6,NAEZ)
         WRITE (6,99002) (SMAT1(IQ1,IQ2),IQ2=1,MIN(6,NAEZ))
      END DO
      WRITE (6,'(4X,A,/)') FMT(1:LFMT)
      WRITE (6,99001) 'B(1,1)',MIN(6,NAEZ),'bvmad.dat'
      WRITE (6,'(4X,A)') FMT(1:LFMT)
      DO IQ1 = 1,MIN(6,NAEZ)
         WRITE (6,99002) (SMAT2(IQ1,IQ2),IQ2=1,MIN(6,NAEZ))
      END DO
      WRITE (6,'(4X,A,/)') FMT(1:LFMT)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      IF ( IPRINT.LT.3 ) RETURN
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      OPEN (78,FILE='avmad.dat',FORM='formatted')
      OPEN (79,FILE='bvmad.dat',FORM='formatted')
      WRITE (78,99003) 'A(IQ1,IQ2,LM1,LM2) '
      WRITE (79,99003) 'B(IQ1,IQ2,LM)'
C
      OPEN (69,ACCESS='direct',RECL=LRECABMAD,FILE='abvmad.unformatted',
     &      FORM='unformatted')
      DO IQ1 = 1,NAEZ
         DO IQ2 = 1,NAEZ
            WRITE (78,99004) IQ1,IQ2
            WRITE (79,99005) IQ1,IQ2
C
            IREC = IQ2 + NAEZ*(IQ1-1)
            READ (69,REC=IREC) AVMAD,BVMAD
            DO LM1 = 1,LMPOTD
               DO LM2 = 1,LMPOTD
                  IF ( ABS(AVMAD(LM1,LM2)).GT.1D-10 ) WRITE (78,99006)
     &                 LM1,LM2,AVMAD(LM1,LM2)
               END DO
               IF ( ABS(BVMAD(LM1)).GT.1D-10 ) WRITE (79,99007)
     &              LM1,BVMAD(LM1)
            END DO
            WRITE (78,'(33(1H-))')
            WRITE (79,'(28(1H-))')
         END DO
      END DO
      CLOSE (69)
      CLOSE (78)
      CLOSE (79)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
99001 FORMAT (5X,'Madelung potential coefficients ',A,' up to NAEZ =',
     &        I2,/,5X,'(full matrix (nonzeros) in file ',A,
     &        ' if IPRINT.GT.2)')
99002 FORMAT (4X,6(D12.4))
99003 FORMAT (' Madelung potential coefficients ',A,/)
99004 FORMAT ('IQ1 =',I3,' IQ2 =',I3,/,17('-'),/,5X,
     &        '  LM1  LM2  A(LM1,LM2)')
99005 FORMAT ('IQ1 =',I3,' IQ2 =',I3,/,17('-'),/,5X,'   LM  B(LM)')
99006 FORMAT (5X,2I5,1P,D18.10)
99007 FORMAT (5X,I5,1P,D18.10)
      END
