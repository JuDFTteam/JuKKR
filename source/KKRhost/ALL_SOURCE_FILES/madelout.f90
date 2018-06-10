! **********************************************************************
SUBROUTINE madel2out(iprint,naez,lrecamad,lmpotd,  &
    nleftoff,nrightoff,nleftall,nrightall)
      IMPLICIT NONE
      INTEGER IPRINT,NAEZ,LMPOTD
      INTEGER LRECAMAD,NLEFTOFF,NRIGHTOFF,NLEFTALL,NRIGHTALL

      INTEGER LFMT,IQ1,IQ2,LM1,LM2
      DOUBLE PRECISION SMAT(6,6)
      DOUBLE PRECISION SMAT1(6,200),SMAT2(6,200)
      DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD)
      CHARACTER*80 FMT
      INTEGER IREC

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,99001) 'A(1,1)',MIN(6,naez),'avmad.dat'

FMT = ' '
lfmt = 0
DO iq1 = 1,MIN(6,naez)
  FMT = FMT(1:lfmt)//'------------'
  lfmt = lfmt + 12
END DO
WRITE (1337,'(4X,A,/,8X," Inside the slab ",/,4X,A)') (FMT(1:lfmt),iq1=1,2)

OPEN (69,ACCESS='direct',RECL=lrecamad,FILE='avmad.unformatted',  &
    FORM='unformatted')
DO iq1 = 1,MIN(6,naez)
  DO iq2 = 1,MIN(6,naez)
    irec = iq2 + naez*(iq1-1)
    READ(69,REC=irec) avmad
    smat(iq1,iq2) = avmad(1,1)
  END DO
  DO iq2 = 1,MIN(200,nleftall)
    irec = iq2 + nleftall*(iq1-1) + nleftoff
    READ(69,REC=irec) avmad
    smat1(iq1,iq2) = avmad(1,1)
  END DO
  DO iq2 = 1,MIN(200,nrightall)
    irec = iq2 + nrightall*(iq1-1) + nrightoff
    READ(69,REC=irec) avmad
    smat2(iq1,iq2) = avmad(1,1)
  END DO
END DO
CLOSE(69)

DO iq1 = 1,MIN(6,naez)
  WRITE (1337,99002) (smat(iq1,iq2),iq2=1,MIN(6,naez))
END DO
WRITE (1337,'(4X,A,/,8X," Slab - left host",/,4X,A)') (FMT(1:lfmt),iq1=1,2)
DO iq2 = 1,MIN(200,nleftall)
  WRITE (1337,99002) (smat1(iq1,iq2),iq1=1,MIN(6,naez))
END DO
WRITE (1337,'(4X,A,/,8X," Slab - right host",/,4X,A)') (FMT(1:lfmt),iq1=1,2)
DO iq2 = 1,MIN(200,nrightall)
  WRITE (1337,99002) (smat2(iq1,iq2),iq1=1,MIN(6,naez))
END DO
WRITE (1337,'(4X,A,/)') FMT(1:lfmt)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

IF ( iprint < 3 ) RETURN

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
OPEN (78,FILE='avmad.dat',FORM='formatted')
WRITE (78,99003) 'A(IQ1,IQ2,LM1,LM2) '

OPEN (69,ACCESS='direct',RECL=lrecamad,FILE='avmad.unformatted',  &
    FORM='unformatted')

WRITE (78,99006) ' Inside the slab '
DO iq1 = 1,naez
  DO iq2 = 1,naez
    WRITE (78,99004) iq1,iq2
    
    irec = iq2 + naez*(iq1-1)
    READ (69,REC=irec) avmad
    
    DO lm1 = 1,lmpotd
      DO lm2 = 1,lmpotd
        IF ( ABS(avmad(lm1,lm2)) > 1D-10 )WRITE (78,99005)  &
            lm1,lm2,avmad(lm1,lm2)
      END DO
    END DO
    WRITE (78,'(33(1H-))')
  END DO
END DO
WRITE (78,99006) ' Slab - Left Host '
DO iq1 = 1,naez
  DO iq2 = 1,nleftall
    WRITE (78,99004) iq1,iq2
    
    irec = iq2 + nleftall*(iq1-1) + nleftoff
    READ (69,REC=irec) avmad
    
    DO lm1 = 1,lmpotd
      DO lm2 = 1,lmpotd
        IF ( ABS(avmad(lm1,lm2)) > 1D-10 )WRITE (78,99005)  &
            lm1,lm2,avmad(lm1,lm2)
      END DO
    END DO
    WRITE (78,'(33(1H-))')
  END DO
END DO
WRITE (78,99006) ' Slab - Right Host '
DO iq1 = 1,naez
  DO iq2 = 1,nrightall
    WRITE (78,99004) iq1,iq2
    
    irec = iq2 + nrightall*(iq1-1) + nrightoff
    READ (69,REC=irec) avmad
    
    DO lm1 = 1,lmpotd
      DO lm2 = 1,lmpotd
        IF ( ABS(avmad(lm1,lm2)) > 1D-10 )WRITE (78,99005)  &
            lm1,lm2,avmad(lm1,lm2)
      END DO
    END DO
    WRITE (78,'(33(1H-))')
  END DO
END DO
CLOSE (69)
CLOSE (78)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

99001 FORMAT (5X,'Madelung potential coefficients ',a,' up to NAEZ =',  &
    i2,/,5X,'(full matrix (nonzeros) in file ',a, ' if IPRINT.GT.2)')
99002 FORMAT (4X,6(d12.4))
99003 FORMAT (' Madelung potential coefficients ',a,/)
99004 FORMAT ('IQ1 =',i3,' IQ2 =',i3,/,17('-'),/,5X, '  LM1  LM2  A(LM1,LM2)')
99005 FORMAT (5X,2I5,1P,d18.10)
99006 FORMAT (70(1H*),/,10X,a,/,70(1H*))
END SUBROUTINE madel2out
! **********************************************************************

! **********************************************************************

SUBROUTINE madel3out(iprint,naez,lrecabmad,smat1,smat2,lmpotd)

      IMPLICIT NONE
      INTEGER IPRINT,NAEZ,LMPOTD
      INTEGER LRECABMAD
      INTEGER LFMT,IQ1,IQ2,LM1,LM2
      DOUBLE PRECISION SMAT1(6,6),SMAT2(6,6)
      DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD),BVMAD(LMPOTD)
      CHARACTER*80 FMT
      INTEGER IREC

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,99001) 'A(1,1)',MIN(6,naez),'avmad.dat'

FMT = ' '
lfmt = 0
DO iq1 = 1,MIN(6,naez)
  FMT = FMT(1:lfmt)//'------------'
  lfmt = lfmt + 12
END DO
WRITE (1337,'(4X,A)') FMT(1:lfmt)
DO iq1 = 1,MIN(6,naez)
  WRITE (1337,99002) (smat1(iq1,iq2),iq2=1,MIN(6,naez))
END DO
WRITE (1337,'(4X,A,/)') FMT(1:lfmt)
WRITE (1337,99001) 'B(1,1)',MIN(6,naez),'bvmad.dat'
WRITE (1337,'(4X,A)') FMT(1:lfmt)
DO iq1 = 1,MIN(6,naez)
  WRITE (1337,99002) (smat2(iq1,iq2),iq2=1,MIN(6,naez))
END DO
WRITE (1337,'(4X,A,/)') FMT(1:lfmt)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

IF ( iprint < 3 ) RETURN

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
OPEN (78,FILE='avmad.dat',FORM='formatted')
OPEN (79,FILE='bvmad.dat',FORM='formatted')
WRITE (78,99003) 'A(IQ1,IQ2,LM1,LM2) '
WRITE (79,99003) 'B(IQ1,IQ2,LM)'

OPEN (69,ACCESS='direct',RECL=lrecabmad,FILE='abvmad.unformatted',  &
    FORM='unformatted')
DO iq1 = 1,naez
  DO iq2 = 1,naez
    WRITE (78,99004) iq1,iq2
    WRITE (79,99005) iq1,iq2
    
    irec = iq2 + naez*(iq1-1)
    READ (69,REC=irec) avmad,bvmad
    DO lm1 = 1,lmpotd
      DO lm2 = 1,lmpotd
        IF ( ABS(avmad(lm1,lm2)) > 1D-10 ) WRITE (78,99006)  &
            lm1,lm2,avmad(lm1,lm2)
      END DO
      IF ( ABS(bvmad(lm1)) > 1D-10 ) WRITE (79,99007) lm1,bvmad(lm1)
    END DO
    WRITE (78,'(33(1H-))')
    WRITE (79,'(28(1H-))')
  END DO
END DO
CLOSE (69)
CLOSE (78)
CLOSE (79)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

99001 FORMAT (5X,'Madelung potential coefficients ',a,' up to NAEZ =',  &
    i2,/,5X,'(full matrix (nonzeros) in file ',a, ' if IPRINT.GT.2)')
99002 FORMAT (4X,6(d12.4))
99003 FORMAT (' Madelung potential coefficients ',a,/)
99004 FORMAT ('IQ1 =',i3,' IQ2 =',i3,/,17('-'),/,5X, '  LM1  LM2  A(LM1,LM2)')
99005 FORMAT ('IQ1 =',i3,' IQ2 =',i3,/,17('-'),/,5X,'   LM  B(LM)')
99006 FORMAT (5X,2I5,1P,d18.10)
99007 FORMAT (5X,i5,1P,d18.10)
END SUBROUTINE madel3out
