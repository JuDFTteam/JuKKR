C*==scfiterang.f    processed by SPAG 6.05Rc at 09:19 on 18 Apr 2003
      SUBROUTINE SCFITERANG(ITRSCF,ITOQ,FACT,MVPHI,MVTET,MVGAM,
     &                     QMPHI,QMTET,QMGAM,NQ,NK,ERRAVANG,
     &                     NQMAX,NTMAX,NMVECMAX,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   * applied to mix the angles specifying                             *
C   * the local frame of reference                                     *
C   *                                                                  *
C   ********************************************************************
      use mod_wunfiles, only: t_params
      IMPLICIT REAL*8(A-H,O-Z)
C
C PARAMETER definitions
C
      INTEGER IXTRMAX
      PARAMETER (IXTRMAX=4)
      DOUBLE PRECISION DANGMAX
      PARAMETER (DANGMAX=3D0)
C
C Dummy arguments
C
      DOUBLE PRECISION ERRAVANG
      INTEGER ITRSCF,NK,NKMMAX,NMVECMAX,NQ,NQMAX,NTMAX
      DOUBLE PRECISION FACT(0:100),MVGAM(NTMAX,NMVECMAX),
     &       MVPHI(NTMAX,NMVECMAX),MVTET(NTMAX,NMVECMAX),
     &       QMGAM(NQMAX),QMPHI(NQMAX),QMTET(NQMAX)
      INTEGER ITOQ(NTMAX,NQMAX)
C
C Local variables
C
      DOUBLE PRECISION A,B,C,PHIXTR,TETXTR
      DOUBLE PRECISION DELPHI,DELTET,LASTERR,MIXING,QMGAMMIX,
     &       QMPHIMIX,QMTETMIX,WN,WO
      INTEGER I,IMV,IPREV,IPRINT,IQ,IT,ITAB,IXTR
      DOUBLE PRECISION QMGAMTAB(NQMAX,3),QMPHITAB(NQMAX,3),
     &                 QMTETTAB(NQMAX,3)
      DOUBLE COMPLEX DROTQ(NKMMAX,NKMMAX,NQMAX)
C
      MIXING = 1.0D0
      WO = 1.0D0 - MIXING
      WN = MIXING
      IPRINT = 0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ITERMDIR I/O
C
      QMTET    = t_params%QMTET
      QMPHI    = t_params%QMPHI
      QMPHITAB = t_params%QMPHITAB
      QMTETTAB = t_params%QMTETTAB
      QMGAMTAB = t_params%QMGAMTAB
      ITAB     = t_params%ITAB
      LASTERR  = t_params%LASTERR
!       READ (67) QMTET,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,ITAB,LASTERR
!       REWIND (67)
C
C      ITERMDIR I/O
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      WRITE (1337,FMT=99002)
C
      IF ( ITRSCF.EQ.0 ) THEN
C
C ----------------------------------------------------- store old angles
C --------------------------------------- dummy for tbkkr, done in main0
         DO IQ = 1,NQ
            QMPHITAB(IQ,1) = QMPHI(IQ)
            QMTETTAB(IQ,1) = QMTET(IQ)
            QMGAMTAB(IQ,1) = QMGAM(IQ)
            DO I = 2,3
               QMPHITAB(IQ,I) = 0D0
               QMTETTAB(IQ,I) = 0D0
               QMGAMTAB(IQ,I) = 0D0
            END DO
         END DO
C
         RETURN
C
      ELSE
C
C------------------------- set new local frame of reference according to
C-------------------------------------------- orientation of spin moment
C
         IMV = 1
C
         DO IQ = 1,NQ
C
            IT = ITOQ(1,IQ)
C
            QMPHI(IQ) = MVPHI(IT,IMV)
            QMTET(IQ) = MVTET(IT,IMV)
            QMGAM(IQ) = MVGAM(IT,IMV)
C
         END DO
C
C=======================================================================
C
         IF ( ITRSCF.LE.2 ) THEN
            ITAB = ITRSCF
            DO IQ = 1,NQ
               QMPHI(IQ) = WO*QMPHITAB(IQ,ITAB) + WN*QMPHI(IQ)
               QMTET(IQ) = WO*QMTETTAB(IQ,ITAB) + WN*QMTET(IQ)
               QMGAM(IQ) = WO*QMGAMTAB(IQ,ITAB) + WN*QMGAM(IQ)
            END DO
C
            ITAB = ITRSCF + 1
            DO IQ = 1,NQ
               
               QMPHITAB(IQ,ITAB) = QMPHI(IQ)
               QMTETTAB(IQ,ITAB) = QMTET(IQ)
               QMGAMTAB(IQ,ITAB) = QMGAM(IQ)
            END DO
C
            IPREV = ITAB - 1
C
         ELSE
            IF ( LASTERR.GT.2D0 ) THEN
               IXTR = 2
            ELSE
               IXTR = IXTRMAX
            END IF
C
            DO IQ = 1,NQ
C
               QMPHIMIX = WO*QMPHITAB(IQ,3) + WN*QMPHI(IQ)
               QMTETMIX = WO*QMTETTAB(IQ,3) + WN*QMTET(IQ)
               QMGAMMIX = WO*QMGAMTAB(IQ,3) + WN*QMGAM(IQ)
C
               A = QMPHITAB(IQ,2)
               B = (QMPHITAB(IQ,3)-QMPHITAB(IQ,1))*0.5D0
               C = (QMPHITAB(IQ,3)+QMPHITAB(IQ,1)-2.D0*A)*0.5D0
               PHIXTR = A + B*DBLE(IXTR) + C*DBLE(IXTR)**2
               IF ( PHIXTR.GE.0.0D0 .AND. PHIXTR.LE.360.0D0 .AND. 
     &              ABS(PHIXTR-QMPHITAB(IQ,3)).LT.DANGMAX ) THEN
                  QMPHI(IQ) = PHIXTR
               ELSE
                  QMPHI(IQ) = QMPHIMIX
               END IF
C
               IF ( QMPHI(IQ).LT.0D0 ) QMPHI(IQ) = QMPHI(IQ) + 360D0
               IF ( QMPHI(IQ).GT.360D0 ) QMPHI(IQ) = QMPHI(IQ) - 360D0
C
               A = QMTETTAB(IQ,2)
               B = (QMTETTAB(IQ,3)-QMTETTAB(IQ,1))*0.5D0
               C = (QMTETTAB(IQ,3)+QMTETTAB(IQ,1)-2.D0*A)*0.5D0
               TETXTR = A + B*DBLE(IXTR) + C*DBLE(IXTR)**2
               D12=QMTETTAB(IQ,1)-QMTETTAB(IQ,2)
               D23=QMTETTAB(IQ,2)-QMTETTAB(IQ,3)
               D3X=QMTETTAB(IQ,3)-TETXTR
               IF ( TETXTR.GE.0.0D0 .AND. TETXTR.LE.180.0D0 .AND. 
     &              ABS(TETXTR-QMTETTAB(IQ,3)).LT.DANGMAX .AND.
     &               (D12*D23.GT.0D0) .AND. (D23*D3X.GT.0D0) ) THEN  
                  QMTET(IQ) = TETXTR
               ELSE
                  QMTET(IQ) = QMTETMIX
               END IF
C
               A = QMGAMTAB(IQ,2)
               B = (QMGAMTAB(IQ,3)-QMGAMTAB(IQ,1))*0.5D0
               C = (QMGAMTAB(IQ,3)+QMGAMTAB(IQ,1)-2.D0*A)*0.5D0
               QMGAM(IQ) = A + B*DBLE(IXTR) + C*DBLE(IXTR)**2
C
               IF ( IPRINT.GT.0 ) WRITE (1337,99001) IQ,IXTR,LASTERR,
     &              ('PHI',I,QMPHITAB(IQ,I),I=1,3),'PHIMIX',QMPHIMIX,
     &              'PHIXTR',PHIXTR,'PHINEW',QMPHI(IQ),
     &              ('TET',I,QMTETTAB(IQ,I),I=1,3),'TETMIX',QMTETMIX,
     &              'TETXTR',TETXTR,'TETNEW',QMTET(IQ)
C
               DO I = 1,2
                  QMPHITAB(IQ,I) = QMPHITAB(IQ,I+1)
                  QMTETTAB(IQ,I) = QMTETTAB(IQ,I+1)
                  QMGAMTAB(IQ,I) = QMGAMTAB(IQ,I+1)
               END DO
               QMPHITAB(IQ,3) = QMPHI(IQ)
               QMTETTAB(IQ,3) = QMTET(IQ)
               QMGAMTAB(IQ,3) = QMGAM(IQ)
C
            END DO
C
            IPREV = 2
         END IF
C
         WRITE (1337,99004)
C
         ERRAVANG = 0D0
         DO IQ = 1,NQ
C
            DELPHI = ABS(QMPHI(IQ)-QMPHITAB(IQ,IPREV))
            DELTET = ABS(QMTET(IQ)-QMTETTAB(IQ,IPREV))
            ERRAVANG = MAX(ERRAVANG,DELPHI,DELTET)
C
            WRITE (1337,99005) IQ,QMPHITAB(IQ,IPREV),QMTETTAB(IQ,IPREV),
     &                      QMPHI(IQ),QMTET(IQ),DELPHI,DELTET
C
C --> update the rotation matrices DROTQ for the new angles
C
            CALL CALCROTMAT(NK,3,QMPHI(IQ),QMTET(IQ),0.0D0,
     &           DROTQ(1,1,IQ),FACT,NKMMAX)
C
         END DO
C
         WRITE (1337,FMT=99003) ITRSCF,ERRAVANG
         WRITE (1337,'(I5,4F10.3,'' #  ANGLES'',/,79(1H+),/)') ITRSCF,
     &          (QMPHI(IQ),QMTET(IQ),IQ=1,MIN(2,NQ))
C
         LASTERR = ERRAVANG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ITERMDIR I/O
C
!          WRITE (67) QMTET,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,ITAB,LASTERR
!          WRITE (67) DROTQ
         t_params%QMTET    = QMTET
         t_params%QMPHI    = QMPHI
         t_params%QMPHITAB = QMPHITAB
         t_params%QMTETTAB = QMTETTAB
         t_params%QMGAMTAB = QMGAMTAB
         t_params%ITAB     = ITAB
         t_params%LASTERR  = LASTERR
         t_params%DROTQ    = DROTQ
C
C      ITERMDIR I/O
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END IF
99001 FORMAT (/,5X,'IQ',I3,'  IXTR',I3,'  LAST ERROR',F12.4,/,
     &        3(3(5X,A,5x,'TAB',I2,2X,F12.6,/),3(5X,A,9X,F12.6,/)))
C   ====================================================================
99002 FORMAT (79('+'),/,28X,'ANGLE - mixing scheme')
99003 FORMAT (5X,'iter.',i4,'     max. CHANGE = ',F12.6)
99004 FORMAT (/,5X,' setting new LOCAL frame of reference ',//,5X,
     &        'IQ  old   phi      tet    --> new   phi      tet',
     &        '     del   phi      tet')
99005 FORMAT (5X,I2,4X,2F9.4,8X,2F9.4,5X,2F9.4)
      END
