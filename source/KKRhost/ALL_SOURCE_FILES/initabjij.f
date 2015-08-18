C*==initabjij.f    processed by SPAG 6.05Rc at 15:51 on 18 Oct 2004
      SUBROUTINE INITABJIJ(IPRINT,NAEZ,NATYP,NATOMIMP,NOFGIJ,NQCALC,
     &                     NSMAX,NSHELL,IQCALC,ATOMIMP,ISH,JSH,
     &                     IJTABCALC,IJTABSH,IJTABSYM,NIJCALC,KIJSH,
     &                     NIJMAX,NSHELL0,NSHELD)
C   ********************************************************************
C   *  subroutine called by < TBXCCPLJIJ > to set up some auxiliary    *
C   *  arrays allowing the indexing of shells, sites, atomic types     *
C   ********************************************************************
C
      IMPLICIT NONE
C ..
C ..  Arguments
      INTEGER IPRINT,NAEZ,NATOMIMP,NATYP,NIJMAX,NOFGIJ,NQCALC,NSHELD,
     &        NSHELL0,NSMAX
      INTEGER ATOMIMP(*),IJTABCALC(*),IJTABSH(*),IJTABSYM(*),IQCALC(*),
     &        ISH(NSHELD,*),JSH(NSHELD,*),KIJSH(NIJMAX,NSHELL0),
     &        NIJCALC(NSHELL0),NSHELL(0:NSHELD)
C ..
C ..  Locals
      INTEGER I1,IA,IDONE(NAEZ),IQTOJQ(NIJMAX),J1,JA,LM1,LM2,NS
      INTEGER NIDONE
C     
C ======================================================================
      DO NS = NSMAX + 1,NSHELL(0)
         DO I1 = 1,NIJMAX
            IQTOJQ(I1) = 0
         END DO
C ----------------------------------------------------------------------
         DO I1 = 1,NSHELL(NS)
            IA = ATOMIMP(ISH(NS,I1))
            JA = 0
            DO J1 = 1,NIJCALC(NS)
               IF ( IA.EQ.IQTOJQ(J1) ) THEN
                  JA = 1
                  GOTO 20
               END IF
            END DO
 20         CONTINUE
            IF ( JA.EQ.0 ) THEN
               NIJCALC(NS) = NIJCALC(NS) + 1
               IF ( NIJCALC(NS).GT.NIJMAX ) THEN
                  WRITE (6,99001) 'local','NIJMAX',NIJCALC(NS)
                  STOP '       in < TBXCCPLJIJ > '
               END IF
               IQTOJQ(NIJCALC(NS)) = IA
               KIJSH(NIJCALC(NS),NS) = I1
            END IF
         END DO
      END DO
C ======================================================================
      IF ( IPRINT.LE.0 ) RETURN
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      NIDONE = 0
      DO IA = 1,NAEZ
         IDONE(IA) = 0
         DO I1 = 1,NQCALC
            IF ( IQCALC(I1).EQ.IA ) THEN
               NIDONE = NIDONE + 1
               IDONE(NIDONE) = IA
            END IF
         END DO
      END DO
C
      LM2 = MIN(25,NATOMIMP)
      WRITE (6,99002) NAEZ,NATYP,NATOMIMP,NOFGIJ,NSHELL(0),LM2
      DO I1 = 1,NIDONE
         DO IA = 1,NATOMIMP
            IF (ATOMIMP(IA).EQ.IDONE(I1)) THEN
               LM1 = (IA-1)*NATOMIMP
               WRITE (6,99003) IA,(IJTABCALC(LM1+JA),JA=1,LM2)
               GOTO 100
            END IF
         END DO
 100     CONTINUE
      END DO
      WRITE (6,99004) LM2
      DO I1 = 1,NIDONE
         DO IA = 1,NATOMIMP
            IF (ATOMIMP(IA).EQ.IDONE(I1)) THEN
               LM1 = (IA-1)*NATOMIMP
               WRITE (6,99003) IA,(IJTABSH(LM1+JA),JA=1,LM2)
               GOTO 110
            END IF
         END DO
 110     CONTINUE
      END DO
      WRITE (6,99005) LM2
      DO I1 = 1,NIDONE
         DO IA = 1,NATOMIMP
            IF (ATOMIMP(IA).EQ.IDONE(I1)) THEN
               LM1 = (IA-1)*NATOMIMP
               WRITE (6,99003) IA,(IJTABSYM(LM1+JA),JA=1,LM2)
               GOTO 120
            END IF
         END DO
 120     CONTINUE
      END DO
      LM2 = 0
      DO NS = NSMAX + 1,NSHELL(0)
         LM2 = MAX(LM2,NIJCALC(NS))
      END DO
      LM2 = MIN(5,LM2)
      WRITE (6,99006)
      DO NS = NSMAX + 1,NSHELL(0)
         WRITE (6,99007) NS,(ISH(NS,KIJSH(I1,NS)),JSH(NS,KIJSH(I1,NS)),
     &                   I1=1,MIN(NIJCALC(NS),LM2))
         WRITE (6,99008) (ATOMIMP(ISH(NS,KIJSH(I1,NS))),ATOMIMP(JSH(NS,
     &                   KIJSH(I1,NS))),
     &                   IJTABSYM((ISH(NS,KIJSH(I1,NS))-1)
     &                   *NATOMIMP+JSH(NS,KIJSH(I1,NS))),I1=1,
     &                   MIN(NIJCALC(NS),LM2))
      END DO
Cccc      WRITE (6,99009) MIN(NATYP,25)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
99001 FORMAT (6X,'Dimension ERROR: please increase the ',A,' parameter',
     &        /,6X,A,' to a value >=',I5,/)
99002 FORMAT (8X,60('-'),/8X,'Data used for J_ij calculation:',//,10X,
     &        'Number of sites/types i        (NAEZ/NATYP) :',2(1X,I3),
     &        /,10X,'Number of atoms in the cluster   (NATOMIMP) :',1X,
     &        I3,/,10X,'Number of ij pairs                 (NOFGIJ) :',
     &        1X,I3,/,10X,
     &        'Number of representative pairs     (NSHELL) :',1X,I3,//,
     &        10X,'ij-pairs calculation table ( 1 = calculated )',/,10X,
     &        'IA   JA = 1 ..',I3)
99003 FORMAT (10X,I3,3X,25(I3))
99004 FORMAT (/,10X,'ij-shells table ',/,10X,'IA   JA = 1 ..',I3)
99005 FORMAT (/,10X,'ij-symmetries table ',/,10X,'IA   JA = 1 ..',I3)
99006 FORMAT (/,10X,'effectively calculated pairs/shells',/,10X,
     &        'SHELL   (IAT,JAT) ',/,10X,
     &        'SHELL   (IQ,JQ - ISYM) ')
99007 FORMAT (10X,I4,3X,5(I3,',',I3,5X))
99008 FORMAT (10X,4X,3X,5(I3,',',I3,' - ',I2))
99009 FORMAT (/,10X,'effectively calculated type-type pairs (shells)',
     &        /,10X,'IT   JT = 1 ..',I3)
      END
