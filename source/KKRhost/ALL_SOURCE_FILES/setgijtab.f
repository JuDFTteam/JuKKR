c      SC*==setgijtab.f    processed by SPAG 6.05Rc at 15:49 on 18 Oct 2004
      SUBROUTINE SETGIJTAB(LINTERFACE,ICC,NAEZ,IQAT,RBASIS,BRAVAIS,
     &                     NATOMIMP,ATOMIMP,RCLSIMP,NOFGIJ,IJTABCALC,
     &                     IOFGIJ,JOFGIJ,NQCALC,IQCALC,NAEZD,NATOMIMPD,
     &                     IJTABCALC_I)
C **********************************************************************
C * Task-specific settings of Gij elements that need to be calculated  *
C * Subroutine (called for ICC=-1) sets up the arrays                  *
C * NATOMIMP    : number of different sites i,j = 1,NATOMIMP           *
C * RCLSIMP     : site coordinates                                     *
C * ATOMIMP     : index of the corresponding site in the unit cell     *
C * IJTABCALC   : flag specifying wehter pair (I,J) needs to be        *
C *               calculated - linear pointer (I-1)*NATOMIMP + J = 1/0 *
C *               for YES/NO                                           *
C * NOFGIJ      : number of all (I,J) pairs - sum of all non-zero I,J  *
C * IOFGIJ      : I index in the list 1..NATOMIMP for pair I,J         *
C * JOFGIJ      : J index                                              *
C **********************************************************************
      IMPLICIT NONE
C ..  
C ..  Scalar arguments
      INTEGER ICC,NAEZ,NAEZD,NATOMIMP,NATOMIMPD,NOFGIJ,NQCALC
      LOGICAL LINTERFACE
C ..   
C ..  Array arguments
      INTEGER ATOMIMP(*),IJTABCALC(*),IJTABCALC_I(*),IOFGIJ(*),IQAT(*),
     &        IQCALC(*),JOFGIJ(*)
      DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,*),RCLSIMP(3,*)
C .. 
C ..  Local scalars
      INTEGER I,IDO,II,J,JJ,NN
      LOGICAL OPT
C ..  
C ..  External subroutines
      EXTERNAL GIJCOND,GIJXCPL,OPT
C     ..
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,'(79(1H=),/,15X,A)') 
     &                      'SETGIJTAB: setting task-specific Gij pairs'
      WRITE (1337,'(79(1H=),/)')
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      IDO = 0
C ======================================================================
      IF ( OPT('CONDUCT ') ) CALL GIJCOND(IDO,NAEZ,RBASIS,IQAT,NATOMIMP,
     &     RCLSIMP,ATOMIMP,IJTABCALC,NAEZD,NATOMIMPD)
C ======================================================================
      IF ( OPT('XCPL    ') ) CALL GIJXCPL(IDO,NAEZ,RBASIS,BRAVAIS,
     &     LINTERFACE,NQCALC,IQCALC,NATOMIMP,RCLSIMP,ATOMIMP,IJTABCALC,
     &     IJTABCALC_I,NATOMIMPD)
C ======================================================================
      IF ( IDO.EQ.0 ) THEN
         ICC = 0
         WRITE (6,99002)
         RETURN
      END IF
C ======================================================================
      NOFGIJ = 0
      DO I = 1,NATOMIMP
         NN = (I-1)*NATOMIMP
         DO J = 1,NATOMIMP
            IF ( IJTABCALC(NN+J).GT.0 ) THEN
               NOFGIJ = NOFGIJ + 1
               IF ( NOFGIJ.GT.NATOMIMPD*NATOMIMPD ) THEN
                  WRITE (6,99001) 'NATOMIMPD',NOFGIJ/NATOMIMP
                  STOP
               END IF
               IOFGIJ(NOFGIJ) = I
               JOFGIJ(NOFGIJ) = J
            END IF
         END DO
      END DO
      IF ( NOFGIJ.EQ.0 ) THEN
         ICC = 0
         WRITE (6,99002)
         RETURN
      END IF
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,99003) NATOMIMP,NOFGIJ
      WRITE (1337,99004)
      WRITE (1337,99005)
      WRITE (1337,99004)
      DO I = 1,NOFGIJ
         II = IOFGIJ(I)
         JJ = JOFGIJ(I)
         WRITE (1337,99006) I,II,ATOMIMP(II),(RCLSIMP(J,II),J=1,3),JJ,
     &                   ATOMIMP(JJ),(RCLSIMP(J,JJ),J=1,3)
      END DO
      WRITE (1337,99004)
      WRITE (1337,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
99001 FORMAT (6X,'brahim ERROR: please increase the global parameter'
     &        ,/,6X,A,' to a value >=',I5,/)
99002 FORMAT (6X,'WARNING: Subroutine entered with invalid task ',
     &        'specification',/,6X,
     &        '         ICC will be set to 0 - no Gij calculated - ',
     &        'input check? ',/)
99003 FORMAT (6X,'Number of different sites (NATOMIMP) :',I4,/,6X,
     &        'Number of pairs set       (NOFGIJ)   :',I4)
99004 FORMAT (8X,71('-'))
99005 FORMAT (9X,'pair|',' I  IQ           position',9X,
     &        'J  JQ           position')
99006 FORMAT (9X,I3,' |',2(I3,1X),3F8.4,1X,2(I3,1X),3F8.4)
99007 FORMAT (I5,2(I5,1X),3F10.6,1X,2(I5,1X),3F10.6)
      END
