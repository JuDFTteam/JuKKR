C*==gijcond.f    processed by SPAG 6.05Rc at 15:45 on 18 Oct 2004
      SUBROUTINE GIJCOND(IDO,NAEZ,RBASIS,IQAT,NATOMIMP,RCLSIMP,ATOMIMP,
     &                   IJTABCALC,NAEZD,NATOMIMPD)
C **********************************************************************
C *                                                                    *
C * In case of tasks requiring Gij blocks calculation, set variables:  *
C *                                                                    *
C * NATOMIMP, RCLSIMP(3,1..NATOMIMP), ATOMIMP(1..NATOMIMP)             *
C * IJTABCALC flag to which pair is needed: I,J --> (I-1)*NATOMIMP + J *
C * IDO takes on the value 1 or 0 if setting up process was OK or not  *
C *                                                                    *
C * CONDUCTANCE calculation case                                       *
C *             still to be implemente the correct read in             *
C **********************************************************************
      IMPLICIT NONE
C ..
C ..  Parameters
      INTEGER NCPAIRD
      PARAMETER (NCPAIRD=10)
C ..
C ..  Arguments
      INTEGER IDO,NAEZ,NAEZD,NATOMIMP,NATOMIMPD
      INTEGER ATOMIMP(*),IJTABCALC(*),IQAT(*)
      DOUBLE PRECISION RBASIS(3,*),RCLSIMP(3,*)
C ..
C ..  Locals
      INTEGER I,IAT,IATCONDL(NCPAIRD),IATCONDR(NCPAIRD),J,JAT,NCONDPAIR,
     &        NN
C ..
      IDO = 0
      DO I = 1,NCPAIRD
         IATCONDL(I) = 0
         IATCONDR(I) = 0
      END DO
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,99001)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C     ---------------------------------------------------- dummy
C     settings so far, need to be replaced by conductance input
C     and some output
      NCONDPAIR = 4
      IF ( NCONDPAIR.GT.NCPAIRD ) THEN
         WRITE (6,99002) 'local','NCPAIRD',NCONDPAIR
         STOP
      END IF
      IATCONDL(1) = 1
      IATCONDR(1) = 2
      IATCONDL(2) = 1
      IATCONDR(2) = 2
      IATCONDL(3) = 2
      IATCONDR(3) = 1
      IATCONDL(4) = 2
      IATCONDR(4) = 1
C
C     ---------------------------------------------------- dummy
      IF ( NCONDPAIR.EQ.0 ) RETURN
      DO I = 1,NCONDPAIR
         IF ( (IATCONDL(I).LE.0) .OR. (IATCONDL(I).GT.NAEZ) ) RETURN
         IF ( (IATCONDR(I).LE.0) .OR. (IATCONDR(I).GT.NAEZ) ) RETURN
      END DO
C
      NATOMIMP = 2*NCONDPAIR
      IF ( NATOMIMP.GT.NATOMIMPD ) THEN
         WRITE (6,99002) 'global','NATOMIMPD',NATOMIMP
         STOP
      END IF
C
      DO I = 1,NATOMIMP
         NN = (I-1)*NATOMIMP
         DO J = 1,NATOMIMP
            IJTABCALC(NN+J) = 0
         END DO
      END DO
C
      NN = 0
      DO I = 1,NCONDPAIR
         IAT = IQAT(IATCONDL(I)) ! left lead
         NN = NN + 1
         DO J = 1,3
            RCLSIMP(J,NN) = RBASIS(J,IAT)
         END DO
         ATOMIMP(NN) = IAT
         IAT = NN
C
         JAT = IQAT(IATCONDR(I)) ! right lead
         NN = NN + 1
         DO J = 1,3
            RCLSIMP(J,NN) = RBASIS(J,JAT)
         END DO
         ATOMIMP(NN) = JAT
         JAT = NN
         IJTABCALC((IAT-1)*NATOMIMP+JAT) = 1
      END DO
      IF ( NATOMIMP.NE.NN ) THEN
         WRITE (6,'(6X,A,/,6X,A,/)') 
     &             'ERROR: Found some inconsistencies in IATCOND arrays'
     &             ,'       Please check your CONDUCTANCE input'
         STOP
      END IF
      IDO = 1
99001 FORMAT (5X,'< GIJCOND > : Conductance/conductivity calculation',/)
99002 FORMAT (6X,'Dimension ERROR: please increase the ',A,' parameter',
     &        /,6X,A,' to a value >=',I5,/)
      END
