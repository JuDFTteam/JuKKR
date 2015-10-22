C*==taustruct.f    processed by SPAG 6.05Rc at 15:50 on 10 Dec 2002
      SUBROUTINE TAUSTRUCT(DROT,NSYM,SYMUNITARY,NKM,NQ,NQMAX,NKMMAX,
     &                     IPRINT,IREL)
C   ********************************************************************
C   *                                                                  *
C   *   find the structure of the site-diagonal TAU - matrices  TAUQ   *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C PARAMETER definitions
C
      COMPLEX*16 C0
      PARAMETER (C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER IPRINT,IREL,NKM,NKMMAX,NQ,NQMAX,NSYM
      COMPLEX*16 DROT(NKMMAX,NKMMAX,48)
      LOGICAL SYMUNITARY(48)
C
C Local variables
C
      REAL*8 ABST,X
      DOUBLE PRECISION DBLE
      INTEGER I,I0,IMWEIGHT,IQ,ISYM,IW,IWR,J,K,L,LIN,NELMT,NKMQ(NQMAX),
     &        NKMTOP,NLIN,NON0(NQMAX)
      COMPLEX*16 ST(NKMMAX,NKMMAX),TAUK(NKMMAX,NKMMAX,NQMAX)
C
      DO IQ = 1,NQ
         NKMQ(IQ) = NKM
      END DO
C
      IMWEIGHT = 0
      NELMT = 0
      NLIN = 0
      IW = 6
C
      DO IQ = 1,NQ
         NON0(IQ) = 0
         NKMTOP = NKMQ(IQ)
C
         IF ( IPRINT.GT.0 ) WRITE (1337,99004) IQ
         DO I = 1,NKMTOP
            DO J = 1,NKMTOP
               ST(I,J) = 0.0D0
C
               CALL CINIT(NKMMAX*NKMMAX*NQMAX,TAUK)
C
               DO ISYM = 1,NSYM
                  I0 = IQ
C
                  IF ( SYMUNITARY(ISYM) ) THEN
                     DO L = 1,NKMTOP
                        DO K = 1,NKMTOP
                           TAUK(K,L,I0) = TAUK(K,L,I0) + DROT(I,K,ISYM)
     &                        *DCONJG(DROT(J,L,ISYM))
                        END DO
                     END DO
                  ELSE
                     DO L = 1,NKMTOP
                        DO K = 1,NKMTOP
                           TAUK(L,K,I0) = TAUK(L,K,I0) + DROT(I,K,ISYM)
     &                        *DCONJG(DROT(J,L,ISYM))
                        END DO
                     END DO
                  END IF
               END DO
C
               LIN = 0
               IWR = 0
               DO K = 1,NKMQ(IQ)
                  DO L = 1,NKMQ(IQ)
                     ABST = ABS(TAUK(K,L,IQ))
                     ST(I,J) = ST(I,J) + ABST
                     IF ( ABST.GT.1D-8 ) THEN
                        IF ( DIMAG(TAUK(K,L,IQ)).GT.1D-5 ) THEN
                           IF ( IPRINT.GT.0 ) WRITE (1337,*)
     &                           ' Im(Weight) > 1D-5 ',I,J,K,L
                           IMWEIGHT = 1
                        END IF
                        X = DREAL(TAUK(K,L,IQ))/DBLE(NSYM)
C
                        IF ( IPRINT.GT.1 ) THEN
                           IF ( IWR.EQ.0 ) THEN
                              IWR = 1
                              WRITE (IW,99002) I,J,IQ,X,K + (IQ-1)*NKM,
     &                               L + (IQ-1)*NKM
                           ELSE
                              WRITE (IW,99003) X,K + (IQ-1)*NKM,
     &                               L + (IQ-1)*NKM
                           END IF
                        END IF
                        LIN = LIN + 1
                     END IF
C
                  END DO
               END DO
C
               IF ( LIN.GT.0 ) THEN
                  NLIN = NLIN + LIN
                  NELMT = NELMT + 1
                  NON0(IQ) = NON0(IQ) + 1
               END IF
C
               IF ( ABS(ST(I,J)).GT.1D-5 ) ST(I,J) = 2
C
            END DO
         END DO
C
         IF ( IPRINT.GT.1 ) CALL CMATSTR('TAU-MAT',7,ST,NKMTOP,NKMMAX,
     &        IREL,IREL,0,1D-8,6)
      END DO
C
      WRITE (1337,99005) NELMT,(NON0(IQ),IQ=1,NQ)
      WRITE (1337,99006) NLIN
C
      IF ( IMWEIGHT.NE.0 ) WRITE (1337,99007)
C
C-----------------------------------------------------------------------
C
      RETURN
99002 FORMAT ('     TAUQ(',I2,',',I2,',',I2,') =  ','   ',F10.4,' * <',
     &        I3,'|T(K)|',I3,'>')
99003 FORMAT (23X,' + ',F10.4,' * <',I3,'|T(K)|',I3,'>')
99004 FORMAT (//,
     &    ' ==========================================================='
     &    ,/,'   structure of  TAU-matrix   INT <i|t(k)|j>     IQ=',I3,
     &    /,
     &    ' ==========================================================='
     &    ,/)
99005 FORMAT (/,5X,'non-0 TAU-elements          ',I5,'   Q:',80I4)
99006 FORMAT (5X,'terms to sum up             ',I5,/)
99007 FORMAT (/,5X,50('#'),/,5X,'WARNING: complex TAU weights found',/,
     &        5X,'this may occur for rotated magnetic moments',/,5X,
     &        'relevant only for tetrahedron BZ-integration',/,5X,
     &        50('#'),/)
      END
