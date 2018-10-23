C*==madelgaunt.f    processed by SPAG 6.05Rc at 11:19 on 17 May 2004
      SUBROUTINE MADELGAUNT(LPOT,YRG,WG,CLEB,ICLEB,IEND,LASSLD,NCLEBD)
      IMPLICIT NONE
C     ..
C     .. Scalar arguments
      INTEGER LPOT,IEND
      INTEGER LASSLD,NCLEBD
C     ..
C     .. Array arguments
C     .. Attention: Dimension NCLEBD appears sometimes as NCLEB1
C     ..            an empirical factor - it has to be optimized
      DOUBLE PRECISION YRG(LASSLD,0:LASSLD,0:LASSLD),WG(LASSLD)
      DOUBLE PRECISION CLEB(NCLEBD)
      INTEGER ICLEB(NCLEBD,3)
C     ..
C     .. Local scalars
      DOUBLE PRECISION CLECG,FACTOR,S
      INTEGER I,J,L1,L2,L3,M1,M1A,M1S,M2,M2A,M2S,M3,M3A,M3S
C     ..
C     .. Intrinsic functions
      INTRINSIC ABS,ATAN,DBLE,SIGN
C
C --> set up of the gaunt coefficients with an index field
C     recognize that they are needed here only for l3=l1+l2
C
      IF ( 2*LPOT.GT.LASSLD ) THEN
         WRITE (6,*) 'Dim ERROR in MADELGAUNT -- 2*LPOT > LASSLD',
     &        2*LPOT,LASSLD
         STOP
      END IF
C
      I = 1
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      DO L1 = 0,LPOT
         DO L2 = 0,LPOT
            L3 = L1 + L2
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO M1 = -L1,L1
               DO M2 = -L2,L2
                  DO M3 = -L3,L3
                     M1S = SIGN(1,M1)
                     M2S = SIGN(1,M2)
                     M3S = SIGN(1,M3)
C **********************************************************************
                     IF ( M1S*M2S*M3S.GE.0 ) THEN
                        M1A = ABS(M1)
                        M2A = ABS(M2)
                        M3A = ABS(M3)
C
                        FACTOR = 0.0D0
                        IF ( M1A+M2A.EQ.M3A ) FACTOR = FACTOR + 
     &                       DBLE(3*M3S+SIGN(1,-M3))/8.0D0
                        IF ( M1A-M2A.EQ.M3A ) FACTOR = FACTOR + 
     &                       DBLE(M1S)/4.0D0
                        IF ( M2A-M1A.EQ.M3A ) FACTOR = FACTOR + 
     &                       DBLE(M2S)/4.0D0
C ======================================================================
                        IF ( FACTOR.NE.0.0D0 ) THEN
                           IF ( M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR. 
     &                          M1S*M3S.NE.1 ) FACTOR = -FACTOR
C
                           S = 0.0D0
                           DO J = 1,LASSLD
                              S = S + WG(J)*YRG(J,L1,M1A)*YRG(J,L2,M2A)
     &                            *YRG(J,L3,M3A)
                           END DO
C
                           CLECG = S*FACTOR
C ----------------------------------------------------------------------
                           IF ( ABS(CLECG).GT.1.D-10 ) THEN
                              CLEB(I) = CLECG
                              ICLEB(I,1) = L1*(L1+1) + M1 + 1
                              ICLEB(I,2) = L2*(L2+1) + M2 + 1
                              ICLEB(I,3) = L3*(L3+1) + M3 + 1
                              I = I + 1
                              IF ( I.GT.NCLEBD ) THEN
                                 WRITE (6,FMT='(2I10)') I,NCLEBD
                                 STOP ' Dim stop in MADELGAUNT '
                              END IF
                           END IF
C ----------------------------------------------------------------------
                        END IF
C ======================================================================
                     END IF
C **********************************************************************
                  END DO
               END DO
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
         END DO
      END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      IEND = I - 1
      END
