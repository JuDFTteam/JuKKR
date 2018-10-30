C*==rrgen.f    processed by SPAG 6.05Rc at 20:37 on 17 May 2004
C 02.08.95 *************************************************************
      SUBROUTINE RRGEN (BV1,LSURF,RR,NR,NRD)
C **********************************************************************
C *                                                                    *
C * generates a number of real space vectors to construct the          *
C * clusters representing the local surrounding of the atoms in        *
C * routine CLSGEN99                                                   *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar arguments ..
      LOGICAL LSURF
      INTEGER NR,NRD
C    ..
C    .. Array arguments ..
      DOUBLE PRECISION BV1(3,3),RR(3,0:NRD)
C    ..
C    .. Local scalars ..
      DOUBLE PRECISION EPSSHL,R,R1,R2,R3,RMAX,RR2,RS
      INTEGER I,J,K,N1,N2,N3,POS,IPRINT
      INTEGER NINT
      DOUBLE PRECISION DBLE
C     ..
C     .. Local arrays
      DOUBLE PRECISION RABS(NRD),RR1(3,NRD),
     +                 V(3),VX(3),VY(3),VZ(3),
     +                 VX0(3),VY0(3),VZ0(3)
      INTEGER IND(NRD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN,SQRT,NINT
C     ..
C     .. External Subroutines ..
      EXTERNAL DSORT,SCALPR,VADD,VEQ
C     ..
C     .. Data Statements ..
      DATA  EPSSHL /1.0D-5/
C     ..................................................................
      WRITE (1337,'(5X,A,/)') 
     &                '< RRGEN > : generation of real space mesh RR(NR)'
C
      IPRINT = 0
C
      CALL SCALPR(BV1(1,1),BV1(1,1),R1)
      CALL SCALPR(BV1(1,2),BV1(1,2),R2)
      CALL SCALPR(BV1(1,3),BV1(1,3),R3)
      RMAX = 5.D0
C
      R1 = SQRT(R1)
      R2 = SQRT(R2)
      R3 = SQRT(R3)
      R = 1.5D0*RMAX + SQRT(R1*R1+R2*R2+R3*R3) + EPSSHL
      RS = R*R
      N1 = NINT(R/R1)
      N2 = NINT(R/R2)
      IF ( .NOT.LSURF ) N3 = NINT(R/R3)
C
      N1 = MIN(12,N1)
      N2 = MIN(12,N2)
      IF ( .NOT.LSURF ) N3 = MIN(12,N3)
C
      N1 = MAX(2,N1)
      N2 = MAX(2,N2)
      IF ( .NOT.LSURF ) N3 = MAX(2,N3)
C
      IF ( LSURF ) N3 = 0
C
      WRITE (1337,99001) R
      WRITE (1337,99002) RS
      IF ( LSURF ) THEN
         WRITE (1337,99003) N1,N2
      ELSE
         WRITE (1337,99004) N1,N2,N3
      END IF
C
      NR = 0
      RR(1,0) = 0.0D0
      RR(2,0) = 0.0D0
      RR(3,0) = 0.0D0
C
      CALL VMUL(BV1(1,1),DBLE(-N1-1),VX0(1))
      CALL VMUL(BV1(1,2),DBLE(-N2-1),VY0(1))
      CALL VMUL(BV1(1,3),DBLE(-N3-1),VZ0(1))
      CALL VEQ(VX0,VX)
C **********************************************************************
      DO I = -N1,N1
         CALL VADD(VX,BV1(1,1),VX)
         CALL VEQ(VY0,VY)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         DO J = -N2,N2
            CALL VADD(VY,BV1(1,2),VY)
            CALL VEQ(VZ0,VZ)
C ----------------------------------------------------------------------
            DO K = -N3,N3
               CALL VADD(VZ,BV1(1,3),VZ)
               CALL VADD(VX,VY,V)
               CALL VADD(V,VZ,V)
               CALL SCALPR(V,V,RR2)
C
               IF ( ((RR2.LE.RS) .OR. (ABS(I)+ABS(J)+ABS(K).LE.6))
     &               .AND. (RR2.GT.EPSSHL) ) THEN
                  NR = NR + 1
C
                  IF ( NR.GT.NRD ) THEN
                     WRITE (6,*) 'Dimension ERROR. Please, change the ',
     &                           'parameter NRD in inc.p to ',NR, NRD
                     STOP
                  END IF
C
                  RR1(1,NR) = V(1)
                  RR1(2,NR) = V(2)
                  RR1(3,NR) = V(3)
                  RABS(NR) = SQRT(RR2)
               END IF
            END DO
C ----------------------------------------------------------------------
         END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END DO
C **********************************************************************
C
      WRITE (1337,99005) NR+1
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      IF ( IPRINT.GT.0 ) THEN
         WRITE (1337,99006)
         WRITE (1337,99008) 0,0.0,0.0,0.0,0.0
      END IF
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      CALL DSORT(RABS,IND,NR,POS)
      DO I = 1,NR
         POS = IND(I)
         RR(1,I) = RR1(1,POS)
         RR(2,I) = RR1(2,POS)
         RR(3,I) = RR1(3,POS)
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      IF ( IPRINT.GT.0 ) WRITE (1337,99008) I,RR(1,I),RR(2,I),
     &                                  RR(3,I),RABS(POS)
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      END DO
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      IF ( IPRINT.GT.0 )  WRITE (1337,99007)
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
99001 FORMAT (10X,'Radius R        : ',F15.6,' (ALAT    units)')
99002 FORMAT (10X,'       R**2     : ',F15.6,' (ALAT**2 units)')
99003 FORMAT (10X,'mesh divisions  : ',5X,2I5)
99004 FORMAT (10X,'mesh divisions  : ',3I5)
99005 FORMAT (10X,'vectors created : ',I15)
99006 FORMAT (/,10X,60('+'),/,18X,
     &        'generated real-space mesh-points (ALAT units)',/,
     &        10X,60('+'),/,13X,
     &        'index      x           y           z          distance  '
     &        ,/,10X,60('-'))
99007 FORMAT (10X,60('+'))
99008 FORMAT (10X,I6,3F12.3,F15.4)
      END
