C*==lattix99.f    processed by SPAG 6.05Rc at 17:56 on 17 May 2004
      SUBROUTINE LATTIX99(ALAT,BRAVAIS,RECBV,VOLUME0,RR,NR,NRD)
C **********************************************************************
C *                                                                    *
C * LATTIX99 generates the real space and reciprocal lattices.         *
C * BRAVAIS(I,J) are basis vectors, with I=X,Y,Z and J=A,B,C           *
C * RECIPROCAL space vectors are in UNITS OF 2*PI/ALATC                *
C * RR are the direct space vectors                                    *
C * NR+1 is the number of direct space vectors created                 *
C * (structure dependent output).                                      *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar arguments ..
      INTEGER NR,NRD            ! number of real space vectors
      INTEGER NAEZD
      DOUBLE PRECISION ALAT,VOLUME0
C     ..
C     .. Array arguments ..
C
C  BRAVAIS(3,3): Real space bravais vectors. Read in normalised to ALAT
C  RECBV(3,3)  : Reciprocal lattice vectors in 2*PI/ALAT
C
      DOUBLE PRECISION BRAVAIS(3,3)
      DOUBLE PRECISION RECBV(3,3),RR(3,0:NRD)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,IER
      DOUBLE PRECISION VOLUC,DET,DDET33,PI,TPIA
      CHARACTER*80 UIO
C     ..
C     .. External declarations ..
      EXTERNAL CROSPR,SPATPR,DDET33,IOINPUT
C     ..
C     .. Intrinsic functions ..
      INTRINSIC ABS,ATAN,DBLE
C     ..
C     ..................................................................
C
C --> initialise
C
      PI = 4D0*ATAN(1D0)
      TPIA = 2D0*PI/ALAT
C
      DO I=1,3
         DO J=1,3
            BRAVAIS(J,I) = 0D0
         END DO
      END DO
C
      DO I=1,3
         DO J=1,3
            RECBV(J,I)=0D0
         ENDDO
      ENDDO
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE(6,'(79(1H=))')
      WRITE (6,'(23X,A)') '  LATTIX99: bulk geometry mode'
      WRITE (6,'(79(1H=))')
      WRITE (6,*)
      WRITE (6,'(5X,A,F12.8,4X,A,F12.8,/)') 
     &     'Lattice constants :  ALAT =',ALAT,' 2*PI/ALAT =',TPIA
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ----------------------------------------------------------------------
C Read in the bravais vectors (normalised to alat)
C Notation: BRAVAIS(J,I) J=x,y,z I=1,2,3
C ----------------------------------------------------------------------
C
      IER = 0
      DO I = 1,3
         CALL IOINPUT('BRAVAIS   ',UIO,I,7,IER)
         READ (UNIT=UIO,FMT=*) (BRAVAIS(J,I),J=1,3)
      END DO
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,'(5X,A,/)') 'Direct lattice cell vectors :'
      WRITE (6,'(9X,A,21X,A)') 'normalised (ALAT)','a.u.'
      WRITE (6,99001)
      DO I = 1,3
         WRITE (6,99002) 'a_',I,(BRAVAIS(J,I),J=1,3),
     &                   (BRAVAIS(J,I)*ALAT,J=1,3)
      END DO
      WRITE (6,99001)
      WRITE (6,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ----------------------------------------------------------------------
C Now generate the reciprocal lattice unit-vectors, 
C and calculate the unit-cell volume in units au**3.
C ----------------------------------------------------------------------
C
         DET = DDET33(BRAVAIS)
         IF ( ABS(DET).LT.1D-8 ) STOP 
     &               ' ERROR: 3D Bravais vectors are linearly dependent'
C
         CALL CROSPR(BRAVAIS(1,2),BRAVAIS(1,3),RECBV(1,1))
         CALL CROSPR(BRAVAIS(1,3),BRAVAIS(1,1),RECBV(1,2))
         CALL CROSPR(BRAVAIS(1,1),BRAVAIS(1,2),RECBV(1,3))
C
         CALL SPATPR(BRAVAIS(1,2),BRAVAIS(1,3),BRAVAIS(1,1),VOLUC)
         VOLUC= ABS(VOLUC)
         DO I=1,3
            DO J=1,3
               RECBV(J,I)=RECBV(J,I)/VOLUC
            ENDDO
         ENDDO   
C ----------------------------------------------------------------------
C
C --> test on volume unit cell:
C
      IF (VOLUC.LT.1.0D-5) STOP 
     &        ' ERROR: Unit-cell volume suspiciously small ( < 1D-5)'

C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,FMT='(5X,A,F8.4,A,F14.8,A,/)') 
     &     'Unit cell volume :  V =',VOLUC,' (ALAT**3) = ',
     &     VOLUC*(ALAT**3),' (a.u.**3)'
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C
C --> check volume of unit cell vs. average WS-radius
C
      VOLUME0 = VOLUC * ALAT**3
C
C ----------------------------------------------------------------------
C  Reciprocal lattice unit-vectors and unit-cell volume calculated
C ----------------------------------------------------------------------
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,'(5X,A,/)') 'Reciprocal lattice cell vectors :'
      WRITE (6,'(9X,A,16X,A)') 'normalised (2*PI/ALAT)','1/a.u.'
      WRITE (6,99001)
      DO I = 1,3
         WRITE (6,99002) 'b_',I,(RECBV(J,I),J=1,3),
     &                      (RECBV(J,I)*TPIA,J=1,3)
      END DO
      WRITE (6,99001)
      WRITE (6,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C --> now generate the real-space lattice vectors for the 
C     cluster generation
C
      CALL RRGEN(BRAVAIS,RR,NR,NRD)
      WRITE(6,*)
C
99001 FORMAT (9X,32(1H-),6X,32(1H-))
99002 FORMAT (5X,A2,I1,':',3F10.6,8X,3F10.6)
      END
