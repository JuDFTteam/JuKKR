C*==lattix99.f    processed by SPAG 6.05Rc at 17:56 on 17 May 2004
      SUBROUTINE LATTIX99(LSURF,ALAT,NATYP,NAEZ,CONC,RWS,BRAVAIS,
     &                    RECBV,VOLUME0,RR,NR,NRD,NATYPD)
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
      LOGICAL LSURF
      INTEGER NR,NRD            ! number of real space vectors
      INTEGER IPRINT,NATYP,NAEZ,NATYPD
      DOUBLE PRECISION ALAT,VOLUME0
C     ..
C     .. Array arguments ..
C
C  BRAVAIS(3,3): Real space bravais vectors normalised to ALAT
C  RECBV(3,3)  : Reciprocal lattice vectors in 2*PI/ALAT
C
      DOUBLE PRECISION BRAVAIS(3,3)
      DOUBLE PRECISION RECBV(3,3),RR(3,0:NRD)
      DOUBLE PRECISION CONC(NATYPD),RWS(NATYPD)
C     ..
C     .. Local Scalars ..
      INTEGER I,J,IER,NDIM
      DOUBLE PRECISION VOLUC,DET,DDET33,PI,TPIA,SWS
      CHARACTER*256 UIO ! NCOLIO=256
C     ..
C     .. External declarations ..
      EXTERNAL CROSPR,SPATPR,DDET33,IDREALS,IOINPUT
C     ..
C     .. Intrinsic functions ..
      INTRINSIC ABS,ATAN,DBLE
C     ..
C     ..................................................................
C
C --> initialise
C
      PI = 4D0*DATAN(1D0)
      TPIA = 2D0*PI/ALAT
      IPRINT = 0
C
      RECBV(1:3,1:3)=0D0
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE(1337,'(79(1H=))')
      IF ( LSURF ) THEN
         NDIM = 2
         WRITE (1337,'(23X,A)') 'LATTIX99: surface geometry mode'
      ELSE
         NDIM = 3
         WRITE (1337,'(23X,A)') '  LATTIX99: bulk geometry mode'
      END IF
      WRITE (1337,'(79(1H=))')
      WRITE (1337,*)
      WRITE (1337,'(5X,A,F12.8,4X,A,F12.8,/)') 
     &     'Lattice constants :  ALAT =',ALAT,' 2*PI/ALAT =',TPIA
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ----------------------------------------------------------------------
C Bravais vectors (normalised to alat)
C Notation: BRAVAIS(J,I) J=x,y,z I=1,2,3
C If LSURF=TRUE (2D geometry) the third Bravais vector and z-components
C of all other vectors are left zero
C ----------------------------------------------------------------------
      CALL IDREALS(BRAVAIS(1,1),9,IPRINT)

C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,'(5X,A,/)') 'Direct lattice cell vectors :'
      WRITE (1337,'(9X,A,21X,A)') 'normalised (ALAT)','a.u.'
      IF ( NDIM.EQ.2 ) THEN
         WRITE (1337,99000)
         DO I = 1,NDIM
            WRITE (1337,99002) 'a_',I,(BRAVAIS(J,I),J=1,NDIM),
     &                      (BRAVAIS(J,I)*ALAT,J=1,NDIM)
         END DO
         WRITE (1337,99000)
      ELSE
         WRITE (1337,99001)
         DO I = 1,NDIM
            WRITE (1337,99003) 'a_',I,(BRAVAIS(J,I),J=1,NDIM),
     &                      (BRAVAIS(J,I)*ALAT,J=1,NDIM)
         END DO
         WRITE (1337,99001)
      END IF
      WRITE (1337,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ----------------------------------------------------------------------
C Now generate the reciprocal lattice unit-vectors, 
C and calculate the unit-cell volume in units au**3.
C ----------------------------------------------------------------------
      IF ( .NOT.LSURF ) THEN
C -------------------------------------------------------------- 3D case
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
      ELSE
C -------------------------------------------------------------- 2D case
C
         DET=BRAVAIS(1,1)*BRAVAIS(2,2)-BRAVAIS(1,2)*BRAVAIS(2,1)
         IF ( ABS(DET).LT.1D-8 ) STOP 
     &               ' ERROR: 2D Bravais vectors are linearly dependent'
C         
         RECBV(1,1)=  BRAVAIS(2,2)/DET
         RECBV(2,1)= -BRAVAIS(1,2)/DET
         RECBV(1,2)= -BRAVAIS(2,1)/DET
         RECBV(2,2)=  BRAVAIS(1,1)/DET
C         
         VOLUC= ABS(DET)
      ENDIF
C
C --> test on volume unit cell:
C
      IF (VOLUC.LT.1.0D-5) STOP 
     &        ' ERROR: Unit-cell volume suspiciously small ( < 1D-5)'

C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,'(5X,A,F14.8,A,I1,A,F14.8,A,I1,A,/)') 
     &     'Unit cell volume :  V =',VOLUC,' (ALAT**',NDIM,
     &     ') = ',VOLUC*(ALAT**NDIM),' (a.u.**',NDIM,')'
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C
C --> check volume of unit cell vs. average WS-radius
C
      VOLUME0 = VOLUC * ALAT**(NDIM)
      IF ( .NOT.LSURF ) THEN
         SWS = 00D0
         DO I=1,NATYP
!            SWS = SWS + CONC(I)*NAT(I)*RWS(I)**3  ! Array NAT removed (was=1) Phivos 13.10.14
            SWS = SWS + CONC(I)*RWS(I)**3
         END DO
         SWS = (SWS/DBLE(NAEZ))**(1D0/3D0)
         SWS = DBLE(NAEZ)*SWS**3*4D0*PI/3D0
         IF( ABS(VOLUME0-SWS).GT.1D-5 ) THEN 
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
            WRITE(1337,'(5X,A,A)') 'WARNING : Unit cell volume',
     &           ' inconsistent with the average WS-radius'
            WRITE(1337,'(15X,A,F14.8)') 'Unit cell volume        =',
     &                                                          VOLUME0
            WRITE(1337,'(15X,A,F14.8)') 'NAEZ * WSRav^3 * 4*PI/3 =',SWS
            WRITE(1337,'(15X,A,F14.8,/)') 'difference              =',
     &                               ABS(VOLUME0-SWS)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
         END IF
      END IF
C
C ----------------------------------------------------------------------
C  Reciprocal lattice unit-vectors and unit-cell volume calculated
C ----------------------------------------------------------------------
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,'(5X,A,/)') 'Reciprocal lattice cell vectors :'
      WRITE (1337,'(9X,A,16X,A)') 'normalised (2*PI/ALAT)','1/a.u.'
      IF ( NDIM.EQ.2 ) THEN
         WRITE (1337,99000)
         DO I = 1,NDIM
            WRITE (1337,99002) 'b_',I,(RECBV(J,I),J=1,NDIM),
     &                      (RECBV(J,I)*TPIA,J=1,NDIM)
         END DO
         WRITE (1337,99000)
      ELSE
         WRITE (1337,99001)
         DO I = 1,NDIM
            WRITE (1337,99003) 'b_',I,(RECBV(J,I),J=1,NDIM),
     &                      (RECBV(J,I)*TPIA,J=1,NDIM)
         END DO
         WRITE (1337,99001)
      END IF
      WRITE (1337,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C --> now generate the real-space lattice vectors for the 
C     cluster generation
C
      CALL RRGEN(BRAVAIS,LSURF,RR,NR,NRD)
      WRITE(1337,*)
C
99000 FORMAT (9X,22(1H-),16X,22(1H-))
99001 FORMAT (9X,32(1H-),6X,32(1H-))
99002 FORMAT (5X,A2,I1,':',2F10.6,18X,2F10.6)
99003 FORMAT (5X,A2,I1,':',3F10.6,8X,3F10.6)
      END
