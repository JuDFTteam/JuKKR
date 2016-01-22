C ************************************************************************
      SUBROUTINE LATTIX12(LSURF,ALATC,BRAVAIS,RECBV,RR,NR,VOLUC)
C ************************************************************************
C LATTIX99 GENERATES THE REAL SPACE AND RECIPROCAL LATTICES.
C BRAVAIS(I,J) ARE BASIS VECTORS, WITH I=X,Y,Z AND J=A,B,C.
C RECIPROCAL SPACE VECTORS ARE IN UNITS OF 2*PI/A0.
C RR ARE THE DIRECT SPACE VECTORS.
C NR+1 IS THE NUMBER OF DIRECT SPACE VECTORS CREATED
C (STRUCTURE DEPENDENT OUTPUT). 
C ************************************************************************
C
      implicit none
      INCLUDE 'inc.geometry'
C
      LOGICAL LSURF
C
      INTEGER
     +     NR,                       ! number of real space vectors
     +     I,J,IER
C
      REAL*8       
     +     VOLUC,ALATC,
     +     DET,PI,TPI,DDET33
C
      REAL*8       
     +     BRAVAIS(3,3),  ! Real space bravais vectors. Read in normalized to
c                         ! altc.
     +     RECBV(3,3),              ! RECIPROCAL BASIS IN 2*PI/A
     +     RR(3,0:NRD)
c
      CHARACTER*200 UIO
c
      EXTERNAL CROSPR,SPATPR,VADD,VSUB,VEQ,DDET33,IOINPUT
c
      PARAMETER (PI   = 3.14159265358979312D0)
c
c ------------------------------------------------------------------------
c
c Read in the bravais vectors (normalized to alatc)
c Notation: BRAVAIS(J,I) J=x,y,z I=1,2,3
c If the third bravais vector is zero, then surface (2-dimentional) geometry
c is implied.
c
      DO I=1,3
         CALL IoInput('BRAVAIS   ',UIO,I,7,IER)
              READ (UNIT=UIO,FMT=*) (BRAVAIS(J,I), J=1,3)
      ENDDO

      IF (BRAVAIS(1,3).EQ.0.D0.AND.BRAVAIS(2,3).EQ.0.D0.
     &                          AND.BRAVAIS(3,3).EQ.0.D0) THEN
         LSURF=.TRUE.
         WRITE(*,*) 'Surface geometry'
         IF (BRAVAIS(3,1).NE.0.D0.OR.BRAVAIS(3,2).NE.0.D0) THEN
            WRITE(6,9010)
            BRAVAIS(3,1) = 0.D0
            BRAVAIS(3,2) = 0.D0
         ENDIF
      ENDIF

      IF (LSURF) THEN
         BRAVAIS(1:3,3) = 0.D0
      ENDIF
c
c
c Now generate the reciprocal lattice unit-vectors, and calculate the unit-cell
c volume in units au**3.
c Initialize:
      DO I=1,3
         DO J=1,3
            RECBV(J,I)=0.D0
         ENDDO
      ENDDO
c Calculate:
      IF (.NOT.(LSURF)) THEN
c 1. 3-dimentional case
         DET=DDET33(BRAVAIS)
         IF (DET.EQ.0.D0) THEN
            WRITE(*,*) '3-d Bravais vectors are linearly dependent:'
            WRITE(*,9000) ((BRAVAIS(J,I),J=1,3),I=1,3)
            STOP 'LATTIX99: 3-d Bravais vectors are 
     +                          linearly dependent'
         END IF

         CALL CROSPR(BRAVAIS(1,2),BRAVAIS(1,3),RECBV(1,1))
         CALL CROSPR(BRAVAIS(1,3),BRAVAIS(1,1),RECBV(1,2))
         CALL CROSPR(BRAVAIS(1,1),BRAVAIS(1,2),RECBV(1,3))

         CALL SPATPR(BRAVAIS(1,2),BRAVAIS(1,3),BRAVAIS(1,1),VOLUC)
         VOLUC= ABS(VOLUC)
         DO I=1,3
            DO J=1,3
               RECBV(J,I)=RECBV(J,I)/VOLUC
            ENDDO
         ENDDO   
      ELSE
c     2. 2-dimentional case
         DET=BRAVAIS(1,1)*BRAVAIS(2,2)-BRAVAIS(1,2)*BRAVAIS(2,1)
         IF (ABS(DET).LT.1.D-8) THEN 
            WRITE(*,*)    '2-d Bravais vectors are linearly dependent:'
            WRITE(*,9000) ((BRAVAIS(J,I),J=1,3),I=1,3)
            STOP 'LATTIX99: 2-d Bravais vectors are  linearly dependent'
         ENDIF
         
         RECBV(1,1)=  BRAVAIS(2,2)/DET
         RECBV(2,1)= -BRAVAIS(1,2)/DET
         RECBV(1,2)= -BRAVAIS(2,1)/DET
         RECBV(2,2)=  BRAVAIS(1,1)/DET
         
         VOLUC= ABS(DET)
      ENDIF
      
c     Reciprocal lattice unit-vectors and unit-cell volume calculated.
c     
      write(6,*) 'True basis vectors (not normalized):'
      write(6,9000) ((BRAVAIS(J,I)*ALATC,J=1,3),I=1,3)
      write(6,*) 'Reciprocal lattice vectors, in units 2pi/a:'
      write(6,9000) ((RECBV(J,I),J=1,3),I=1,3)
c
c ---> now generate the real-space lattice vectors for the cluster generation:
c      The parameter LATT is not used in RRGEN
      CALL RRGEN(BRAVAIS,ALATC,LSURF,RR,NR)
c
c ---> test on volume unit cell:
c
      IF (VOLUC.LT.1.0D-5) THEN
         WRITE(6,909) VOLUC
  909    FORMAT(//' STOP. VOLUC IN LATTIX99 IS TOO SMALL:',
     +   ' VOLUC= ',D13.4)
         STOP
      END IF
C
 1000 FORMAT(i4,3F15.10)
 1001 FORMAT( ' vol,voluc: ',2F10.4)
 9000 FORMAT(3F14.8/3F14.8/3F14.8)
 9010 FORMAT('Caution: the two-dimentional Bravais vectors have also',/,
     &     'a z-component. It will be set to zero.')
      RETURN
      END

