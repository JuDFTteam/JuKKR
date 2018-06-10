SUBROUTINE lattix99(lsurf,alat,natyp,naez,conc,rws,bravais,  &
    recbv,volume0,rr,nr,nrd,natypd)
! **********************************************************************
! *                                                                    *
! * LATTIX99 generates the real space and reciprocal lattices.         *
! * BRAVAIS(I,J) are basis vectors, with I=X,Y,Z and J=A,B,C..         *
! * RECIPROCAL space vectors are in UNITS OF 2*PI/ALATC..              *
! * RR are the direct space vectors                                    *
! * NR+1 is the number of direct space vectors created                 *
! * (structure dependent output).                                      *
! *                                                                    *
! **********************************************************************
      IMPLICIT NONE
!..
!.. Scalar arguments ..
      LOGICAL LSURF
      INTEGER NR,NRD            ! number of real space vectors
      INTEGER IPRINT,NATYP,NAEZ,NATYPD
      DOUBLE PRECISION ALAT,VOLUME0
!..
!.. Array arguments ..
!
!  BRAVAIS(3,3): Real space bravais vectors normalised to ALAT
!  RECBV(3,3)  : Reciprocal lattice vectors in 2*PI/ALAT
!
      DOUBLE PRECISION BRAVAIS(3,3)
      DOUBLE PRECISION RECBV(3,3),RR(3,0:NRD)
      DOUBLE PRECISION CONC(NATYPD),RWS(NATYPD)
!..
!.. Local Scalars ..
      INTEGER I,J,NDIM
      DOUBLE PRECISION VOLUC,DET,DDET33,PI,TPIA,SWS
!..
!.. External declarations ..
      EXTERNAL CROSPR,SPATPR,DDET33,IDREALS,IOINPUT
!..
!.. Intrinsic functions ..
      INTRINSIC ABS,ATAN,DBLE
!     ..
!     ..................................................................

! --> initialise

pi = 4D0*DATAN(1D0)
tpia = 2D0*pi/alat
iprint = 0

recbv(1:3,1:3)=0D0

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE(1337,'(79(1H=))')
IF ( lsurf ) THEN
  ndim = 2
  WRITE (1337,'(23X,A)') 'LATTIX99: surface geometry mode'
ELSE
  ndim = 3
  WRITE (1337,'(23X,A)') '  LATTIX99: bulk geometry mode'
END IF
WRITE (1337,'(79(1H=))')
WRITE (1337,*)
WRITE (1337,'(5X,A,F12.8,4X,A,F12.8,/)')  &
    'Lattice constants :  ALAT =',alat,' 2*PI/ALAT =',tpia
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ----------------------------------------------------------------------
! Bravais vectors (normalised to alat)
! Notation: BRAVAIS(J,I) J=x,y,z I=1,2,3
! If LSURF=TRUE (2D geometry) the third Bravais vector and z-components
! of all other vectors are left zero
! ----------------------------------------------------------------------
CALL idreals(bravais(1,1),9,iprint)

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,'(5X,A,/)') 'Direct lattice cell vectors :'
WRITE (1337,'(9X,A,21X,A)') 'normalised (ALAT)','a.u.'
IF ( ndim == 2 ) THEN
  WRITE (1337,99000)
  DO i = 1,ndim
    WRITE (1337,99002) 'a_',i,(bravais(j,i),j=1,ndim),  &
        (bravais(j,i)*alat,j=1,ndim)
  END DO
  WRITE (1337,99000)
ELSE
  WRITE (1337,99001)
  DO i = 1,ndim
    WRITE (1337,99003) 'a_',i,(bravais(j,i),j=1,ndim),  &
        (bravais(j,i)*alat,j=1,ndim)
  END DO
  WRITE (1337,99001)
END IF
WRITE (1337,*)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ----------------------------------------------------------------------
! Now generate the reciprocal lattice unit-vectors,
! and calculate the unit-cell volume in units au**3.
! ----------------------------------------------------------------------
IF ( .NOT.lsurf ) THEN
! -------------------------------------------------------------- 3D case
  
  det = ddet33(bravais)
  IF ( ABS(det) < 1D-8 ) STOP  &
      ' ERROR: 3D Bravais vectors are linearly dependent'
  
  CALL crospr(bravais(1,2),bravais(1,3),recbv(1,1))
  CALL crospr(bravais(1,3),bravais(1,1),recbv(1,2))
  CALL crospr(bravais(1,1),bravais(1,2),recbv(1,3))
  
  CALL spatpr(bravais(1,2),bravais(1,3),bravais(1,1),voluc)
  voluc= ABS(voluc)
  DO i=1,3
    DO j=1,3
      recbv(j,i)=recbv(j,i)/voluc
    END DO
  END DO
! ----------------------------------------------------------------------
ELSE
! -------------------------------------------------------------- 2D case
  
  det=bravais(1,1)*bravais(2,2)-bravais(1,2)*bravais(2,1)
  IF ( ABS(det) < 1D-8 ) STOP  &
      ' ERROR: 2D Bravais vectors are linearly dependent'
  
  recbv(1,1)=  bravais(2,2)/det
  recbv(2,1)= -bravais(1,2)/det
  recbv(1,2)= -bravais(2,1)/det
  recbv(2,2)=  bravais(1,1)/det
  
  voluc= ABS(det)
END IF

! --> test on volume unit cell:

IF (voluc < 1.0D-5) STOP  &
    ' ERROR: Unit-cell volume suspiciously small ( < 1D-5)'


! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,'(5X,A,F14.8,A,I1,A,F14.8,A,I1,A,/)')  &
    'Unit cell volume :  V =',voluc,' (ALAT**',ndim,  &
    ') = ',voluc*(alat**ndim),' (a.u.**',ndim,')'
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


! --> check volume of unit cell vs. average WS-radius

volume0 = voluc * alat**(ndim)
IF ( .NOT.lsurf ) THEN
  sws = 00D0
  DO i=1,natyp
!            SWS = SWS + CONC(I)*NAT(I)*RWS(I)**3  ! Array NAT removed (was=1) Phivos 13.10.14
    sws = sws + conc(i)*rws(i)**3
  END DO
  sws = (sws/DBLE(naez))**(1D0/3D0)
  sws = DBLE(naez)*sws**3*4D0*pi/3D0
  IF( ABS(volume0-sws) > 1D-5 ) THEN
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    WRITE(1337,'(5X,A,A)') 'WARNING : Unit cell volume',  &
        ' inconsistent with the average WS-radius'
    WRITE(1337,'(15X,A,F14.8)') 'Unit cell volume        =', volume0
    WRITE(1337,'(15X,A,F14.8)') 'NAEZ * WSRav^3 * 4*PI/3 =',sws
    WRITE(1337,'(15X,A,F14.8,/)') 'difference              =',  &
        ABS(volume0-sws)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  END IF
END IF

! ----------------------------------------------------------------------
!  Reciprocal lattice unit-vectors and unit-cell volume calculated
! ----------------------------------------------------------------------

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,'(5X,A,/)') 'Reciprocal lattice cell vectors :'
WRITE (1337,'(9X,A,16X,A)') 'normalised (2*PI/ALAT)','1/a.u.'
IF ( ndim == 2 ) THEN
  WRITE (1337,99000)
  DO i = 1,ndim
    WRITE (1337,99002) 'b_',i,(recbv(j,i),j=1,ndim),  &
        (recbv(j,i)*tpia,j=1,ndim)
  END DO
  WRITE (1337,99000)
ELSE
  WRITE (1337,99001)
  DO i = 1,ndim
    WRITE (1337,99003) 'b_',i,(recbv(j,i),j=1,ndim),  &
        (recbv(j,i)*tpia,j=1,ndim)
  END DO
  WRITE (1337,99001)
END IF
WRITE (1337,*)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! --> now generate the real-space lattice vectors for the
!     cluster generation

CALL rrgen(bravais,lsurf,rr,nr,nrd)
WRITE(1337,*)

99000 FORMAT (9X,22(1H-),16X,22(1H-))
99001 FORMAT (9X,32(1H-),6X,32(1H-))
99002 FORMAT (5X,a2,i1,':',2F10.6,18X,2F10.6)
99003 FORMAT (5X,a2,i1,':',3F10.6,8X,3F10.6)
END SUBROUTINE lattix99
