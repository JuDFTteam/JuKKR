!*==lattix99.f    processed by SPAG 6.05Rc at 17:56 on 17 May 2004
subroutine LATTIX99(ALAT,BRAVAIS,RECBV,VOLUME0, output)
  ! **********************************************************************
  ! *                                                                    *
  ! * LATTIX99 generates the real space and reciprocal lattices.         *
  ! * BRAVAIS(I,J) are basis vectors, with I=X,Y,Z and J=A,B,C           *
  ! *              input - normalised vectors                            *
  ! * RECIPROCAL space vectors are in UNITS OF 2*PI/ALATC - output       *
  ! * RR are the direct space vectors - output                           *
  ! * NR+1 is the number of direct space vectors created - output        *
  ! * (structure dependent output).                                      *
  ! *                                                                    *
  ! **********************************************************************
  implicit none
  !     ..
  !     .. Scalar arguments ..
  logical output

  double precision ALAT,VOLUME0
  !     ..
  !     .. Array arguments ..

  !  BRAVAIS(3,3): Real space bravais vectors. Read in normalised to ALAT
  !  RECBV(3,3)  : Reciprocal lattice vectors in 2*PI/ALAT

  double precision BRAVAIS(3,3)
  double precision RECBV(3,3)
  !     ..
  !     .. Local Scalars ..
  integer I,J
  double precision VOLUC,DET,DDET33,PI,TPIA
  !     ..
  !     .. External declarations ..
  external CROSPR,SPATPR,DDET33,IOINPUT
  !     ..
  !     .. Intrinsic functions ..
  intrinsic ABS,ATAN,DBLE
  !     ..
  !     ..................................................................

  ! --> initialise

  PI = 4D0*ATAN(1D0)
  TPIA = 2D0*PI/ALAT

  do I=1,3
    do J=1,3
      RECBV(J,I)=0D0
    enddo
  enddo

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  if (output) then
    write(6,'(79(1H=))')
    write (6,'(23X,A)') '  LATTIX99: bulk geometry mode'
    write (6,'(79(1H=))')
    write (6,*)
    write (6,'(5X,A,F12.8,4X,A,F12.8,/)') &
    'Lattice constants :  ALAT =',ALAT,' 2*PI/ALAT =',TPIA
  end if
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  if (output) then
    write (6,'(5X,A,/)') 'Direct lattice cell vectors :'
    write (6,'(9X,A,21X,A)') 'normalised (ALAT)','a.u.'
    write (6,99001)
    do I = 1,3
      write (6,99002) 'a_',I,(BRAVAIS(J,I),J=1,3), &
      (BRAVAIS(J,I)*ALAT,J=1,3)
    end do
    write (6,99001)
    write (6,*)
  end if
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  ! ----------------------------------------------------------------------
  ! Now generate the reciprocal lattice unit-vectors,
  ! and calculate the unit-cell volume in units au**3.
  ! ----------------------------------------------------------------------

  DET = DDET33(BRAVAIS)
  if ( ABS(DET).lt.1D-8 ) stop &
  ' ERROR: 3D Bravais vectors are linearly dependent'

  call CROSPR(BRAVAIS(1,2),BRAVAIS(1,3),RECBV(1,1))
  call CROSPR(BRAVAIS(1,3),BRAVAIS(1,1),RECBV(1,2))
  call CROSPR(BRAVAIS(1,1),BRAVAIS(1,2),RECBV(1,3))

  call SPATPR(BRAVAIS(1,2),BRAVAIS(1,3),BRAVAIS(1,1),VOLUC)
  VOLUC= ABS(VOLUC)
  do I=1,3
    do J=1,3
      RECBV(J,I)=RECBV(J,I)/VOLUC
    enddo
  enddo
  ! ----------------------------------------------------------------------

  ! --> test on volume unit cell:

  if (VOLUC.lt.1.0D-5) stop &
  ' ERROR: Unit-cell volume suspiciously small ( < 1D-5)'


  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  if (output) then
    write (6,fmt='(5X,A,F8.4,A,F14.8,A,/)') &
    'Unit cell volume :  V =',VOLUC,' (ALAT**3) = ', &
    VOLUC*(ALAT**3),' (a.u.**3)'
  end if
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


  ! --> check volume of unit cell vs. average WS-radius

  VOLUME0 = VOLUC * ALAT**3

  ! ----------------------------------------------------------------------
  !  Reciprocal lattice unit-vectors and unit-cell volume calculated
  ! ----------------------------------------------------------------------

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  if (output) then
    write (6,'(5X,A,/)') 'Reciprocal lattice cell vectors :'
    write (6,'(9X,A,16X,A)') 'normalised (2*PI/ALAT)','1/a.u.'
    write (6,99001)
    do I = 1,3
      write (6,99002) 'b_',I,(RECBV(J,I),J=1,3), &
      (RECBV(J,I)*TPIA,J=1,3)
    end do
    write (6,99001)
    write (6,*)
  end if
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


99001 format (9X,32(1H-),6X,32(1H-))
99002 format (5X,A2,I1,':',3F10.6,8X,3F10.6)
end
