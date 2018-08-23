module mod_lattix99

contains

subroutine lattix99(lsurf, alat, natyp, naez, conc, rws, bravais, recbv, &
  volume0, rr, nrd, natypd)
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
  use :: mod_datatypes, only: dp
   use mod_rrgen
   use mod_spatpr
   use mod_crospr
   use mod_ddet33
   use mod_idreals
  implicit none
  ! ..
  ! .. Scalar arguments ..
  logical :: lsurf
  integer :: nrd                   ! number of real space vectors
  integer :: iprint, natyp, naez, natypd
  real (kind=dp) :: alat, volume0
  ! ..
  ! .. Array arguments ..

  ! BRAVAIS(3,3): Real space bravais vectors normalised to ALAT
  ! RECBV(3,3)  : Reciprocal lattice vectors in 2*PI/ALAT

  real (kind=dp) :: bravais(3, 3)
  real (kind=dp) :: recbv(3, 3), rr(3, 0:nrd)
  real (kind=dp) :: conc(natypd), rws(natypd)
  ! ..
  ! .. Local Scalars ..
  integer :: i, j, ndim
  real (kind=dp) :: voluc, det, ddet33, pi, tpia, sws
  ! ..
  ! .. External declarations ..
  external :: crospr, spatpr, ddet33, idreals, ioinput
  ! ..
  ! .. Intrinsic functions ..
  intrinsic :: abs, atan, real
  ! ..
  ! ..................................................................

  ! --> initialise

  pi = 4e0_dp*atan(1e0_dp)
  tpia = 2e0_dp*pi/alat
  iprint = 0

  recbv(1:3, 1:3) = 0e0_dp

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  write (1337, '(79("="))')
  if (lsurf) then
    ndim = 2
    write (1337, '(23X,A)') 'LATTIX99: surface geometry mode'
  else
    ndim = 3
    write (1337, '(23X,A)') '  LATTIX99: bulk geometry mode'
  end if
  write (1337, '(79("="))')
  write (1337, *)
  write (1337, '(5X,A,F12.8,4X,A,F12.8,/)') 'Lattice constants :  ALAT =', &
    alat, ' 2*PI/ALAT =', tpia
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  ! ----------------------------------------------------------------------
  ! Bravais vectors (normalised to alat)
  ! Notation: BRAVAIS(J,I) J=x,y,z I=1,2,3
  ! If LSURF=TRUE (2D geometry) the third Bravais vector and z-components
  ! of all other vectors are left zero
  ! ----------------------------------------------------------------------
  call idreals(bravais(1,1), 9, iprint)

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  write (1337, '(5X,A,/)') 'Direct lattice cell vectors :'
  write (1337, '(9X,A,21X,A)') 'normalised (ALAT)', 'a.u.'
  if (ndim==2) then
    write (1337, 100)
    do i = 1, ndim
      write (1337, 120) 'a_', i, (bravais(j,i), j=1, ndim), &
        (bravais(j,i)*alat, j=1, ndim)
    end do
    write (1337, 100)
  else
    write (1337, 110)
    do i = 1, ndim
      write (1337, 130) 'a_', i, (bravais(j,i), j=1, ndim), &
        (bravais(j,i)*alat, j=1, ndim)
    end do
    write (1337, 110)
  end if
  write (1337, *)
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  ! ----------------------------------------------------------------------
  ! Now generate the reciprocal lattice unit-vectors,
  ! and calculate the unit-cell volume in units au**3.
  ! ----------------------------------------------------------------------
  if (.not. lsurf) then
    ! -------------------------------------------------------------- 3D case

    det = ddet33(bravais)
    if (abs(det)<1e-8_dp) stop &
      ' ERROR: 3D Bravais vectors are linearly dependent'

    call crospr(bravais(1,2), bravais(1,3), recbv(1,1))
    call crospr(bravais(1,3), bravais(1,1), recbv(1,2))
    call crospr(bravais(1,1), bravais(1,2), recbv(1,3))

    call spatpr(bravais(1,2), bravais(1,3), bravais(1,1), voluc)
    voluc = abs(voluc)
    do i = 1, 3
      do j = 1, 3
        recbv(j, i) = recbv(j, i)/voluc
      end do
    end do
    ! ----------------------------------------------------------------------
  else
    ! -------------------------------------------------------------- 2D case

    det = bravais(1, 1)*bravais(2, 2) - bravais(1, 2)*bravais(2, 1)
    if (abs(det)<1e-8_dp) stop &
      ' ERROR: 2D Bravais vectors are linearly dependent'

    recbv(1, 1) = bravais(2, 2)/det
    recbv(2, 1) = -bravais(1, 2)/det
    recbv(1, 2) = -bravais(2, 1)/det
    recbv(2, 2) = bravais(1, 1)/det

    voluc = abs(det)
  end if

  ! --> test on volume unit cell:

  if (voluc<1.0e-5_dp) stop &
    ' ERROR: Unit-cell volume suspiciously small ( < 1D-5)'


  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  write (1337, '(5X,A,F14.8,A,I1,A,F14.8,A,I1,A,/)') &
    'Unit cell volume :  V =', voluc, ' (ALAT**', ndim, ') = ', &
    voluc*(alat**ndim), ' (a.u.**', ndim, ')'
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


  ! --> check volume of unit cell vs. average WS-radius

  volume0 = voluc*alat**(ndim)
  if (.not. lsurf) then
    sws = 00e0_dp
    do i = 1, natyp
      ! SWS = SWS + CONC(I)*NAT(I)*RWS(I)**3  ! Array NAT removed (was=1)
      ! Phivos 13.10.14
      sws = sws + conc(i)*rws(i)**3
    end do
    sws = (sws/real(naez,kind=dp))**(1e0_dp/3e0_dp)
    sws = real(naez, kind=dp)*sws**3*4e0_dp*pi/3e0_dp
    if (abs(volume0-sws)>1e-5_dp) then
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      write (1337, '(5X,A,A)') 'WARNING : Unit cell volume', &
        ' inconsistent with the average WS-radius'
      write (1337, '(15X,A,F14.8)') 'Unit cell volume        =', volume0
      write (1337, '(15X,A,F14.8)') 'NAEZ * WSRav^3 * 4*PI/3 =', sws
      write (1337, '(15X,A,F14.8,/)') 'difference              =', &
        abs(volume0-sws)
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    end if
  end if

  ! ----------------------------------------------------------------------
  ! Reciprocal lattice unit-vectors and unit-cell volume calculated
  ! ----------------------------------------------------------------------

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  write (1337, '(5X,A,/)') 'Reciprocal lattice cell vectors :'
  write (1337, '(9X,A,16X,A)') 'normalised (2*PI/ALAT)', '1/a.u.'
  if (ndim==2) then
    write (1337, 100)
    do i = 1, ndim
      write (1337, 120) 'b_', i, (recbv(j,i), j=1, ndim), &
        (recbv(j,i)*tpia, j=1, ndim)
    end do
    write (1337, 100)
  else
    write (1337, 110)
    do i = 1, ndim
      write (1337, 130) 'b_', i, (recbv(j,i), j=1, ndim), &
        (recbv(j,i)*tpia, j=1, ndim)
    end do
    write (1337, 110)
  end if
  write (1337, *)
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  ! --> now generate the real-space lattice vectors for the
  ! cluster generation

  call rrgen(bravais, lsurf, rr, nrd)
  write (1337, *)

100 format (9x, 22('-'), 16x, 22('-'))
110 format (9x, 32('-'), 6x, 32('-'))
120 format (5x, a2, i1, ':', 2f10.6, 18x, 2f10.6)
130 format (5x, a2, i1, ':', 3f10.6, 8x, 3f10.6)
end subroutine lattix99

end module mod_lattix99
