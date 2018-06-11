    Subroutine lattix99(lsurf, alat, natyp, naez, conc, rws, bravais, recbv, &
      volume0, rr, nr, nrd, natypd)
      Use mod_datatypes, Only: dp
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
      Implicit None
!..
!.. Scalar arguments ..
      Logical :: lsurf
      Integer :: nr, nrd ! number of real space vectors
      Integer :: iprint, natyp, naez, natypd
      Real (Kind=dp) :: alat, volume0
!..
!.. Array arguments ..
!
!  BRAVAIS(3,3): Real space bravais vectors normalised to ALAT
!  RECBV(3,3)  : Reciprocal lattice vectors in 2*PI/ALAT
!
      Real (Kind=dp) :: bravais(3, 3)
      Real (Kind=dp) :: recbv(3, 3), rr(3, 0:nrd)
      Real (Kind=dp) :: conc(natypd), rws(natypd)
!..
!.. Local Scalars ..
      Integer :: i, j, ndim
      Real (Kind=dp) :: voluc, det, ddet33, pi, tpia, sws
!..
!.. External declarations ..
      External :: crospr, spatpr, ddet33, idreals, ioinput
!..
!.. Intrinsic functions ..
      Intrinsic :: abs, atan, real
!     ..
!     ..................................................................

! --> initialise

      pi = 4E0_dp*atan(1E0_dp)
      tpia = 2E0_dp*pi/alat
      iprint = 0

      recbv(1:3, 1:3) = 0E0_dp

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, '(79("="))')
      If (lsurf) Then
        ndim = 2
        Write (1337, '(23X,A)') 'LATTIX99: surface geometry mode'
      Else
        ndim = 3
        Write (1337, '(23X,A)') '  LATTIX99: bulk geometry mode'
      End If
      Write (1337, '(79("="))')
      Write (1337, *)
      Write (1337, '(5X,A,F12.8,4X,A,F12.8,/)') 'Lattice constants :  ALAT =', &
        alat, ' 2*PI/ALAT =', tpia
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ----------------------------------------------------------------------
! Bravais vectors (normalised to alat)
! Notation: BRAVAIS(J,I) J=x,y,z I=1,2,3
! If LSURF=TRUE (2D geometry) the third Bravais vector and z-components
! of all other vectors are left zero
! ----------------------------------------------------------------------
      Call idreals(bravais(1,1), 9, iprint)

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, '(5X,A,/)') 'Direct lattice cell vectors :'
      Write (1337, '(9X,A,21X,A)') 'normalised (ALAT)', 'a.u.'
      If (ndim==2) Then
        Write (1337, 100)
        Do i = 1, ndim
          Write (1337, 120) 'a_', i, (bravais(j,i), j=1, ndim), &
            (bravais(j,i)*alat, j=1, ndim)
        End Do
        Write (1337, 100)
      Else
        Write (1337, 110)
        Do i = 1, ndim
          Write (1337, 130) 'a_', i, (bravais(j,i), j=1, ndim), &
            (bravais(j,i)*alat, j=1, ndim)
        End Do
        Write (1337, 110)
      End If
      Write (1337, *)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ----------------------------------------------------------------------
! Now generate the reciprocal lattice unit-vectors,
! and calculate the unit-cell volume in units au**3.
! ----------------------------------------------------------------------
      If (.Not. lsurf) Then
! -------------------------------------------------------------- 3D case

        det = ddet33(bravais)
        If (abs(det)<1E-8_dp) Stop &
          ' ERROR: 3D Bravais vectors are linearly dependent'

        Call crospr(bravais(1,2), bravais(1,3), recbv(1,1))
        Call crospr(bravais(1,3), bravais(1,1), recbv(1,2))
        Call crospr(bravais(1,1), bravais(1,2), recbv(1,3))

        Call spatpr(bravais(1,2), bravais(1,3), bravais(1,1), voluc)
        voluc = abs(voluc)
        Do i = 1, 3
          Do j = 1, 3
            recbv(j, i) = recbv(j, i)/voluc
          End Do
        End Do
! ----------------------------------------------------------------------
      Else
! -------------------------------------------------------------- 2D case

        det = bravais(1, 1)*bravais(2, 2) - bravais(1, 2)*bravais(2, 1)
        If (abs(det)<1E-8_dp) Stop &
          ' ERROR: 2D Bravais vectors are linearly dependent'

        recbv(1, 1) = bravais(2, 2)/det
        recbv(2, 1) = -bravais(1, 2)/det
        recbv(1, 2) = -bravais(2, 1)/det
        recbv(2, 2) = bravais(1, 1)/det

        voluc = abs(det)
      End If

! --> test on volume unit cell:

      If (voluc<1.0E-5_dp) Stop &
        ' ERROR: Unit-cell volume suspiciously small ( < 1D-5)'


! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, '(5X,A,F14.8,A,I1,A,F14.8,A,I1,A,/)') &
        'Unit cell volume :  V =', voluc, ' (ALAT**', ndim, ') = ', &
        voluc*(alat**ndim), ' (a.u.**', ndim, ')'
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


! --> check volume of unit cell vs. average WS-radius

      volume0 = voluc*alat**(ndim)
      If (.Not. lsurf) Then
        sws = 00E0_dp
        Do i = 1, natyp
!            SWS = SWS + CONC(I)*NAT(I)*RWS(I)**3  ! Array NAT removed (was=1) Phivos 13.10.14
          sws = sws + conc(i)*rws(i)**3
        End Do
        sws = (sws/real(naez,kind=dp))**(1E0_dp/3E0_dp)
        sws = real(naez, kind=dp)*sws**3*4E0_dp*pi/3E0_dp
        If (abs(volume0-sws)>1E-5_dp) Then
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
          Write (1337, '(5X,A,A)') 'WARNING : Unit cell volume', &
            ' inconsistent with the average WS-radius'
          Write (1337, '(15X,A,F14.8)') 'Unit cell volume        =', volume0
          Write (1337, '(15X,A,F14.8)') 'NAEZ * WSRav^3 * 4*PI/3 =', sws
          Write (1337, '(15X,A,F14.8,/)') 'difference              =', &
            abs(volume0-sws)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        End If
      End If

! ----------------------------------------------------------------------
!  Reciprocal lattice unit-vectors and unit-cell volume calculated
! ----------------------------------------------------------------------

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, '(5X,A,/)') 'Reciprocal lattice cell vectors :'
      Write (1337, '(9X,A,16X,A)') 'normalised (2*PI/ALAT)', '1/a.u.'
      If (ndim==2) Then
        Write (1337, 100)
        Do i = 1, ndim
          Write (1337, 120) 'b_', i, (recbv(j,i), j=1, ndim), &
            (recbv(j,i)*tpia, j=1, ndim)
        End Do
        Write (1337, 100)
      Else
        Write (1337, 110)
        Do i = 1, ndim
          Write (1337, 130) 'b_', i, (recbv(j,i), j=1, ndim), &
            (recbv(j,i)*tpia, j=1, ndim)
        End Do
        Write (1337, 110)
      End If
      Write (1337, *)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! --> now generate the real-space lattice vectors for the
!     cluster generation

      Call rrgen(bravais, lsurf, rr, nr, nrd)
      Write (1337, *)

100   Format (9X, 22('-'), 16X, 22('-'))
110   Format (9X, 32('-'), 6X, 32('-'))
120   Format (5X, A2, I1, ':', 2F10.6, 18X, 2F10.6)
130   Format (5X, A2, I1, ':', 3F10.6, 8X, 3F10.6)
    End Subroutine
