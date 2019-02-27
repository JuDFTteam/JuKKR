!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Generates the real space and reciprocal lattices
!> Author: 
!> Generates the real space and reciprocal lattices, this lattice is **not** the one
!> used for the Eawld summation
!> @note this lattice is **not** the one used for the Ewald summation
!> @endnote
!------------------------------------------------------------------------------------
module mod_lattix99
  use :: mod_datatypes, only: dp
  use :: mod_constants, only: pi
  private :: dp

contains
    !-------------------------------------------------------------------------------
    !> Summary: Generates the real space and reciprocal lattices
    !> Author: 
    !> Category: geometry, k-points, KKRhost
    !> Deprecated: False 
    !> generates the real space and reciprocal lattices.
    !> `BRAVAIS(I,J)` are basis vectors, with `I=X,Y,Z` and `J=A,B,C`
    !> **reciprocal** space vectors are in **units** of \(2\pi/a_{lat}^c\)
    !> `RR` are the direct space vectors `NR+1` is the number of direct space vectors created
    !> (structure dependent output)
    !-------------------------------------------------------------------------------
    !> @note this lattice is **not** the one used for the Ewald summation
    !> @endnote
    !-------------------------------------------------------------------------------
  subroutine lattix99(lsurf,alat,natyp,naez,conc,rws,bravais,recbv,volume0,rr,nrd,  &
    natypd)

    use :: mod_rrgen
    use :: mod_spatpr
    use :: mod_crospr
    use :: mod_ddet33
    use :: mod_idreals
    implicit none
    ! ..
    ! .. Scalar arguments ..
    logical, intent(in) :: lsurf  !! If True a matching with semi-inifinite surfaces must be performed
    integer, intent(inout) :: nrd !! Number of real space vectors, modified in rrgen
    integer, intent(in) :: naez   !! Number of atoms in unit cell
    integer, intent(in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent(in) :: natypd !! Auxiliary number of kinds of atoms in the unit cell
    real (kind=dp), intent(in) :: alat !! Lattice constant in a.u.
    real (kind=dp), intent(out) :: volume0 !! Unit cell volume
    ! ..
    ! .. Array arguments ..

    ! BRAVAIS(3,3): Real space bravais vectors normalised to ALAT
    ! RECBV(3,3)  : Reciprocal lattice vectors in 2*PI/ALAT
    real (kind=dp), dimension(natypd), intent(in) :: conc !! Concentration of a given atom
    real (kind=dp), dimension(natypd), intent(in) :: rws !! Wigner Seitz radius
    real (kind=dp), dimension(3,0:nrd), intent(inout):: rr !! Set of real space vectors (in a.u.), modified in rrgen
    real (kind=dp), dimension(3,3), intent(inout) :: bravais !! Bravais lattice vectors, modified in idreals
    real (kind=dp), dimension(3,3), intent(out) :: recbv !! Reciprocal basis vectors
    ! ..
    ! .. Local Scalars ..
    integer :: i, j, ndim
    integer :: iprint
    real (kind=dp) :: voluc, det, tpia, sws
    ! ..................................................................

    ! --> initialise

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
    write (1337, '(5X,A,F12.8,4X,A,F12.8,/)') 'Lattice constants :  ALAT =', alat, ' 2*PI/ALAT =', tpia
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
        write (1337, 120) 'a_', i, (bravais(j,i), j=1, ndim), (bravais(j,i)*alat, j=1, ndim)
      end do
      write (1337, 100)
    else
      write (1337, 110)
      do i = 1, ndim
        write (1337, 130) 'a_', i, (bravais(j,i), j=1, ndim), (bravais(j,i)*alat, j=1, ndim)
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
      if (abs(det)<1e-8_dp) stop ' ERROR: 3D Bravais vectors are linearly dependent'

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
      if (abs(det)<1e-8_dp) stop ' ERROR: 2D Bravais vectors are linearly dependent'

      recbv(1, 1) = bravais(2, 2)/det
      recbv(2, 1) = -bravais(1, 2)/det
      recbv(1, 2) = -bravais(2, 1)/det
      recbv(2, 2) = bravais(1, 1)/det

      voluc = abs(det)
    end if

    ! --> test on volume unit cell:

    if (voluc<1.0e-5_dp) stop ' ERROR: Unit-cell volume suspiciously small ( < 1D-5)'


    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, '(5X,A,F14.8,A,I1,A,F14.8,A,I1,A,/)') 'Unit cell volume :  V =', voluc, ' (ALAT**', ndim, ') = ', voluc*(alat**ndim), ' (a.u.**', ndim, ')'
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
        write (1337, '(5X,A,A)') 'WARNING : Unit cell volume', ' inconsistent with the average WS-radius'
        write (1337, '(15X,A,F14.8)') 'Unit cell volume        =', volume0
        write (1337, '(15X,A,F14.8)') 'NAEZ * WSRav^3 * 4*PI/3 =', sws
        write (1337, '(15X,A,F14.8,/)') 'difference              =', abs(volume0-sws)
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
        write (1337, 120) 'b_', i, (recbv(j,i), j=1, ndim), (recbv(j,i)*tpia, j=1, ndim)
      end do
      write (1337, 100)
    else
      write (1337, 110)
      do i = 1, ndim
        write (1337, 130) 'b_', i, (recbv(j,i), j=1, ndim), (recbv(j,i)*tpia, j=1, ndim)
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
