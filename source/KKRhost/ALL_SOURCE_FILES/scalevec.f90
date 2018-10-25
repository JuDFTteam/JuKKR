!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Transforms all the basis positions into the cartesian reference system
!> Author: 
!> Transforms all the basis positions into the cartesian reference system
!------------------------------------------------------------------------------------
module mod_scalevec
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Transforms all the basis positions into the cartesian reference system
  !> Author: 
  !> Category: geometry, KKRhost
  !> Deprecated: False
  !> Transforms all the basis positions into the cartesian reference system
  !-------------------------------------------------------------------------------
  subroutine scalevec(lcartesian,rbasis,abasis,bbasis,cbasis,nlbasis,nrbasis,nleft, &
    nright,zperleft,zperight,tleft,tright,linterface,naez,nemb,bravais,kaoez,noq,   &
    naezd,natypd,nembd)

    implicit none
    real (kind=dp), parameter :: eps = 1e-14_dp
    ! ..
    ! .. Arguments ..
    integer :: naezd, natypd, nembd
    integer :: naez, nemb, nlbasis, nleft, nrbasis, nright
    real (kind=dp) :: abasis, bbasis, cbasis
    logical :: linterface
    integer :: kaoez(natypd, *), noq(naezd)
    real (kind=dp) :: bravais(3, 3), rbasis(3, *), tleft(3, *), tright(3, *), zperight(3), zperleft(3)
    ! ..
    ! .. Locals ..
    integer :: i, i1, j
    logical :: lcartesian
    real (kind=dp) :: rbasis1(3, naezd+nembd), temp(3), tx, ty, tz

    write (1337, '(79("="))')
    write (1337, '(23X,A)') 'SCALEVEC: scale site coordinates'
    write (1337, '(23X,A)') '          bring all to CARTESIAN system'
    write (1337, '(79("="))')
    write (1337, *)

    ! -->   normalization of basis vectors
    ! multiplication instead of division 04/2004

    do i = 1, naez + nemb
      rbasis1(1, i) = rbasis(1, i)*abasis
      rbasis1(2, i) = rbasis(2, i)*bbasis
      rbasis1(3, i) = rbasis(3, i)*cbasis
    end do

    if (linterface) then
      do i = 1, nlbasis
        tleft(1, i) = tleft(1, i)*abasis
        tleft(2, i) = tleft(2, i)*bbasis
        tleft(3, i) = tleft(3, i)*cbasis
      end do
      zperleft(1) = zperleft(1)*abasis
      zperleft(2) = zperleft(2)*bbasis
      zperleft(3) = zperleft(3)*cbasis

      do i = 1, nrbasis
        tright(1, i) = tright(1, i)*abasis
        tright(2, i) = tright(2, i)*bbasis
        tright(3, i) = tright(3, i)*cbasis
      end do
      zperight(1) = zperight(1)*abasis
      zperight(2) = zperight(2)*bbasis
      zperight(3) = zperight(3)*cbasis
    end if

    if (abs(abasis)<eps .or. abs(bbasis)<eps .or. abs(cbasis)<eps) then
      write (1337, '(5X,A,2(/,34X,F12.8,A))') 'Scaling site coordinates with:', abasis, '  x', bbasis, '  y'
      if (.not. linterface) write (1337, '(34X,F12.8,A)') cbasis, '  z'
      write (1337, '(5X,44("-"))')
    else
      write (1337, '(5X,A)') 'Site coordinates will not be scaled'
    end if

    ! ---> normalization of atomic positions in the unit cell

    ! if lcartesian is true cartesian coordinates are used
    ! else the basis atoms are in units of the lattice vectors

    if (lcartesian) then
      write (1337, '(A)') ' CARTESIAN coordinates'
    else
      write (1337, '(A)') ' LATTICE VECTOR coordinates will be', ' changed to CARTESIAN coordinates'
    end if

    ! **********************************************************************
    ! Change to cartesian coordinates
    if (linterface) then
      ! ======================================================================
      if (.not. lcartesian) then
        ! ----------------------------------------------------------------------
        write (1337, *)
        write (1337, '(12X,49("-"))')
        write (1337, '(13X,A)') 'Input positions transformed to CARTESIAN system'
        write (1337, '(12X,49("-"),/,13X,A,/,12X,49("-"))') 'IQ        x             y             z        IT'
        do i = 1, naez + nemb
          do j = 1, 2
            rbasis(j, i) = (rbasis1(1,i)*bravais(j,1)+rbasis1(2,i)*bravais(j,2))
          end do
          rbasis(3, i) = rbasis1(3, i)

          if (i<=naez) then
            write (1337, 140) i, (rbasis(j,i), j=1, 3), (kaoez(j,i), j=1, noq(i))
          else
            write (1337, 140) i, (rbasis(j,i), j=1, 3), kaoez(1, i)
          end if
          if (i==naez) write (1337, '(12X,49("."))')
        end do
        write (1337, '(12X,49("-"),/)')
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! -->  Do the same for the boundary vectors

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! left side

        do i = 1, nlbasis
          do i1 = 1, 2
            temp(i1) = tleft(i1, i)
          end do
          do j = 1, 2
            tleft(j, i) = (temp(1)*bravais(j,1)+temp(2)*bravais(j,2))
          end do
        end do

        do i1 = 1, 2
          temp(i1) = zperleft(i1)
        end do
        do j = 1, 2
          zperleft(j) = (temp(1)*bravais(j,1)+temp(2)*bravais(j,2))
        end do
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! right side

        do i = 1, nrbasis
          do i1 = 1, 2
            temp(i1) = tright(i1, i)
          end do
          do j = 1, 2
            tright(j, i) = (temp(1)*bravais(j,1)+temp(2)*bravais(j,2))
          end do
        end do

        do i1 = 1, 2
          temp(i1) = zperight(i1)
        end do
        do j = 1, 2
          zperight(j) = (temp(1)*bravais(j,1)+temp(2)*bravais(j,2))
        end do
        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! ----------------------------------------------------------------------
      else
        ! ----------------------------------------------------------------------
        write (1337, '(42X,A)') '---> No transformation required'
        do i = 1, 3
          do j = 1, naez + nemb
            rbasis(i, j) = rbasis1(i, j)
          end do
        end do
        ! ----------------------------------------------------------------------
      end if                       ! IF (.NOT.LCARTESIAN)
      ! ======================================================================

      write (1337, 110)
      do i = nleft, 1, -1
        do i1 = nlbasis, 1, -1
          tx = tleft(1, i1) + (i-1)*zperleft(1)
          ty = tleft(2, i1) + (i-1)*zperleft(2)
          tz = tleft(3, i1) + (i-1)*zperleft(3)
          write (1337, 100)(i-1)*nlbasis + i1, tx, ty, tz, kaoez(1, naez+i1)
        end do
      end do

      write (1337, 120)
      do i = 1, naez
        write (1337, 100) i, (rbasis(i1,i), i1=1, 3), (kaoez(i1,i), i1=1, noq(i))
      end do

      write (1337, 130)
      do i = 1, nright
        do i1 = 1, nrbasis
          tx = tright(1, i1) + (i-1)*zperight(1)
          ty = tright(2, i1) + (i-1)*zperight(2)
          tz = tright(3, i1) + (i-1)*zperight(3)
          write (1337, 100)(i-1)*nrbasis + i1, tx, ty, tz, kaoez(1, naez+nlbasis+i1)
        end do
      end do
      write (1337, '(14X,45("-"),/)')
      ! ======================================================================
    else if (.not. lcartesian) then ! Rescale lattice
      ! ----------------------------------------------------------------------
      write (1337, *)
      write (1337, '(12X,49("-"))')
      write (1337, '(13X,A)') 'Input positions transformed to CARTESIAN system'
      write (1337, '(12X,49("-"),/,13X,A,/,12X,49("-"))') 'IQ        x             y             z        IT'
      do i = 1, naez + nemb
        do j = 1, 3
          rbasis(j, i) = (rbasis1(1,i)*bravais(j,1)+rbasis1(2,i)*bravais(j,2)+rbasis1(3,i)*bravais(j,3))
        end do

        if (i<=naez) then
          write (1337, 140) i, (rbasis(j,i), j=1, 3), (kaoez(j,i), j=1, noq(i))
        else
          write (1337, 140) i, (rbasis(j,i), j=1, 3), kaoez(1, i)
        end if
        if (i==naez .and. nemb>0) write (1337, '(12X,49("."))')
      end do
      write (1337, '(12X,49("-"),/)')
      ! ----------------------------------------------------------------------
    else
      ! ----------------------------------------------------------------------
      write (1337, '(42X,A,/)') '---> No transformation required'
      write (1337, 150)
      ! changed by v.Bellini 21/10/99
      do j = 1, naez + nemb
        do i = 1, 3
          rbasis(i, j) = rbasis1(i, j)
        end do

        if (j<=naez) then
          write (1337, 100) j, (rbasis(i,j), i=1, 3), (kaoez(i,j), i=1, noq(j))
        else
          write (1337, 100) j, (rbasis(i,j), i=1, 3), kaoez(1, j)
        end if
        if (i==naez .and. nemb>0) write (1337, '(12X,51("."))')
      end do
      ! end of the change
      write (1337, '(12X,51("-"),/)')
      ! ----------------------------------------------------------------------
      ! ======================================================================
    end if                         ! IF (.NOT.LINTERFACE )
    ! **********************************************************************

    ! FROM NOW ON after < SCALEVEC > RBASIS are the basis vectors
    ! in units of au/alat in (xyz) reference

    ! **********************************************************************
100 format (13x, i5, 3f12.6, 10i3)
110 format (14x, 45('-'), /, 15x, '     Positions of ALL generated sites ', /, 15x, '   in CARTESIAN coordinates (ALAT units)', /, 14x, 45('-'), /, 15x, &
      'IQ       x           y           z       IT', /, 15x, '**************** Left  Host ***************')
120 format (15x, '****************   S L A B  ***************')
130 format (15x, '**************** Right Host ***************')
140 format (12x, i3, 3f14.8, 10i3)
150 format (12x, 51('-'), /, 16x, '    Positions of (ALL) generated sites', /, 16x, '   in CARTESIAN coordinates (ALAT units)', /, 12x, 51('-'), /, 15x, &
      'IQ       x           y           z       IT', /, 12x, 51('-'))
  end subroutine scalevec

end module mod_scalevec
