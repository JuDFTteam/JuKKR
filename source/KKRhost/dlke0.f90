!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_dlke0

contains

  !-------------------------------------------------------------------------------
  !> Summary: Driver for lattice fourier transform
  !> Author: 
  !> Category: KKRhost, k-points, structural-greensfunction
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Function, set up in the spin-independent non-relativstic
  !> (l,m_l)-representation
  !-------------------------------------------------------------------------------
  subroutine dlke0(gllke, alat, naez, cls, nacls, naclsmax, rr, ezoa, atom, bzkp, rcls, ginp)
use :: mod_runoptions, only: calc_complex_bandstructure --manopt-- 

    use :: global_variables, only: lmgf0d, almgf0, naclsd, nrd
    use :: mod_datatypes, only: dp
    use :: mod_dlke1, only: dlke1
    use :: mod_cinit, only: cinit
    implicit none

    real (kind=dp) :: alat
    integer :: naez, naclsmax
    complex (kind=dp) :: ginp(lmgf0d*naclsmax, lmgf0d, *), gllke(almgf0, *)
    real (kind=dp) :: bzkp(*), rcls(3, naclsd, *), rr(3, 0:nrd)
    integer :: atom(naclsd, *), cls(*), ezoa(naclsd, *), nacls(*)
    ! ..
    ! .. External Subroutines ..
    integer :: i, ic, im, j, jn, m, n
    ! ..
    ! .. Save statement ..
    complex (kind=dp) :: gllke1(almgf0, lmgf0d)
    real (kind=dp) :: kp(6)

    logical, external :: opt

    ! write(6,*) '>>> DLKE0 : Fourier-transforms the ',
    ! +           'GF of reference system'

    call cinit(almgf0*almgf0, gllke(1,1))

    do i = 1, naez


      kp(1) = bzkp(1)
      kp(2) = bzkp(2)
      kp(3) = bzkp(3)
      if (calc_complex_bandstructure) then
        kp(4) = bzkp(4)
        kp(5) = bzkp(5)
        kp(6) = bzkp(6)
      end if

      ic = cls(i)
      call dlke1(gllke1, alat, nacls, naclsmax, rr, ezoa(1,i), atom(1,i), kp, ic, ginp(1,1,ic), rcls(1,1,ic))

      do m = 1, lmgf0d
        im = (i-1)*lmgf0d + m
        do jn = 1, lmgf0d*naez
          gllke(jn, im) = gllke(jn, im) + gllke1(jn, m)
        end do
      end do
      ! ----------------------------------------------------------------------

    end do


    if (opt('symG(k) ')) then

      ! -->   symmetrization

      do i = 1, naez

        kp(1) = -bzkp(1)
        kp(2) = -bzkp(2)
        kp(3) = -bzkp(3)
        if (calc_complex_bandstructure) then
          kp(4) = -bzkp(4)
          kp(5) = -bzkp(5)
          kp(6) = -bzkp(6)
        end if
        ! ----------------------------------------------------------------------
        ic = cls(i)
        call dlke1(gllke1, alat, nacls, naclsmax, rr, ezoa(1,i), atom(1,i), kp, ic, ginp(1,1,ic), rcls(1,1,ic))

        do j = 1, naez
          do m = 1, lmgf0d
            im = (i-1)*lmgf0d + m
            do n = 1, lmgf0d
              jn = (j-1)*lmgf0d + n
              gllke(im, jn) = (gllke(im,jn)+gllke1(jn,m))/2.0e0_dp
            end do
          end do
        end do

      end do

    end if

    return

  end subroutine dlke0

end module mod_dlke0
