!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Routine for the mixing of the potential matrix in the case of LDA+U
!> Author:
!> Routine for the mixing of the potential matrix in the case of LDA+U 
!------------------------------------------------------------------------------------
module mod_mixldau
  use :: mod_datatypes, only: dp
  private :: dp

contains
   !-------------------------------------------------------------------------------  
   !> Summary: Routine for the mixing of the potential matrix in the case of LDA+U
   !> Author: 
   !> Category: potential, mixing, lda+u, KKRhost 
   !> Deprecated: False 
   !> Routine for the mixing of the potential matrix in the case of LDA+U as well as 
   !> the rms error in the interaction matrix
   !-------------------------------------------------------------------------------  
  subroutine mixldau(mmaxd, nspind, natypd, natyp, nspin, lopt, wldauold, wldau)
    use :: mod_ioinput
    implicit none
    ! Input:
    integer, intent(in) :: mmaxd   !! 2*lmax + 1
    integer, intent(in) :: nspind  !! krel + (1-krel)*2
    integer, intent(in) :: natypd  !! Number of kinds of atoms in unit cell
    integer, intent(in) :: nspin   !! Counter for spin directions 
    integer, intent(in) :: natyp   !! Number of kinds of atoms in unit cell
    integer, dimension(natypd), intent(in) :: lopt !! angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
    real (kind=dp), dimension(mmaxd, mmaxd, nspind, natypd), intent(in) :: wldauold
    ! Input/Output:
    real (kind=dp), dimension(mmaxd, mmaxd, nspind, natypd), intent(inout) :: wldau !! potential matrix
    ! Inside:
    integer :: iat, is, m1, m2, mmax
    integer :: ier
    real (kind=dp) :: xmix, xmix2, rmserr
    character (len=:), allocatable :: uio                             ! NCOLIO=256

    ! First calculate rms error in interaction matrix
    do iat = 1, natyp
      rmserr = 0.e0_dp
      if (lopt(iat)>=0) then
        mmax = 2*lopt(iat) + 1
        do is = 1, nspin
          do m2 = 1, mmax
            do m1 = 1, mmax
              rmserr = rmserr + (wldau(m1,m2,is,iat)-wldauold(m1,m2,is,iat))**2
            end do
          end do
        end do
        rmserr = sqrt(rmserr)
        write (1337, 100) iat, rmserr
100     format ('LDA+U interaction matrix rms error for atom', i6, ' = ', e10.2)
      end if
    end do

    ! Now mix old/new interaction matrices
    ier = 0
    call ioinput('MIXLDAU         ', uio, 1, 7, ier)
    if (ier/=0) then
      write (*, *) 'MIXLDAU not found, setting to 1.'
      return
    else
      read (unit=uio, fmt=*) xmix
      write (1337, *) 'Using MIXLDAU = ', xmix
    end if

    xmix2 = 1.e0_dp - xmix

    do iat = 1, natyp
      if (lopt(iat)>=0) then
        do is = 1, nspin
          do m2 = 1, mmaxd
            do m1 = 1, mmaxd
              wldau(m1, m2, is, iat) = xmix*wldau(m1, m2, is, iat) + xmix2*wldauold(m1, m2, is, iat)
            end do
          end do
        end do
      end if
    end do

  end subroutine mixldau

end module mod_mixldau
