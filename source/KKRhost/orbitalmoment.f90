!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper module for the calculation of the orbital moment 
!> Author:
!> Wrapper module for the calculation of the orbital moment 
!------------------------------------------------------------------------------------
module mod_orbitalmoment

contains

!-------------------------------------------------------------------------------
  !> Summary: Calculation of the orbital moment
  !> Author: 
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Calculation of the orbital moment 
  !-------------------------------------------------------------------------------
  subroutine calc_orbitalmoment(lmax, lmmaxd, loperator)

    use :: mod_constants, only: czero
    use :: mod_profiling
    use :: mod_datatypes, only: dp

    implicit none
    integer, intent (in) :: lmax    !! Maximum l component in wave function expansion
    integer, intent (in) :: lmmaxd  !! (KREL+KORBIT+1)*(LMAX+1)**2
    complex (kind=dp), dimension (lmmaxd, lmmaxd, 3), intent (out) :: loperator
    integer :: lval
    complex (kind=dp), dimension (:, :, :), allocatable :: lorbit_onel
    integer :: lmmax0d, lstart, lstop, i_stat, i_all

    lmmax0d = (lmax+1)**2
    loperator = czero

    loperator(1, 1, 1) = czero
    loperator(1, 1, 2) = czero
    loperator(1, 1, 3) = czero

    do lval = 1, lmax
      allocate (lorbit_onel(2*lval+1,2*lval+1,3), stat=i_stat)
      call memocc(i_stat, product(shape(lorbit_onel))*kind(lorbit_onel), 'lorbit_onel', 'calc_orbitalmoment')
      lorbit_onel = czero
      call calc_orbit_onel(lval, lorbit_onel)
      lstart = ((lval-1)+1)**2 + 1
      lstop = ((lval)+1)**2
      loperator(lstart:lstop, lstart:lstop, 1) = lorbit_onel(:, :, 1)
      loperator(lstart:lstop, lstart:lstop, 2) = lorbit_onel(:, :, 2)
      loperator(lstart:lstop, lstart:lstop, 3) = lorbit_onel(:, :, 3)

      i_all = -product(shape(lorbit_onel))*kind(lorbit_onel)
      deallocate (lorbit_onel, stat=i_stat)
      call memocc(i_stat, i_all, 'lorbit_onel', 'calc_orbitalmoment')
    end do
    if (lmmaxd/=lmmax0d) then
      loperator(lmmax0d+1:lmmaxd, lmmax0d+1:lmmaxd, :) = loperator(:lmmax0d, :lmmax0d, :)
    end if

  end subroutine calc_orbitalmoment

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of the contibution of a given orbital to the orbital moment 
  !> Author: 
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Calculation of the contibution of a given orbital to the orbital moment 
  !-------------------------------------------------------------------------------
  subroutine calc_orbit_onel(lval, lorbit_onel)

    use :: mod_constants
    use :: mod_datatypes, only: dp

    implicit none
    integer, intent (in) :: lval
    complex (kind=dp), dimension (2*lval+1, 2*lval+1, 3), intent (out) :: lorbit_onel
    complex (kind=dp), dimension (-lval:lval, -lval:lval) :: l_z
    complex (kind=dp), dimension (-lval:lval, -lval:lval) :: l_x
    complex (kind=dp), dimension (-lval:lval, -lval:lval) :: l_y
    complex (kind=dp), dimension (-lval:lval, -lval:lval) :: l_dn
    complex (kind=dp), dimension (-lval:lval, -lval:lval) :: l_up
    integer :: i1
    real (kind=dp) :: lfac

    l_z = czero
    l_x = czero
    l_y = czero
    l_dn = czero
    l_up = czero

    lorbit_onel = czero
    do i1 = -lval, lval
      l_z(-i1, i1) = -ci*i1
    end do

    if (lval>0) then

      lfac = sqrt(lval*(lval+1e0_dp))/sqrt(2e0_dp)
      l_dn(0, -1) = -ci*lfac
      l_dn(0, 1) = lfac
      l_dn(-1, 0) = ci*lfac
      l_dn(1, 0) = -lfac

      if (lval>1) then
        do i1 = 2, lval
          lfac = 0.5e0_dp*sqrt(lval*(lval+1e0_dp)-i1*(i1-1e0_dp))
          l_dn(-i1, -i1+1) = -lfac
          l_dn(-i1, i1-1) = ci*lfac
          l_dn(i1, -i1+1) = -ci*lfac
          l_dn(i1, i1-1) = -lfac

          lfac = 0.5e0_dp*sqrt(lval*(lval+1e0_dp)-(i1-1e0_dp)*i1)
          l_dn(-i1+1, -i1) = lfac
          l_dn(-i1+1, i1) = ci*lfac
          l_dn(i1-1, -i1) = -ci*lfac
          l_dn(i1-1, i1) = lfac
        end do
      end if
    end if

    if (lval>0) then

      lfac = sqrt(lval*(lval+1e0_dp))/sqrt(2e0_dp)
      l_up(0, -1) = -ci*lfac
      l_up(0, 1) = -lfac
      l_up(-1, 0) = ci*lfac
      l_up(1, 0) = lfac

      if (lval>1) then
        do i1 = 2, lval
          lfac = 0.5e0_dp*sqrt(lval*(lval+1e0_dp)-i1*(i1-1e0_dp))
          l_up(-i1, -i1+1) = lfac
          l_up(-i1, i1-1) = ci*lfac
          l_up(i1, -i1+1) = -ci*lfac
          l_up(i1, i1-1) = lfac

          lfac = 0.5e0_dp*sqrt(lval*(lval+1e0_dp)-(i1-1e0_dp)*i1)
          l_up(-i1+1, -i1) = -lfac
          l_up(-i1+1, i1) = ci*lfac
          l_up(i1-1, -i1) = -ci*lfac
          l_up(i1-1, i1) = -lfac
        end do
      end if
    end if

    l_x = 0.5e0_dp*(l_up+l_dn)
    l_y = -0.5e0_dp*ci*(l_up-l_dn)

    lorbit_onel(:, :, 1) = l_x
    lorbit_onel(:, :, 2) = l_y
    lorbit_onel(:, :, 3) = l_z

  end subroutine calc_orbit_onel

end module mod_orbitalmoment
