!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_mathtools

  implicit none

  private
  public :: crossprod, bubblesort, bubblesort_int, findminindex, simple_integration, simple_integration_general, simpson2D_integration
  !obsolete functions: machpi, mach_tpiimag

  double precision, parameter, public :: pi=4d0*atan(1d0), tpi=8d0*atan(1d0)
  double complex,   parameter, public :: tpiimag=(0d0, tpi)
! logical, save :: piinit=.true., tpiimag_init=.true.

contains

  subroutine findminindex(nin, Xin, imin)

    implicit none
    integer,          intent(in)  :: nin
    double precision, intent(in)  :: Xin(nin)
    integer,          intent(out) :: imin

    integer :: ii
    double precision :: dmin

    dmin=1d38
    imin=0
    do ii=1,nin
      if(Xin(ii)<dmin)then
        imin=ii
        dmin = Xin(ii)
      end if
    end do
    if(imin<1 .or. imin>nin) stop 'error in findminindex'

  end subroutine findminindex





  !-------------------------------------------------------------------------------
  !> Summary: Sorts an array in ascending order
  !> Author: 
  !> Category: PKKprime, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine bubblesort(n, Xin, Iout)

    implicit none
    integer,          intent(in)  :: n
    double precision, intent(in)  :: Xin(n)
    integer,          intent(out) :: Iout(n)

    logical          :: swaped
    integer          :: ii, ntmp, itmp

    swaped = .true.
    ntmp   = n

    do ii=1,n
      Iout(ii) = ii
    end do

    do while ( swaped .and. ntmp>1 )
      swaped = .false.
      do ii=1,ntmp-1
        if( Xin(Iout(ii))>Xin(Iout(ii+1)) ) then
          itmp       = Iout(ii)
          Iout(ii)   = Iout(ii+1)
          Iout(ii+1) = itmp
          swaped = .true.
        end if
      end do!ii
      ntmp = ntmp-1
    end do!while

  end subroutine bubblesort





  subroutine bubblesort_int(n, Iin, Iout)

    implicit none
    integer, intent(in)  :: n
    integer, intent(in)  :: Iin(n)
    integer, intent(out) :: Iout(n)

    logical          :: swaped
    integer          :: ii, ntmp, itmp

    swaped = .true.
    ntmp   = n

    do ii=1,n
      Iout(ii) = ii
    end do

    do while ( swaped .and. ntmp>1 )
      swaped = .false.
      do ii=1,ntmp-1
        if( Iin(Iout(ii))>Iin(Iout(ii+1)) ) then
          itmp       = Iout(ii)
          Iout(ii)   = Iout(ii+1)
          Iout(ii+1) = itmp
          swaped = .true.
        end if
      end do!ii
      ntmp = ntmp-1
    end do!while

  end subroutine bubblesort_int

  subroutine crossprod(a,b,c)
  ! calculates the cross product of two vectors, c = a x b
    implicit none

    double precision, intent(in)  :: a(3), b(3)
    double precision, intent(out) :: c(3)

    c(:) = 0d0
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end subroutine crossprod


!===========================================================!
!===========================================================!
! the functions 'machpi' and 'mach_tpiimag' have been made  !
! obsolete through defining 'pi' and 'tpiimag' as parameter !
!===========================================================!
!===========================================================!
 
! double precision function machpi()
!   implicit none

!   if(piinit) then
!     pi=4d0*atan(1d0)
!     piinit=.false.
!   end if

!   machpi = pi

! end function machpi


! double complex function mach_tpiimag()
!   implicit none
!   double precision :: pitmp
!   double complex, parameter :: CI=(0d0,1d0)

!   if(tpiimag_init) then
!     pitmp = machpi()
!     tpiimag=2d0*pitmp*CI
!     tpiimag_init=.false.
!   end if

!   mach_tpiimag = tpiimag

! end function mach_tpiimag

  subroutine simple_integration_general(nBZdim, area, fermi_velocity, function_value, integratedvalue, densityofstates)
  implicit none

      integer,          intent(in)  :: nBZdim
      double precision, intent(in)  :: area, fermi_velocity(3,nBZdim), function_value(nBZdim)
      double precision, intent(out) :: integratedvalue, densityofstates

      double precision :: fermi_velocity_mean(3), function_value_mean, fermi_velocity_mean_abs

      integratedvalue = 0d0
      densityofstates = 0d0

      if(nBZdim==3)then
        fermi_velocity_mean = ( fermi_velocity(:,1) + fermi_velocity(:,2) + fermi_velocity(:,3)  ) /3d0
      elseif(nBZdim==2)then
        fermi_velocity_mean = ( fermi_velocity(:,1) + fermi_velocity(:,2))/2d0
      else!nBZdim
      end if!nBZdim
      function_value_mean = sum(function_value)/nBZdim

      fermi_velocity_mean_abs = sqrt(sum(fermi_velocity_mean**2))

      integratedvalue = function_value_mean/fermi_velocity_mean_abs*area
      densityofstates = 1d0/fermi_velocity_mean_abs * area

  end subroutine simple_integration_general

  subroutine simpson2D_integration(d, fermi_velocity, function_value, integratedvalue, densityofstates)
  implicit none
      double precision, intent(in)  :: d(2), fermi_velocity(3,3), function_value(3)
      double precision, intent(out) :: integratedvalue, densityofstates

      double precision :: fermi_velocity_abs(3)

      integratedvalue = 0d0
      densityofstates = 0d0

      fermi_velocity_abs(1) = sqrt(sum(fermi_velocity(:,1)**2))
      fermi_velocity_abs(2) = sqrt(sum(fermi_velocity(:,2)**2))
      fermi_velocity_abs(3) = sqrt(sum(fermi_velocity(:,3)**2))

      integratedvalue  =   function_value(1)/fermi_velocity_abs(1) * ( 1d0/3d0*d(1)  + 1/6d0*d(2)            - 1d0/6d0*d(2)**2/d(1)) &
                         + function_value(2)/fermi_velocity_abs(2) * ( 1d0/6d0*(d(1) + d(2))**3/(d(1)*d(2)))                         &
                         + function_value(3)/fermi_velocity_abs(3) * ( 1d0/6d0*d(1)  + 1/3d0*d(2)            - 1d0/6d0*d(1)**2/d(2))

      densityofstates  =   1/fermi_velocity_abs(1) * ( 1d0/3d0*d(1)  + 1d0/6d0*d(2)            - 1d0/6d0*d(2)**2/d(1)) &
                         + 1/fermi_velocity_abs(2) * ( 1d0/6d0*(d(1) + d(2))**3/(d(1)*d(2)))                           &
                         + 1/fermi_velocity_abs(3) * ( 1d0/6d0*d(1)  + 1d0/3d0*d(2)            - 1d0/6d0*d(1)**2/d(2))
  end subroutine


  subroutine simple_integration(area, fermi_velocity, function_value, integratedvalue, densityofstates)
  implicit none

      double precision, intent(in)  :: area, fermi_velocity(3,3), function_value(3)
      double precision, intent(out) :: integratedvalue, densityofstates

      double precision :: fermi_velocity_mean(3), function_value_mean, fermi_velocity_mean_abs

      integratedvalue = 0d0
      densityofstates = 0d0

      fermi_velocity_mean = (  fermi_velocity(:,1) + fermi_velocity(:,2) + fermi_velocity(:,3)  ) /3d0
      function_value_mean = (  function_value(1)   + function_value(2)   + function_value(3)    ) /3d0

      fermi_velocity_mean_abs = sqrt(sum(fermi_velocity_mean**2))

      integratedvalue = function_value_mean/fermi_velocity_mean_abs*area
      densityofstates = 1d0/fermi_velocity_mean_abs * area

  end subroutine simple_integration


end module mod_mathtools
