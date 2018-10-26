!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_intcheb_cell
  
  private
  public :: intcheb_cell

contains

  !-------------------------------------------------------------------------------
  !> Summary: Panel-wise Chebychev integration of complex density
  !> Author: 
  !> Category: KKRhost, radial-grid, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Integrate the complex density of states for LM=1
  !> gives the total complex charge which is then
  !> transformed to the xyz component of the magnetic
  !> moment
  !-------------------------------------------------------------------------------
  subroutine intcheb_cell(cden, den, rpan_intervall, ipan_intervall, npan_tot, ncheb, irmdnew)

    use :: mod_datatypes, only: dp
    use :: mod_constants, only: czero
    implicit none

    ! arguments
    integer, intent(in) :: ncheb !! number of mesh points in Chebychev mesh
    integer, intent(in) :: npan_tot !! total number of panels (each having `ncheb` radial points)
    integer, intent(in) :: irmdnew !! total number of radial points
    integer, intent(in) :: ipan_intervall(0:npan_tot) !! indices of radial points for boundaries of panels
    real (kind=dp), intent(in) :: rpan_intervall(0:npan_tot) !! radial values at panel boundaries
    complex (kind=dp), intent(in)  :: cden(irmdnew) !! complex input density
    complex (kind=dp), intent(out)  :: den !! integrated density
    ! locals
    integer :: irstart, irstop, ipan
    real (kind=dp) :: widthfac
    complex (kind=dp) :: int1

    den = czero

    do ipan = 1, npan_tot
      irstart = ipan_intervall(ipan-1) + 1
      irstop = ipan_intervall(ipan)
      widthfac = 0.5e0_dp*(rpan_intervall(ipan)-rpan_intervall(ipan-1))
      call intcheb_complex(ncheb, cden(irstart:irstop), int1)
      den = den + int1*widthfac
    end do

  end subroutine intcheb_cell


  !-------------------------------------------------------------------------------
  !> Summary: Chebychev integration of complex array
  !> Author: 
  !> Category: KKRhost, radial-grid, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine intcheb_complex(ncheb, arr1, result1)

    use :: mod_datatypes, only: dp
    use :: mod_constants, only: pi, czero
    implicit none

    ! arguments
    integer, intent (in) :: ncheb !! number of radial points in Chebychev mesh 
    complex (kind=dp), intent (in) :: arr1(0:ncheb) !! values to be integrated
    complex (kind=dp), intent (out) :: result1 !! result of integral
    ! local
    real (kind=dp) :: intweight(0:ncheb)
    integer :: icheb1, icheb2

    intweight = 1.0e0_dp
    do icheb1 = 0, ncheb
      do icheb2 = 2, ncheb, 2
        intweight(icheb1) = intweight(icheb1) + (-2.0e0_dp/(icheb2**2-1.0e0_dp))*cos(icheb2*pi*(icheb1+0.5e0_dp)/(ncheb+1))
      end do
      intweight(icheb1) = intweight(icheb1)*2.0e0_dp/(ncheb+1)
    end do

    result1 = czero
    do icheb1 = 0, ncheb
      result1 = result1 + intweight(icheb1)*arr1(icheb1)
    end do

  end subroutine intcheb_complex

end module mod_intcheb_cell
