module mod_intcheb_cell
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine intcheb_cell(cden, den, rpan_intervall, ipan_intervall, npan_tot, ncheb, irmdnew)
    ! ***********************************************************************
    ! integrate the complex density of states for LM=1
    ! gives the total complex charge which is then
    ! transformed to the xyz component of the magnetic
    ! moment
    ! ***********************************************************************
    implicit none

    integer :: ncheb, npan_tot, irmdnew
    integer :: ipan_intervall(0:npan_tot)
    real (kind=dp) :: rpan_intervall(0:npan_tot)
    complex (kind=dp) :: cden(irmdnew), den
    integer :: irstart, irstop, ipan
    real (kind=dp) :: widthfac
    complex (kind=dp) :: int1

    den = (0.0e0_dp, 0.0e0_dp)

    do ipan = 1, npan_tot
      irstart = ipan_intervall(ipan-1) + 1
      irstop = ipan_intervall(ipan)
      widthfac = 0.5e0_dp*(rpan_intervall(ipan)-rpan_intervall(ipan-1))
      call intcheb_complex(ncheb, cden(irstart:irstop), int1)
      den = den + int1*widthfac
    end do

  end subroutine intcheb_cell

  subroutine intcheb_complex(ncheb, arr1, result1)
    implicit none
    integer, intent (in) :: ncheb
    complex (kind=dp), intent (in) :: arr1(0:ncheb)
    complex (kind=dp), intent (out) :: result1
    real (kind=dp) :: pi
    real (kind=dp) :: intweight(0:ncheb)
    integer :: icheb1, icheb2

    pi = 4e0_dp*atan(1e0_dp)
    intweight = 1.0e0_dp
    do icheb1 = 0, ncheb
      do icheb2 = 2, ncheb, 2
        intweight(icheb1) = intweight(icheb1) + (-2.0e0_dp/(icheb2**2-1.0e0_dp))*cos(icheb2*pi*(icheb1+0.5e0_dp)/(ncheb+1))
      end do
      intweight(icheb1) = intweight(icheb1)*2.0e0_dp/(ncheb+1)
    end do

    result1 = (0.0e0_dp, 0.0e0_dp)
    do icheb1 = 0, ncheb
      result1 = result1 + intweight(icheb1)*arr1(icheb1)
    end do

  end subroutine intcheb_complex

end module mod_intcheb_cell
