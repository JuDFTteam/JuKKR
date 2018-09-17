module mod_rintsimp
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine rintsimp(fx, jtop, cint)
    ! ********************************************************************
    ! *                                                                  *
    ! *  SIMPSON - INTERGRATION FOR  REAL   INTEGRAND  FX FROM 1 TO JTOP *
    ! *  AND EQUIDISTANT MESH    I                                       *
    ! *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
    ! *                                                                  *
    ! ********************************************************************

    implicit none

    ! Dummy arguments
    real (kind=dp) :: cint
    integer :: jtop
    real (kind=dp) :: fx(jtop)

    ! Local variables
    integer :: i
    real (kind=dp) :: simp

    cint = fx(1)
    simp = -1.0e0_dp

    do i = 2, jtop - 1
      simp = -simp
      cint = cint + (3.0e0_dp+simp)*fx(i)
    end do

    cint = (cint+fx(jtop))/3.0e0_dp

  end subroutine rintsimp

end module mod_rintsimp
