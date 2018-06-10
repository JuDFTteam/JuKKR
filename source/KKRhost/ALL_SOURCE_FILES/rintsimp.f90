subroutine rintsimp(fx, jtop, cint)
!   ********************************************************************
!   *                                                                  *
!   *  SIMPSON - INTERGRATION FOR  REAL   INTEGRAND  FX FROM 1 TO JTOP *
!   *  AND EQUIDISTANT MESH    I                                       *
!   *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
!   *                                                                  *
!   ********************************************************************

  implicit none

! Dummy arguments
  real *8 :: cint
  integer :: jtop
  real *8 :: fx(jtop)

! Local variables
  integer :: i
  real *8 :: simp

  cint = fx(1)
  simp = -1.0d0

  do i = 2, jtop - 1
    simp = -simp
    cint = cint + (3.0d0+simp)*fx(i)
  end do

  cint = (cint+fx(jtop))/3.0d0

end subroutine
