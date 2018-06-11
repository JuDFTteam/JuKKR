    Subroutine rintsimp(fx, jtop, cint)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *  SIMPSON - INTERGRATION FOR  REAL   INTEGRAND  FX FROM 1 TO JTOP *
!   *  AND EQUIDISTANT MESH    I                                       *
!   *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
!   *                                                                  *
!   ********************************************************************

      Implicit None

! Dummy arguments
      Real (Kind=dp) :: cint
      Integer :: jtop
      Real (Kind=dp) :: fx(jtop)

! Local variables
      Integer :: i
      Real (Kind=dp) :: simp

      cint = fx(1)
      simp = -1.0E0_dp

      Do i = 2, jtop - 1
        simp = -simp
        cint = cint + (3.0E0_dp+simp)*fx(i)
      End Do

      cint = (cint+fx(jtop))/3.0E0_dp

    End Subroutine
