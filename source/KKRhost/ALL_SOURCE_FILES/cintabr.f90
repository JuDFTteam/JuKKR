    Subroutine cintabr(ag, bg, agbg, af, bf, afbf, rpw, nka, nkb, jtop, nrmax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *  SIMPSON - INTERGRATION FOR COMPLEX INTEGRAND  FX FROM 1 TO JTOP *
!   *  AND EQUIDISTANT MESH    I                                       *
!   *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
!   *                                                                  *
!   *            FX = AG*AG*RPW   and   FX = AF*AF*RPW                 *
!   *                                                                  *
!   ********************************************************************

      Implicit None


!Dummy arguments
      Integer :: jtop, nka, nkb, nrmax
      Complex (Kind=dp) :: af(nrmax, 2), afbf(2, 2), ag(nrmax, 2), agbg(2, 2), &
        bf(nrmax, 2), bg(nrmax, 2)
      Real (Kind=dp) :: rpw(nrmax)

!Local variables
      Real (Kind=dp) :: f, simp
      Integer :: i, ka, kb

      Do kb = 1, nkb
        Do ka = 1, nka
          agbg(ka, kb) = ag(1, ka)*bg(1, kb)*rpw(1)
          afbf(ka, kb) = af(1, ka)*bf(1, kb)*rpw(1)
        End Do
      End Do

      If (mod(jtop,2)==0) Stop '<CINTABR>  JTOP is even !!!'

      simp = -1.0E0_dp
      Do i = 2, jtop - 1
        simp = -simp
        f = (3.0E0_dp+simp)*rpw(i)
        Do kb = 1, nkb
          Do ka = 1, nka
            agbg(ka, kb) = agbg(ka, kb) + ag(i, ka)*bg(i, kb)*f
            afbf(ka, kb) = afbf(ka, kb) + af(i, ka)*bf(i, kb)*f
          End Do
        End Do
      End Do

      Do kb = 1, nkb
        Do ka = 1, nka
          agbg(ka, kb) = (agbg(ka,kb)+ag(jtop,ka)*bg(jtop,kb)*rpw(jtop))/ &
            3.0E0_dp
          afbf(ka, kb) = (afbf(ka,kb)+af(jtop,ka)*bf(jtop,kb)*rpw(jtop))/ &
            3.0E0_dp
        End Do
      End Do

    End Subroutine
