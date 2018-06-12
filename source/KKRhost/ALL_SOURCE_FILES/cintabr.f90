subroutine cintabr(ag, bg, agbg, af, bf, afbf, rpw, nka, nkb, jtop, nrmax)
  use :: mod_datatypes, only: dp
  ! ********************************************************************
  ! *                                                                  *
  ! *  SIMPSON - INTERGRATION FOR COMPLEX INTEGRAND  FX FROM 1 TO JTOP *
  ! *  AND EQUIDISTANT MESH    I                                       *
  ! *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
  ! *                                                                  *
  ! *            FX = AG*AG*RPW   and   FX = AF*AF*RPW                 *
  ! *                                                                  *
  ! ********************************************************************

  implicit none


  ! Dummy arguments
  integer :: jtop, nka, nkb, nrmax
  complex (kind=dp) :: af(nrmax, 2), afbf(2, 2), ag(nrmax, 2), agbg(2, 2), &
    bf(nrmax, 2), bg(nrmax, 2)
  real (kind=dp) :: rpw(nrmax)

  ! Local variables
  real (kind=dp) :: f, simp
  integer :: i, ka, kb

  do kb = 1, nkb
    do ka = 1, nka
      agbg(ka, kb) = ag(1, ka)*bg(1, kb)*rpw(1)
      afbf(ka, kb) = af(1, ka)*bf(1, kb)*rpw(1)
    end do
  end do

  if (mod(jtop,2)==0) stop '<CINTABR>  JTOP is even !!!'

  simp = -1.0e0_dp
  do i = 2, jtop - 1
    simp = -simp
    f = (3.0e0_dp+simp)*rpw(i)
    do kb = 1, nkb
      do ka = 1, nka
        agbg(ka, kb) = agbg(ka, kb) + ag(i, ka)*bg(i, kb)*f
        afbf(ka, kb) = afbf(ka, kb) + af(i, ka)*bf(i, kb)*f
      end do
    end do
  end do

  do kb = 1, nkb
    do ka = 1, nka
      agbg(ka, kb) = (agbg(ka,kb)+ag(jtop,ka)*bg(jtop,kb)*rpw(jtop))/3.0e0_dp
      afbf(ka, kb) = (afbf(ka,kb)+af(jtop,ka)*bf(jtop,kb)*rpw(jtop))/3.0e0_dp
    end do
  end do

end subroutine cintabr
