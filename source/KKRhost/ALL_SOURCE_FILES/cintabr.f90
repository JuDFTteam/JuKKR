SUBROUTINE cintabr(ag,bg,agbg,af,bf,afbf,rpw,nka,nkb,jtop,nrmax)
!   ********************************************************************
!   *                                                                  *
!   *  SIMPSON - INTERGRATION FOR COMPLEX INTEGRAND  FX FROM 1 TO JTOP *
!   *  AND EQUIDISTANT MESH    I                                       *
!   *   INT = [ F1 + 4*F2 + 2*F3 + .... + 4F(N-1) + FN ]/3     N:ODD   *
!   *                                                                  *
!   *            FX = AG*AG*RPW   and   FX = AF*AF*RPW                 *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE


!Dummy arguments
INTEGER JTOP,NKA,NKB,NRMAX
COMPLEX*16 AF(NRMAX,2),AFBF(2,2),AG(NRMAX,2),AGBG(2,2),BF(NRMAX,2) &
           ,BG(NRMAX,2)
REAL*8 RPW(NRMAX)

!Local variables
REAL*8 F,SIMP
INTEGER I,KA,KB

DO kb = 1,nkb
  DO ka = 1,nka
    agbg(ka,kb) = ag(1,ka)*bg(1,kb)*rpw(1)
    afbf(ka,kb) = af(1,ka)*bf(1,kb)*rpw(1)
  END DO
END DO

IF ( MOD(jtop,2) == 0 ) STOP '<CINTABR>  JTOP is even !!!'

simp = -1.0D0
DO i = 2,jtop - 1
  simp = -simp
  f = (3.0D0+simp)*rpw(i)
  DO kb = 1,nkb
    DO ka = 1,nka
      agbg(ka,kb) = agbg(ka,kb) + ag(i,ka)*bg(i,kb)*f
      afbf(ka,kb) = afbf(ka,kb) + af(i,ka)*bf(i,kb)*f
    END DO
  END DO
END DO

DO kb = 1,nkb
  DO ka = 1,nka
    agbg(ka,kb) = (agbg(ka,kb)+ag(jtop,ka)*bg(jtop,kb)*rpw(jtop) )/3.0D0
    afbf(ka,kb) = (afbf(ka,kb)+af(jtop,ka)*bf(jtop,kb)*rpw(jtop) )/3.0D0
  END DO
END DO

END SUBROUTINE cintabr
