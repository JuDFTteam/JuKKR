      subroutine beshan(hl, jl, nl, z, lmax)
!-----------------------------------------------------------------------
!  calculates spherical bessel, hankel and neumann functions
!  for the orders lmin .le. l .le. lmax.
!  for |z|  <   1 the taylor expansions of jl and nl are used.
!  for |z|  >=  1 the explicit expressions for hl(+), hl(-) are used.
!
!                            r. zeller   jan. 1990
!-----------------------------------------------------------------------
      double complex, intent(in) :: z
      integer, intent(in) :: lmax
      double complex, intent(out) :: hl(0:lmax), jl(0:lmax), nl(0:lmax)
      
      double complex, parameter :: ci=(0.d0,1.d0)
      double complex :: termj, termn, z2, zj, zn
      double precision :: rl, rn, rnm
      integer :: l, m, n

      zj = 1.d0
      zn = 1.d0
      z2 = z*z
      if (abs(z) < lmax+1.d0) then
        do l = 0, lmax
          rl = l + l
          termj = -0.5d0/(rl+3.d0)*z2
          termn =  0.5d0/(rl-1.d0)*z2
          jl(l) = 1.d0
          nl(l) = 1.d0
          do n = 2, 25
            jl(l) = jl(l) + termj
            nl(l) = nl(l) + termn
            rn = n + n
            termj = -termj/(rl+rn+1.d0)/rn*z2
            termn =  termn/(rl-rn+1.d0)/rn*z2
          enddo ! n
          jl(l) = jl(l)*zj
          nl(l) = -nl(l)*zn/z
          hl(l) = jl(l) + nl(l)*ci

          zj = zj*z/(rl+3.d0)
          zn = zn/z*(rl+1.d0)
        enddo ! l
      endif

      do l = 0, lmax
        if (abs(z) >= l+1.d0) then
          hl(l) = 0.d0
          nl(l) = 0.d0
          rnm = 1.d0
          do m = 0, l
            hl(l) = hl(l) + rnm/(-ci*(z+z))**m
            nl(l) = nl(l) + rnm/( ci*(z+z))**m
            rnm = rnm*(l*l+l-m*m-m)/(m+1.d0)
          enddo ! m
          hl(l) = hl(l)*(-ci)**l*exp( ci*z)/( ci*z)
          nl(l) = nl(l)*( ci)**l*exp(-ci*z)/(-ci*z)
          jl(l) = (hl(l)+nl(l))*0.5d0
          nl(l) = (hl(l)-jl(l))/ci
        endif
      enddo ! l

      endsubroutine
