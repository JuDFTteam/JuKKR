      subroutine bessel(jl, nl, hl, z, lmx)
!**********************************************************************
!    attention : contrary to abramowitz and stegun and
!                contrary to subroutine beshan
!
!                the bessel functions of third kind ( hankel functions)
!                are defined as:      hl(l) = nl(l) - i * jl(l)
!**********************************************************************
      double complex, intent(in) :: z
      integer, intent(in) :: lmx
      double complex, intent(out) :: hl(0:lmx), jl(0:lmx), nl(0:lmx)
      
      external :: beshan
      double complex, parameter :: ci=(0.d0,1.d0)

      call beshan(hl, jl, nl, z, lmx)
      hl(0:lmx) = -ci*hl(0:lmx) ! scale with -i

      endsubroutine
