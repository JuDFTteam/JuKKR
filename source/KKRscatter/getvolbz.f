      SUBROUTINE GETVOLBZ(recbv, bravais, VOLBZ)
      IMPLICIT NONE

      DOUBLE PRECISION recbv(3,3), bravais(3,3), volbz, v1
      LOGICAL lsurf
      double precision ddet33
      external ddet33

      lsurf = .false.
      if (bravais(1,3).eq.0.d0.and.bravais(2,3).eq.0.d0.and.
     &     bravais(3,3).eq.0.d0) lsurf = .true.

      if (lsurf) then
         v1 = dabs(recbv(1,1)*recbv(2,2)-recbv(1,2)*recbv(2,1))
      else
         v1 = ddet33(recbv) 
      endif

      VOLBZ = v1

      END SUBROUTINE
