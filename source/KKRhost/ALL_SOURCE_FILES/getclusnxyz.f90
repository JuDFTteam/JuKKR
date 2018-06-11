    Subroutine getclusnxyz(clurad, bravais, ndim, cluradsq, nbr)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * Given a spherical cluster of radius CLURAD it determines the three *
! * integers N1,N2,N3 such that any vector                             *
! *                                                                    *
! *    R_i = r_i + SUM_j  N_j * a_j                                    *
! *                                                                    *
! *  with i = 1,NAEZ and a_j the primitive Bravais vectors, is inside  *
! *  the cluster. Subroutine also returns the CLURAD**2 value          *
! *                                                                    *
! **********************************************************************
      Implicit None
! ..
! ..  Arguments
      Real (Kind=dp) :: clurad, cluradsq
      Real (Kind=dp) :: bravais(3, 3)
      Integer :: ndim, nbr(3)
! .. 
! ..  Locals
      Real (Kind=dp) :: dr(3)
      Integer :: i, j
      Integer :: int
! ..
! ..
      Do i = 1, ndim
        dr(i) = 0E0_dp
        Do j = 1, ndim
          dr(i) = dr(i) + bravais(j, i)*bravais(j, i)
        End Do
        dr(i) = sqrt(dr(i))
      End Do

      If (abs(clurad)<1E-6_dp) Then
        Do i = 1, ndim
          nbr(i) = 0
        End Do
        cluradsq = 1E10_dp
      Else
        Do i = 1, ndim
          nbr(i) = int(clurad/dr(i)) + 2
        End Do
        cluradsq = clurad*clurad
      End If

    End Subroutine
