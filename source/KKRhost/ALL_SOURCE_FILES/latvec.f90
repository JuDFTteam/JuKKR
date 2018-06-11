    Logical Function latvec(n, qlat, vec)
      Use mod_datatypes, Only: dp
!- Checks if a set of vectors are lattice vectors
! ----------------------------------------------------------------------
!i Inputs:
!i   n     :number of vectors
!i   qlat  :primitive translation vectors in reciprocal space
!i   vec   :double-precision vector
!o Outputs:
!o   latvec:.true. if all vectors are lattice vectors
!r Remarks:
! ----------------------------------------------------------------------
      Implicit None
! Passed parameters:                                                    
      Integer :: n
      Real (Kind=dp) :: qlat(3, *), vec(3, *)
! Local parameters:                                                     
      Integer :: i, m
      Real (Kind=dp) :: tol, vdiff
      Parameter (tol=1.E-3_dp)
! Common block:                                                         
! Intrinsic functions:                                                  
      Intrinsic :: abs, anint

      latvec = .False.
      Do i = 1, n
        Do m = 1, 3
          vdiff = vec(1, i)*qlat(1, m) + vec(2, i)*qlat(2, m) + &
            vec(3, i)*qlat(3, m)
          vdiff = abs(vdiff-anint(vdiff))
          If (vdiff>tol) Return
        End Do
      End Do
      latvec = .True.
    End Function
