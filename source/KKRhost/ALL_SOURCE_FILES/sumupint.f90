    Subroutine sumupint(sum, vg, g, wg, vf, f, wf, n)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   ********************************************************************
      Implicit None


! Dummy arguments
      Integer :: n
      Complex (Kind=dp) :: sum
      Real (Kind=dp) :: vf, vg
      Complex (Kind=dp) :: f(2, 2), g(2, 2)
      Real (Kind=dp) :: wf(2, 2), wg(2, 2)

! Local variables
      Integer :: i, j

      sum = 0.0E0_dp
      Do j = 1, n
        Do i = 1, n
          sum = sum + vg*g(i, j)*wg(i, j) + vf*f(i, j)*wf(i, j)
        End Do
      End Do

    End Subroutine
