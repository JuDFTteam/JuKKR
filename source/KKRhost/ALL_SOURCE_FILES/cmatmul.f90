    Subroutine cmatmul(n, m, a, b, c)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   perform  the matrix-matrix operation           C = A * B       *
!   *                                                                  *
!   *   A,B,C   complex  SQUARE  N x N - matrices                      *
!   *   N       dimension of A, B and C                                *
!   *   M       array size of A, B, C with M >= N                      *
!   *                                                                  *
!   ********************************************************************

      Implicit None

! PARAMETER definitions
      Complex (Kind=dp) :: c0
      Parameter (c0=(0.0E0_dp,0.0E0_dp))

! Dummy arguments
      Integer :: m, n
      Complex (Kind=dp) :: a(m, m), b(m, m), c(m, m)

! Local variables
      Complex (Kind=dp) :: blj
      Integer :: i, j, l

      Do j = 1, n
        Do i = 1, n
          c(i, j) = c0
        End Do
      End Do

      Do j = 1, n
        Do l = 1, n
          blj = b(l, j)
          If (blj/=c0) Then
            Do i = 1, n
              c(i, j) = c(i, j) + a(i, l)*blj
            End Do
          End If
        End Do
      End Do

    End Subroutine
