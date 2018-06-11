! ************************************************************************
    Subroutine veq(a, b)
      Use mod_datatypes, Only: dp
! ************************************************************************
      Implicit None

      Real (Kind=dp) :: a(*)
      Real (Kind=dp) :: b(*)

      Integer :: i

      Do i = 1, 3
        b(i) = a(i)
      End Do
    End Subroutine
