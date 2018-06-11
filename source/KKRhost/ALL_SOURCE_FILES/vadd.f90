! ************************************************************************
    Subroutine vadd(a, b, c)
      Use mod_datatypes, Only: dp
! ************************************************************************

      Real (Kind=dp), Intent (In) :: a(*)
      Real (Kind=dp), Intent (In) :: b(*)
      Real (Kind=dp), Intent (Out) :: c(*)

      Integer :: i

      Do i = 1, 3
        c(i) = a(i) + b(i)
      End Do
      Return
    End Subroutine
