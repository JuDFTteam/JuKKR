! ************************************************************************
    Integer Function length(s, max)
! ************************************************************************

      Character (Len=1), Intent (Inout) :: s(*)
      Integer, Intent (In) :: max

      Integer :: i
! ------------------------------------------------------------------------
      i = max

      Do While (s(i)==' ')
        i = i - 1
      End Do

      length = i

      Return
    End Function
