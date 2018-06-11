    Integer Function lngstring(string, lstrmax)
!   ********************************************************************
!   *                                                                  *
!   *  find position of last non-blank character in STRING(1:LSTRMAX)  *
!   *                                                                  *
!   ********************************************************************
      Implicit None

! Dummy arguments
      Integer :: lstrmax
      Character (Len=*) :: string

! Local variables
      Character :: c
      Integer :: i, ichar

      lngstring = 0
      Do i = lstrmax, 1, -1
        c = string(i:i)
        If (c/=' ' .And. ichar(c)>0) Then
          lngstring = i
          Return
        End If
      End Do
    End Function
