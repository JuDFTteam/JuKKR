    Function cdnlzdz(l, z, mode)
!   ********************************************************************
!   *                                                                  *
!   *     d n(L,Z) / dz    analytically                                *
!   *                                                                  *
!   ********************************************************************
      Use mod_datatypes, Only: dp
      Implicit None

! Dummy arguments
      Integer :: l, mode
      Complex (Kind=dp) :: z
      Complex (Kind=dp) :: cdnlzdz

! Local variables
      Complex (Kind=dp) :: cnlz

      If (mode==1) Then

        If (l==0) Then

          cdnlzdz = l*cnlz(l, z)/z - cnlz(l+1, z)
        Else
          cdnlzdz = (l*cnlz(l-1,z)-(l+1)*cnlz(l+1,z))/real(2*l+1, kind=dp)
          Return
        End If
      Else If (mode==2) Then

        If (l==0) Then
          cdnlzdz = l*cnlz(l, z)/z - cnlz(l+1, z)
        Else
          cdnlzdz = cnlz(l-1, z) - (l+1)*cnlz(l, z)/z
          Return
        End If
      Else
        cdnlzdz = l*cnlz(l, z)/z - cnlz(l+1, z)
      End If
    End Function
