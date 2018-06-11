    Function cdjlzdz(l, z, mode)
!   ********************************************************************
!   *                                                                  *
!   *     d j(L,Z) / dz    analytically                                *
!   *                                                                  *
!   ********************************************************************
      Use mod_datatypes, Only: dp
      Implicit None

! Dummy arguments
      Integer :: l, mode
      Complex (Kind=dp) :: z
      Complex (Kind=dp) :: cdjlzdz

! Local variables
      Complex (Kind=dp) :: cjlz

      If (mode==1) Then

        If (l==0) Then

          cdjlzdz = l*cjlz(l, z)/z - cjlz(l+1, z)
        Else
          cdjlzdz = (l*cjlz(l-1,z)-(l+1)*cjlz(l+1,z))/real(2*l+1, kind=dp)
          Return
        End If
      Else If (mode==2) Then

        If (l==0) Then
          cdjlzdz = l*cjlz(l, z)/z - cjlz(l+1, z)
        Else
          cdjlzdz = cjlz(l-1, z) - (l+1)*cjlz(l, z)/z
          Return
        End If
      Else
        cdjlzdz = l*cjlz(l, z)/z - cjlz(l+1, z)
      End If
    End Function
