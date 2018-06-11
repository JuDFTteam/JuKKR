    Function ikapmue(kappa, muem05)
!   ********************************************************************
!   *                                                                  *
!   *  INDEXING OF MATRIX-ELEMENTS:                                    *
!   *                                                                  *
!   *  I = 2*L*(J+1/2) + J + MUE + 1                                   *
!   *                                                                  *
!   ********************************************************************
      Implicit None

! Dummy arguments
      Integer :: kappa, muem05
      Integer :: ikapmue

! Local variables
      Integer :: iabs
      Integer :: jp05, l

      jp05 = iabs(kappa)

      If (kappa<0) Then
        l = -kappa - 1
      Else
        l = kappa
      End If

      ikapmue = 2*l*jp05 + jp05 + muem05 + 1

    End Function
