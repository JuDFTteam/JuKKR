    Function ssum(n, v, iv)
      Use mod_datatypes, Only: dp
! **********************************************************************
!        sum up the first N elements of the double precision
!        array V(*) with a stepwidth of IV
! ----------------------------------------------------------------------
      Implicit None
      Real (Kind=dp) :: ssum
!.. Scalar Arguments ..
      Integer :: iv, n
!..
!.. Array Arguments ..
      Real (Kind=dp) :: v(*)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: vsum
      Integer :: i, ibot, itop
!..
      If (iv>=0) Then
        ibot = 1
        itop = 1 + (n-1)*iv

      Else
        ibot = 1 - (n-1)*iv
        itop = 1
      End If

      vsum = 0.0E0_dp
      Do i = ibot, itop, iv
        vsum = vsum + v(i)
      End Do
      ssum = vsum
    End Function
