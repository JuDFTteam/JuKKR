! 19.10.95 *************************************************************
    Complex *16 Function csum(n, v, iv)
      Use mod_datatypes, Only: dp
! **********************************************************************
!        sum up the first N elements of the double complex
!        array V(*) with a stepwidth of IV
! ----------------------------------------------------------------------
!.. Scalar Arguments ..
      Integer :: iv, n
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: v(*)
!..
!.. Local Scalars ..
      Complex (Kind=dp) :: vsum
      Integer :: i, ibot, itop
!..
      If (iv>=0) Then
        ibot = 1
        itop = 1 + (n-1)*iv

      Else
        ibot = 1 - (n-1)*iv
        itop = 1
      End If

      vsum = (0E0_dp, 0E0_dp)
      Do i = ibot, itop, iv
        vsum = vsum + v(i)
      End Do
      csum = vsum
      Return
    End Function
