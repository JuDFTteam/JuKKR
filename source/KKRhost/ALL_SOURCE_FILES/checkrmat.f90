    Function checkrmat(rmat, co1, si1, co2, si2, co3, si3, i, j)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *  check whether the values of the cosinus and sinus found for the *
!   *  Euler angles TET1, TET2, TET3 are consistent with the           *
!   *  rotation matrix   RMAT                                          *
!   *                                                                  *
!   ********************************************************************

      Implicit None

! Dummy arguments
      Real (Kind=dp) :: co1, co2, co3, si1, si2, si3
      Integer :: i, j
      Logical :: checkrmat
      Real (Kind=dp) :: rmat(3, 3)

! Local variables
      Real (Kind=dp) :: a, b
      Logical :: equal
      Logical :: result

      equal(a, b) = (abs(a-b)<1E-7_dp)

      result = .False.

      If (i==1) Then
        If (j==1) Then
          result = equal(rmat(1,1), co3*co2*co1-si3*si1)
        Else If (j==2) Then
          result = equal(rmat(1,2), co3*co2*si1+si3*co1)
        Else If (j==3) Then
          result = equal(rmat(1,3), -co3*si2)
        End If
      Else If (i==2) Then
        If (j==1) Then
          result = equal(rmat(2,1), -si3*co2*co1-co3*si1)
        Else If (j==2) Then
          result = equal(rmat(2,2), -si3*co2*si1+co3*co1)
        Else If (j==3) Then
          result = equal(rmat(2,3), si3*si2)
        End If
      Else If (j==1) Then
        result = equal(rmat(3,1), si2*co1)
      Else If (j==2) Then
        result = equal(rmat(3,2), si2*si1)
      Else If (j==3) Then
        result = equal(rmat(3,3), co2)
      End If
      checkrmat = result
    End Function
