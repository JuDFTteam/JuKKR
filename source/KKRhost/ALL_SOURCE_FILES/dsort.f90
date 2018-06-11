! ************************************************************************
    Subroutine dsort(w, ind, max, pos)
      Use mod_datatypes, Only: dp
! ************************************************************************
!     p.zahn, april 96
!     W   is the original array returned unchanged
!     IND is an array that holds the new positions
!     max number of ellements to be sorted
!     pos the position where the first element is found
! ------------------------------------------------------------------------
      Implicit None
      Integer :: max, pos
      Real (Kind=dp) :: w(*)
      Integer :: ind(*)

      Integer :: i, ii, j, jj, k
      Real (Kind=dp) :: bound, diff
      Data bound/1.0E-12_dp/
! ------------------------------------------------------------------------
      Do i = 1, max
        ind(i) = i
      End Do

      j = max
      j = 1
      Do While (j<max/3)
        j = 3*j + 1
      End Do

      Do While (j>1)
        j = j/3
        jj = 1
        Do While (jj==1)
          jj = 0
          Do k = 1, max - j
            diff = abs(w(ind(k))-w(ind(k+j)))
            If (w(ind(k))>w(ind(k+j)) .And. diff>bound) Then
              ii = ind(k)
              ind(k) = ind(k+j)
              ind(k+j) = ii
              jj = 1
            End If
          End Do ! K=1,MAX-J
        End Do
      End Do ! WHILE (JJ.EQ.1)

      Do i = 1, max
        If (ind(i)==1) pos = i
      End Do

      Return
    End Subroutine
