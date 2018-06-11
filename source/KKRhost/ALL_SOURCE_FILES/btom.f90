! **********************************************************************
    Subroutine btom(pl1, pl2, block, nsize, gin, almd, lsub)
      Use mod_datatypes, Only: dp
!     This subroutine copies or subtracts a block to a matrix
! **********************************************************************
      Implicit None
!.. Scalar Arguments ..
      Integer :: almd, nsize, pl1, pl2
      Logical :: lsub
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: block(nsize, nsize), gin(almd, almd)
!..
!.. Local Scalars ..
      Integer :: i1, i1s, i2, i2s
!     ..
      i1s = (pl1-1)*nsize
      i2s = (pl2-1)*nsize
      If (lsub) Then
        Do i1 = 1, nsize
          Do i2 = 1, nsize
            gin(i1s+i1, i2s+i2) = gin(i1s+i1, i2s+i2) - block(i1, i2)
          End Do
        End Do
      Else
        Do i1 = 1, nsize
          Do i2 = 1, nsize
            gin(i1s+i1, i2s+i2) = block(i1, i2)
          End Do
        End Do
      End If

      Return

    End Subroutine
