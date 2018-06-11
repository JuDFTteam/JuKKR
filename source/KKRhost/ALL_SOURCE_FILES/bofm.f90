! **********************************************************************
    Subroutine bofm(pl1, pl2, block, nsize, gin, almd)
! **********************************************************************

      Use mod_datatypes, Only: dp
      Implicit None
!.. Scalar Arguments ..
      Integer :: almd, nsize, pl1, pl2
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: block(nsize, nsize), gin(almd, almd)
!..
!.. Local Scalars ..
      Integer :: i1, i1s, i2, i2s
!..
      i1s = (pl1-1)*nsize
      i2s = (pl2-1)*nsize
      Do i1 = 1, nsize
        Do i2 = 1, nsize
          block(i1, i2) = gin(i1s+i1, i2s+i2)
        End Do
      End Do

      Return

    End Subroutine
