! **********************************************************************
subroutine btom(pl1, pl2, block, nsize, gin, almd, lsub)
!     This subroutine copies or subtracts a block to a matrix
! **********************************************************************
  implicit none
!.. Scalar Arguments ..
  integer :: almd, nsize, pl1, pl2
  logical :: lsub
!..
!.. Array Arguments ..
  double complex :: block(nsize, nsize), gin(almd, almd)
!..
!.. Local Scalars ..
  integer :: i1, i1s, i2, i2s
!     ..
  i1s = (pl1-1)*nsize
  i2s = (pl2-1)*nsize
  if (lsub) then
    do i1 = 1, nsize
      do i2 = 1, nsize
        gin(i1s+i1, i2s+i2) = gin(i1s+i1, i2s+i2) - block(i1, i2)
      end do
    end do
  else
    do i1 = 1, nsize
      do i2 = 1, nsize
        gin(i1s+i1, i2s+i2) = block(i1, i2)
      end do
    end do
  end if

  return

end subroutine
