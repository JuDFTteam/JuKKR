module mod_bofm

contains

! **********************************************************************
subroutine bofm(pl1, pl2, block, nsize, gin, almd)
  ! **********************************************************************

  use :: mod_datatypes, only: dp
  implicit none
  ! .. Scalar Arguments ..
  integer :: almd, nsize, pl1, pl2
  ! ..
  ! .. Array Arguments ..
  complex (kind=dp) :: block(nsize, nsize), gin(almd, almd)
  ! ..
  ! .. Local Scalars ..
  integer :: i1, i1s, i2, i2s
  ! ..
  i1s = (pl1-1)*nsize
  i2s = (pl2-1)*nsize
  do i1 = 1, nsize
    do i2 = 1, nsize
      block(i1, i2) = gin(i1s+i1, i2s+i2)
    end do
  end do

  return

end subroutine bofm

end module mod_bofm
