! ************************************************************************
subroutine dsort(w, ind, max, pos)
! ************************************************************************
!     p.zahn, april 96
!     W   is the original array returned unchanged
!     IND is an array that holds the new positions
!     max number of ellements to be sorted
!     pos the position where the first element is found
! ------------------------------------------------------------------------
  implicit none
  integer :: max, pos
  double precision :: w(*)
  integer :: ind(*)

  integer :: i, ii, j, jj, k
  double precision :: bound, diff
  data bound/1.0d-12/
! ------------------------------------------------------------------------
  do i = 1, max
    ind(i) = i
  end do

  j = max
  j = 1
  do while (j<max/3)
    j = 3*j + 1
  end do

  do while (j>1)
    j = j/3
    jj = 1
    do while (jj==1)
      jj = 0
      do k = 1, max - j
        diff = abs(w(ind(k))-w(ind(k+j)))
        if (w(ind(k))>w(ind(k+j)) .and. diff>bound) then
          ii = ind(k)
          ind(k) = ind(k+j)
          ind(k+j) = ii
          jj = 1
        end if
      end do ! K=1,MAX-J
    end do
  end do ! WHILE (JJ.EQ.1)

  do i = 1, max
    if (ind(i)==1) pos = i
  end do

  return
end subroutine


