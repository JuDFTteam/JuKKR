! 20.07.96 ***************************************************************
subroutine crtstar(ratom, nshell, nd, irot, isymindex, rrot)
! ************************************************************************
!  THE SYMMETRY OPERATIONS OF THE SYMMETRY GROUP ARE APPLIED TO THE
!  INPUT VECTOR RATOM
! ------------------------------------------------------------------------
  implicit none
  integer :: irot, nshell
  double precision :: nd(64, 3, *), ratom(3, *), rrot(48, 3, *)
  integer :: isymindex(*)

  integer :: i, id, ns, k, j, isym
  logical :: test
  external :: test
! ------------------------------------------------------------------------

  do ns = 1, nshell
    do id = 1, irot
      isym = isymindex(id)
      do i = 1, 3
        rrot(id, i, ns) = nd(isym, i, 1)*ratom(1, ns) + &
          nd(isym, i, 2)*ratom(2, ns) + nd(isym, i, 3)*ratom(3, ns)
      end do
    end do
  end do

  if (test('ND      ')) write (1337, fmt='((I3,3(/,3f6.2)))')(k, ((nd(k,i,j), &
    j=1,3),i=1,3), k=1, irot)

  return

end subroutine
