! 20.07.96 ***************************************************************
    Subroutine crtstar(ratom, nshell, nd, irot, isymindex, rrot)
      Use mod_datatypes, Only: dp
! ************************************************************************
!  THE SYMMETRY OPERATIONS OF THE SYMMETRY GROUP ARE APPLIED TO THE
!  INPUT VECTOR RATOM
! ------------------------------------------------------------------------
      Implicit None
      Integer :: irot, nshell
      Real (Kind=dp) :: nd(64, 3, *), ratom(3, *), rrot(48, 3, *)
      Integer :: isymindex(*)

      Integer :: i, id, ns, k, j, isym
      Logical :: test
      External :: test
! ------------------------------------------------------------------------

      Do ns = 1, nshell
        Do id = 1, irot
          isym = isymindex(id)
          Do i = 1, 3
            rrot(id, i, ns) = nd(isym, i, 1)*ratom(1, ns) + &
              nd(isym, i, 2)*ratom(2, ns) + nd(isym, i, 3)*ratom(3, ns)
          End Do
        End Do
      End Do

      If (test('ND      ')) Write (1337, Fmt='((I3,3(/,3f6.2)))')(k, (( &
        nd(k,i,j),j=1,3),i=1,3), k=1, irot)

      Return

    End Subroutine
