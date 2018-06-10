subroutine rotate(t1, mode, t2, n, rot, nkmmax)
!   ********************************************************************
!   *                                                                  *
!   *   performs the rotation of the matrix  T1  using the rotation-   *
!   *   matrix  ROT, set up by <CALCROTMAT>                            *
!   *                                                                  *
!   *          T2 = ROT  * T1 * ROT+     IF  MODE = 'L->G'             *
!   *          T2 = ROT+ * T1 * ROT      IF  MODE = 'G->L'             *
!   *                                                                  *
!   *   see:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
!   *                                                                  *
!   * 01/11/00                                                         *
!   ********************************************************************
  implicit none

! PARAMETER definitions
  complex *16 :: c0, c1
  parameter (c0=(0.0d0,0.0d0), c1=(1.0d0,0.0d0))

! Dummy arguments
  character (len=4) :: mode
  integer :: n, nkmmax
  complex *16 :: rot(nkmmax, nkmmax), t1(nkmmax, nkmmax), t2(nkmmax, nkmmax)

! Local variables
  character (len=1) :: fl1, fl2
  complex *16 :: w1(nkmmax, nkmmax)


  if (mode=='L->G') then
    fl1 = 'N'
    fl2 = 'C'
  else if (mode=='G->L') then
    fl1 = 'C'
    fl2 = 'N'
  else
    write (*, *) ' MODE = ', mode
    stop 'in <ROTATE>  MODE not allowed'
  end if

  call zgemm(fl1, 'N', n, n, n, c1, rot, nkmmax, t1, nkmmax, c0, w1, nkmmax)

  call zgemm('N', fl2, n, n, n, c1, w1, nkmmax, rot, nkmmax, c0, t2, nkmmax)

end subroutine
