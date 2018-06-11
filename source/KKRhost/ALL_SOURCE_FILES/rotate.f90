    Subroutine rotate(t1, mode, t2, n, rot, nkmmax)
      Use mod_datatypes, Only: dp
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
      Implicit None

! PARAMETER definitions
      Complex (Kind=dp) :: c0, c1
      Parameter (c0=(0.0E0_dp,0.0E0_dp), c1=(1.0E0_dp,0.0E0_dp))

! Dummy arguments
      Character (Len=4) :: mode
      Integer :: n, nkmmax
      Complex (Kind=dp) :: rot(nkmmax, nkmmax), t1(nkmmax, nkmmax), &
        t2(nkmmax, nkmmax)

! Local variables
      Character (Len=1) :: fl1, fl2
      Complex (Kind=dp) :: w1(nkmmax, nkmmax)


      If (mode=='L->G') Then
        fl1 = 'N'
        fl2 = 'C'
      Else If (mode=='G->L') Then
        fl1 = 'C'
        fl2 = 'N'
      Else
        Write (*, *) ' MODE = ', mode
        Stop 'in <ROTATE>  MODE not allowed'
      End If

      Call zgemm(fl1, 'N', n, n, n, c1, rot, nkmmax, t1, nkmmax, c0, w1, &
        nkmmax)

      Call zgemm('N', fl2, n, n, n, c1, w1, nkmmax, rot, nkmmax, c0, t2, &
        nkmmax)

    End Subroutine
