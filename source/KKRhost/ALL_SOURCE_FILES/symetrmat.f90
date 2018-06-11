    Subroutine symetrmat(nsym, cpref, dsymll, symunitary, matq, iqs, matsym, &
      lmmaxd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! *  Symmetrising the t/G matrix (or their inverses):                  *
! *                                                                    *
! *                       nsym                                         *
! *     MATSYM = CPREF *  SUM  [ DLL(i) * MATQ(iqs(i)) * DLL(i)^T ]    *
! *                      i = 1                                         *
! *                                                                    *
! *  IQS - set outside the routine - has either the same value         *
! *        regardless of i (e.g. in case of the single-site matrices)  *
! *        or is taking on the value of i => MATSYM is a cummulative   *
! *        sum over MATQ[1...nsym] (e.g. in case of the BZ-integration *
! *        of G)                                                       *
! *                                                                    *
! * CPREF = 1/NSYM  or 1/VBZ                                           *
! *                                                                    *
! *                                    v.popescu, munich nov. 2004     *
! **********************************************************************
      Implicit None
!..
!.. Parameters ..
      Complex (Kind=dp) :: czero, cone
      Parameter (czero=(0E0_dp,0E0_dp), cone=(1E0_dp,0E0_dp))
!..
!.. Arguments ..
      Integer :: lmmaxd, nsym
      Integer :: iqs(*)
      Complex (Kind=dp) :: cpref
      Logical :: symunitary(*)
      Complex (Kind=dp) :: dsymll(lmmaxd, lmmaxd, *), matsym(lmmaxd, lmmaxd), &
        matq(lmmaxd, lmmaxd, *)
!..
!.. Locals ..
      Integer :: isym, l1, l2
      Character (Len=1) :: cnt
      Complex (Kind=dp) :: w1(lmmaxd, lmmaxd), ts(lmmaxd, lmmaxd)
!..
      External :: zcopy, zgemm
!..

      Call zcopy(lmmaxd*lmmaxd, matq(1,1,iqs(1)), 1, ts, 1)

! ----------------------------------------------------------------------
      Do isym = 2, nsym

! --> look if the symmetry operation is unitary / anti-unitary

        cnt = 'N'
        If (.Not. symunitary(isym)) cnt = 'T'

        Call zgemm('N', cnt, lmmaxd, lmmaxd, lmmaxd, cone, dsymll(1,1,isym), &
          lmmaxd, matq(1,1,iqs(isym)), lmmaxd, czero, w1, lmmaxd)

        Call zgemm('N', 'C', lmmaxd, lmmaxd, lmmaxd, cone, w1, lmmaxd, &
          dsymll(1,1,isym), lmmaxd, cone, ts, lmmaxd)
      End Do
! ----------------------------------------------------------------------

      Do l2 = 1, lmmaxd
        Do l1 = 1, lmmaxd
          matsym(l1, l2) = cpref*ts(l1, l2)
        End Do
      End Do

    End Subroutine
