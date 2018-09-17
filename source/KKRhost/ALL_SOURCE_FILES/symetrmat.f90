module mod_symetrmat
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine symetrmat(nsym, cpref, dsymll, symunitary, matq, iqs, matsym, lmmaxd, nsymaxd)
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
    implicit none
    ! ..
    ! .. Parameters ..
    complex (kind=dp) :: czero, cone
    parameter (czero=(0e0_dp,0e0_dp), cone=(1e0_dp,0e0_dp))
    ! ..
    ! .. Arguments ..
    integer :: lmmaxd, nsym, nsymaxd
    integer :: iqs(nsymaxd)
    complex (kind=dp) :: cpref
    logical :: symunitary(nsymaxd)
    complex (kind=dp) :: dsymll(lmmaxd, lmmaxd, nsymaxd), matsym(lmmaxd, lmmaxd), matq(lmmaxd, lmmaxd, *)
    ! ..
    ! .. Locals ..
    integer :: isym, l1, l2
    character (len=1) :: cnt
    complex (kind=dp) :: w1(lmmaxd, lmmaxd), ts(lmmaxd, lmmaxd)
    ! ..

    call zcopy(lmmaxd*lmmaxd, matq(1,1,iqs(1)), 1, ts, 1)

    ! ----------------------------------------------------------------------
    do isym = 2, nsym

      ! --> look if the symmetry operation is unitary / anti-unitary

      cnt = 'N'
      if (.not. symunitary(isym)) cnt = 'T'

      call zgemm('N', cnt, lmmaxd, lmmaxd, lmmaxd, cone, dsymll(1,1,isym), lmmaxd, matq(1,1,iqs(isym)), lmmaxd, czero, w1, lmmaxd)

      call zgemm('N', 'C', lmmaxd, lmmaxd, lmmaxd, cone, w1, lmmaxd, dsymll(1,1,isym), lmmaxd, cone, ts, lmmaxd)
    end do
    ! ----------------------------------------------------------------------

    do l2 = 1, lmmaxd
      do l1 = 1, lmmaxd
        matsym(l1, l2) = cpref*ts(l1, l2)
      end do
    end do

  end subroutine symetrmat

end module mod_symetrmat
