subroutine gijdmat(tauq, tsst, mssq, dmat, dtil, cfctorinv, iprint, ie, it, &
  krel, lmmaxd)
  use :: mod_datatypes, only: dp
  ! **********************************************************************
  ! * Subroutine to get the projection matrices                          *
  ! *                                                                    *
  ! *                ii     i       i    (-1)                            *
  ! *    D  = [ 1 + G   * (t    -  t  ) ]                                *
  ! *     a          CPA    CPA     a                                    *
  ! *                                                                    *
  ! *    _            i       i      ii   (-1)                           *
  ! *    D  = [ 1 + (t    -  t  ) * G    ]                               *
  ! *     a           CPA     a      CPA                                 *
  ! *                                                                    *
  ! * Used is made of the same routine GETDMAT as in the case of TAU     *
  ! * projection. Note, however, that there the matrices are defined for *
  ! * example as                                                         *
  ! *                                                                    *
  ! *                ij     i       i      (-1)                          *
  ! *    D  = [ 1 + tau * (m    -  m    ) ]                              *
  ! *     a          CPA    a       CPA                                  *
  ! *                                                                    *
  ! * with m = (t)**(-1)                                                 *
  ! *                                                                    *
  ! *                                      v.popescu Oct. 2004           *
  ! **********************************************************************

  implicit none
  ! ..
  ! .. Arguments ..
  integer :: iprint, lmmaxd
  integer :: ie, it, krel
  complex (kind=dp) :: cfctorinv
  complex (kind=dp) :: tauq(lmmaxd, lmmaxd), tsst(lmmaxd, lmmaxd), &
    mssq(lmmaxd, lmmaxd)
  complex (kind=dp) :: dmat(lmmaxd, lmmaxd), dtil(lmmaxd, lmmaxd)
  ! ..
  ! .. Locals ..
  integer :: ik, info, i1, i2
  integer :: ipvt(lmmaxd)
  complex (kind=dp) :: cone, czero
  complex (kind=dp) :: gll(lmmaxd, lmmaxd), tssq(lmmaxd, lmmaxd), &
    tpg(lmmaxd, lmmaxd), xc(lmmaxd, lmmaxd)
  character (len=18) :: banner
  ! ..
  ! .. Externals
  external :: cmatstr, getdmat, zcopy, zgemm, zgetrf, zgetri, zscal
  ! ..
  ! .. Data
  data cone/(1e0_dp, 0e0_dp)/
  data czero/(0e0_dp, 0e0_dp)/
  ! --> get G(CPA) using the same procedure as for GMATLL in < KLOOPZ >
  ! G(CPA) = -MSSQ - MSSQ * TAUQ * MSSQ

  do i2 = 1, lmmaxd
    do i1 = 1, lmmaxd
      tpg(i1, i2) = -cone*tauq(i1, i2)
    end do
  end do

  call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, mssq, lmmaxd, tpg, &
    lmmaxd, czero, xc, lmmaxd)

  call zcopy(lmmaxd*lmmaxd, mssq, 1, gll, 1)

  call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, -cone, xc, lmmaxd, mssq, &
    lmmaxd, -cone, gll, lmmaxd)

  call zscal(lmmaxd*lmmaxd, cfctorinv, gll, 1)

  ! --> invert MSSQ to get TSSQ

  do i2 = 1, lmmaxd
    do i1 = 1, lmmaxd
      tssq(i1, i2) = mssq(i1, i2)*cfctorinv
    end do
  end do
  call zgetrf(lmmaxd, lmmaxd, tssq, lmmaxd, ipvt, info)
  call zgetri(lmmaxd, tssq, lmmaxd, ipvt, xc, lmmaxd*lmmaxd, info)

  ! --> get projection matrices from G(CPA)
  ! call getdmat with t,c for the local arg. c,t

  call getdmat(gll, dmat, dtil, xc, lmmaxd, tsst, tssq, lmmaxd)

  if (iprint<=1) return

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
  ik = 2*krel + 1
  write (banner, '("DMAT IE=",I3," IT=",I3)') ie, it
  call cmatstr(banner, 18, dmat, lmmaxd, lmmaxd, ik, ik, 0, 1e-10_dp, 6)

  write (banner, '("DTIL IE=",I3," IT=",I3)') ie, it
  call cmatstr(banner, 18, dtil, lmmaxd, lmmaxd, ik, ik, 0, 1e-10_dp, 6)
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

end subroutine gijdmat
