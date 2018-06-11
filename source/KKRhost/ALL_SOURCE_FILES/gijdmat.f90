    Subroutine gijdmat(tauq, tsst, mssq, dmat, dtil, cfctorinv, iprint, ie, &
      it, krel, lmmaxd)
      Use mod_datatypes, Only: dp
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

      Implicit None
!..
!.. Arguments ..
      Integer :: iprint, lmmaxd
      Integer :: ie, it, krel
      Complex (Kind=dp) :: cfctorinv
      Complex (Kind=dp) :: tauq(lmmaxd, lmmaxd), tsst(lmmaxd, lmmaxd), &
        mssq(lmmaxd, lmmaxd)
      Complex (Kind=dp) :: dmat(lmmaxd, lmmaxd), dtil(lmmaxd, lmmaxd)
!..
!.. Locals ..
      Integer :: ik, info, i1, i2
      Integer :: ipvt(lmmaxd)
      Complex (Kind=dp) :: cone, czero
      Complex (Kind=dp) :: gll(lmmaxd, lmmaxd), tssq(lmmaxd, lmmaxd), &
        tpg(lmmaxd, lmmaxd), xc(lmmaxd, lmmaxd)
      Character (Len=18) :: banner
!..
!.. Externals
      External :: cmatstr, getdmat, zcopy, zgemm, zgetrf, zgetri, zscal
!..
!.. Data
      Data cone/(1E0_dp, 0E0_dp)/
      Data czero/(0E0_dp, 0E0_dp)/
! --> get G(CPA) using the same procedure as for GMATLL in < KLOOPZ >
!           G(CPA) = -MSSQ - MSSQ * TAUQ * MSSQ

      Do i2 = 1, lmmaxd
        Do i1 = 1, lmmaxd
          tpg(i1, i2) = -cone*tauq(i1, i2)
        End Do
      End Do

      Call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, mssq, lmmaxd, tpg, &
        lmmaxd, czero, xc, lmmaxd)

      Call zcopy(lmmaxd*lmmaxd, mssq, 1, gll, 1)

      Call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, -cone, xc, lmmaxd, mssq, &
        lmmaxd, -cone, gll, lmmaxd)

      Call zscal(lmmaxd*lmmaxd, cfctorinv, gll, 1)

! --> invert MSSQ to get TSSQ

      Do i2 = 1, lmmaxd
        Do i1 = 1, lmmaxd
          tssq(i1, i2) = mssq(i1, i2)*cfctorinv
        End Do
      End Do
      Call zgetrf(lmmaxd, lmmaxd, tssq, lmmaxd, ipvt, info)
      Call zgetri(lmmaxd, tssq, lmmaxd, ipvt, xc, lmmaxd*lmmaxd, info)

! --> get projection matrices from G(CPA)
!     call getdmat with t,c for the local arg. c,t

      Call getdmat(gll, dmat, dtil, xc, lmmaxd, tsst, tssq, lmmaxd)

      If (iprint<=1) Return

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      ik = 2*krel + 1
      Write (banner, '("DMAT IE=",I3," IT=",I3)') ie, it
      Call cmatstr(banner, 18, dmat, lmmaxd, lmmaxd, ik, ik, 0, 1E-10_dp, 6)

      Write (banner, '("DTIL IE=",I3," IT=",I3)') ie, it
      Call cmatstr(banner, 18, dtil, lmmaxd, lmmaxd, ik, ik, 0, 1E-10_dp, 6)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    End Subroutine
