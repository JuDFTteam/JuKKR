module mod_overlap
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine overlap(result, phi, pz, qz, pqns, acr, dr, lirreg, ipan, ircut, drdi, irmin, lphi, ipand, lmaxd, lmmaxd, mmaxd, lmpotd, irmind, irmd)
    ! **********************************************************************
    ! *                                                                    *
    ! * Calculates the overlap integral of test function PHI with regular  *
    ! * or irregular wavefunction (PZ, QZ: spherical part of the large     *
    ! * component, PQNS: nonspherical part of the large component) for     *
    ! * LDA+U.                                                             *
    ! *                                                                    *
    ! * THETAS are the shape functions (needed only in case of cell        *
    ! *        integration)                                                *
    ! * LPHI is the angular momentum of PHI.                               *
    ! * LIRREG is true if the irrregular wavefunction is to be used.       *
    ! * The overlap is given out in array RESULT                           *
    ! *                                                                    *
    ! * In the inner region the non-spherical wavefunctions are            *
    ! * approximated as:                                                   *
    ! *                                                                    *
    ! *  * the regular one (ir < irmin = irws-irns) :                      *
    ! *                                                                    *
    ! *          pns(ir,lm1,lm2) = pz(ir,l1) * ar(lm1,lm2)                 *
    ! *                                                                    *
    ! *   where pz is the regular wavefct of the spherically symmetric     *
    ! *   part of the potential and ar the alpha matrix (see sub. regns)   *
    ! *                                                                    *
    ! *  * the irregular one (ir < irmin) :                                *
    ! *                                                                    *
    ! *          qns(ir,lm1,lm2) = pz(ir,l1) * cr(lm1,lm2)                 *
    ! *                                + qz(ir,l1) * dr(lm1,lm2)           *
    ! *                                                                    *
    ! *   where pz is the regular and qz is the irregular                  *
    ! *   wavefunction of the spherically symmetric part of the            *
    ! *   potential and cr, dr the matrices calculated at the point irmin  *
    ! *   (see sub. irwns)                                                 *
    ! *                                                                    *
    ! *  Attention: last index of pns,qns is not the spin index but        *
    ! *             sra index.                                             *
    ! *                                                                    *
    ! *  The integrand is convoluted with the shape functions:             *
    ! *  (commented out now)                                               *
    ! *    Result = Int conjg(phi_l)                                       *
    ! *                 sum_{L''L'''} R_L''L' Theta_L''' GAUNT_LL''L'''    *
    ! *                                                                    *
    ! *  Finally, the result should then be transformed to complex         *
    ! *  spherical harmonics basis.                                        *
    ! *                                                                    *
    ! *  Querry: Here only large component is used, but was PHI normalised *
    ! *  according to large + small component? (qldau)                     *
    ! *                                                                    *
    ! *                             ph. mavropoulos, juelich, 2002         *
    ! *                                                                    *
    ! **********************************************************************
    use :: mod_csimpk
    use :: mod_rinit
    use :: mod_cinit
    implicit none
    ! ..
    ! .. Scalar Arguments ..
    integer :: ipand, lmaxd, lmmaxd, mmaxd, lmpotd, irmind, irmd
    integer :: irmin, ipan, lphi   ! l-value for LDA+U
    logical :: lirreg
    ! ..
    ! .. Array arguments ..
    complex (kind=dp) :: phi(irmd), pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd)
    complex (kind=dp) :: pqns(lmmaxd, lmmaxd, irmind:irmd, 2)
    complex (kind=dp) :: acr(lmmaxd, lmmaxd), dr(lmmaxd, lmmaxd)
    complex (kind=dp) :: result(mmaxd, mmaxd)
    real (kind=dp) :: drdi(irmd)
    integer :: ircut(0:ipand)
    ! ..
    ! .. Locals ..
    integer :: lphisq, mmax, irs1, irc1
    integer :: lm1, lm2, lm3, mm1, mm2, mm3, ir
    real (kind=dp) :: gaunt(lmmaxd, lmmaxd, lmpotd)
    complex (kind=dp) :: wint(irmd)
    ! ..
    mmax = 2*lphi + 1
    lphisq = lphi*lphi
    irs1 = ircut(1)
    irc1 = ircut(ipan)

    call rinit(lmmaxd*lmmaxd*lmpotd, gaunt)

    do lm1 = 1, lmmaxd
      gaunt(lm1, lm1, 1) = 1.e0_dp
    end do
    call cinit(mmaxd*mmaxd, result)

    ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    ! Loop over mm1,mm2:
    do mm1 = 1, mmax
      lm1 = lphisq + mm1
      do mm2 = 1, mmax
        lm2 = lphisq + mm2
        ! Set up integrand
        call cinit(irmd, wint)
        ! ----------------------------------------------------------------------
        do mm3 = 1, mmax
          lm3 = lphisq + mm3

          ! -> First inner part (up to irmin, no thetas)
          do ir = 2, irmin
            wint(ir) = wint(ir) + phi(ir)*pz(ir, lphi)*acr(lm3, lm2)*gaunt(lm1, lm3, 1)
          end do

          if (lirreg) then
            do ir = 2, irmin
              wint(ir) = wint(ir) + phi(ir)*qz(ir, lphi)*dr(lm3, lm2)*gaunt(lm1, lm3, 1)
            end do
          end if

          ! -> Next middle part (from irmin+1 to irc1, no thetas)
          do ir = irmin + 1, irs1
            wint(ir) = wint(ir) + phi(ir)*pqns(lm3, lm2, ir, 1)*gaunt(lm1, lm3, 1)
          end do

          ! -> Finally last part - from irc1+1 to irs1 - with THETAS and proper
          ! GAUNTS if we integrate in cell:
          ! ccc               DO  LM4 = 1,LMPOTD
          ! ccc                  IF (LMSP(LM4).GT.0) THEN
          ! ccc                     DO IR = IRS1+1,IRC1
          ! ccc                        WINT(IR) = WINT(IR)
          ! ccc     &                           + PHI(IR) * PQNS(LM3,LM2,IR,1)
          ! ccc     &                           * GAUNT(LM1,LM3,LM4)
          ! ccc     &                           * THETAS(IR-IRS1,LM4)
          ! ccc                     ENDDO
          ! ccc                  END IF
          ! ccc               ENDDO

          ! or still without THETAS if we integrate in sphere

          do ir = irs1 + 1, irc1
            wint(ir) = wint(ir) + phi(ir)*pqns(lm3, lm2, ir, 1)*gaunt(lm1, lm3, 1)
          end do
        end do
        ! ----------------------------------------------------------------------

        call csimpk(wint, result(mm1,mm2), ipan, ircut, drdi(1))
      end do
    end do
    ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  end subroutine overlap

end module mod_overlap
