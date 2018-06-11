    Subroutine densitymat(df, pz, qz, pns, qns, ar, cr, dr, gmatll, ipan, &
      ircut, drdi, ek, irmin, lopt, mmax, lmstart, lmend, phi, denmatc, den, &
      ie) ! test fivos 19.9.08
! **********************************************************************
! *                                                                    *
! * Calculation of density matrix needed in evaluating the Coulomb     *
! * interaction potential in LDA+U                                     *
! * non-relativisti! case -- otherwise matrices DENMAT and VLDAU must  *
! *                          have double dimension                     *
! *                                                                    *
! * The density matrix n is calculated using the Green function and    *
! * the reference functions Phi as                                     *
! *                                                                    *
! *  n_{m,s,m',s'} = -1/pi Int dE                                      *
! *    [ Sum_{L''L'''} (Phi_L,R_L'') G_{L''L'''}(E) (R_L''',Phi_L') +  *
! *                                                                    *
! *      + Sum_L'' (Phi_L,R_L'')(H_L'',Phi_L') ]                       *
! *                                                                    *
! * This is expressed in complex spherical harmonics basis. It is then *
! * transformed into the real spherical harmonics basis - subroutine   *
! * WMATLDAU                                                           *
! *                                                                    *
! * Here, only the l-th subblock of G is used. (Is this correct? qldau)*
! *                                                                    *
! * The integration is implicit: this routine is called within an      *
! * energy loop and the result is summed up.                           *
! *                                                                    *
! *                  ph. mavropoulos, h.ebert munich/juelich 2002-2004 *
! **********************************************************************
      Use mod_datatypes
      Implicit None
      Include 'inc.p'
!
! Dummy arguments

      Integer :: lmmaxd
      Parameter (lmmaxd=(krel+1)*(lmaxd+1)**2)
      Integer :: mmaxd
      Parameter (mmaxd=2*lmaxd+1)
      Integer :: irmind
      Parameter (irmind=irmd-irnsd)
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
!
      Complex (Kind=dp) :: czero, cone
      Parameter (czero=(0.0E0_dp,0.0E0_dp), cone=(1.E0_dp,0.E0_dp))
! Local variables

! test fivos 19.9.08
      Complex (Kind=dp) :: df, ek
      Integer :: irmin, lopt
      Integer :: lmstart, lmend, mmax
      Complex (Kind=dp) :: ar(lmmaxd, lmmaxd), cr(lmmaxd, lmmaxd), &
        denmatc(mmaxd, mmaxd), dr(lmmaxd, lmmaxd), gmatll(lmmaxd, lmmaxd), &
        phi(irmd), pns(lmmaxd, lmmaxd, irmind:irmd, 2), pz(irmd, 0:lmaxd), &
        qns(lmmaxd, lmmaxd, irmind:irmd, 2), qz(irmd, 0:lmaxd)
      Real (Kind=dp) :: drdi(irmd)
      Integer :: ipan, ircut(0:ipand)
! test fivos 19.9.08
! test fivos 19.9.08

      Complex (Kind=dp) :: denmatc2(mmaxd, mmaxd), gtemp(mmaxd, mmaxd), &
        phiq(mmaxd, mmaxd), phir(mmaxd, mmaxd)
      Integer :: lm1, lm2, m1, m2
      Integer :: lmaxd1, ie ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      Parameter (lmaxd1=lmaxd+1) ! 1.  Within implicit energy loop:
      Complex (Kind=dp) :: den(0:lmaxd1, iemxd*(1+krel)) ! Calculate density matrix.
      External :: cinit, overlap, rinit
!
      Call cinit(mmaxd*mmaxd, denmatc2(1,1))
! 1a. Calculate overlap integral (inner product) with wavefunctions of
! current energy (Phi,R) (result: PHIR) and (Phi,Q) (result:PHIQ)
! Result is in real Ylm basis.
!

!
! 1b. Use overlap integrals and Green function matrix to find density 
! matrix (implicit integration)
      Call cinit(mmaxd*mmaxd, phir)
      Call overlap(phir, phi, pz, qz, pns, ar, dr, .False., ipan, ircut, drdi, &
        irmin, lopt, ipand, lmaxd, lmmaxd, mmaxd, lmpotd, irmind, irmd)
! Result is in real Ylm basis.
      Call cinit(mmaxd*mmaxd, phiq)
      Call overlap(phiq, phi, pz, qz, qns, cr, dr, .True., ipan, ircut, drdi, &
        irmin, lopt, ipand, lmaxd, lmmaxd, mmaxd, lmpotd, irmind, irmd)
!
! Copy l-th subblock of G into Gtemp (Is this correct? qldau)
!
! First step: PHIQ = G*PHIR+EK*PHIQ.
! (If phi=pz, the trace of this should be similar to the dos).
!
      Do lm2 = lmstart, lmend
        m2 = lm2 - lmstart + 1
        Do lm1 = lmstart, lmend
          m1 = lm1 - lmstart + 1
          gtemp(m1, m2) = gmatll(lm1, lm2)
        End Do
      End Do
! Second step: DENMATC2 = PHIR*PHIQ
!
! Third step: Integration step: DENMAT! = DF*DENMATC2 + DENMATC
!
      Call zgemm('N', 'N', mmax, mmax, mmax, cone, gtemp, mmaxd, phir, mmaxd, &
        ek, phiq, mmaxd)
! test fivos 19.9.08
      Call zgemm('T', 'N', mmax, mmax, mmax, cone, phir, mmaxd, phiq, mmaxd, &
        czero, denmatc2, mmaxd)
! test fivos 19.9.08
! test fivos 19.9.08
! test fivos 19.9.08
      Call cinit(mmaxd*mmaxd, denmatc2(1,1)) ! test fivos
      Do m1 = 1, mmax !         write(*,9001) denmatc(1,1),denmatc(2,2),denmatc(3,3),
        denmatc2(m1, m1) = den(lopt, ie)/real(mmax, kind=dp) !    &        denmatc(4,4),denmatc(5,5),denmatc(6,6),denmatc(7,7)
      End Do !        write(*,9001) -denmatc2(1,1)/3.14159,
      Do m2 = 1, mmax
        Do m1 = 1, mmax
          denmatc(m1, m2) = denmatc(m1, m2) + df*denmatc2(m1, m2)
        End Do
      End Do
!    &                 -denmatc2(2,2)/3.14159,
!    &                 -denmatc2(3,3)/3.14159,
!    &                 -denmatc2(4,4)/3.14159,
!    &                 -denmatc2(5,5)/3.14159,
!    &                 -denmatc2(6,6)/3.14159,
!    &                 -denmatc2(7,7)/3.14159
!9001 format(14e12.4,' test fivos')
!
! Energy loop ends
! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! test fivos 19.9.08
! **********************************************************************
! *                                                                    *
! * Calculation of density matrix needed in evaluating the Coulomb     *
    End Subroutine
