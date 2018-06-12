subroutine densitymat(df, pz, qz, pns, qns, ar, cr, dr, gmatll, ipan, ircut, &
  drdi, ek, irmin, lopt, mmax, lmstart, lmend, phi, denmatc, den, ie) ! test
                                                                      ! fivos
                                                                      ! 19.9.08
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
  use :: mod_datatypes
  implicit none
  include 'inc.p'

  ! Dummy arguments

  integer :: lmmaxd
  parameter (lmmaxd=(krel+1)*(lmaxd+1)**2)
  integer :: mmaxd
  parameter (mmaxd=2*lmaxd+1)
  integer :: irmind
  parameter (irmind=irmd-irnsd)
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)

  complex (kind=dp) :: czero, cone
  parameter (czero=(0.0e0_dp,0.0e0_dp), cone=(1.e0_dp,0.e0_dp))
  ! Local variables

  ! test fivos 19.9.08
  complex (kind=dp) :: df, ek
  integer :: irmin, lopt
  integer :: lmstart, lmend, mmax
  complex (kind=dp) :: ar(lmmaxd, lmmaxd), cr(lmmaxd, lmmaxd), &
    denmatc(mmaxd, mmaxd), dr(lmmaxd, lmmaxd), gmatll(lmmaxd, lmmaxd), &
    phi(irmd), pns(lmmaxd, lmmaxd, irmind:irmd, 2), pz(irmd, 0:lmaxd), &
    qns(lmmaxd, lmmaxd, irmind:irmd, 2), qz(irmd, 0:lmaxd)
  real (kind=dp) :: drdi(irmd)
  integer :: ipan, ircut(0:ipand)
  ! test fivos 19.9.08
  ! test fivos 19.9.08

  complex (kind=dp) :: denmatc2(mmaxd, mmaxd), gtemp(mmaxd, mmaxd), &
    phiq(mmaxd, mmaxd), phir(mmaxd, mmaxd)
  integer :: lm1, lm2, m1, m2
  integer :: lmaxd1, ie            ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  parameter (lmaxd1=lmaxd+1)       ! 1.  Within implicit energy loop:
  complex (kind=dp) :: den(0:lmaxd1, iemxd*(1+krel)) ! Calculate density
                                                     ! matrix.
  external :: cinit, overlap, rinit

  call cinit(mmaxd*mmaxd, denmatc2(1,1))
  ! 1a. Calculate overlap integral (inner product) with wavefunctions of
  ! current energy (Phi,R) (result: PHIR) and (Phi,Q) (result:PHIQ)
  ! Result is in real Ylm basis.



  ! 1b. Use overlap integrals and Green function matrix to find density
  ! matrix (implicit integration)
  call cinit(mmaxd*mmaxd, phir)
  call overlap(phir, phi, pz, qz, pns, ar, dr, .false., ipan, ircut, drdi, &
    irmin, lopt, ipand, lmaxd, lmmaxd, mmaxd, lmpotd, irmind, irmd)
  ! Result is in real Ylm basis.
  call cinit(mmaxd*mmaxd, phiq)
  call overlap(phiq, phi, pz, qz, qns, cr, dr, .true., ipan, ircut, drdi, &
    irmin, lopt, ipand, lmaxd, lmmaxd, mmaxd, lmpotd, irmind, irmd)

  ! Copy l-th subblock of G into Gtemp (Is this correct? qldau)

  ! First step: PHIQ = G*PHIR+EK*PHIQ.
  ! (If phi=pz, the trace of this should be similar to the dos).

  do lm2 = lmstart, lmend
    m2 = lm2 - lmstart + 1
    do lm1 = lmstart, lmend
      m1 = lm1 - lmstart + 1
      gtemp(m1, m2) = gmatll(lm1, lm2)
    end do
  end do
  ! Second step: DENMATC2 = PHIR*PHIQ

  ! Third step: Integration step: DENMAT! = DF*DENMATC2 + DENMATC

  call zgemm('N', 'N', mmax, mmax, mmax, cone, gtemp, mmaxd, phir, mmaxd, ek, &
    phiq, mmaxd)
  ! test fivos 19.9.08
  call zgemm('T', 'N', mmax, mmax, mmax, cone, phir, mmaxd, phiq, mmaxd, &
    czero, denmatc2, mmaxd)
  ! test fivos 19.9.08
  ! test fivos 19.9.08
  ! test fivos 19.9.08
  call cinit(mmaxd*mmaxd, denmatc2(1,1)) ! test fivos
  do m1 = 1, mmax                  ! write(*,9001)
                                   ! denmatc(1,1),denmatc(2,2),denmatc(3,3),
    denmatc2(m1, m1) = den(lopt, ie)/real(mmax, kind=dp) ! &
                                                         ! denmatc(4,4),denmatc(5,5),denmatc(6,6),denmatc(7,7)
  end do                           ! write(*,9001) -denmatc2(1,1)/3.14159,
  do m2 = 1, mmax
    do m1 = 1, mmax
      denmatc(m1, m2) = denmatc(m1, m2) + df*denmatc2(m1, m2)
    end do
  end do
  ! &                 -denmatc2(2,2)/3.14159,
  ! &                 -denmatc2(3,3)/3.14159,
  ! &                 -denmatc2(4,4)/3.14159,
  ! &                 -denmatc2(5,5)/3.14159,
  ! &                 -denmatc2(6,6)/3.14159,
  ! &                 -denmatc2(7,7)/3.14159
  ! 9001 format(14e12.4,' test fivos')

  ! Energy loop ends
  ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  ! test fivos 19.9.08
  ! **********************************************************************
  ! *                                                                    *
  ! * Calculation of density matrix needed in evaluating the Coulomb     *
end subroutine densitymat
