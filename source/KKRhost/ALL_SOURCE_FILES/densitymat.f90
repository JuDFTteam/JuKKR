module mod_densitymat

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of density matrix for LDA+U
  !> Author: Phivos Mavropoulos, H. Ebert
  !> Date: 2002
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculation of density matrix needed in evaluating the Coulomb     
  !> interaction potential in LDA+U
  !> non-relativistic case -- otherwise matrices DENMAT and VLDAU must  
  !>                          have double dimension
  !>
  !> The density matrix n is calculated using the Green function and    
  !> the reference functions Phi as
  !>
  !>  n_{m,s,m',s'} = -1/pi Int dE
  !>    [ Sum_{L''L'''} (Phi_L,R_L'') G_{L''L'''}(E) (R_L''',Phi_L') +  
  !>
  !>      + Sum_L'' (Phi_L,R_L'')(H_L'',Phi_L') ]
  !>
  !> This is expressed in complex spherical harmonics basis. It is then 
  !> transformed into the real spherical harmonics basis - subroutine   
  !> WMATLDAU
  !>
  !> Here, only the l-th subblock of G is used. (Is this correct? qldau)
  !>
  !> The integration is implicit: this routine is called within an      
  !> energy loop and the result is summed up.
  !>
  !>                  Ph. Mavropoulos, H. Ebert Munich/Juelich 2002-2004
  !-------------------------------------------------------------------------------
  subroutine densitymat(df, pz, qz, pns, qns, ar, cr, dr, gmatll, ipan, ircut, drdi, ek, irmin, lopt, mmax, lmstart, lmend, phi, denmatc, den, ie)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: lmmaxd, mmaxd, irmd, irmind, lmaxd, ipand, krel, iemxd, lmpotd
    use :: mod_overlap, only: overlap
    use :: mod_cinit, only: cinit
    use :: mod_constants, only: czero, cone
    implicit none

    ! Local variables
    complex (kind=dp) :: df, ek
    integer :: irmin, lopt
    integer :: lmstart, lmend, mmax
    complex (kind=dp) :: ar(lmmaxd, lmmaxd), cr(lmmaxd, lmmaxd), denmatc(mmaxd, mmaxd), dr(lmmaxd, lmmaxd), gmatll(lmmaxd, lmmaxd), phi(irmd), pns(lmmaxd, lmmaxd, irmind:irmd, 2), &
      pz(irmd, 0:lmaxd), qns(lmmaxd, lmmaxd, irmind:irmd, 2), qz(irmd, 0:lmaxd)
    real (kind=dp) :: drdi(irmd)
    integer :: ipan, ircut(0:ipand)

    complex (kind=dp) :: denmatc2(mmaxd, mmaxd), gtemp(mmaxd, mmaxd), phiq(mmaxd, mmaxd), phir(mmaxd, mmaxd)
    integer :: lm1, lm2, m1, m2

    integer :: ie
    complex (kind=dp) :: den(0:(lmaxd+1), iemxd*(1+krel))


    call cinit(mmaxd*mmaxd, denmatc2(1,1))

    ! 1.  Within implicit energy loop:
    ! Calculate density matrix.

    ! 1a. Calculate overlap integral (inner product) with wavefunctions of
    ! current energy (Phi,R) (result: PHIR) and (Phi,Q) (result:PHIQ)
    ! Result is in real Ylm basis.
    call cinit(mmaxd*mmaxd, phir)
    call overlap(phir, phi, pz, qz, pns, ar, dr, .false., ipan, ircut, drdi, irmin, lopt, ipand, lmaxd, lmmaxd, mmaxd, lmpotd, irmind, irmd)
    
    call cinit(mmaxd*mmaxd, phiq)
    call overlap(phiq, phi, pz, qz, qns, cr, dr, .true., ipan, ircut, drdi, irmin, lopt, ipand, lmaxd, lmmaxd, mmaxd, lmpotd, irmind, irmd)


    ! 1b. Use overlap integrals and Green function matrix to find density
    ! matrix (implicit integration)
    ! Result is in real Ylm basis.

    ! Copy l-th subblock of G into Gtemp (Is this correct? qldau)
    do lm2 = lmstart, lmend
      m2 = lm2 - lmstart + 1
      do lm1 = lmstart, lmend
        m1 = lm1 - lmstart + 1
        gtemp(m1, m2) = gmatll(lm1, lm2)
      end do
    end do

    ! First step: PHIQ = G*PHIR+EK*PHIQ.
    ! (If phi=pz, the trace of this should be similar to the dos).
    call zgemm('N', 'N', mmax, mmax, mmax, cone, gtemp, mmaxd, phir, mmaxd, ek, phiq, mmaxd)
    ! Second step: DENMATC2 = PHIR*PHIQ
    call zgemm('T', 'N', mmax, mmax, mmax, cone, phir, mmaxd, phiq, mmaxd, czero, denmatc2, mmaxd)

    ! Third step: Integration step: DENMAT! = DF*DENMATC2 + DENMATC
    call cinit(mmaxd*mmaxd, denmatc2(1,1)) 
    do m1 = 1, mmax 
      denmatc2(m1, m1) = den(lopt, ie)/real(mmax, kind=dp)
    end do
    do m2 = 1, mmax
      do m1 = 1, mmax
        denmatc(m1, m2) = denmatc(m1, m2) + df*denmatc2(m1, m2)
      end do
    end do

    ! Energy loop ends
    
  end subroutine densitymat

end module mod_densitymat
