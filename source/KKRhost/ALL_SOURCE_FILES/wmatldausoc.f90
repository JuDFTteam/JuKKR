! -------------------------------------------------------------------------------
! SUBROUTINE: WMATLDAUSOC
! > @brief Calculation of Coulomb interaction potential in LDA+U relativistic
! + SOC (new solver)
! > @details The expression evaluated (array VLDAU) is:
! > \f$V_{m1,s,m2,s'} =\delta_{ss'} \sum_{s'',m3,m4} U_{m1,m2,m3,m4}
! n_{m3,s'',m4,s''}-\sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s} - \left[Ueff
! (dentot-1/2) - Jeff (n_s - 1/2)\right] \delta_{ss'} \delta_{m1,m2}\f$
! >
! > For details see H. Ebert at al., Sol. Stat. Comm. 127 (2003) 443
! > @author N. long
! > @date 04.2016
! -------------------------------------------------------------------------------
subroutine wmatldausoc(ntldau, itldau, nspin, denmatn, lopt, ueff, jeff, &
  uldau, wldau, eu, edc, mmaxd, natyp, nspind, lmax)
  ! **********************************************************************
  ! *                                                                    *
  ! * Calculation of Coulomb interaction potential in LDA+U              *
  ! * relativistic + SOC (new solver)                                    *
  ! *                                                                    *
  ! * The expression evaluated (array VLDAU) is                          *
  ! *                                                                    *
  ! *       V_{m1,s,m2,s'} =                                             *
  ! * delta_{ss'} Sum_{s'',m3,m4} U_{m1,m2,m3,m4} n_{m3,s'',m4,s''}      *
  ! * - Sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s}                       *
  ! * - [Ueff (dentot-1/2) - Jeff (n_s - 1/2)] delta_{ss'} delta_{m1,m2} *
  ! *                                                                    *
  ! * details see H. Ebert at al., Sol. Stat. Comm. 127 (2003) 443       *
  ! *                                                                    *
  ! *                  n.long,  April 2016, Juelich                      *
  ! **********************************************************************
  use :: constants
  use :: mod_datatypes

  implicit none

  ! .. Input variables
  integer, intent (in) :: lmax     ! < Maximum l component in wave function
                                   ! expansion
  integer, intent (in) :: natyp    ! < Number of kinds of atoms in unit cell
  integer, intent (in) :: nspin    ! < Counter for spin directions
  integer, intent (in) :: mmaxd    ! < 2*LMAX+1
  integer, intent (in) :: nspind   ! < KREL+(1-KREL)*(NSPIN+1)
  integer, intent (in) :: ntldau   ! < number of atoms on which LDA+U is
                                   ! applied
  integer, dimension (natyp), intent (in) :: lopt ! < angular momentum QNUM
                                                  ! for the atoms on which
                                                  ! LDA+U should be applied
                                                  ! (-1 to switch it OFF)
  integer, dimension (natyp), intent (in) :: itldau ! < integer pointer
                                                    ! connecting the NTLDAU
                                                    ! atoms to heir
                                                    ! corresponding index in
                                                    ! the unit cell
  real (kind=dp), dimension (natyp), intent (in) :: ueff ! < input U parameter
                                                         ! for each atom
  real (kind=dp), dimension (natyp), intent (in) :: jeff ! < input J parameter
                                                         ! for each atom
  ! .. Input/Output variables
  real (kind=dp), dimension (natyp), intent (inout) :: eu ! < Total energy
                                                          ! corrections
  real (kind=dp), dimension (natyp), intent (inout) :: edc ! < Double-counting
                                                           ! correction
  real (kind=dp), dimension (mmaxd, mmaxd, nspind, natyp), &
    intent (inout) :: wldau        ! < potential matrix
  real (kind=dp), dimension (mmaxd, mmaxd, mmaxd, mmaxd, natyp), &
    intent (in) :: uldau           ! < calculated Coulomb matrix elements
                                   ! (EREFLDAU)
  complex (kind=dp), dimension (mmaxd, mmaxd, 2, 2, natyp), &
    intent (inout) :: denmatn

  ! .. Local variables
  integer :: iprint
  integer :: i1, it, is, js, m1, m2, m3, m4, mm, mmax
  real (kind=dp) :: dentot
  real (kind=dp) :: factor
  complex (kind=dp) :: csum, csum2
  character (len=15) :: str15
  real (kind=dp), dimension (nspind) :: dentots
  real (kind=dp), dimension (mmaxd, mmaxd, 2, 2) :: denmat
  complex (kind=dp), dimension (mmaxd, mmaxd, 2, 2) :: vldau
  ! ..
  data iprint/1/
  data factor/1.e0_dp/             ! if this is 1. then: n*(n-1) in Edc and
                                   ! potential
  ! if this is 0. then: n**2 in Edc and potential

  write (1337, '(/,79("#"),/,16X,A,/,79("#"))') &
    'LDA+U: Calculating interaction potential VLDAU'
  ! ----------------------------------------------------------------------------
  do it = 1, ntldau
    i1 = itldau(it)
    ! -------------------------------------------------------------------------
    if (lopt(i1)>=0) then
      call rinit(mmaxd*mmaxd*2*2, denmat(1,1,1,1))
      mmax = 2*lopt(i1) + 1
      write (1337, 100) i1, lopt(i1)
      ! ----------------------------------------------------------------------
      ! Result is in real Ylm basis.
      ! It must be converted to complex Ylm basis:
      ! ----------------------------------------------------------------------
      if (iprint>1) write (1337, 110) 'Occupation matrix in REAL basis:'
      ! ----------------------------------------------------------------------
      do is = 1, nspin
        if (iprint>1) then
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, denmatn(1,1,is,is,i1), mmaxd, mmax, 0, 0, 0, &
            1e-8_dp, 1337)
        end if
        ! -------------------------------------------------------------------
        ! Convert DENMATC and DENMAT to complex spherical harmonics.
        ! -------------------------------------------------------------------
        do js = 1, nspin
          call rclm(1, lopt(i1), lmax, denmatn(1,1,js,is,i1))
        end do                     ! js
      end do                       ! is
      ! ----------------------------------------------------------------------
      if (iprint>1) write (1337, 110) 'Occupation matrix in COMPLEX basis:'
      dentot = 0.e0_dp
      ! ----------------------------------------------------------------------
      do is = 1, nspin
        if (iprint>1) then
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, denmatn(1,1,is,is,i1), mmaxd, mmax, 0, 0, 0, &
            1e-8_dp, 1337)
        end if
        ! -------------------------------------------------------------------
        ! DENMAT is real: (imag(denmatc))
        ! -------------------------------------------------------------------
        do js = 1, nspin
          do m2 = 1, mmax
            do m1 = 1, mmax
              denmat(m1, m2, js, is) = (denmatn(m1,m2,js,is,i1))
            end do
          end do
        end do                     ! js
      end do                       ! is
      ! ----------------------------------------------------------------------
      ! 2.  Calculate total occupation numbers:
      ! ntot_s = Sum_m n_{m,s,m,s}, ntot = n_1 + n_2
      ! ----------------------------------------------------------------------
      do is = 1, nspin
        dentots(is) = 0.e0_dp
        do js = 1, nspin
          do mm = 1, mmax
            dentots(is) = dentots(is) + denmat(mm, mm, js, is)
          end do
        end do                     ! JS
        dentot = dentot + dentots(is)
      end do                       ! IS
      ! ----------------------------------------------------------------------
      if (iprint>0) then
        write (1337, 110) 'Occupation matrix (real):'
        do is = 1, nspin
          write (1337, 120) is
          call rwrite(denmat(1,1,is,is), mmaxd, mmax, 1337)
          write (1337, 130) 'Trace     =', dentots(is)
        end do
        write (1337, 140) 'Spins sum =', dentot
      end if
      ! ----------------------------------------------------------------------
      call cinit(mmaxd*mmaxd*2*2, vldau(1,1,1,1))
      do is = 1, nspin
        ! -------------------------------------------------------------------
        ! 3.  Use density matrix and Coulomb matrix ULDAU to calculate the
        ! interaction potential VLDAU
        ! 3a. First part (always diagonal in spin).
        ! -------------------------------------------------------------------
        do m2 = 1, mmax
          do m1 = 1, mmax
            csum = czero
            do m4 = 1, mmax
              do m3 = 1, mmax
                csum2 = czero
                do js = 1, nspin
                  csum2 = csum2 + denmat(m3, m4, js, js)
                end do
                csum = csum + uldau(m1, m2, m3, m4, i1)*csum2
              end do
            end do
            vldau(m1, m2, is, is) = vldau(m1, m2, is, is) + csum
          end do
        end do
        ! -------------------------------------------------------------------
        ! 3b. Second part
        ! -------------------------------------------------------------------
        do js = 1, nspin
          do m2 = 1, mmax
            do m1 = 1, mmax
              csum = czero
              do m4 = 1, mmax
                do m3 = 1, mmax
                  csum = csum - uldau(m1, m4, m3, m2, i1)*denmat(m3, m4, js, &
                    is)
                end do
              end do
              vldau(m1, m2, js, is) = vldau(m1, m2, js, is) + csum
            end do
          end do
        end do                     ! js
        ! -------------------------------------------------------------------
        ! 3c. Third part (always spin- and m-diagonal).
        ! -------------------------------------------------------------------
        do m1 = 1, mmax
          vldau(m1, m1, is, is) = vldau(m1, m1, is, is) - &
            ueff(i1)*(dentot-0.5e0_dp*factor) + jeff(i1)*(dentots(is)-0.5e0_dp &
            *factor)
        end do
      end do                       ! IS
      ! ----------------------------------------------------------------------
      ! 4. Calculate total-energy corrections EU and EDC (double-counting).
      ! Then the correction is EU - EDC.
      ! L[LDA+U]=E[LDA]+E[U]-E[DC]
      ! > @note: EU,EDC initialised outside the routine
      ! ----------------------------------------------------------------------
      ! Calculate EDC
      do is = 1, nspin
        edc(i1) = edc(i1) + jeff(i1)*dentots(is)*(dentots(is)-factor)
      end do

      edc(i1) = 0.5e0_dp*(ueff(i1)*dentot*(dentot-1.e0_dp)-edc(i1))

      ! Calculate EU
      do is = 1, nspin
        do js = 1, nspin
          do m4 = 1, mmax
            do m3 = 1, mmax
              do m2 = 1, mmax
                do m1 = 1, mmax
                  eu(i1) = eu(i1) + denmat(m1, m2, is, is)*uldau(m1, m2, m3, &
                    m4, i1)*denmat(m3, m4, js, js)
                end do
              end do
            end do
          end do
        end do
      end do

      do is = 1, nspin
        do js = 1, nspin
          do m4 = 1, mmax
            do m3 = 1, mmax
              do m2 = 1, mmax
                do m1 = 1, mmax
                  eu(i1) = eu(i1) - denmat(m1, m2, is, js)*uldau(m1, m4, m3, &
                    m2, i1)*denmat(m3, m4, js, is)
                end do
              end do
            end do
          end do
        end do
      end do

      eu(i1) = 0.5e0_dp*eu(i1)
      ! ----------------------------------------------------------------------
      if (iprint>0) write (1337, 110) &
        'Interaction potential in COMPLEX basis:'
      ! ----------------------------------------------------------------------
      do is = 1, nspin
        wldau(:, :, is, i1) = 0e0_dp
        if (iprint>0) then
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, vldau(1,1,is,is), mmaxd, mmax, 0, 0, 0, &
            1e-8_dp, 1337)
        end if
        ! -------------------------------------------------------------------
        ! 5.  Transform VLDAU into real spherical harmonics basis
        ! -------------------------------------------------------------------
        do js = 1, nspin
          call rclm(2, lopt(i1), lmax, vldau(1,1,js,is))
          ! ----------------------------------------------------------------
          ! Copy transformed VLDAU to real WLDAU
          ! Apply damping to the interaction matrix WLDAU ? Here not.
          ! ----------------------------------------------------------------
          do m2 = 1, mmax
            do m1 = 1, mmax
              wldau(m1, m2, is, i1) = wldau(m1, m2, is, i1) + &
                real(vldau(m1,m2,js,is))
            end do
          end do
        end do                     ! js
      end do                       ! is
      ! ----------------------------------------------------------------------
      if (iprint>0) then
        write (1337, 110) 'Interaction potential in REAL basis:'
        do is = 1, nspin
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, vldau(1,1,is,is), mmaxd, mmax, 0, 0, 0, &
            1e-8_dp, 1337)
        end do
      end if
      ! ----------------------------------------------------------------------
      write (1337, 110) 'Interaction potential (real):'
      do is = 1, nspin
        write (1337, 120) is
        call rwrite(wldau(1,1,is,i1), mmaxd, mmax, 1337)
      end do
      write (1337, *)
      ! ----------------------------------------------------------------------
      ! Corrections in total energy:
      ! Write out corrections on energy:
      ! E[LDA+U] = E[LDA] + EU - EDC
      ! ----------------------------------------------------------------------
      write (1337, 110) 'Corrections to the total energy:'
      write (1337, *)
      write (1337, 130) 'EU  =', eu(i1)
      write (1337, 130) 'Edc =', edc(i1)
      write (1337, 150) 'E[LDA+U] = E[LDA] + EU - Edc'
    end if
    ! -------------------------------------------------------------------------
  end do                           ! I1 = 1,NTLDAU
  ! ----------------------------------------------------------------------------
100 format (/, 6x, 65('='), /, 6x, 'Atom :', i3, ' (l =', i2, ')', /, 6x, &
    18('='))
110 format (8x, '* ', a)
120 format (/, 15x, '> ISPIN =', i1)
130 format (10x, a, f10.6)
140 format (10x, 21('-'), /, 10x, a, f10.6, /, 10x, 60('-'), /)
150 format (27x, a, /)
end subroutine wmatldausoc
