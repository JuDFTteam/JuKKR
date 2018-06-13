subroutine wmatldau(ntldau, itldau, nspin, denmatc, lopt, ueff, jeff, uldau, &
  wldau, eu, edc, mmaxd, npotd)
  ! **********************************************************************
  ! *                                                                    *
  ! * Calculation of Coulomb interaction potential in LDA+U              *
  ! * non-relativistic case -- otherwise matrices DENMAT and VLDAU must  *
  ! *                          have double dimension                     *
  ! *                                                                    *
  ! * Uses the Coulomb matrix U (array ULDAU), the density matrix n      *
  ! * (array DENMAT) and the occupation numbers dentot (total) and n_s   *
  ! * (array DENTOTS) (per spin).                                        *
  ! *                                                                    *
  ! * The expression evaluated (array VLDAU) is                          *
  ! *                                                                    *
  ! *       V_{m1,s,m2,s'} =                                             *
  ! * delta_{ss'} Sum_{s'',m3,m4} U_{m1,m2,m3,m4} n_{m3,s'',m4,s''}      *
  ! * - Sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s}                       *
  ! * - [Ueff (dentot-1/2) - Jeff (n_s - 1/2)] delta_{ss'} delta_{m1,m2} *
  ! *                                                                    *
  ! *                  ph. mavropoulos, h.ebert munich/juelich 2002-2004 *
  ! **********************************************************************
  use :: mod_datatypes, only: dp
  use global_variables
  implicit none

  complex (kind=dp) :: czero
  parameter (czero=(0.0d0,0.0d0))
  ! Local variables
  ! ..
  ! ..
  integer :: ntldau, nspin, mmaxd, npotd
  integer :: itldau(natypd), lopt(natypd)
  real (kind=dp) :: ueff(natypd), jeff(natypd), edc(natypd), eu(natypd), &
    wldau(mmaxd, mmaxd, nspind, natypd)
  real (kind=dp), allocatable :: uldau(:, :, :, :, :)
  complex (kind=dp) :: denmatc(mmaxd, mmaxd, npotd)


  complex (kind=dp) :: csum, csum2, vldau(mmaxd, mmaxd, nspind)
  real (kind=dp) :: denmat(mmaxd, mmaxd, nspind), dentot, dentots(nspind)
  integer :: i1, it, ipot, is, js, m1, m2, m3, m4, mm, mmax
  integer :: iprint
  character (len=15) :: str15

  data iprint/1/
  ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
  write (1337, '(/,79(1H#),/,16X,A,/,79(1H#))') &
    'LDA+U: Calculating interaction potential VLDAU'

  allocate (uldau(mmaxd,mmaxd,mmaxd,mmaxd,natypd))
  ! Result is in real Ylm basis.
  ! It must be converted to complex Ylm basis:
  do it = 1, ntldau
    i1 = itldau(it)

    if (lopt(i1)>=0) then
      call rinit(mmaxd*mmaxd*nspind, denmat(1,1,1))
      mmax = 2*lopt(i1) + 1
      write (1337, 100) i1, lopt(i1)
      ! ----------------------------------------------------------------------

      ! -> Convert DENMATC and DENMAT to complex spherical harmonics.

      if (iprint>1) write (1337, 110) 'Occupation matrix in REAL basis:'
      ! ----------------------------------------------------------------------
      do is = 1, nspin
        ipot = (i1-1)*nspin + is
        if (iprint>1) then
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, denmatc(1,1,ipot), mmaxd, mmax, 0, 0, 0, &
            1d-8, 6)
        end if
        ! ----------------------------------------------------------------------

        ! -> DENMAT is real: (imag(denmatc))
        call rclm(1, lopt(i1), lmaxd, denmatc(1,1,ipot))
      end do

      if (iprint>1) write (1337, 110) 'Occupation matrix in COMPLEX basis:'
      dentot = 0.d0

      do is = 1, nspin
        ipot = (i1-1)*nspin + is
        if (iprint>1) then
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, denmatc(1,1,ipot), mmaxd, mmax, 0, 0, 0, &
            1d-8, 6)
        end if
        ! 2.  Calculate total occupation numbers:
        ! ntot_s = Sum_m n_{m,s,m,s}, ntot = n_1 + n_2

        do m2 = 1, mmax
          do m1 = 1, mmax
            denmat(m1, m2, is) = dimag(denmatc(m1,m2,ipot))
          end do
        end do
        ! ----------------------------------------------------------------------

        ! In paramagnetic case the spin degeneracy has been accounted
        ! for by the weight DF in tmatrho.
        dentots(is) = 0.d0
        do mm = 1, mmax
          dentots(is) = dentots(is) + denmat(mm, mm, is)
        end do
        dentot = dentot + dentots(is)
      end do

      if (iprint>0) then
        write (1337, 110) 'Occupation matrix (real):'
        do is = 1, nspin
          write (1337, 120) is
          call rwrite(denmat(1,1,is), mmaxd, mmax, 1337)
          write (1337, 130) 'Trace     =', dentots(is)
        end do
        write (1337, 140) 'Spins sum =', dentot
      end if
      ! ----------------------------------------------------------------------

      ! 3.  Use density matrix and Coulomb matrix ULDAU to calculate the
      ! interaction potential VLDAU
      ! 3a. First part (always diagonal in spin).
      call cinit(mmaxd*mmaxd*nspind, vldau(1,1,1))
      do is = 1, nspin


        ! 3b. Second part (in fully rel. case not diagonal in spin; then this
        ! loop must be changed accordingly).

        do m2 = 1, mmax
          do m1 = 1, mmax
            csum = czero
            do m4 = 1, mmax
              do m3 = 1, mmax
                csum2 = czero
                do js = 1, nspin
                  csum2 = csum2 + denmat(m3, m4, js)
                end do
                csum = csum + uldau(m1, m2, m3, m4, i1)*csum2
              end do
            end do
            vldau(m1, m2, is) = vldau(m1, m2, is) + csum
          end do
        end do

        ! 3c. Third part (always spin- and m-diagonal).


        do m2 = 1, mmax
          do m1 = 1, mmax
            csum = czero
            do m4 = 1, mmax
              do m3 = 1, mmax
                csum = csum - uldau(m1, m4, m3, m2, i1)*denmat(m3, m4, is)
              end do
            end do
            vldau(m1, m2, is) = vldau(m1, m2, is) + csum
          end do
        end do
        ! 4. Calculate total-energy corrections EU and EDC (double-counting).
        ! Then the correction is EU-EDC.
        ! Note: EU,EDC initialised outside the routine
        do m1 = 1, mmax
          vldau(m1, m1, is) = vldau(m1, m1, is) - ueff(i1)*(dentot-0.5d0) + &
            jeff(i1)*(dentots(is)-0.5d0)
        end do

        ! Here VLDAU is assumed spin-diagonal (contrary to the spin-orbit
        ! case).

        ! ----------------------------------------------------------------------
        ! ----------------------------------------------------------------------

        ! 5.  Transform VLDAU into real spherical harmonics basis
        do m2 = 1, mmax
          do m1 = 1, mmax
            eu(i1) = eu(i1) + denmat(m1, m2, is)*dreal(vldau(m1,m2,is))
          end do
        end do
        edc(i1) = edc(i1) + jeff(i1)*dentots(is)*(dentots(is)-1.d0)
      end do

      if (iprint>0) write (1337, 110) &
        'Interaction potential in COMPLEX basis:'

      do is = 1, nspin
        if (iprint>0) then
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, vldau(1,1,is), mmaxd, mmax, 0, 0, 0, 1d-8, &
            6)
        end if
        ! Copy transformed VLDAU to real WLDAU

        ! Apply damping to the interaction matrix WLDAU ? Here not.
        call rclm(2, lopt(i1), lmaxd, vldau(1,1,is))



        ! ----------------------------------------------------------------------
        ! ----------------------------------------------------------------------
        ! ----------------------------------------------------------------------
        do m2 = 1, mmax
          do m1 = 1, mmax
            wldau(m1, m2, is, i1) = dreal(vldau(m1,m2,is))
          end do
        end do

      end do
      ! Corrections in total energy:
      if (iprint>0) then
        write (1337, 110) 'Interaction potential in REAL basis:'
        do is = 1, nspin
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, vldau(1,1,is), mmaxd, mmax, 0, 0, 0, 1d-8, &
            6)
        end do
      end if

      write (1337, 110) 'Interaction potential (real):'
      do is = 1, nspin
        write (1337, 120) is
        call rwrite(wldau(1,1,is,i1), mmaxd, mmax, 1337)
      end do
      write (1337, *)

      ! -> Write out corrections on energy:
      ! E[LDA+U] = E[LDA] + EU - EDC

      eu(i1) = 0.5d0*eu(i1)
      edc(i1) = 0.5d0*(ueff(i1)*dentot*(dentot-1.d0)-edc(i1))
      ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      ! I1 = 1,NTLDAU
      ! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      ! **********************************************************************
      write (1337, 110) 'Corrections to the total energy:'
      write (1337, *)
      write (1337, 130) 'EU  =', eu(i1)
      write (1337, 130) 'Edc =', edc(i1)
      write (1337, 150) 'E[LDA+U] = E[LDA] + EU - Edc'
    end if
    ! *                                                                    *
  end do                           ! * Calculation of Coulomb interaction
                                   ! potential in LDA+U              *
  ! * non-relativistic case -- otherwise matrices DENMAT and VLDAU must  *
100 format (/, 6x, 65('='), /, 6x, 'Atom :', i3, ' (l =', i2, ')', /, 6x, &
    18('='))
110 format (8x, '* ', a)
120 format (/, 15x, '> ISPIN =', i1)
130 format (10x, a, f10.6)
140 format (10x, 21('-'), /, 10x, a, f10.6, /, 10x, 60('-'), /)
150 format (27x, a, /)
end subroutine wmatldau


subroutine rwrite(z, mmaxd, mmax, ifile)

  real (kind=dp), intent (inout) :: z(mmaxd, mmaxd)
  integer, intent (inout) :: mmaxd
  integer, intent (in) :: mmax
  integer, intent (inout) :: ifile
  implicit none

  integer :: m1, m2

  write (ifile, 100)

  do m2 = 1, mmax
    write (ifile, 110)(z(m1,m2), m1=1, min(mmax,7))
  end do
  write (ifile, 100)
100 format (10x, 60('-'))
110 format (10x, 7f10.6)
end subroutine rwrite
