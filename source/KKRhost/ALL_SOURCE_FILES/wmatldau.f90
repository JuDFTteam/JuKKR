!-------------------------------------------------------------------------------
! SUBROUTINE: WMATLDAU
!> @brief Calculation of Coulomb interaction potential in LDA+U non-relativistic case
!> otherwise matrices DENMAT and VLDAU must have double dimension
!> @details Uses the Coulomb matrix U (array ULDAU), the density matrix \f$n\f$
!> (array DENMAT) and the occupation numbers dentot (total) and \f$n_s\f$  (array DENTOTS) (per spin).
!> The expression evaluated (array VLDAU) is
!>
!> \f$V_{m1,s,m2,s'} = \delta_{ss'} \sum_{s^{''},m3,m4} U_{m1,m2,m3,m4} n_{m3,s'',m4,s^{''}} - \sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s}-\left[Ueff (dentot-1/2) - Jeff (n_s - 1/2)\right]\delta_{ss'} \delta_{m1,m2} \f$
!> @author Ph. Mavropoulos, H. Ebert (Munich)
!> @date 2002-2004
!> @note Modifications by N. Long Xmas Juelich 2015
!-------------------------------------------------------------------------------
subroutine wmatldau(ntldau, itldau, nspin, denmatc, lopt, ueff, jeff, uldau, &
  wldau, eu, edc, mmaxd, npotd, natyp, nspind, lmax)

  use :: constants

  implicit none
! .. Input variables
  integer, intent (in) :: lmax !< Maximum l component in wave function expansion
  integer, intent (in) :: nspin !< Counter for spin directions
  integer, intent (in) :: mmaxd !< 2*LMAX+1
  integer, intent (in) :: npotd !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
  integer, intent (in) :: natyp !< Number of kinds of atoms in unit cell
  integer, intent (in) :: ntldau !< number of atoms on which LDA+U is applied
  integer, intent (in) :: nspind !< KREL+(1-KREL)*(NSPIN+1)
  integer, dimension (natyp), intent (in) :: lopt !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
  integer, dimension (natyp), intent (in) :: itldau !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
  double precision, dimension (natyp), intent (in) :: ueff !< input U parameter for each atom
  double precision, dimension (natyp), intent (in) :: jeff !< input J parameter for each atom
! .. Input/Output variables
  double precision, dimension (natyp), intent (inout) :: edc !< Double-counting correction
  double precision, dimension (natyp), intent (inout) :: eu !< Total energy corrections
  double precision, dimension (mmaxd, mmaxd, nspind, natyp), &
    intent (inout) :: wldau !< potential matrix
  double precision, dimension (mmaxd, mmaxd, mmaxd, mmaxd, natyp), &
    intent (inout) :: uldau !< calculated Coulomb matrix elements (EREFLDAU)
  double complex, dimension (mmaxd, mmaxd, npotd), intent (inout) :: denmatc
! .. Local variables
  integer :: iprint
  integer :: i1, it, ipot, is, js, m1, m2, m3, m4, mm, mmax
  double precision :: factor
  double precision :: dentot
  character (len=15) :: str15
  double complex :: csum, csum2
  double precision, dimension (nspind) :: dentots
  double precision, dimension (mmaxd, mmaxd, nspind) :: denmat
  double complex, dimension (mmaxd, mmaxd, nspind) :: vldau

!  ..
  data iprint/1/
  data factor/1.d0/ ! if this is 1. then: n*(n-1) in Edc and potential
! if this is 0. then: n**2    in Edc and potential


  write (1337, '(/,79(1H#),/,16X,A,/,79(1H#))') &
    'LDA+U: Calculating interaction potential VLDAU'
!----------------------------------------------------------------------------
  do it = 1, ntldau
    i1 = itldau(it)
!-------------------------------------------------------------------------
    if (lopt(i1)>=0) then
      call rinit(mmaxd*mmaxd*nspind, denmat(1,1,1))
      mmax = 2*lopt(i1) + 1
      write (1337, 100) i1, lopt(i1)
!----------------------------------------------------------------------
! Result is in real Ylm basis.
! It must be converted to complex Ylm basis:
!----------------------------------------------------------------------
      if (iprint>1) write (1337, 110) 'Occupation matrix in REAL basis:'
!----------------------------------------------------------------------
      do is = 1, nspin
        ipot = (i1-1)*nspin + is
        if (iprint>1) then
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, denmatc(1,1,ipot), mmaxd, mmax, 0, 0, 0, &
            1d-8, 1337)
        end if
!-------------------------------------------------------------------
! Convert DENMATC and DENMAT to complex spherical harmonics.
!-------------------------------------------------------------------
        call rclm(1, lopt(i1), lmax, denmatc(1,1,ipot))
      end do
!----------------------------------------------------------------------
      if (iprint>1) then
        write (1337, 110) 'Occupation matrix in COMPLEX basis:'
      end if
      dentot = 0.d0
!----------------------------------------------------------------------
      do is = 1, nspin
        ipot = (i1-1)*nspin + is
        if (iprint>1) then
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, denmatc(1,1,ipot), mmaxd, mmax, 0, 0, 0, &
            1d-8, 1337)
        end if
!-------------------------------------------------------------------
! DENMAT is real: (imag(denmatc))
!-------------------------------------------------------------------
        do m2 = 1, mmax
          do m1 = 1, mmax
! DENMAT(M1,M2,IS) = DIMAG(DENMATC(M1,M2,IPOT))
! in the new medthod, taking imaginary part included
            denmat(m1, m2, is) = (denmatc(m1,m2,ipot))
          end do
        end do
!-------------------------------------------------------------------
! 2.  Calculate total occupation numbers:
! ntot_s = Sum_m n_{m,s,m,s}, ntot = n_1 + n_2
!-------------------------------------------------------------------
        dentots(is) = 0.d0
        do mm = 1, mmax
          dentots(is) = dentots(is) + denmat(mm, mm, is)
        end do
        dentot = dentot + dentots(is)
        dentots(is) = dentots(is)/dfloat(3-nspin)
      end do
!----------------------------------------------------------------------
      if (iprint>0) then
        write (1337, 110) 'Occupation matrix (real):'
        do is = 1, nspin
          write (1337, 120) is
          call rwrite(denmat(1,1,is), mmaxd, mmax, 1337)
          write (1337, 130) 'Trace     =', dentots(is)
        end do
        write (1337, 140) 'Spins sum =', dentot
      end if
!----------------------------------------------------------------------
! In paramagnetic case the spin degeneracy has been accounted
! for by the weight DF in tmatrho.
!----------------------------------------------------------------------
      call cinit(mmaxd*mmaxd*nspind, vldau(1,1,1))
      do is = 1, nspin
!-------------------------------------------------------------------
! 3.  Use density matrix and Coulomb matrix ULDAU to calculate the
! interaction potential VLDAU
! 3a. First part (always diagonal in spin).
!-------------------------------------------------------------------
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
!-------------------------------------------------------------------
! 3b. Second part (in fully rel. case not diagonal in spin; then this
! loop must be changed accordingly).
!-------------------------------------------------------------------
        do m2 = 1, mmax
          do m1 = 1, mmax
            csum = czero
            do m4 = 1, mmax
              do m3 = 1, mmax
                csum = csum - uldau(m1, m4, m3, m2, i1)*denmat(m3, m4, is)/ &
                  dfloat(3-nspin)
              end do
            end do
            vldau(m1, m2, is) = vldau(m1, m2, is) + csum
          end do
        end do
!-------------------------------------------------------------------
! 3c. Third part (always spin- and m-diagonal).
!-------------------------------------------------------------------
        do m1 = 1, mmax
          vldau(m1, m1, is) = vldau(m1, m1, is) - ueff(i1)*(dentot-0.5d0* &
            factor) + jeff(i1)*(dentots(is)-0.5d0*factor)
        end do
      end do
!----------------------------------------------------------------------
! 4. Calculate total-energy corrections EU and EDC (double-counting).
! Then the correction is EU - EDC.
! L[LDA+U]=E[LDA]+E[U]-E[DC]
!> @note EU,EDC initialised outside the routine
!----------------------------------------------------------------------

! Here VLDAU is assumed spin-diagonal (contrary to the spin-orbit case).
!c              DO M2 = 1,MMAX
!c                 DO M1 = 1,MMAX
!c                    EU(I1) = EU(I1) +
!c    &                    DENMAT(M1,M2,IS) * DREAL(VLDAU(M1,M2,IS))
!c                 END DO
!c              END DO

! Relativistic case, see the paper H. Ebert et al., Sol. Stat. Comm. 127 (2003) 443
! Calculate EDC
      do is = 1, nspin
        edc(i1) = edc(i1) + jeff(i1)*dentots(is)*(dentots(is)-factor)
      end do

      edc(i1) = 0.5d0*(ueff(i1)*dentot*(dentot-1.d0)-edc(i1))

! Calculate EU
      do is = 1, nspin
        do js = 1, nspin
          do m4 = 1, mmax
            do m3 = 1, mmax
              do m2 = 1, mmax
                do m1 = 1, mmax
                  eu(i1) = eu(i1) + denmat(m1, m2, is)*uldau(m1, m2, m3, m4, &
                    i1)*denmat(m3, m4, js)
                end do
              end do
            end do
          end do
        end do
      end do

      do is = 1, nspin
        do m4 = 1, mmax
          do m3 = 1, mmax
            do m2 = 1, mmax
              do m1 = 1, mmax
                eu(i1) = eu(i1) - denmat(m1, m2, is)*uldau(m1, m4, m3, m2, i1) &
                  *denmat(m3, m4, is)
              end do
            end do
          end do
        end do
      end do

      eu(i1) = 0.5d0*eu(i1)

!----------------------------------------------------------------------
      if (iprint>0) write (1337, 110) &
        'Interaction potential in COMPLEX basis:'
!----------------------------------------------------------------------
      do is = 1, nspin
        if (iprint>0) then
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, vldau(1,1,is), mmaxd, mmax, 0, 0, 0, 1d-8, &
            1337)
        end if
!-------------------------------------------------------------------
! 5.  Transform VLDAU into real spherical harmonics basis
!-------------------------------------------------------------------
        call rclm(2, lopt(i1), lmax, vldau(1,1,is))
!-------------------------------------------------------------------
! Copy transformed VLDAU to real WLDAU
!
! Apply damping to the interaction matrix WLDAU ? Here not.
!-------------------------------------------------------------------
        do m2 = 1, mmax
          do m1 = 1, mmax
            wldau(m1, m2, is, i1) = dreal(vldau(m1,m2,is))
          end do
        end do
      end do
!----------------------------------------------------------------------
      if (iprint>0) then
        write (1337, 110) 'Interaction potential in REAL basis:'
        do is = 1, nspin
          write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
          call cmatstr(str15, 15, vldau(1,1,is), mmaxd, mmax, 0, 0, 0, 1d-8, &
            1337)
        end do
      end if
!----------------------------------------------------------------------
      write (1337, 110) 'Interaction potential (real):'
      do is = 1, nspin
        write (1337, 120) is
        call rwrite(wldau(1,1,is,i1), mmaxd, mmax, 1337)
      end do
      write (1337, *)
!----------------------------------------------------------------------
! Corrections in total energy:
! Write out corrections on energy:
! E[LDA+U] = E[LDA] + EU - EDC
!----------------------------------------------------------------------
      write (1337, 110) 'Corrections to the total energy:'
      write (1337, *)
      write (1337, 130) 'EU  =', eu(i1)
      write (1337, 130) 'Edc =', edc(i1)
      write (1337, 150) 'E[LDA+U] = E[LDA] + EU - Edc'
    end if
!-------------------------------------------------------------------------
  end do ! I1 = 1,NTLDAU
!----------------------------------------------------------------------------
100 format (/, 6x, 65('='), /, 6x, 'Atom :', i3, ' (l =', i2, ')', /, 6x, &
    18('='))
110 format (8x, '* ', a)
120 format (/, 15x, '> ISPIN =', i1)
130 format (10x, a, f10.6)
140 format (10x, 21('-'), /, 10x, a, f10.6, /, 10x, 60('-'), /)
150 format (27x, a, /)
end subroutine

!-------------------------------------------------------------------------------
! SUBROUTINE: RWRITE
!> @brief Auxiliary subroutine to write the entries of the different potentials
!-------------------------------------------------------------------------------
subroutine rwrite(z, mmaxd, mmax, ifile)

  implicit none
!.. Input variables
  integer, intent (in) :: ifile
  integer, intent (in) :: mmax
  integer, intent (in) :: mmaxd
  double precision, dimension (mmaxd, mmaxd), intent (in) :: z
!.. Local variables
  integer :: m1, m2
!
  write (ifile, 100)
  do m2 = 1, mmax
    write (ifile, 110)(z(m1,m2), m1=1, min(mmax,7))
  end do
  write (ifile, 100)
100 format (10x, 60('-'))
110 format (10x, 7f10.6)
end subroutine
