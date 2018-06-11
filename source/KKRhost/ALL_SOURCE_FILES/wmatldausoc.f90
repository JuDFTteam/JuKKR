!-------------------------------------------------------------------------------
! SUBROUTINE: WMATLDAUSOC
!> @brief Calculation of Coulomb interaction potential in LDA+U relativistic + SOC (new solver)
!> @details The expression evaluated (array VLDAU) is:
!> \f$V_{m1,s,m2,s'} =\delta_{ss'} \sum_{s'',m3,m4} U_{m1,m2,m3,m4} n_{m3,s'',m4,s''}-\sum_{m3,m4} U_{m1,m4,m3,m2} n_{m3,s',m4,s} - \left[Ueff (dentot-1/2) - Jeff (n_s - 1/2)\right] \delta_{ss'} \delta_{m1,m2}\f$
!>
!> For details see H. Ebert at al., Sol. Stat. Comm. 127 (2003) 443
!> @author N. long
!> @date 04.2016
!-------------------------------------------------------------------------------
    Subroutine wmatldausoc(ntldau, itldau, nspin, denmatn, lopt, ueff, jeff, &
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
      Use constants
      Use mod_datatypes

      Implicit None
!
! .. Input variables
      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: mmaxd !< 2*LMAX+1
      Integer, Intent (In) :: nspind !< KREL+(1-KREL)*(NSPIN+1)
      Integer, Intent (In) :: ntldau !< number of atoms on which LDA+U is applied
      Integer, Dimension (natyp), Intent (In) :: lopt !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
      Integer, Dimension (natyp), Intent (In) :: itldau !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
      Real (Kind=dp), Dimension (natyp), Intent (In) :: ueff !< input U parameter for each atom
      Real (Kind=dp), Dimension (natyp), Intent (In) :: jeff !< input J parameter for each atom
! .. Input/Output variables
      Real (Kind=dp), Dimension (natyp), Intent (Inout) :: eu !< Total energy corrections
      Real (Kind=dp), Dimension (natyp), Intent (Inout) :: edc !< Double-counting correction
      Real (Kind=dp), Dimension (mmaxd, mmaxd, nspind, natyp), &
        Intent (Inout) :: wldau !< potential matrix
      Real (Kind=dp), Dimension (mmaxd, mmaxd, mmaxd, mmaxd, natyp), &
        Intent (In) :: uldau !< calculated Coulomb matrix elements (EREFLDAU)
      Complex (Kind=dp), Dimension (mmaxd, mmaxd, 2, 2, natyp), &
        Intent (Inout) :: denmatn

! .. Local variables
      Integer :: iprint
      Integer :: i1, it, is, js, m1, m2, m3, m4, mm, mmax
      Real (Kind=dp) :: dentot
      Real (Kind=dp) :: factor
      Complex (Kind=dp) :: csum, csum2
      Character (Len=15) :: str15
      Real (Kind=dp), Dimension (nspind) :: dentots
      Real (Kind=dp), Dimension (mmaxd, mmaxd, 2, 2) :: denmat
      Complex (Kind=dp), Dimension (mmaxd, mmaxd, 2, 2) :: vldau
!     ..
      Data iprint/1/
      Data factor/1.E0_dp/ ! if this is 1. then: n*(n-1) in Edc and potential
! if this is 0. then: n**2 in Edc and potential

      Write (1337, '(/,79("#"),/,16X,A,/,79("#"))') &
        'LDA+U: Calculating interaction potential VLDAU'
!----------------------------------------------------------------------------
      Do it = 1, ntldau
        i1 = itldau(it)
!-------------------------------------------------------------------------
        If (lopt(i1)>=0) Then
          Call rinit(mmaxd*mmaxd*2*2, denmat(1,1,1,1))
          mmax = 2*lopt(i1) + 1
          Write (1337, 100) i1, lopt(i1)
!----------------------------------------------------------------------
! Result is in real Ylm basis.
! It must be converted to complex Ylm basis:
!----------------------------------------------------------------------
          If (iprint>1) Write (1337, 110) 'Occupation matrix in REAL basis:'
!----------------------------------------------------------------------
          Do is = 1, nspin
            If (iprint>1) Then
              Write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
              Call cmatstr(str15, 15, denmatn(1,1,is,is,i1), mmaxd, mmax, 0, &
                0, 0, 1E-8_dp, 1337)
            End If
!-------------------------------------------------------------------
! Convert DENMATC and DENMAT to complex spherical harmonics.
!-------------------------------------------------------------------
            Do js = 1, nspin
              Call rclm(1, lopt(i1), lmax, denmatn(1,1,js,is,i1))
            End Do ! js
          End Do ! is
!----------------------------------------------------------------------
          If (iprint>1) Write (1337, 110) &
            'Occupation matrix in COMPLEX basis:'
          dentot = 0.E0_dp
!----------------------------------------------------------------------
          Do is = 1, nspin
            If (iprint>1) Then
              Write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
              Call cmatstr(str15, 15, denmatn(1,1,is,is,i1), mmaxd, mmax, 0, &
                0, 0, 1E-8_dp, 1337)
            End If
!-------------------------------------------------------------------
! DENMAT is real: (imag(denmatc))
!-------------------------------------------------------------------
            Do js = 1, nspin
              Do m2 = 1, mmax
                Do m1 = 1, mmax
                  denmat(m1, m2, js, is) = (denmatn(m1,m2,js,is,i1))
                End Do
              End Do
            End Do ! js
          End Do ! is
!----------------------------------------------------------------------
! 2.  Calculate total occupation numbers:
! ntot_s = Sum_m n_{m,s,m,s}, ntot = n_1 + n_2
!----------------------------------------------------------------------
          Do is = 1, nspin
            dentots(is) = 0.E0_dp
            Do js = 1, nspin
              Do mm = 1, mmax
                dentots(is) = dentots(is) + denmat(mm, mm, js, is)
              End Do
            End Do ! JS
            dentot = dentot + dentots(is)
          End Do ! IS
!----------------------------------------------------------------------
          If (iprint>0) Then
            Write (1337, 110) 'Occupation matrix (real):'
            Do is = 1, nspin
              Write (1337, 120) is
              Call rwrite(denmat(1,1,is,is), mmaxd, mmax, 1337)
              Write (1337, 130) 'Trace     =', dentots(is)
            End Do
            Write (1337, 140) 'Spins sum =', dentot
          End If
!----------------------------------------------------------------------
          Call cinit(mmaxd*mmaxd*2*2, vldau(1,1,1,1))
          Do is = 1, nspin
!-------------------------------------------------------------------
! 3.  Use density matrix and Coulomb matrix ULDAU to calculate the
! interaction potential VLDAU
! 3a. First part (always diagonal in spin).
!-------------------------------------------------------------------
            Do m2 = 1, mmax
              Do m1 = 1, mmax
                csum = czero
                Do m4 = 1, mmax
                  Do m3 = 1, mmax
                    csum2 = czero
                    Do js = 1, nspin
                      csum2 = csum2 + denmat(m3, m4, js, js)
                    End Do
                    csum = csum + uldau(m1, m2, m3, m4, i1)*csum2
                  End Do
                End Do
                vldau(m1, m2, is, is) = vldau(m1, m2, is, is) + csum
              End Do
            End Do
!-------------------------------------------------------------------
! 3b. Second part
!-------------------------------------------------------------------
            Do js = 1, nspin
              Do m2 = 1, mmax
                Do m1 = 1, mmax
                  csum = czero
                  Do m4 = 1, mmax
                    Do m3 = 1, mmax
                      csum = csum - uldau(m1, m4, m3, m2, i1)*denmat(m3, m4, &
                        js, is)
                    End Do
                  End Do
                  vldau(m1, m2, js, is) = vldau(m1, m2, js, is) + csum
                End Do
              End Do
            End Do ! js
!-------------------------------------------------------------------
! 3c. Third part (always spin- and m-diagonal).
!-------------------------------------------------------------------
            Do m1 = 1, mmax
              vldau(m1, m1, is, is) = vldau(m1, m1, is, is) - &
                ueff(i1)*(dentot-0.5E0_dp*factor) + jeff(i1)*(dentots(is)- &
                0.5E0_dp*factor)
            End Do
          End Do ! IS
!----------------------------------------------------------------------
! 4. Calculate total-energy corrections EU and EDC (double-counting).
! Then the correction is EU - EDC.
! L[LDA+U]=E[LDA]+E[U]-E[DC]
!> @note: EU,EDC initialised outside the routine
!----------------------------------------------------------------------
! Calculate EDC
          Do is = 1, nspin
            edc(i1) = edc(i1) + jeff(i1)*dentots(is)*(dentots(is)-factor)
          End Do

          edc(i1) = 0.5E0_dp*(ueff(i1)*dentot*(dentot-1.E0_dp)-edc(i1))

! Calculate EU
          Do is = 1, nspin
            Do js = 1, nspin
              Do m4 = 1, mmax
                Do m3 = 1, mmax
                  Do m2 = 1, mmax
                    Do m1 = 1, mmax
                      eu(i1) = eu(i1) + denmat(m1, m2, is, is)*uldau(m1, m2, &
                        m3, m4, i1)*denmat(m3, m4, js, js)
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do

          Do is = 1, nspin
            Do js = 1, nspin
              Do m4 = 1, mmax
                Do m3 = 1, mmax
                  Do m2 = 1, mmax
                    Do m1 = 1, mmax
                      eu(i1) = eu(i1) - denmat(m1, m2, is, js)*uldau(m1, m4, &
                        m3, m2, i1)*denmat(m3, m4, js, is)
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do

          eu(i1) = 0.5E0_dp*eu(i1)
!----------------------------------------------------------------------
          If (iprint>0) Write (1337, 110) &
            'Interaction potential in COMPLEX basis:'
!----------------------------------------------------------------------
          Do is = 1, nspin
            wldau(:, :, is, i1) = 0E0_dp
            If (iprint>0) Then
              Write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
              Call cmatstr(str15, 15, vldau(1,1,is,is), mmaxd, mmax, 0, 0, 0, &
                1E-8_dp, 1337)
            End If
!-------------------------------------------------------------------
! 5.  Transform VLDAU into real spherical harmonics basis
!-------------------------------------------------------------------
            Do js = 1, nspin
              Call rclm(2, lopt(i1), lmax, vldau(1,1,js,is))
!----------------------------------------------------------------
! Copy transformed VLDAU to real WLDAU
! Apply damping to the interaction matrix WLDAU ? Here not.
!----------------------------------------------------------------
              Do m2 = 1, mmax
                Do m1 = 1, mmax
                  wldau(m1, m2, is, i1) = wldau(m1, m2, is, i1) + &
                    real(vldau(m1,m2,js,is))
                End Do
              End Do
            End Do ! js
          End Do ! is
!----------------------------------------------------------------------
          If (iprint>0) Then
            Write (1337, 110) 'Interaction potential in REAL basis:'
            Do is = 1, nspin
              Write (str15, '(4X,"> ",A,I1)') 'ISPIN = ', is
              Call cmatstr(str15, 15, vldau(1,1,is,is), mmaxd, mmax, 0, 0, 0, &
                1E-8_dp, 1337)
            End Do
          End If
!----------------------------------------------------------------------
          Write (1337, 110) 'Interaction potential (real):'
          Do is = 1, nspin
            Write (1337, 120) is
            Call rwrite(wldau(1,1,is,i1), mmaxd, mmax, 1337)
          End Do
          Write (1337, *)
!----------------------------------------------------------------------
! Corrections in total energy:
! Write out corrections on energy:
!    E[LDA+U] = E[LDA] + EU - EDC
!----------------------------------------------------------------------
          Write (1337, 110) 'Corrections to the total energy:'
          Write (1337, *)
          Write (1337, 130) 'EU  =', eu(i1)
          Write (1337, 130) 'Edc =', edc(i1)
          Write (1337, 150) 'E[LDA+U] = E[LDA] + EU - Edc'
        End If
!-------------------------------------------------------------------------
      End Do ! I1 = 1,NTLDAU
!----------------------------------------------------------------------------
100   Format (/, 6X, 65('='), /, 6X, 'Atom :', I3, ' (l =', I2, ')', /, 6X, &
        18('='))
110   Format (8X, '* ', A)
120   Format (/, 15X, '> ISPIN =', I1)
130   Format (10X, A, F10.6)
140   Format (10X, 21('-'), /, 10X, A, F10.6, /, 10X, 60('-'), /)
150   Format (27X, A, /)
    End Subroutine
