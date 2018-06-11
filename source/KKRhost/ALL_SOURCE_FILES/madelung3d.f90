    Subroutine madelung3d(lpot, yrg, wg, naez, alat, volume0, bravais, recbv, &
      rbasis, rmax, gmax, naezd, lmxspd, lassld, lpotd, lmpotd, nmaxd, ishld, &
      nembd, wlength)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * This subroutine calculates the Madelung potential coefficients     *
! * in the 3D case and stores them in the DA-file abvmad.unformatted   *
! * The record index is simply (IQ1-1)*NAEZ + IQ2 for record (IQ1,IQ2) *
! *                                                                    *
! **********************************************************************

      Implicit None
!..
!.. Scalar Arguments ..
      Integer :: lpot, naez, wlength
      Integer :: naezd, lmxspd, lassld, lpotd, lmpotd, nmaxd, ishld, nembd
      Real (Kind=dp) :: alat, volume0, rmax, gmax
!..
!.. Array Arguments ..
      Real (Kind=dp) :: yrg(lassld, 0:lassld, 0:lassld), wg(lassld)
      Real (Kind=dp) :: bravais(3, 3), recbv(3, 3)
      Real (Kind=dp) :: rbasis(3, naezd+nembd)
!..
!.. Local Scalars ..
      Integer :: iend, iprint, iq1, iq2, nclebd
      Integer :: ngmax, nrmax, nshlg, nshlr
      Integer :: lrecabmad, irec
!..
!.. Local Arrays ..
!.. Attention: Dimension LMXSPD*LMPOTD appears sometimes as NCLEB1
      Real (Kind=dp) :: avmad(lmpotd, lmpotd), bvmad(lmpotd)
      Real (Kind=dp) :: cleb(lmxspd*lmpotd)
      Real (Kind=dp) :: madelsmat(lmxspd, naezd, naezd)
      Real (Kind=dp) :: smat1(6, 6), smat2(6, 6)
      Real (Kind=dp) :: gn(3, nmaxd), rm(3, nmaxd)
      Integer :: icleb(lmxspd*lmpotd, 3)
      Integer :: nsg(ishld), nsr(ishld)
!..
!.. External subroutines
      External :: lattice3d, strmat, madelgaunt, madelcoef
! ......................................................................
      iprint = 0
      nclebd = lmxspd*lmpotd

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, '(79("="))')
      Write (1337, '(18X,A)') 'MADELUNG3D: setting bulk Madelung coefficients'
      Write (1337, '(79("="))')
      Write (1337, *)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ======================================================================
      Call lattice3d(alat, bravais, recbv, ngmax, nrmax, nshlg, nshlr, nsg, &
        nsr, gn, rm, rmax, gmax, iprint, nmaxd, ishld)

      Call strmat(alat, lpot, naez, ngmax, nrmax, nsg, nsr, nshlg, nshlr, gn, &
        rm, rbasis, madelsmat, volume0, iprint, lassld, lmxspd, naezd)
! ======================================================================

      lrecabmad = wlength*2*lmpotd*lmpotd + wlength*2*lmpotd
      write(*,*) lrecabmad, 2*lmpotd*lmpotd, 2*lmpotd
      Open (69, Access='direct', Recl=lrecabmad, File='abvmad.unformatted', &
        Form='unformatted')

! --> calculate the gaunt coefficients

      Call madelgaunt(lpot, yrg, wg, cleb, icleb, iend, lassld, nclebd)

! --> calculate the madelung coefficients to be used for VMAD
!     call MADELCOEF with first arg. .FALSE. = 3D case

      Do iq1 = 1, naez
        Do iq2 = 1, naez
          Call madelcoef(.False., lpot, avmad, bvmad, madelsmat(1,iq1,iq2), &
            cleb, icleb, iend, lpotd, lmpotd, lmxspd, nclebd)

          irec = iq2 + naez*(iq1-1)
          Write (69, Rec=irec) avmad, bvmad
!-----------------------------------------------------------------------
          If ((iq1<=6) .And. (iq2<=6)) Then
            smat1(iq1, iq2) = avmad(1, 1)
            smat2(iq1, iq2) = bvmad(1)
          End If
!-----------------------------------------------------------------------
        End Do
      End Do
      Close (69)

      If (iprint<1) Return
! ======================================================================

      Call madel3out(iprint, naez, lrecabmad, smat1, smat2, lmpotd)

    End Subroutine
