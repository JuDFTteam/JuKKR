    Subroutine madelung2d(lpot, yrg, wg, naez, alat, vol, bravais, recbv, &
      rbasis, rmax, gmax, nlbasis, nleft, zperleft, tleft, nrbasis, nright, &
      zperight, tright, lmxspd, lassld, lpotd, lmpotd, nmaxd, ishld, nembd1, &
      wlength)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * This subroutine calculates the Madelung potential coefficients     *
! * in the 2D case and stores them in the DA-file abmad.unformatted    *
! * For each layer in the slab, the summation is split into three      *
! * parts (see also VINTERFACE):                                       *
! * within the slab, over the NLEFT*NLBASIS left host sites and over   *
! * the NRIGHT*NRBASIS right host sites, the last two steps only in    *
! * case of decimation run                                             *
! *                                                                    *
! * all positions must be scaled with ALAT to get them correct         *
! * (done in EWALD2D)                                                  *
! *                                                                    *
! * The record index is:                                               *
! *   (IQ1-1)*NAEZ + IQ2                 for (IQ1,IQ2) within the slab *
! *   NAEZ*NAEZ + (IQ1-1)*NLEFT*NLBASIS  for (IQ1,(IL,IBL)), IQ1 in    *
! *                  + (IL-1)*NLEFT+IBL  slab, (IL,IBL) in the left    *
! *   NAEZ*NAEZ + NAEZ*NLEFT*NLBASIS                                   *
! *             + (IQ1-1)*NRIGHT*NRBASIS for (IQ1,(IR,IBR)), IQ1 in    *
! *             + (IR-1)*NRIGHT+IBR      slab, (IR,IBR) in the right   *
! *                                                                    *
! **********************************************************************
      Implicit None
!..
!.. Scalar Arguments ..
      Integer :: lpot, naez, wlength
      Integer :: nlbasis, nleft, nrbasis, nright
      Integer :: lassld, lpotd, lmpotd, lmxspd, nmaxd, ishld, nembd1
      Real (Kind=dp) :: alat, vol, rmax, gmax
!..
!.. Array Arguments ..
      Real (Kind=dp) :: bravais(3, 3), recbv(3, 3)
      Real (Kind=dp) :: yrg(lassld, 0:lassld, 0:lassld), wg(lassld)
      Real (Kind=dp) :: rbasis(3, *)
      Real (Kind=dp) :: zperight(3), zperleft(3)
      Real (Kind=dp) :: tleft(3, nembd1), tright(3, nembd1)
!..
!.. Local Scalars ..
      Integer :: iq1, iq2, iend, nclebd, iprint
      Integer :: i, ib, ih, ileft, iright
      Integer :: lrecamad, irec, nleftoff, nrightoff, nleftall, nrightall
      Integer :: ngmax, nrmax, nshlg, nshlr
      Logical :: opt
!..
!.. Local Arrays ..
!.. Attention: LMXSPD*LMPOTD appears as NCLEB1 in other routines
      Real (Kind=dp) :: cleb(lmxspd*lmpotd)
      Real (Kind=dp) :: bm(lmpotd), vec2(3), sum(lmxspd)
      Real (Kind=dp) :: gn2(2, nmaxd), rm2(2, nmaxd)
      Real (Kind=dp) :: avmad(lmpotd, lmpotd)
      Integer :: nsg(ishld), nsr(ishld)
      Integer :: icleb(lmxspd*lmpotd, 3)
!..
!.. External Functions/Subroutines
! ......................................................................
      iprint = 0
      nclebd = lmxspd*lmpotd

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      Write (1337, '(79(1H=))')
      Write (1337, '(18X,A)') 'MADELUNG2D: setting 2D Madelung coefficients'
      Write (1337, '(79(1H=))')
      Write (1337, *)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ======================================================================
      Call lattice2d(alat, bravais, recbv, ngmax, nrmax, nshlg, nshlr, nsg, &
        nsr, gn2, rm2, rmax, gmax, iprint, nmaxd, ishld)
! ======================================================================

      lrecamad = wlength*2*lmpotd*lmpotd
      Open (69, Access='direct', Recl=lrecamad, File='avmad.unformatted', &
        Form='unformatted')

! --> calculate the gaunt coefs

      Call madelgaunt(lpot, yrg, wg, cleb, icleb, iend, lassld, nclebd)

! --> calculate the madelung coefficients to be used for VMAD

! **********************************************************************
! ********************************************** loop over atoms in slab
      Do iq1 = 1, naez

!        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                     1.  Summation in all layers in the slab
!        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! ++++++++++++++++++++++++++++++++ loop over all other sites in the slab
        Do iq2 = 1, naez

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
          If (iq1==1 .And. iq2==1) Then
            Write (1337, '(5X,2A,/)') &
              '< EWALD2D > : calculating 2D-lattice sums ', 'inside the slab'
            If (iprint>=2) Write (1337, 100)
          End If
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


! make ewald sumation in plane and inverse space
! sum if rz<>0 (out of plane)

!             WRITE(99,*) 'Layer pair:',IQ1,IQ2
          Call ewald2d(lpot, alat, rbasis(1,iq1), rbasis(1,iq2), iq1, iq2, &
            rm2, nrmax, nshlr, nsr, gn2, ngmax, nshlg, nsg, sum, vol, lassld, &
            lmxspd)
!             WRITE(99,*) 'SUM: ',SUM

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
          If (iprint>=2) Then
            Write (1337, 110) iq1, iq2, sum(1)
            If (iq2==naez .And. iq1/=naez) Write (1337, '(20X,20(1H-))')
          End If
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

          Call madelcoef(.True., lpot, avmad, bm, sum, cleb, icleb, iend, &
            lpotd, lmpotd, lmxspd, nclebd)

          irec = iq2 + naez*(iq1-1)
          Write (69, Rec=irec) avmad
        End Do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      End Do
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      If (iprint>=2) Write (1337, '(18X,22(1H-),/)')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
! ********************************************** loop over atoms in slab

! ######################################################################
      If (opt('DECIMATE')) Then

        nleftoff = naez*naez ! record offsets
        nrightoff = nleftoff + naez*nleft*nlbasis ! left and right
        nleftall = nleft*nlbasis
        nrightall = nright*nrbasis

! ********************************************** loop over atoms in slab
        Do iq1 = 1, naez
!     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                     2.  Summation in the LEFT bulk side
!     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
          ileft = 0
! ++++++++++++++++++++++++++++++++ loop over all sites in the left host
          Do ih = 1, nleft
            Do ib = 1, nlbasis
              Do i = 1, 3
                vec2(i) = (tleft(i,ib)+(ih-1)*zperleft(i))
              End Do
              ileft = ileft + 1

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
              If (iq1==1 .And. ileft==1) Then
                Write (1337, '(5X,2A,/)') &
                  '< EWALD2D > : calculating 2D-lattice sums ', &
                  'slab - left host'
                If (iprint>=2) Write (1337, 100)
              End If
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


!-->  make ewald sumation for m= 0 l<5 rz=0 (in plane) and
!     Inverse space sum if rz<>0 (out of plane)

              Call ewald2d(lpot, alat, rbasis(1,iq1), vec2, iq1, ih, rm2, &
                nrmax, nshlr, nsr, gn2, ngmax, nshlg, nsg, sum, vol, lassld, &
                lmxspd)

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
              If (iprint>=2) Then
                Write (1337, 110) iq1, ileft, sum(1)
                If (ileft==nleftall .And. iq1/=naez) Write (1337, &
                  '(20X,20(1H-))')
              End If
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

              Call madelcoef(.True., lpot, avmad, bm, sum, cleb, icleb, iend, &
                lpotd, lmpotd, lmxspd, nclebd)

              irec = ileft + nleftall*(iq1-1) + nleftoff
              Write (69, Rec=irec) avmad
            End Do ! ib loop in left host basis
          End Do ! ih loop in layers to get convergence
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          If (ileft/=nleftall) Then
            Write (6, *) ' < MADELUNG2D > : index error ', &
              'ILEFT <> NLEFT*NLBASIS'
            Stop
          End If
        End Do ! ILAY1 loop
! **********************************************************************

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        If (iprint>=2) Write (1337, '(18X,22(1H-),/)')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

! ********************************************** loop over atoms in slab
        Do iq1 = 1, naez
!        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                     3.  Summation in the RIGHT bulk side
!        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
! ++++++++++++++++++++++++++++++++ loop over all sites in the right host
          iright = 0
          Do ih = 1, nright
            Do ib = 1, nrbasis
              Do i = 1, 3
                vec2(i) = (tright(i,ib)+(ih-1)*zperight(i))
              End Do
              iright = iright + 1

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
              If (iq1==1 .And. iright==1) Then
                Write (1337, '(5X,2A,/)') &
                  '< EWALD2D > : calculating 2D-lattice sums ', &
                  'slab - right host'
                If (iprint>=2) Write (1337, 100)
              End If
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

!-->  make ewald sumation (in plane) and
!     Inverse space sum if rz<>0 (out of plane)

              Call ewald2d(lpot, alat, rbasis(1,iq1), vec2, iq1, ih, rm2, &
                nrmax, nshlr, nsr, gn2, ngmax, nshlg, nsg, sum, vol, lassld, &
                lmxspd)

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
              If (iprint>=2) Then
                Write (1337, 110) iq1, iright, sum(1)
                If (iright==nrightall .And. iq1/=naez) Write (1337, &
                  '(20X,20(1H-))')
              End If
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

              Call madelcoef(.True., lpot, avmad, bm, sum, cleb, icleb, iend, &
                lpotd, lmpotd, lmxspd, nclebd)

              irec = iright + nrightall*(iq1-1) + nrightoff
              Write (69, Rec=irec) avmad
            End Do ! ib loop in right host basis
          End Do ! ih loop in layers to get convergence
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          If (iright/=nrightall) Then
            Write (6, *) ' < MADELUNG2D > : index error ', &
              'IRIGHT <> NRIGHT*NRBASIS'
            Stop
          End If
        End Do ! ILAY1 loop
! **********************************************************************

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        If (iprint>=2) Write (1337, '(18X,22(1H-),/)')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      End If
! ######################################################################
      Close (69)

      If (iprint<1) Return
! ======================================================================

      Call madel2out(iprint, naez, lrecamad, lmpotd, nleftoff, nrightoff, &
        nleftall, nrightall)

100   Format (8X, '2D Lattice sum (LMXSP = 1)', /, 18X, '  IQ1  IQ2  SUM', /, &
        18X, 23('-'))
110   Format (18X, 2I5, D12.4)
    End Subroutine
