    Subroutine deciopt(alat, ins, krel, kvrel, kmrot, nspin, naez, lmmax, &
      bravais, tk, npol, npnt1, npnt2, npnt3, ez, ielast, kaoez, lefttinvll, &
      righttinvll, vacflag, nlbasis, nrbasis, cmomhost, vref, rmtref, nref, &
      refpot, lmaxd, lmgf0d, lmmaxd, lm2d, nembd1, iemxd, nspind, lmpotd, &
      natypd, irmd, ipand)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * This routine treats the DECIMATION case setting up the single-site *
! * (Delta t)^(-1) matrices and the charge moments of the host(s).     *
! *                                                                    *
! * This is realised in two ways:                                      *
! *      - either reading in the matrices (and moments - if SCFSTEPS   *
! *        is greater than 1) as written out in a previous (bulk) run  *
! *        -- DECIFILES token points to the files containing the nece- *
! *        ssary information                                           *
! *      - or reading in the self-consistent potential for each host   *
! *        and effectively calculating the matrices; the potential     *
! *        must have the specific format set in < OUTPOTHOST > routine *
! *        and the DECIPOTS token should point to the corresponding    *
! *        potential file(s)                                           *
! *                                                                    *
! * Notes:                                                             *
! *        - DECIFILES token is considered by default and sought first *
! *                          is requiring the same energy mesh for the *
! *                          system as for the host                    *
! *        - DECIPOTS token is not restrictive in this sense           *
! *                         however, it does not calculate charge mo-  *
! *                         ments -- does not work in SCF mode         *
! *                         is not dealing with CPA host so far        *
! *                         is not dealing with NON-SPHERICAL poten-   *
! *                         tials so far                               *
! *                                                                    *
! *                                     v.popescu - munich, Dec 04     *
! *                                                                    *
! **********************************************************************

      Implicit None
!..
!.. Scalar arguments
      Integer :: lmmaxd, nembd1, iemxd, nspind, lmpotd, natypd, ipand, irmd, &
        lmaxd
      Integer :: lm2d, nref, lmgf0d
      Integer :: ins, krel, kmrot, nspin, naez, lmmax, npol, npnt1, npnt2, &
        npnt3
      Integer :: ielast, nlbasis, nrbasis, kvrel
      Real (Kind=dp) :: alat, tk
!..
!.. Array arguments
      Integer :: kaoez(natypd, *), refpot(nembd1)
      Real (Kind=dp) :: bravais(3, 3), cmomhost(lmpotd, *)
      Real (Kind=dp) :: vref(*), rmtref(*)
      Complex (Kind=dp) :: lefttinvll(lmmaxd, lmmaxd, nembd1, nspind, iemxd), &
        righttinvll(lmmaxd, lmmaxd, nembd1, nspind, iemxd)
      Complex (Kind=dp) :: ez(iemxd)
      Logical :: vacflag(2)
!..
!.. Local scalars
      Integer :: ierror, il, ie, ispin, nspinso ! ruess: for tmat new solver
      Complex (Kind=dp) :: cfctor
      Character (Len=40) :: fileleft, fileright
      Character (Len=256) :: uio ! NCOLIO=256

!..                                  ! ruess: for NEWSOSOL running option
!.. External Functions ..
      Logical :: opt
      External :: opt

! ======================================================================
      Write (1337, '(79("="))')
      Write (1337, '(15X,A,/,79("="),/)') &
        'DECIOPT: reading left/right host decimation files'
      il = 1
      ierror = 0
      Call ioinput('DECIFILES       ', uio, il, 7, ierror)
! :::::::::::::::::::::::::::::::::::::::::::::::: decifiles (tmatrices)
      If (ierror==0) Then
        Read (Unit=uio, Fmt='(A40)') fileleft
        Call ioinput('DECIFILES       ', uio, il+1, 7, ierror)
        Read (Unit=uio, Fmt='(A40)') fileright
! ----------------------------------------------------------------------

! --> first call to read the header ( IE = 0 )

        ie = 0
        Call decimaread(ez, tk, npnt1, npnt2, npnt3, npol, nspin, &
          lefttinvll(1,1,1,1,1), righttinvll(1,1,1,1,1), vacflag, ie, nlbasis, &
          nrbasis, naez, kaoez, kmrot, ins, nspin, lmmax, ielast, fileleft, &
          fileright, krel, natypd, lmmaxd, nembd1)

! --> get the left and right host Delta_t matrices

        cfctor = alat/(8.E0_dp*atan(1.0E0_dp)) ! = ALAT/(2*PI)
        nspinso = nspin
        If (opt('NEWSOSOL')) nspinso = 1 ! ruess: only combined l-s index for newsolver
        Do ispin = 1, nspinso
          Do ie = 1, ielast
            Call decimaread(ez, tk, npnt1, npnt2, npnt3, npol, ispin, &
              lefttinvll(1,1,1,ispin,ie), righttinvll(1,1,1,ispin,ie), &
              vacflag, ie, nlbasis, nrbasis, naez, kaoez, kmrot, ins, nspin, &
              lmmax, ielast, fileleft, fileright, krel, natypd, lmmaxd, &
              nembd1)

! --> host matrices have been written out in true units
!     they are used in p.u. units (see kloopz) --> convert them here

            Call zscal(lmmaxd*lmmaxd*nembd1, cfctor, &
              lefttinvll(1,1,1,ispin,ie), 1)
            Call zscal(lmmaxd*lmmaxd*nembd1, cfctor, &
              righttinvll(1,1,1,ispin,ie), 1)
          End Do
        End Do

! --> get the left and right host charge moments
!     ( not needed in single-iteration mode calculations )

!fivos        IF ( SCFSTEPS.GT.1 )
        Call cmomsread(nlbasis, nrbasis, naez, cmomhost, vacflag, kaoez, &
          natypd, nembd1, lmpotd)
        Close (37)
        Close (38)
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      Else
! :::::::::::::::::::::::::::::::::::::::::::::::: decipots (calc tmats)
        ierror = 0
        Call ioinput('DECIPOTS        ', uio, il, 7, ierror)
        If (ierror/=0) Then
          Write (6, 100)
          Stop
        End If
        Read (Unit=uio, Fmt='(A40)') fileleft
        Call ioinput('DECIPOTS        ', uio, il+1, 7, ierror)
        Read (Unit=uio, Fmt='(A40)') fileright
        Call decitset(alat, bravais, ez, ielast, nlbasis, nrbasis, fileleft, &
          fileright, ins, kvrel, krel, nspin, kmrot, vref, rmtref, nref, &
          refpot, lefttinvll, righttinvll, vacflag, nembd1, iemxd, irmd, &
          ipand, lmaxd, lmgf0d, lmmaxd, lm2d, nspind)
      End If
! ======================================================================

100   Format (/, 6X, 'ERROR : Missing decimation files (t-mat or pot)', /, &
        14X, 'Please use one of the tokens DECIFILES/DECIPOTS', /)
    End Subroutine
