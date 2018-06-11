    Subroutine startldau(itrunldau, idoldau, kreadldau, lopt, ueff, jeff, &
      erefldau, natyp, nspin, wldau, uldau, phildau, irws, ntldau, itldau, &
      irmd, natypd, nspind, mmaxd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * Reads in LDA+U arrays from formatted file 'ldaupot'                *
! *                                                                    *
! **********************************************************************
      Implicit None
!..
      Integer :: irmd, mmaxd, natypd, nspind, irws(natypd)
!..
!.. Arguments ..
      Integer :: itrunldau, idoldau, kreadldau, natyp, nspin, ntldau
      Integer :: lopt(natypd), itldau(natypd)
      Real (Kind=dp) :: ueff(natypd), jeff(natypd), erefldau(natypd)
      Real (Kind=dp) :: wldau(mmaxd, mmaxd, nspind, natypd)
      Real (Kind=dp) :: uldau(mmaxd, mmaxd, mmaxd, mmaxd, natypd)
!      DOUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
      Complex (Kind=dp) :: phildau(irmd, natypd)
!.. 
!.. Locals ..
      Integer :: i1, im1, im3, is, it, ll
!     ..
! ----------------------------------------------------------------------


!      ALLOCATE( ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) )

      itrunldau = 0
      idoldau = 1
      ntldau = 0
      Do it = 1, natyp
        If (lopt(it)+1/=0) Then
          ntldau = ntldau + 1
          itldau(ntldau) = it
        End If
      End Do

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      Write (1337, '(79(1H=),/,27X,A,/, 79(1H=),/)') &
        'LDA+U: starting parameters'
      Write (1337, 100) natyp, ntldau
      Write (1337, 110)
      Write (1337, 120)
      Do it = 1, ntldau
        i1 = itldau(it)
        Write (1337, 130) i1, ueff(i1), jeff(i1), erefldau(i1)
      End Do
      Write (1337, 120)
      Write (1337, *)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

! -> read in LDA+U from file if available (KREADLDAU=1)

      Call rinit(mmaxd*mmaxd*nspind*natypd, wldau)
      Call cinit(irmd*natypd, phildau)
      If (kreadldau==1) Then
        Write (1337, 140)
        Call readldaupot(itrunldau, lopt, ueff, jeff, erefldau, natyp, wldau, &
          uldau, phildau, irws, ntldau, itldau, irmd, natypd, nspind, mmaxd)
      Else
        Write (1337, 150)
      End If

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      If (itrunldau/=0) Then
        Write (1337, 160) 'Coulomb matrix U(m1,m1,m3,m3)'
        Do it = 1, ntldau
          i1 = itldau(it)
          ll = lopt(i1)
          ll = min(3, ll)
          Write (1337, 170) i1
          Do im1 = 1, 2*ll + 1
            Write (1337, 180)(uldau(im1,im1,im3,im3,i1), im3=1, 2*ll+1)
          End Do
          Write (1337, *)
          If (it<ntldau) Write (1337, 190)
        End Do
        Write (1337, 160) 'Interaction potential W(m1,m2)'
        Do it = 1, ntldau
          i1 = itldau(it)
          ll = lopt(i1)
          ll = min(3, ll)
          Do is = 1, nspin
            Write (1337, 200) i1, is
            Do im1 = 1, 2*ll + 1
              Write (1337, 180)(wldau(im1,im3,is,i1), im3=1, 2*ll+1)
            End Do
            Write (1337, *)
          End Do
          If (it<ntldau) Write (1337, 190)
        End Do
        Write (1337, '(9X,60(1H-))')
      Else
        Call rinit(mmaxd*mmaxd*mmaxd*mmaxd*natypd, uldau)
        Call rinit(mmaxd*mmaxd*nspind*natypd, wldau)
        Call cinit(irmd*natypd, phildau)
      End If
      Write (1337, *)

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

100   Format (6X, 'Number of atoms ', '  in the u.c. :', I4, /, 24X, &
        'using LDA+U :', I4, /)
110   Format (9X, ' IT ', '   Ueff   ', '   Jeff   ', '   Eref   ', ' (Ry)')
120   Format (9X, 40('-'))
130   Format (9X, I3, 1X, 3F10.6)
140   Format (9X, 'Reading in LDA+U potential information (file ldaupot)', /)
150   Format (9X, 'LDA+U potential initialised (set to zero)')
160   Format (9X, 60('-'), /, 9X, A, /, 9X, 60('-'), /)
170   Format (9X, 'IT =', I3, /)
180   Format (9X, 7F10.6)
190   Format (11X, 58('~'))
200   Format (9X, 'IT =', I3, ' ISPIN =', I2, /)
    End Subroutine
