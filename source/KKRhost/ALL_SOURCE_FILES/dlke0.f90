! 04.10.95 *************************************************************
    Subroutine dlke0(gllke, alat, naez, cls, nacls, naclsmax, rr, ezoa, atom, &
      bzkp, rcls, ginp)
      Use mod_datatypes, Only: dp
! **********************************************************************
      Implicit None
!     .. Parameters ..
      Include 'inc.p'
! *          function, set up in the spin-independent non-relativstic *
! *          (l,m_l)-representation                                   *
! *                                                                   *
! *********************************************************************
!     ..
!..
!.. Scalar Arguments ..
!..
!.. Array Arguments ..
      Integer :: lmax
      Parameter (lmax=lmaxd)
      Integer :: lmgf0d
      Parameter (lmgf0d=(lmax+1)**2)
      Integer :: almgf0
      Parameter (almgf0=lmgf0d*naezd)


      Real (Kind=dp) :: alat
      Integer :: naez, naclsmax
!..
!.. Local Scalars ..
!..
!.. Local Arrays ..
      Complex (Kind=dp) :: ginp(lmgf0d*naclsmax, lmgf0d, *), gllke(almgf0, *)
      Real (Kind=dp) :: bzkp(*), rcls(3, naclsd, *), rr(3, 0:nrd)
      Integer :: atom(naclsd, *), cls(*), ezoa(naclsd, *), nacls(*)
!..
!.. External Subroutines ..
      Integer :: i, ic, im, j, jn, m, n
!..
!.. Save statement ..
      Complex (Kind=dp) :: gllke1(almgf0, lmgf0d)
      Real (Kind=dp) :: kp(6)
!      write(6,*) '>>> DLKE0 : Fourier-transforms the ',
!     +           'GF of reference system'
      External :: cinit, dlke1
! ----------------------------------------------------------------------

      Save
!     .. External Functions ..
!     ..



      Logical :: opt
      External :: opt

      Call cinit(almgf0*almgf0, gllke(1,1))

      Do i = 1, naez


        kp(1) = bzkp(1)
        kp(2) = bzkp(2)
        kp(3) = bzkp(3)
        If (opt('COMPLEX ')) Then
          kp(4) = bzkp(4)
          kp(5) = bzkp(5)
          kp(6) = bzkp(6)
        End If

        ic = cls(i)
        Call dlke1(gllke1, alat, nacls, naclsmax, rr, ezoa(1,i), atom(1,i), &
          kp, ic, ginp(1,1,ic), rcls(1,1,ic))

        Do m = 1, lmgf0d
          im = (i-1)*lmgf0d + m
          Do jn = 1, lmgf0d*naez
            gllke(jn, im) = gllke(jn, im) + gllke1(jn, m)
          End Do
        End Do
! ----------------------------------------------------------------------

! -->   symmetrization
      End Do


      If (opt('symG(k) ')) Then



        Do i = 1, naez

          kp(1) = -bzkp(1)
          kp(2) = -bzkp(2)
          kp(3) = -bzkp(3)
          If (opt('COMPLEX ')) Then
            kp(4) = -bzkp(4)
            kp(5) = -bzkp(5)
            kp(6) = -bzkp(6)
          End If
! ----------------------------------------------------------------------
          ic = cls(i)
          Call dlke1(gllke1, alat, nacls, naclsmax, rr, ezoa(1,i), atom(1,i), &
            kp, ic, ginp(1,1,ic), rcls(1,1,ic))

          Do j = 1, naez
            Do m = 1, lmgf0d
              im = (i-1)*lmgf0d + m
              Do n = 1, lmgf0d
                jn = (j-1)*lmgf0d + n
                gllke(im, jn) = (gllke(im,jn)+gllke1(jn,m))/2.0E0_dp
              End Do
            End Do
          End Do

        End Do
! 04.10.95 *************************************************************
      End If
! **********************************************************************
!     .. Parameters ..
      Return
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
    End Subroutine
