    Subroutine impcoefs(natomimp, naez, atomimp, rclsimp, nshell, nsh1, nsh2, &
      ratom, nsymat, isymindex, rotname, hostimp, natypd, lmaxd, nsheld, &
      nsize)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! * Writes out the auxiliary file impurity.coefs which is needed for   *
! * impurity calculations                                              *
! * Sets up the array HOSTIMP -- also needed for impurity case         *
! *                                adopted from N. Papanikolaou        *
! **********************************************************************

      Implicit None
!.. 
!.. Scalar arguments
      Integer :: lmaxd, naez, natomimp, natypd, nsheld, nsize, nsymat
!..
!.. Array arguments
      Integer :: atomimp(*), hostimp(0:natypd), isymindex(*), nsh1(*), &
        nsh2(*), nshell(0:nsheld)
      Real (Kind=dp) :: ratom(3, nsheld), rclsimp(3, *)
      Character (Len=10) :: rotname(*)
!..
!.. Local scalars
      Integer :: ai, i, ii, j, nb, ndim, nhost, nrep, ns
      Real (Kind=dp) :: r1
!..
!.. Local arrays
      Integer :: imphost(naez), nshout(natomimp)
      Logical :: exist(naez)
!     ..

! -->  shells around atom icc are prepared for storing the
!      cluster-gf in subroutine kkrmat in GMATLL(LMMAXD,LMMAXD,*)

      Do i = 1, naez
        exist(i) = .False.
      End Do
      Do i = 1, natomimp
        exist(atomimp(i)) = .True.
      End Do

      nhost = 0
      Do i = 1, naez
        imphost(i) = 0
        If (exist(i)) Then
          nhost = nhost + 1
          imphost(i) = nhost
          hostimp(nhost) = i
        End If
      End Do
      hostimp(0) = nhost
      If (nhost/=naez) Write (6, 100)

      nrep = 1
      ndim = 1

      Do i = 1, natomimp
        nshout(i) = 1
      End Do

      Open (58, File='impurity.coefs', Form='FORMATTED')
      Write (58, 110) nrep, natomimp, lmaxd, natomimp, &
        (nshout(i), i=1, natomimp)
      Write (58, 120)
!-----------------------------------------------------------------------
      Do i = 1, natomimp

        r1 = sqrt(rclsimp(1,i)**2+rclsimp(2,i)**2+rclsimp(3,i)**2)

        If (naez==nhost) Then
          ai = atomimp(i)
        Else
          ai = imphost(atomimp(i))
        End If

        Write (58, 130)(rclsimp(j,i), j=1, 3), ai, i, i, r1, atomimp(i)
      End Do
!-----------------------------------------------------------------------

      nb = 0
      Do ns = 1, nshell(0)
        nb = nb + nshell(ns)
      End Do

      Write (58, 200) nsize, nb
      Write (58, 110) ndim
      Write (58, 140) nhost
      Write (58, 150)(hostimp(i), i=1, nhost)
      Write (58, 160) nsymat
      Write (58, 170)(rotname(isymindex(i)), i=1, nsymat)
      Write (58, 180)
      Write (58, 110) nshell(0)
      Write (58, 190)(ns, nsh1(ns), nsh2(ns), (ratom(ii, &
        ns),ii=1,3), nshell(ns), sqrt(ratom(1,ns)**2+ratom(2,ns)**2+ratom(3, &
        ns)**2), ns=1, nshell(0))
      Close (58)
! ======================================================================
100   Format (8X, 'WARNING: Some host atoms are missing in the ', &
        'impurity cluster', /, 8X, '         Indexing will be changed. Check ' &
        , 'impurity.coefs file?', /)
110   Format (11I5)
120   Format ('     Position of Impurity            Host Imp Shell', &
        '   Dist     Host id in Bulk')
130   Format (3F12.8, I4, I4, I5, F10.6, I5)
140   Format ('Host order, no of host sites: ', I5)
150   Format (12I4)
160   Format (I5, '    Symmetries for the Bulk')
170   Format (5A10)
180   Format ('Shells in the reduced format')
190   Format (3I5, 3F12.8, I8, F10.5)
200   Format (11I20)
    End Subroutine
