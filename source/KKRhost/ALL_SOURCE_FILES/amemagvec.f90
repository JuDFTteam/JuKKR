    Subroutine amemagvec(irel, iprint, nkm, amemvec, ikmllim1, ikmllim2, &
      imkmtab, cgc, nlmax, nkmmax, nkmpmax, nmvecmax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *   calculate the angular matrix elements connected with           *
!   *                                                                  *
!   *   1:       < LAM | sigma(ipol) | LAM' >    spin moment           *
!   *   2:       < LAM |     l(ipol) | LAM' >    orbital moment        *
!   *   3:       < LAM |     T(ipol) | LAM' >    spin dipole moment    *
!   *   4:       < LAM |  B_hf(ipol) | LAM' >    hyperfine field       *
!   *                                                                  *
!   *   ipol= 1,2,3  ==  (+),(-),(z):                                  *
!   *                                                                  *
!   ********************************************************************
      Implicit Real *8(A-H, O-Z)

! Dummy arguments

      Integer :: iprint, irel, nkm, nkmmax, nkmpmax, nlmax, nmvecmax
      Real (Kind=dp) :: amemvec(nkmmax, nkmmax, 3, nmvecmax), cgc(nkmpmax, 2)
      Integer :: ikmllim1(nkmmax), ikmllim2(nkmmax), imkmtab(nkmmax)

! Local variables

      Character (Len=1) :: chpol(3)
      Real (Kind=dp) :: dble
      Integer :: i, ikm1, ikm2, imkm, imv, imvec, ipol, j1p05, j2p05, k, k1, &
        k2, kap1, kap2, l, l1, l2, lb1, lb2, m2, msm05, mue1m05, mue2m05, nk, &
        nmvec, ixm
      Integer :: iabs, nint
      Character (Len=20) :: str20
      Real (Kind=dp) :: sum, xj, xjm, xjp, xm, xynorm
      Character (Len=4) :: txtmvec(4)

      Data chpol/'+', '-', 'z'/
      Data txtmvec/'spin', 'orb ', 'T_z ', 'B_hf'/

      nk = 2*nlmax - 1
!     XYNORM = DSQRT(2.0D0)
      xynorm = 2.0E0_dp

      Call rinit(nkmmax*nkmmax*3*nmvecmax, amemvec)

      If (irel<=1) Return

! ----------------------------------------------------------------------
!     find the bounding indices  IKMLLIM1  and  IKMLLIM2  for IKM-loops
!     assuming that the matrix elements are diagonal with respect to l
!     this does not hold for B_hf for which there are l-(l+/-2)-terms
! ----------------------------------------------------------------------

      i = 0
      Do k = 1, nk
        l = k/2
        xjm = l - 0.5E0_dp
        xjp = l + 0.5E0_dp
        If (mod(k,2)==1) Then
          xj = l + 0.5E0_dp
        Else
          xj = l - 0.5E0_dp
        End If
!DO XM = -XJ, + XJ
        Do ixm = 1, 2*nint(xj) + 1
          xm = -xj + dble(ixm-1)
          i = i + 1
          ikmllim1(i) = nint(l*2*(xjm+0.5E0_dp)+1)
          ikmllim2(i) = nint(l*2*(xjp+0.5E0_dp)+2*xjp+1)
        End Do
      End Do

! ----------------------------------------------------------------------

      ikm1 = 0
      Do k1 = 1, nk
        l1 = k1/2
        If (mod(k1,2)==0) Then
          kap1 = l1
          lb1 = l1 - 1
        Else
          kap1 = -l1 - 1
          lb1 = l1 + 1
        End If
        j1p05 = iabs(kap1)

        Do mue1m05 = -j1p05, j1p05 - 1
          ikm1 = ikm1 + 1
          imkm = lb1*2*j1p05 + j1p05 + mue1m05 + 1
          imkmtab(ikm1) = imkm

          ikm2 = 0
          Do k2 = 1, nk
            l2 = k2/2
            If (mod(k2,2)==0) Then
              kap2 = l2
              lb2 = l2 - 1
            Else
              kap2 = -l2 - 1
              lb2 = l2 + 1
            End If
            j2p05 = iabs(kap2)

            Do mue2m05 = -j2p05, j2p05 - 1
              ikm2 = ikm2 + 1
! ----------------------------------------------------------------------
              If ((mue1m05-mue2m05)==+1) Then
                amemvec(ikm1, ikm2, 1, 1) = xynorm*cgc(ikm1, 2)*cgc(ikm2, 1)

                sum = 0E0_dp
                Do msm05 = -1, 0
                  m2 = mue2m05 - msm05
                  If (abs(m2)<=l2) sum = sum + cgc(ikm1, msm05+2)*cgc(ikm2, &
                    msm05+2)*sqrt(dble((l2-m2)*(l2+m2+1)))
                End Do
                amemvec(ikm1, ikm2, 1, 2) = sum
              End If

              If ((mue1m05-mue2m05)==-1) Then
                amemvec(ikm1, ikm2, 2, 1) = xynorm*cgc(ikm1, 1)*cgc(ikm2, 2)

                sum = 0E0_dp
                Do msm05 = -1, 0
                  m2 = mue2m05 - msm05
                  If (abs(m2)<=l2) sum = sum + cgc(ikm1, msm05+2)*cgc(ikm2, &
                    msm05+2)*sqrt(dble((l2+m2)*(l2-m2+1)))

                End Do
                amemvec(ikm1, ikm2, 2, 2) = sum
              End If

              If ((mue1m05-mue2m05)==0) Then
                amemvec(ikm1, ikm2, 3, 1) = cgc(ikm1, 2)*cgc(ikm2, 2) - &
                  cgc(ikm1, 1)*cgc(ikm2, 1)

                sum = 0E0_dp
                Do msm05 = -1, 0
                  m2 = mue2m05 - msm05
                  sum = sum + cgc(ikm1, msm05+2)*cgc(ikm2, msm05+2)*m2
                End Do
                amemvec(ikm1, ikm2, 3, 2) = sum
              End If

! ----------------------------------------------------------------------
            End Do
          End Do
        End Do
      End Do

      nmvec = 2
! ----------------------------------------------------------------------
      If (iprint<=90) Return

      Do imvec = 1, nmvec
        imv = min(imvec, 2)
        Do ipol = 1, 3
          str20 = 'A  ' // txtmvec(imv) // '  (' // chpol(ipol) // ')'
          Call rmatstr(str20, 12, amemvec(1,1,ipol,imvec), nkm, nkmmax, 3, 3, &
            1E-8_dp, 6)
        End Do
      End Do
    End Subroutine
