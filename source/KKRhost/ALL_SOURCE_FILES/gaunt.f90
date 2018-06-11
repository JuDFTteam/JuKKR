! ************************************************************************
    Subroutine gaunt(lmax, lpot, w, yr, cleb, loflm, icleb, iend, jend, ncleb, &
      lmaxd, lmgf0d, lmpotd)
      Use mod_datatypes, Only: dp
! ************************************************************************

!   - fills the array cleb with the gaunt coeffients ,i.e.
!      the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)
!      but only for lm2.le.lm1 and lm3>1
!   - calculate the pointer array jend  to project the indices
!      array cleb with the same lm3,l1,l2 values - because of
!      the special ordering of array cleb only the last index
!      has to be determined .
!     (the parameter n has to be chosen that l1+l2+l3 .lt. 2*n)
!     using gaussian quadrature as given by
!     m. abramowitz and i.a. stegun, handbook of mathematical functions,
!     nbs applied mathematics series 55 (1968), pages 887 and 916
!     m. weinert and e. wimmer
!     northwestern university march 1980

!     an index array -icleb- is used to save storage place .
!     fills the array loflm which is used to determine the
!     l-value of a given lm-value .
!     this subroutine has to be called only once !

!                               b.drittler   november 1987

!     modified gaunt coefficients are als calculated defined by
!     the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)*i**(l2-l1+l3)
!-----------------------------------------------------------------------

!---> attention : ncleb is an empirical factor - it has to be optimized

      Implicit None
!..
      Complex (Kind=dp) :: ci
      Parameter (ci=(0.0E0_dp,1.0E0_dp))
!..
      Integer :: lmpotd, lmgf0d, lmaxd, ncleb
!..
!.. Scalar Arguments ..
      Integer :: iend, lmax, lpot
!..
!.. Array Arguments ..
      Real (Kind=dp) :: cleb(ncleb, 2), w(*), yr(4*lmaxd, 0:4*lmaxd, 0:4*lmaxd &
        )
      Integer :: icleb(ncleb, 4), jend(lmpotd, 0:lmaxd, 0:lmaxd), loflm(*)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: clecg, factor, fci, s
      Integer :: i, j, l, l1, l1p, l2, l2p, l3, lm1, lm2, lm3, lm3p, lmpot, m, &
        m1, m1a, m1s, m2, m2a, m2s, m3, m3a, m3s
!..
!.. Intrinsic Functions ..
      Intrinsic :: abs, mod, real, sign
!..
!.. External Subroutines ..
      External :: rcstop
!..

      i = 1
      Do l = 0, 2*lmax
        Do m = -l, l
          loflm(i) = l
          i = i + 1
        End Do
      End Do

      icleb = 0
      cleb = 0E0_dp
      If (lpot==0) Then
        iend = 1
        icleb(1, 1) = (lmax+1)**2
        icleb(1, 3) = 1
      End If

      If (lpot/=0) Then

!---> set up of the gaunt coefficients with an index field

        i = 1
        Do l3 = 1, lpot
          Do m3 = -l3, l3

            Do l1 = 0, lmax
              Do l2 = 0, l1

                If (mod((l1+l2+l3),2)/=1 .And. (l1+l2-l3)>=0 .And. (l1-l2+l3) &
                  >=0 .And. (l2-l1+l3)>=0) Then

                  fci = real(ci**(l2-l1+l3))
                  Do m1 = -l1, l1
                    Do m2 = -l2, l2

!---> store only gaunt coeffients for lm2.le.lm1

                      lm1 = l1*(l1+1) + m1 + 1
                      lm2 = l2*(l2+1) + m2 + 1
                      If (lm2<=lm1) Then

                        m1s = sign(1, m1)
                        m2s = sign(1, m2)
                        m3s = sign(1, m3)

                        If (m1s*m2s*m3s>=0) Then

                          m1a = abs(m1)
                          m2a = abs(m2)
                          m3a = abs(m3)

                          factor = 0.0_dp

                          If (m1a+m2a==m3a) factor = factor + &
                            real(3*m3s+sign(1,-m3), kind=dp)/8.0E0_dp
                          If (m1a-m2a==m3a) factor = factor + &
                            real(m1s, kind=dp)/4.0E0_dp
                          If (m2a-m1a==m3a) factor = factor + &
                            real(m2s, kind=dp)/4.0E0_dp

                          If (factor/=0.0_dp) Then

                            If (m1s*m2s/=1 .Or. m2s*m3s/=1 .Or. m1s*m3s/=1) &
                              factor = -factor

                            s = 0.0_dp
                            Do j = 1, 4*lmaxd
                              s = s + w(j)*yr(j, l1, m1a)*yr(j, l2, m2a)*yr(j, &
                                l3, m3a)
                            End Do
                            clecg = s*factor
                            If (abs(clecg)>1.E-10_dp) Then
                              cleb(i, 1) = clecg
                              cleb(i, 2) = fci*clecg
                              icleb(i, 1) = lm1
                              icleb(i, 2) = lm2
                              icleb(i, 3) = l3*(l3+1) + m3 + 1
                              icleb(i, 4) = lm2*lmgf0d - (lm2*lm2-lm2)/2 + &
                                lm1 - lmgf0d
                              i = i + 1
                            End If

                          End If

                        End If

                      End If

                    End Do
                  End Do
                End If

              End Do
            End Do
          End Do
        End Do
        iend = i - 1
        If (ncleb<iend) Then
          Write (6, Fmt=100) ncleb, iend
          Call rcstop('33      ')

        Else

!---> set up of the pointer array jend,use explicitly
!     the ordering of the gaunt coeffients

          lmpot = (lpot+1)*(lpot+1)
          Do l1 = 0, lmax
            Do l2 = 0, l1
              Do lm3 = 2, lmpot
                jend(lm3, l1, l2) = 0
              End Do
            End Do
          End Do

          lm3 = icleb(1, 3)
          l1 = loflm(icleb(1,1))
          l2 = loflm(icleb(1,2))

          Do j = 2, iend
            lm3p = icleb(j, 3)
            l1p = loflm(icleb(j,1))
            l2p = loflm(icleb(j,2))

            If (lm3/=lm3p .Or. l1/=l1p .Or. l2/=l2p) Then
              jend(lm3, l1, l2) = j - 1
              lm3 = lm3p
              l1 = l1p
              l2 = l2p
            End If

          End Do
          jend(lm3, l1, l2) = iend


        End If

      End If



100   Format (13X, 'error stop in gaunt : dimension of NCLEB = ', I10, &
        ' too small ', /, 13X, 'change NCLEB to ', I6)
    End Subroutine
