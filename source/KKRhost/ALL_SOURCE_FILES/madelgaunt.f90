    Subroutine madelgaunt(lpot, yrg, wg, cleb, icleb, iend, lassld, nclebd)
      Use mod_datatypes, Only: dp
      Implicit None
!..
!.. Scalar arguments
      Integer :: lpot, iend
      Integer :: lassld, nclebd
!..
!.. Array arguments
!.. Attention: Dimension NCLEBD appears sometimes as NCLEB1
!..            an empirical factor - it has to be optimized
      Real (Kind=dp) :: yrg(lassld, 0:lassld, 0:lassld), wg(lassld)
      Real (Kind=dp) :: cleb(nclebd)
      Integer :: icleb(nclebd, 3)
!..
!.. Local scalars
      Real (Kind=dp) :: clecg, factor, s
      Integer :: i, j, l1, l2, l3, m1, m1a, m1s, m2, m2a, m2s, m3, m3a, m3s
!..
!.. Intrinsic functions
      Intrinsic :: abs, atan, real, sign

! --> set up of the gaunt coefficients with an index field
!     recognize that they are needed here only for l3=l1+l2

      If (2*lpot>lassld) Then
        Write (6, *) 'Dim ERROR in MADELGAUNT -- 2*LPOT > LASSLD', 2*lpot, &
          lassld
        Stop
      End If

      i = 1
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      Do l1 = 0, lpot
        Do l2 = 0, lpot
          l3 = l1 + l2
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
          Do m1 = -l1, l1
            Do m2 = -l2, l2
              Do m3 = -l3, l3
                m1s = sign(1, m1)
                m2s = sign(1, m2)
                m3s = sign(1, m3)
! **********************************************************************
                If (m1s*m2s*m3s>=0) Then
                  m1a = abs(m1)
                  m2a = abs(m2)
                  m3a = abs(m3)

                  factor = 0.0E0_dp
                  If (m1a+m2a==m3a) factor = factor + &
                    real(3*m3s+sign(1,-m3), kind=dp)/8.0E0_dp
                  If (m1a-m2a==m3a) factor = factor + &
                    real(m1s, kind=dp)/4.0E0_dp
                  If (m2a-m1a==m3a) factor = factor + &
                    real(m2s, kind=dp)/4.0E0_dp
! ======================================================================
                  If (factor/=0.0E0_dp) Then
                    If (m1s*m2s/=1 .Or. m2s*m3s/=1 .Or. m1s*m3s/=1) &
                      factor = -factor

                    s = 0.0E0_dp
                    Do j = 1, lassld
                      s = s + wg(j)*yrg(j, l1, m1a)*yrg(j, l2, m2a)*yrg(j, l3, &
                        m3a)
                    End Do

                    clecg = s*factor
! ----------------------------------------------------------------------
                    If (abs(clecg)>1.E-10_dp) Then
                      cleb(i) = clecg
                      icleb(i, 1) = l1*(l1+1) + m1 + 1
                      icleb(i, 2) = l2*(l2+1) + m2 + 1
                      icleb(i, 3) = l3*(l3+1) + m3 + 1
                      i = i + 1
                      If (i>nclebd) Then
                        Write (6, Fmt='(2I10)') i, nclebd
                        Stop ' Dim stop in MADELGAUNT '
                      End If
                    End If
! ----------------------------------------------------------------------
                  End If
! ======================================================================
                End If
! **********************************************************************
              End Do
            End Do
          End Do
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        End Do
      End Do
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      iend = i - 1
    End Subroutine
