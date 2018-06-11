    Subroutine shape_corr(lpot, natyp, gsh, ilm_map, imaxsh, lmsp, ntcell, w, &
      yr, lassld, lmpotd, natypd, ngshd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *  Prepares shape corrections using gaussian quadrature as given by  *
! *  m. abramowitz and i.a. stegun, handbook of mathematical functions *
! *  nbs applied mathematics series 55 (1968), pages 887 and 916       *
! *                                                                    *
! *  the parameter LASSLD has to be chosen such that                   *
! *                        l1+l2+l3 .le. 2*LASSLD                      *
! *                                                                    *
! **********************************************************************

      Implicit None
!..
!.. Scalar Arguments ..
      Integer :: lassld, lmpotd, natypd, ngshd
      Integer :: lpot, natyp
!..
!.. Array Arguments ..
      Real (Kind=dp) :: gsh(*), w(lassld), yr(lassld, 0:lassld, 0:lassld)
      Integer :: ilm_map(ngshd, 3), imaxsh(0:lmpotd), lmsp(natypd, *), &
        ntcell(*)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: factor, gaunt, s
      Integer :: i, iat, icell, isum, j, l1, l2, l3, lm1, lm2, lm3, m1, m1a, &
        m1s, m2, m2a, m2s, m3, m3a, m3s
      Logical :: triangle
!..
!.. Intrinsic Functions ..
      Intrinsic :: abs, real, sign
!..
!.. External Subroutines ..
      External :: rcstop, triangle
!..

! -> set up of the gaunt coefficients with an index field
!    so that  c(lm,lm',lm'') is mapped to c(i)
      i = 1
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      Do l1 = 0, lpot
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        Do m1 = -l1, l1

          lm1 = l1*l1 + l1 + m1 + 1
          imaxsh(lm1-1) = i - 1
! llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
          Do l3 = 0, lpot*2
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
            Do m3 = -l3, l3

              lm3 = l3*l3 + l3 + m3 + 1
              isum = 0

              Do iat = 1, natyp
                icell = ntcell(iat)
!                      write(*,*) 'test icell=ntcell(iat) in shape_corr.f',
!      +                icell,iat
                isum = isum + lmsp(icell, lm3)
              End Do

! ======================================================================
              If (isum>0) Then
                Do l2 = 0, lpot
! ----------------------------------------------------------------------
                  If (triangle(l1,l2,l3)) Then
                    Do m2 = -l2, l2

                      lm2 = l2*l2 + l2 + m2 + 1

! -> use the m-conditions for the gaunt coefficients not to be 0

                      m1s = sign(1, m1)
                      m2s = sign(1, m2)
                      m3s = sign(1, m3)
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
! ......................................................................
                        If (factor/=0.0E0_dp) Then

                          If (m1s*m2s/=1 .Or. m2s*m3s/=1 .Or. m1s*m3s/=1) &
                            factor = -factor

                          s = 0.0E0_dp
                          Do j = 1, lassld
                            s = s + w(j)*yr(j, l1, m1a)*yr(j, l2, m2a)*yr(j, &
                              l3, m3a)
                          End Do

                          gaunt = s*factor
                          If (abs(gaunt)>1E-10_dp) Then
                            gsh(i) = gaunt
                            ilm_map(i, 1) = lm1
                            ilm_map(i, 2) = lm2
                            ilm_map(i, 3) = lm3
                            i = i + 1
                          End If
                        End If
! ......................................................................
                      End If
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                    End Do
                  End If
! ----------------------------------------------------------------------
                End Do
              End If
! ======================================================================
            End Do
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
          End Do
! llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
        End Do
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      End Do
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

      imaxsh(lm1) = i - 1
      Write (1337, Fmt=100) imaxsh(lm1), ngshd
      If (imaxsh(lm1)>ngshd) Call rcstop('SHAPE   ')

100   Format (' >>> SHAPE : IMAXSH(', I4, '),NGSHD :', 2I6)

    End Subroutine

    Function triangle(l1, l2, l3)
      Implicit None
      Integer :: l1, l2, l3
      Logical :: triangle
      Intrinsic :: mod
!     ..
      triangle = (l1>=abs(l3-l2)) .And. (l1<=(l3+l2)) .And. &
        (mod((l1+l2+l3),2)==0)
    End Function
