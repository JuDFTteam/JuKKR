    Subroutine sphere_gga(lmax, yr, wtyr, rij, ijd, lmmaxd, thet, ylm, dylmt1, &
      dylmt2, dylmf1, dylmf2, dylmtf)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!     generate an angular mesh and spherical harmonics at those
!     mesh points. For an angular integration the weights are ge-
!     rated .

!     R. Zeller      Feb. 1996
!     Small change for GGA implementation
!     Nikos          Dec. 1996
!-----------------------------------------------------------------------
      Implicit None
!.. Scalar Arguments ..
      Integer :: ijd, lmax, lmmaxd
!..
!.. Local Scalars ..
      Real (Kind=dp) :: dx1, dx2, dx3, f0, pi, r, r1, r2, r3
      Integer :: ij, lm1
!..
!.. External Subroutines ..
      External :: cylm02, ymy
!..
!.. Array Arguments ..
      Real (Kind=dp) :: dylmf1(ijd, lmmaxd), dylmf2(ijd, lmmaxd), &
        dylmt1(ijd, lmmaxd), dylmt2(ijd, lmmaxd), dylmtf(ijd, lmmaxd), &
        rij(ijd, 3), thet(ijd), wtyr(ijd, *), ylm(ijd, lmmaxd), yr(ijd, *)
!..
!.. Local Arrays ..

      Real (Kind=dp) :: cosfi(ijd), cosx(ijd), fai(ijd), nd(3, 3), sinfi(ijd), &
        wght, y(1000)
!..
!.. Intrinsic Functions ..
      Intrinsic :: abs, acos, atan, cos, sin, sqrt
!     ..
      pi = 4.E0_dp*atan(1.E0_dp)
      Write (6, *) 'SPHERE for GGA: read LEBEDEV mesh'
      If (ijd>1000) Stop 'SPHERE'


      Do ij = 1, ijd
        Call lebedev(ij, r1, r2, r3, wght)

!      make a small rotation

        f0 = 0.08E0_dp
        nd(1, 1) = cos(f0)
        nd(1, 2) = 0E0_dp
        nd(1, 3) = sin(f0)
        nd(2, 1) = 0E0_dp
        nd(2, 2) = 1E0_dp
        nd(2, 3) = 0E0_dp
        nd(3, 1) = -sin(f0)
        nd(3, 2) = 0E0_dp
        nd(3, 3) = cos(f0)

        dx1 = nd(1, 1)*r1 + nd(2, 1)*r2 + nd(3, 1)*r3
        dx2 = nd(1, 2)*r1 + nd(2, 2)*r2 + nd(3, 2)*r3
        dx3 = nd(1, 3)*r1 + nd(2, 3)*r2 + nd(3, 3)*r3

        r1 = dx1
        r2 = dx2
        r3 = dx3

        rij(ij, 1) = r1
        rij(ij, 2) = r2
        rij(ij, 3) = r3

        Call ymy(r1, r2, r3, r, y, lmax)
        Do lm1 = 1, (lmax+1)**2
          yr(ij, lm1) = y(lm1)
        End Do

!---> multiply the spherical harmonics with the weights

        Do lm1 = 1, (lmax+1)**2
          wtyr(ij, lm1) = yr(ij, lm1)*wght*pi*4.E0_dp
        End Do

!---> produce what is needed for GGA

        cosx(ij) = r3
        If (abs(r3)/=1.E0_dp) Then
          cosfi(ij) = r1/sqrt(1.E0_dp-r3*r3)
          sinfi(ij) = r2/sqrt(1.E0_dp-r3*r3)
          If (abs(cosfi(ij))>1.E0_dp) cosfi(ij) = cosfi(ij)/abs(cosfi(ij))
          If (abs(sinfi(ij))>1.E0_dp) sinfi(ij) = sinfi(ij)/abs(sinfi(ij))
          fai(ij) = acos(cosfi(ij))
        Else If (sinfi(ij)==0.E0_dp) Then
          fai(ij) = pi/2.E0_dp
        Else
          cosfi(ij) = r1
          sinfi(ij) = r2
          If (abs(cosfi(ij))>1.E0_dp) cosfi(ij) = cosfi(ij)/abs(cosfi(ij))
        End If
        fai(ij) = acos(cosfi(ij))
      End Do

      Call cylm02(lmax, cosx, fai, 2*lmax+1, lmmaxd, thet, ylm, dylmt1, &
        dylmt2, dylmf1, dylmf2, dylmtf)

    End Subroutine
