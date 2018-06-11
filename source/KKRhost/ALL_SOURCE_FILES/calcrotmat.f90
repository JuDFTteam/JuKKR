    Subroutine calcrotmat(nk, irel, alfdeg, betdeg, gamdeg, rot, fact, nkmmax)
!   ********************************************************************
!   *                                                                  *
!   *   SETS UP THE ROTATION-MATRICES FOR THE EULER ANGLES             *
!   *           ( ALFDEG, BETDEG, GAMDEG )                             *
!   *                                                                  *
!   *   SEE:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
!   *            EQS. (4.8), (4.12) AND (4.13)                         *
!   *                                                                  *
!   *   for IREL=0,1   NK == NL           non-relativistic (l,m_l)     *
!   *       IREL=3     NK == odd          relativistic (kappa,mue)     *
!   *                                                                  *
!   *   12/11/96  HE  deal with beta = 0                               *
!   ********************************************************************

      Use mod_datatypes, Only: dp
      Implicit None

      Complex (Kind=dp) :: ci, c0
      Parameter (ci=(0.0E0_dp,1.0E0_dp), c0=(0.0E0_dp,0.0E0_dp))
      Real (Kind=dp) :: pi
      Parameter (pi=3.141592653589793238462643E0_dp)

      Real (Kind=dp) :: num, msb05, msb05sq, msb05pw, j, m1, m2, rfac, dom, x, &
        cb05, cb05sq, alfdeg, betdeg, gamdeg, sum, cb05pw
      Real (Kind=dp) :: fact(0:100)

      Integer :: s, slow, shigh, off, nkmmax, irel, nk, i1, i2, k, l, im1, &
        im2, nmue
      Complex (Kind=dp) :: emim2a, emim1g, rot(nkmmax, nkmmax)

! INLINE FUNCTION    FACTORIAL FOR REAL ARGUMENT
      rfac(x) = fact(nint(x))

      If (irel==2) Call errortrap('calcrotmat', 12, 1)
      If (irel==3 .And. mod(nk,2)==0) Call errortrap('CALCROTMAT', 13, 1)

      Do i2 = 1, nkmmax
        Do i1 = 1, nkmmax
          rot(i1, i2) = c0
        End Do
      End Do

      cb05 = cos(betdeg*0.5E0_dp*pi/180.0E0_dp)
      cb05sq = cb05*cb05
      msb05 = -sin(betdeg*0.5E0_dp*pi/180.0E0_dp)
      msb05sq = msb05*msb05

      off = 0
      Do k = 1, nk
        If (irel<2) Then
          l = k - 1
          j = l
        Else
          l = k/2
          If (l*2==k) Then
            j = l - 0.5E0_dp
          Else
            j = l + 0.5E0_dp
          End If
        End If

        nmue = nint(2*j+1)

        Do im2 = 1, nmue
          m2 = -j + (im2-1.0E0_dp)
          emim2a = exp(-ci*m2*alfdeg*pi/180.0E0_dp)

          Do im1 = 1, nmue
            m1 = -j + (im1-1.0E0_dp)
            emim1g = exp(-ci*m1*gamdeg*pi/180.0E0_dp)

            If (abs(betdeg)<1E-8_dp) Then
              If (im1==im2) Then
                sum = 1.0E0_dp
              Else
                sum = 0.0E0_dp
              End If
            Else
              slow = max(0, nint(m1-m2))
              shigh = min(nint(j-m2), nint(j+m1))
              cb05pw = cb05**nint(2*j+m1-m2-2*slow+2)
              msb05pw = msb05**nint(m2-m1+2*slow-2)
              dom = (-1.0E0_dp)**(slow-1)*sqrt(rfac(j+m1)*rfac(j-m1)*rfac(j+m2 &
                )*rfac(j-m2))
              sum = 0.0E0_dp

              Do s = slow, shigh
                dom = -dom
                num = fact(s)*rfac(j-m2-s)*rfac(j+m1-s)*rfac(m2-m1+s)
                cb05pw = cb05pw/cb05sq
                msb05pw = msb05pw*msb05sq
                sum = sum + (dom/num)*cb05pw*msb05pw
              End Do
            End If

            rot(off+im2, off+im1) = emim1g*sum*emim2a
          End Do

        End Do

        off = off + nmue
      End Do

      Return
    End Subroutine
