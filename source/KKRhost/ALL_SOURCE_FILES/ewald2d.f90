    Subroutine ewald2d(lpot, alat, vec1, vec2, iq1, iq2, rm2, nrmax, nshlr, &
      nsr, gn2, ngmax, nshlg, nsg, sum2d, vol, lassld, lmxspd)
      Use mod_datatypes, Only: dp
! **********************************************************************
! *                                                                    *
! *   calculation of lattice sums for l .le. 2*lpot :                  *
! *                                                                    *
! *                    ylm( q(i) - q(j) + rm )                         *
! *         sum      ===========================                       *
! *                  | q(i) - q(j) + rm |**(l+1)                       *
! *                                                                    *
! *          - summed over all 2D lattice vectors rm  -                *
! *                                                                    *
! *   ylm       : real spherical harmic to given l,m                   *
! *                                                                    *
! *   The sum is done different in the plane (qi-qj)z = 0              *
! *   and out of the plane. In plane an Ewald procedure similar        *
! *   to the 3d is used and we perform 2 sums (real and reciprocal)    *
! *   the l= 2,4 m=0 terms are calculated with a different method      *
! *                                                                    *
! *   The l=0 term is calculated with a extra factror sqrt(4*pi) this  *
! *   is for transparency reasons (so that the correction terms        *
! *   *r=0,g=0* can be followed in the program)                        *
! *   Literature : lm = (0,0), (1,0) terms :PRB 40, 12164 (1989)       *
! *                                         PRB 49, 2721 (1994)        *
! *                                         PRB 47, 16525 (1993)       *
! *                                         Ziman , p.39-40            *
! *       l=2,4 (m=0) terms are done with recursive diferentiation     *
! *                                                  v. 16.8.99        *
! *       The l multipoles are treated using the expansion             *
! *       for a complex plane wave.                                    *
! *       eq.(33) , M. Weinert, J. Math Phys. 22, 2439 (1981)          *
! *                                                                    *
! *   Final version : 11.01.2000   (No direct sum needed everything    *
! *                                 is done with Ewald method)         *
! *   Programmed by N. Papanikolaou                                    *
! *                                                                    *
! **********************************************************************
      Use mod_datatypes, Only: dp
      Implicit None
!..
!.. Scalar arguments ..
      Integer :: lassld, lmxspd
! parameter (lassld replaces old l2maxd=2*LPOTD=4*lmaxd)
      Integer :: lpot, nrmax, nshlr, ngmax, nshlg, iq1, iq2
      Real (Kind=dp) :: alat, vol
!..
!.. Array arguments ..
      Real (Kind=dp) :: vec1(3), vec2(3)
      Real (Kind=dp) :: rm2(2, *), gn2(2, *), sum2d(lmxspd)
      Integer :: nsr(*), nsg(*)
!..
!.. Local scalars ..
      Real (Kind=dp) :: alpha, beta, bound, con, con1, dot1, dq1, dq2, dq3, &
        dqdotg, expbsq, fpi, g1, g2, g3, ga, lamda, pi, r, r0, r1, r2, r3, &
        rfac, s, signrz, stest0, tpi
      Complex (Kind=dp) :: aprefmm, aprefpp, bfac, ci, factexp, simag
      Real (Kind=dp) :: derfc, exponent, crit
      Integer :: i, im, ir, l, lm, lmax, lmmax, m, icall, ngmax1
!..
!.. Local arrays ..
      Real (Kind=dp) :: dfac(0:2*lassld+1), g(0:lassld), gal(0:lassld), &
        gi(0:4), gr(0:4), pref0(0:lassld), signrzl(0:lassld), ylm(lmxspd), &
        ylmpref(0:lassld), ylmpref1(0:lassld, 0:lassld)
      Complex (Kind=dp) :: cim(0:lassld), exponl(0:lassld), pref2(lassld), &
        s0(lmxspd), stest(lmxspd), stestnew(lmxspd)
!..
!.. External subroutines ..
      External :: fplaneg, fplaner
!..
!.. Intrinsic functions ..
      Intrinsic :: abs, atan, real, exp, sqrt
!..
!.. Data statements
      Data icall/0/
!..
!.. Save statements
      Save :: icall, ci, bound, pi, fpi, tpi
!..................................................................
!
      icall = 1
!ICALL = ICALL + 1
      If (icall==1) Then
        ci = (0.0E0_dp, 1.0E0_dp)
        bound = 1.0E-9_dp
        pi = 4.0E0_dp*atan(1.0E0_dp)
        fpi = 4.0E0_dp*pi
        tpi = 2.0E0_dp*pi
      End If
!
! Factorial
!
      dfac(0) = 1
      Do l = 1, 2*lassld + 1
        dfac(l) = dfac(l-1)*l
      End Do
      Do l = 0, lassld
        pref0(l) = 0E0_dp
      End Do

      lmax = 2*lpot
      lmmax = (lmax+1)**2
!
      pref0(2) = sqrt(5E0_dp/pi)/2E0_dp/2E0_dp
      pref0(4) = 3E0_dp*sqrt(9E0_dp/pi)/16E0_dp/9E0_dp
!
!choose proper splitting parameter
!
      lamda = sqrt(pi)/alat
!
      dq1 = (vec2(1)-vec1(1))*alat ! SCALE WITH ALAT
      dq2 = (vec2(2)-vec1(2))*alat
      dq3 = (vec2(3)-vec1(3))*alat
!
!Initialise
!
      Do lm = 1, lmmax
        stest(lm) = 0E0_dp
        stestnew(lm) = 0E0_dp
      End Do
!
!Add correction if rz = 0
!
      If (abs(dq3)<1E-6_dp) Then
        stest(1) = stest(1) - 2E0_dp*lamda/sqrt(pi) - &
          2E0_dp*sqrt(pi)/lamda/vol
        stestnew(1) = stestnew(1) - 2E0_dp*lamda/sqrt(pi) - &
          2E0_dp*sqrt(pi)/lamda/vol
!
        If ((dq1*dq1+dq2*dq2)>1E-6_dp) Then
          stest(1) = stest(1) + 2E0_dp*lamda/sqrt(pi)
          stestnew(1) = stestnew(1) + 2E0_dp*lamda/sqrt(pi)
        End If
      Else
!
!Add correction if rz<> 0
!
        stest(1) = stest(1) - abs(dq3)*fpi/2E0_dp/vol
        stest(3) = stest(3) - abs(dq3)/dq3*sqrt(3E0_dp*fpi)/2E0_dp/vol ! -d/dz
! the correction for higher l vanishes...
      End If
!
! **********************************************************************
! ******************    I N-P L A N E      M = 0  **********************
! **********************************************************************
      If (abs(dq3)<1E-6_dp) Then
!
!Real space sum
!
!==================================================================
        Do ir = 1, nrmax
          r1 = dq1 - rm2(1, ir)
          r2 = dq2 - rm2(2, ir)
          r3 = 0E0_dp
          r = sqrt(r1*r1+r2*r2)
!------------------------------------------------------------------
          If (r>1E-8_dp) Then
            alpha = lamda*r
            Call fplaner(alpha, gr, r)
            Do l = 0, 4
              lm = l*(l+1) + 1 ! m =0
              stest(lm) = stest(lm) + gr(l)
            End Do
            Call ymy(r1, r2, r3, r0, ylm, lassld)
            Call gamfc(alpha, g, lassld, r)
            ylm(1) = 1E0_dp ! just definition matter
!
            Do l = 0, lmax
              rfac = g(l)/sqrt(pi)
              Do m = -l, l
                lm = l*(l+1) + m + 1
                stestnew(lm) = stestnew(lm) + ylm(lm)*rfac
              End Do
            End Do
!
            If (ir==(nrmax-nsr(nshlr))) Then
!keep the value before the last shell to test convergence
              Do l = 0, lmax
                Do m = -l, l
                  lm = l*(l+1) + m + 1
                  s0(lm) = stestnew(lm)
                End Do
              End Do
            End If
          End If ! r <> 0
!------------------------------------------------------------------
        End Do ! ir loop
!==================================================================
!
!Check convergence
!
        s = 0E0_dp
!==================================================================
        Do l = 0, lmax
          Do m = -l, l
            lm = l*(l+1) + m + 1
            stest0 = abs(s0(lm)-stestnew(lm))
            If (s<stest0) s = stest0
          End Do
        End Do
!==================================================================
        If (s>bound) Write (1337, Fmt=100) abs(s), bound, iq1, iq2
!
!Sum in reciprocal lattice
!
        con = fpi/2E0_dp/vol
! Find an upper cutoff for G-vectors
        i = 1
        exponent = 1.E0_dp
        crit = 8.E0_dp
        Do While (exponent<crit .And. i<ngmax)
          i = i + 1
          g1 = gn2(1, i) ! G vectors are assumed to be sorted with increasing length
          g2 = gn2(2, i)
          ga = sqrt(g1*g1+g2*g2)
          exponent = ga/2.E0_dp/lamda ! If EXPONENT>8., then ERFC(EXPONENT) and EXP(-EXPONENT**2) 
! (see sub. fplaneg) are below is 1E-27, considered negligible.
        End Do
        ngmax1 = i
!     WRITE(99,FMT='(A7,I8,2E10.2)') 'NGMAX1:',NGMAX1,
! &                EXP(-EXPONENT**2),ERFC(EXPONENT)

        If (ngmax1>ngmax) Stop 'ewald2d: 1: NGMAX1.GT.NGMAX' ! should never occur
!==================================================================
        Do im = 1, ngmax1
          g1 = gn2(1, im)
          g2 = gn2(2, im)
          g3 = 0E0_dp
          ga = sqrt(g1*g1+g2*g2)
!------------------------------------------------------------------
          dot1 = dq1*g1 + dq2*g2
          Call fplaneg(lamda, gi, pref0, lassld, ga, vol)
          simag = exp(ci*dot1)
          Do l = 0, 4
            lm = l*(l+1) + 1
            stest(lm) = stest(lm) + gi(l)*simag
          End Do
!------------------------------------------------------------------
          If (ga>1E-6_dp) Then
            Call ymy(g1, g2, g3, ga, ylm, lassld)
            beta = ga/lamda
            expbsq = derfc(beta/2E0_dp)
!
            bfac = con*simag*expbsq
            stestnew(1) = stestnew(1) + bfac/ga
!
            Do l = 0, lmax
              If (l/=0) Then
                Do m = -l, l
                  lm = l*(l+1) + m + 1
                  stestnew(lm) = stestnew(lm) + ylm(lm)*bfac*ga**(l-1)
                End Do
              End If
              bfac = bfac/ci/real(2*l+1, kind=dp)
            End Do
          End If
!------------------------------------------------------------------
          If (im==(ngmax1-nsg(nshlg))) Then
!keep the value before the last shell to test convergence
            Do lm = 1, lmmax
              s0(lm) = stestnew(lm)
            End Do
          End If
        End Do
!==================================================================
!
!Check convergence
!
        Do lm = 1, lmmax
          stest0 = abs(s0(lm)-stestnew(lm))
          If (s<stest0) s = stest0
        End Do
!
!Correction due to r=0 term only for DRn = 0
!
!------------------------------------------------------------------
        If ((dq1*dq1+dq2*dq2)<1E-6_dp) Then
          Do l = 2, 4, 2
            pref0(l) = pref0(l)*dfac(l)/dfac(l/2)/(l+1)
          End Do
!
          i = 1
          Do l = 2, 4, 2
            i = i + 1
            lm = l*(l+1) + 1
            stest(lm) = stest(lm) + (-1)**i*2E0_dp/sqrt(pi)*pref0(l)*lamda**(l &
              +1)
          End Do
        End If
!------------------------------------------------------------------
        Do l = 2, 4, 2
          lm = l*(l+1) + 1
          stestnew(lm) = stest(lm)
        End Do
!end of correction
!------------------------------------------------------------------
!
        If (s>bound .And. ngmax1==ngmax) Write (1337, Fmt=110) abs(s), bound, &
          iq1, iq2
!
        Do lm = 1, lmmax
          stest(lm) = stestnew(lm)
        End Do
!******************************************************************
      Else
!******************************************************************
!********************* OUT OF THE PLANE ***************************
!******************************************************************
!
!Prepare arrays for speed up
!
        signrz = dq3/abs(dq3)
        con1 = tpi/vol
        Do l = 0, lmax
          ylmpref(l) = sqrt((2E0_dp*l+1E0_dp)/fpi)/dfac(l)*con1
          signrzl(l) = (-signrz)**l
          cim(l) = (-ci)**l
          Do m = 1, l
            ylmpref1(l, m) = con1*sqrt(real(2*l+1,kind=dp)/2E0_dp/fpi/dfac(l+m &
              )/dfac(l-m))
          End Do
        End Do
        ylmpref(0) = 1E0_dp*con1
!
!Sum in reciprocal space
!
! Find an upper cutoff for G-vectors
        i = 1
        exponent = 1.E0_dp
        crit = 60.E0_dp
        Do While (exponent<crit .And. i<ngmax)
          i = i + 1
          g1 = gn2(1, i) ! G vectors are assumed to be sorted with increasing length
          g2 = gn2(2, i)
          ga = sqrt(g1*g1+g2*g2)
          exponent = ga*abs(dq3) ! If EXPONENT>60., then EXPBSQ below is 8.7E-27, considered negligible.
        End Do
        ngmax1 = i
        If (ngmax1>ngmax) Stop 'ewald2d: 2: NGMAX1.GT.NGMAX' ! should never occur
!     WRITE(99,FMT='(A7,I8,2E10.2)') 'NGMAX1:',NGMAX1,EXP(-EXPONENT)
!==================================================================
        Do i = 2, ngmax1
!
!   Exclude the origin all terms vanish for g=0
!   except the l = 0 components which are treated
!   sepparately. (look at the begining of the sub)
!
          g1 = gn2(1, i)
          g2 = gn2(2, i)
          ga = sqrt(g1*g1+g2*g2)
          Do l = 0, lmax
            gal(l) = ga**l
          End Do
!
          expbsq = exp(-ga*abs(dq3))
          dqdotg = (dq1*g1+dq2*g2)
          factexp = exp(ci*dqdotg)*expbsq/ga
          exponl(0) = 1E0_dp
          Do l = 1, lmax
            exponl(l) = (g1+ci*g2)/ga*exponl(l-1) ! exp(i*m*fi)
          End Do
!
!In case rz < 0 then multiply by (-1)**(l-m)
!(M. Weinert J. Math Phys. 22, 2439 (1981) formula 33
! compare also formula A9)
!
          Do l = 0, lmax
!m = 0
            lm = l*(l+1) + 1
            stest(lm) = stest(lm) + ylmpref(l)*gal(l)*signrzl(l)*factexp
!m <> 0
            Do m = 1, l
              pref2(m) = cim(m)*ylmpref1(l, m)*gal(l)*signrz**m*signrzl(l)
!
! Go from the <usual> Jackson Ylm to the ones we use
!
              aprefpp = (1E0_dp/exponl(m)+exponl(m))
              aprefmm = (1E0_dp/exponl(m)-exponl(m))*ci
!     m > 0
              lm = l*(l+1) + m + 1
              stest(lm) = stest(lm) + pref2(m)*aprefpp*factexp
!     m < 0
              lm = l*(l+1) - m + 1
              stest(lm) = stest(lm) + pref2(m)*aprefmm*factexp
            End Do
          End Do
! ----------------------------------------------------------------------
          If (i==(ngmax1-nsg(nshlg))) Then
!     keep the value before the last shell to test convergence
            Do lm = 1, lmmax
              s0(lm) = stest(lm)
            End Do
          End If
! ----------------------------------------------------------------------
        End Do
! ======================================================================
!
! --> Check convergence
!
        s = 0E0_dp
        Do lm = 2, lmmax
          stest0 = abs(s0(lm)-stest(lm))
          If (s<stest0) s = stest0
        End Do
        If (s>bound .And. ngmax1==ngmax) Write (1337, Fmt=120) abs(s), bound, &
          iq1, iq2
      End If
! **********************************************************************
!
      Do lm = 1, lmmax
        If (abs(aimag(stest(lm)))>bound) Then
          Write (6, *) ' ERROR: Imaginary contribution', &
            ' to REAL lattice sum', aimag(stest(lm)), bound
          Stop
        End If
        sum2d(lm) = real(stest(lm))
        stest(lm) = 0.0E0_dp
      End Do
!
      sum2d(1) = sum2d(1)/sqrt(fpi)
!
100   Format (5X, 'WARNING 1 : Convergence of 2D-sum is ', 1P, E9.2, ' > ', &
        E9.2, 'LAYER PAIR', 2I6, /, 15X, &
        'You should use more lattice vectors (RMAX)')
110   Format (5X, 'WARNING 2 : Convergence of 2D-sum is ', 1P, E9.2, ' > ', &
        E9.2, 'LAYER PAIR', 2I6, /, 15X, &
        'You should use more lattice vectors (GMAX)')
120   Format (5X, 'WARNING 3 : Convergence of 2D-sum is ', 1P, E9.2, ' > ', &
        E9.2, 'LAYER PAIR', 2I6, /, 15X, &
        'You should use more lattice vectors (GMAX)')
    End Subroutine
