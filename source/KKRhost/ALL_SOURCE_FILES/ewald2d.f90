subroutine ewald2d(lpot, alat, vec1, vec2, iq1, iq2, rm2, nrmax, nshlr, nsr, &
  gn2, ngmax, nshlg, nsg, sum2d, vol, lassld, lmxspd)
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
  implicit none
!..
!.. Scalar arguments ..
  integer :: lassld, lmxspd
! parameter (lassld replaces old l2maxd=2*LPOTD=4*lmaxd)
  integer :: lpot, nrmax, nshlr, ngmax, nshlg, iq1, iq2
  double precision :: alat, vol
!..
!.. Array arguments ..
  double precision :: vec1(3), vec2(3)
  double precision :: rm2(2, *), gn2(2, *), sum2d(lmxspd)
  integer :: nsr(*), nsg(*)
!..
!.. Local scalars ..
  double precision :: alpha, beta, bound, con, con1, dot1, dq1, dq2, dq3, &
    dqdotg, expbsq, fpi, g1, g2, g3, ga, lamda, pi, r, r0, r1, r2, r3, rfac, &
    s, signrz, stest0, tpi
  double complex :: aprefmm, aprefpp, bfac, ci, factexp, simag
  double precision :: derfc, exponent, crit
  integer :: i, im, ir, l, lm, lmax, lmmax, m, icall, ngmax1
!..
!.. Local arrays ..
  double precision :: dfac(0:2*lassld+1), g(0:lassld), gal(0:lassld), gi(0:4), &
    gr(0:4), pref0(0:lassld), signrzl(0:lassld), ylm(lmxspd), &
    ylmpref(0:lassld), ylmpref1(0:lassld, 0:lassld)
  double complex :: cim(0:lassld), exponl(0:lassld), pref2(lassld), &
    s0(lmxspd), stest(lmxspd), stestnew(lmxspd)
!..
!.. External subroutines ..
  external :: fplaneg, fplaner
!..
!.. Intrinsic functions ..
  intrinsic :: dabs, atan, dble, exp, sqrt
!..
!.. Data statements
  data icall/0/
!..
!.. Save statements
  save :: icall, ci, bound, pi, fpi, tpi
!..................................................................
!
  icall = 1
!ICALL = ICALL + 1
  if (icall==1) then
    ci = (0.0d0, 1.0d0)
    bound = 1.0d-9
    pi = 4.0d0*atan(1.0d0)
    fpi = 4.0d0*pi
    tpi = 2.0d0*pi
  end if
!
! Factorial
!
  dfac(0) = 1
  do l = 1, 2*lassld + 1
    dfac(l) = dfac(l-1)*l
  end do
  do l = 0, lassld
    pref0(l) = 0d0
  end do

  lmax = 2*lpot
  lmmax = (lmax+1)**2
!
  pref0(2) = sqrt(5d0/pi)/2d0/2d0
  pref0(4) = 3d0*sqrt(9d0/pi)/16d0/9d0
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
  do lm = 1, lmmax
    stest(lm) = 0d0
    stestnew(lm) = 0d0
  end do
!
!Add correction if rz = 0
!
  if (dabs(dq3)<1d-6) then
    stest(1) = stest(1) - 2d0*lamda/sqrt(pi) - 2d0*sqrt(pi)/lamda/vol
    stestnew(1) = stestnew(1) - 2d0*lamda/sqrt(pi) - 2d0*sqrt(pi)/lamda/vol
!
    if ((dq1*dq1+dq2*dq2)>1d-6) then
      stest(1) = stest(1) + 2d0*lamda/sqrt(pi)
      stestnew(1) = stestnew(1) + 2d0*lamda/sqrt(pi)
    end if
  else
!
!Add correction if rz<> 0
!
    stest(1) = stest(1) - dabs(dq3)*fpi/2d0/vol
    stest(3) = stest(3) - dabs(dq3)/dq3*sqrt(3d0*fpi)/2d0/vol ! -d/dz
! the correction for higher l vanishes...
  end if
!
! **********************************************************************
! ******************    I N-P L A N E      M = 0  **********************
! **********************************************************************
  if (dabs(dq3)<1d-6) then
!
!Real space sum
!
!==================================================================
    do ir = 1, nrmax
      r1 = dq1 - rm2(1, ir)
      r2 = dq2 - rm2(2, ir)
      r3 = 0d0
      r = sqrt(r1*r1+r2*r2)
!------------------------------------------------------------------
      if (r>1d-8) then
        alpha = lamda*r
        call fplaner(alpha, gr, r)
        do l = 0, 4
          lm = l*(l+1) + 1 ! m =0
          stest(lm) = stest(lm) + gr(l)
        end do
        call ymy(r1, r2, r3, r0, ylm, lassld)
        call gamfc(alpha, g, lassld, r)
        ylm(1) = 1d0 ! just definition matter
!
        do l = 0, lmax
          rfac = g(l)/sqrt(pi)
          do m = -l, l
            lm = l*(l+1) + m + 1
            stestnew(lm) = stestnew(lm) + ylm(lm)*rfac
          end do
        end do
!
        if (ir==(nrmax-nsr(nshlr))) then
!keep the value before the last shell to test convergence
          do l = 0, lmax
            do m = -l, l
              lm = l*(l+1) + m + 1
              s0(lm) = stestnew(lm)
            end do
          end do
        end if
      end if ! r <> 0
!------------------------------------------------------------------
    end do ! ir loop
!==================================================================
!
!Check convergence
!
    s = 0d0
!==================================================================
    do l = 0, lmax
      do m = -l, l
        lm = l*(l+1) + m + 1
        stest0 = abs(s0(lm)-stestnew(lm))
        if (s<stest0) s = stest0
      end do
    end do
!==================================================================
    if (s>bound) write (1337, fmt=100) abs(s), bound, iq1, iq2
!
!Sum in reciprocal lattice
!
    con = fpi/2d0/vol
! Find an upper cutoff for G-vectors
    i = 1
    exponent = 1.d0
    crit = 8.d0
    do while (exponent<crit .and. i<ngmax)
      i = i + 1
      g1 = gn2(1, i) ! G vectors are assumed to be sorted with increasing length
      g2 = gn2(2, i)
      ga = sqrt(g1*g1+g2*g2)
      exponent = ga/2.d0/lamda ! If EXPONENT>8., then ERFC(EXPONENT) and EXP(-EXPONENT**2) 
! (see sub. fplaneg) are below is 1E-27, considered negligible.
    end do
    ngmax1 = i
!     WRITE(99,FMT='(A7,I8,2E10.2)') 'NGMAX1:',NGMAX1,
! &                EXP(-EXPONENT**2),ERFC(EXPONENT)

    if (ngmax1>ngmax) stop 'ewald2d: 1: NGMAX1.GT.NGMAX' ! should never occur
!==================================================================
    do im = 1, ngmax1
      g1 = gn2(1, im)
      g2 = gn2(2, im)
      g3 = 0d0
      ga = sqrt(g1*g1+g2*g2)
!------------------------------------------------------------------
      dot1 = dq1*g1 + dq2*g2
      call fplaneg(lamda, gi, pref0, lassld, ga, vol)
      simag = exp(ci*dot1)
      do l = 0, 4
        lm = l*(l+1) + 1
        stest(lm) = stest(lm) + gi(l)*simag
      end do
!------------------------------------------------------------------
      if (ga>1d-6) then
        call ymy(g1, g2, g3, ga, ylm, lassld)
        beta = ga/lamda
        expbsq = derfc(beta/2d0)
!
        bfac = con*simag*expbsq
        stestnew(1) = stestnew(1) + bfac/ga
!
        do l = 0, lmax
          if (l/=0) then
            do m = -l, l
              lm = l*(l+1) + m + 1
              stestnew(lm) = stestnew(lm) + ylm(lm)*bfac*ga**(l-1)
            end do
          end if
          bfac = bfac/ci/dble(2*l+1)
        end do
      end if
!------------------------------------------------------------------
      if (im==(ngmax1-nsg(nshlg))) then
!keep the value before the last shell to test convergence
        do lm = 1, lmmax
          s0(lm) = stestnew(lm)
        end do
      end if
    end do
!==================================================================
!
!Check convergence
!
    do lm = 1, lmmax
      stest0 = abs(s0(lm)-stestnew(lm))
      if (s<stest0) s = stest0
    end do
!
!Correction due to r=0 term only for DRn = 0
!
!------------------------------------------------------------------
    if ((dq1*dq1+dq2*dq2)<1d-6) then
      do l = 2, 4, 2
        pref0(l) = pref0(l)*dfac(l)/dfac(l/2)/(l+1)
      end do
!
      i = 1
      do l = 2, 4, 2
        i = i + 1
        lm = l*(l+1) + 1
        stest(lm) = stest(lm) + (-1)**i*2d0/sqrt(pi)*pref0(l)*lamda**(l+1)
      end do
    end if
!------------------------------------------------------------------
    do l = 2, 4, 2
      lm = l*(l+1) + 1
      stestnew(lm) = stest(lm)
    end do
!end of correction
!------------------------------------------------------------------
!
    if (s>bound .and. ngmax1==ngmax) write (1337, fmt=110) abs(s), bound, iq1, &
      iq2
!
    do lm = 1, lmmax
      stest(lm) = stestnew(lm)
    end do
!******************************************************************
  else
!******************************************************************
!********************* OUT OF THE PLANE ***************************
!******************************************************************
!
!Prepare arrays for speed up
!
    signrz = dq3/dabs(dq3)
    con1 = tpi/vol
    do l = 0, lmax
      ylmpref(l) = sqrt((2d0*l+1d0)/fpi)/dfac(l)*con1
      signrzl(l) = (-signrz)**l
      cim(l) = (-ci)**l
      do m = 1, l
        ylmpref1(l, m) = con1*sqrt(dble(2*l+1)/2d0/fpi/dfac(l+m)/dfac(l-m))
      end do
    end do
    ylmpref(0) = 1d0*con1
!
!Sum in reciprocal space
!
! Find an upper cutoff for G-vectors
    i = 1
    exponent = 1.d0
    crit = 60.d0
    do while (exponent<crit .and. i<ngmax)
      i = i + 1
      g1 = gn2(1, i) ! G vectors are assumed to be sorted with increasing length
      g2 = gn2(2, i)
      ga = sqrt(g1*g1+g2*g2)
      exponent = ga*dabs(dq3) ! If EXPONENT>60., then EXPBSQ below is 8.7E-27, considered negligible.
    end do
    ngmax1 = i
    if (ngmax1>ngmax) stop 'ewald2d: 2: NGMAX1.GT.NGMAX' ! should never occur
!     WRITE(99,FMT='(A7,I8,2E10.2)') 'NGMAX1:',NGMAX1,EXP(-EXPONENT)
!==================================================================
    do i = 2, ngmax1
!
!   Exclude the origin all terms vanish for g=0
!   except the l = 0 components which are treated
!   sepparately. (look at the begining of the sub)
!
      g1 = gn2(1, i)
      g2 = gn2(2, i)
      ga = sqrt(g1*g1+g2*g2)
      do l = 0, lmax
        gal(l) = ga**l
      end do
!
      expbsq = exp(-ga*dabs(dq3))
      dqdotg = (dq1*g1+dq2*g2)
      factexp = exp(ci*dqdotg)*expbsq/ga
      exponl(0) = 1d0
      do l = 1, lmax
        exponl(l) = (g1+ci*g2)/ga*exponl(l-1) ! exp(i*m*fi)
      end do
!
!In case rz < 0 then multiply by (-1)**(l-m)
!(M. Weinert J. Math Phys. 22, 2439 (1981) formula 33
! compare also formula A9)
!
      do l = 0, lmax
!m = 0
        lm = l*(l+1) + 1
        stest(lm) = stest(lm) + ylmpref(l)*gal(l)*signrzl(l)*factexp
!m <> 0
        do m = 1, l
          pref2(m) = cim(m)*ylmpref1(l, m)*gal(l)*signrz**m*signrzl(l)
!
! Go from the <usual> Jackson Ylm to the ones we use
!
          aprefpp = (1d0/exponl(m)+exponl(m))
          aprefmm = (1d0/exponl(m)-exponl(m))*ci
!     m > 0
          lm = l*(l+1) + m + 1
          stest(lm) = stest(lm) + pref2(m)*aprefpp*factexp
!     m < 0
          lm = l*(l+1) - m + 1
          stest(lm) = stest(lm) + pref2(m)*aprefmm*factexp
        end do
      end do
! ----------------------------------------------------------------------
      if (i==(ngmax1-nsg(nshlg))) then
!     keep the value before the last shell to test convergence
        do lm = 1, lmmax
          s0(lm) = stest(lm)
        end do
      end if
! ----------------------------------------------------------------------
    end do
! ======================================================================
!
! --> Check convergence
!
    s = 0d0
    do lm = 2, lmmax
      stest0 = abs(s0(lm)-stest(lm))
      if (s<stest0) s = stest0
    end do
    if (s>bound .and. ngmax1==ngmax) write (1337, fmt=120) abs(s), bound, iq1, &
      iq2
  end if
! **********************************************************************
!
  do lm = 1, lmmax
    if (abs(aimag(stest(lm)))>bound) then
      write (6, *) ' ERROR: Imaginary contribution', ' to REAL lattice sum', &
        aimag(stest(lm)), bound
      stop
    end if
    sum2d(lm) = dble(stest(lm))
    stest(lm) = 0.0d0
  end do
!
  sum2d(1) = sum2d(1)/sqrt(fpi)
!
100 format (5x, 'WARNING 1 : Convergence of 2D-sum is ', 1p, e9.2, ' > ', &
    e9.2, 'LAYER PAIR', 2i6, /, 15x, &
    'You should use more lattice vectors (RMAX)')
110 format (5x, 'WARNING 2 : Convergence of 2D-sum is ', 1p, e9.2, ' > ', &
    e9.2, 'LAYER PAIR', 2i6, /, 15x, &
    'You should use more lattice vectors (GMAX)')
120 format (5x, 'WARNING 3 : Convergence of 2D-sum is ', 1p, e9.2, ' > ', &
    e9.2, 'LAYER PAIR', 2i6, /, 15x, &
    'You should use more lattice vectors (GMAX)')
end subroutine
