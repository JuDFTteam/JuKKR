    Subroutine mkxcpe2(ir, np, rv, rholm, vxcp, excp, ylm, dylmt1, dylmf1, &
      dylmf2, dylmtf, drrl, ddrrl, drrul, ddrrul, irmd, lmpotd, lmmax, &
      use_sol)
      Use mod_datatypes, Only: dp
!     ------------------------------------------------------------------
!     Calculation of the exchange-correlation potential.
!     coded by M. Ogura, Apr. 2015, Munich
!     ------------------------------------------------------------------
      Implicit None
      Integer :: ijd
      Parameter (ijd=434)

      Real (Kind=dp) :: rv
      Integer :: ir, irmd, lmmax, lmpotd, np
      Real (Kind=dp) :: ddrrl(irmd, lmpotd), ddrrul(irmd, lmpotd), &
        drrl(irmd, lmpotd), drrul(irmd, lmpotd), dylmf1(ijd, lmpotd), &
        dylmf2(ijd, lmpotd), dylmt1(ijd, lmpotd), dylmtf(ijd, lmpotd), &
        excp(ijd), rholm(lmpotd, 2), vxcp(ijd, 2), ylm(ijd, lmpotd)

      Real (Kind=dp) :: c, pi, s
      Integer :: ip, ispin, l1, lm, n
      Real (Kind=dp) :: d(2), d1(3, 2), d2(5, 2), dl(2)
!use_sol=0 -> PBE, use_sol=1 -> PBEsol
      Logical :: use_sol
      Real (Kind=dp) :: um, bet

      pi = 4E0_dp*atan(1E0_dp)

! Set 'UM' for subroutine 'excpbex' and 'BET' for subroutine 'excpbec' for PBE or PBEsol
      If (use_sol) Then
        um = 0.123456790123456E0_dp
        bet = 0.046E0_dp
      Else
        um = 0.2195149727645171E0_dp
        bet = 0.06672455060314922E0_dp
      End If


!     --- surface integration
      Do ip = 1, np
        Do ispin = 1, 2
          d(ispin) = 0E0_dp
          dl(ispin) = 0E0_dp
          Do n = 1, 3
            d1(n, ispin) = 0E0_dp
          End Do
          Do n = 1, 5
            d2(n, ispin) = 0E0_dp
          End Do
        End Do
        Do lm = 1, lmmax
          l1 = sqrt(real(lm,kind=dp)-5E-1_dp)
          d(1) = d(1) + rholm(lm, 1)*ylm(ip, lm)
          d(2) = d(2) + rholm(lm, 2)*ylm(ip, lm)
          dl(1) = dl(1) + real(l1*(l1+1), kind=dp)*rholm(lm, 1)*ylm(ip, lm)
          dl(2) = dl(2) + real(l1*(l1+1), kind=dp)*rholm(lm, 2)*ylm(ip, lm)
          d1(1, 2) = d1(1, 2) + drrul(ir, lm)*ylm(ip, lm)
          d1(1, 1) = d1(1, 1) + (drrl(ir,lm)-drrul(ir,lm))*ylm(ip, lm)
          d1(2, 1) = d1(2, 1) + rholm(lm, 1)*dylmt1(ip, lm)
          d1(2, 2) = d1(2, 2) + rholm(lm, 2)*dylmt1(ip, lm)
          d1(3, 1) = d1(3, 1) + rholm(lm, 1)*dylmf1(ip, lm)
          d1(3, 2) = d1(3, 2) + rholm(lm, 2)*dylmf1(ip, lm)
          d2(1, 2) = d2(1, 2) + ddrrul(ir, lm)*ylm(ip, lm)
          d2(1, 1) = d2(1, 1) + (ddrrl(ir,lm)-ddrrul(ir,lm))*ylm(ip, lm)
          d2(2, 2) = d2(2, 2) + drrul(ir, lm)*dylmt1(ip, lm)
          d2(2, 1) = d2(2, 1) + (drrl(ir,lm)-drrul(ir,lm))*dylmt1(ip, lm)
          d2(3, 2) = d2(3, 2) + drrul(ir, lm)*dylmf1(ip, lm)
          d2(3, 1) = d2(3, 1) + (drrl(ir,lm)-drrul(ir,lm))*dylmf1(ip, lm)
          d2(4, 1) = d2(4, 1) + rholm(lm, 1)*dylmtf(ip, lm)
          d2(4, 2) = d2(4, 2) + rholm(lm, 2)*dylmtf(ip, lm)
          d2(5, 1) = d2(5, 1) + rholm(lm, 1)*dylmf2(ip, lm)
          d2(5, 2) = d2(5, 2) + rholm(lm, 2)*dylmf2(ip, lm)
        End Do
        c = ylm(ip, 3)*sqrt(4E0_dp*pi/3E0_dp)
        s = sqrt(1E0_dp-c**2)
        Do ispin = 1, 2
          dl(ispin) = dl(ispin)/rv**2
          d1(2, ispin) = d1(2, ispin)/rv
          d2(2, ispin) = d2(2, ispin)/rv
          d2(4, ispin) = d2(4, ispin)/rv**2
          If (s>1E-8_dp) Then
            d1(3, ispin) = d1(3, ispin)/rv/s
            d2(3, ispin) = d2(3, ispin)/rv/s
            d2(5, ispin) = d2(5, ispin)/rv**2/s
            Call fpexcpbe(d, dl, d1, d2, rv, s, c, vxcp(ip,1), vxcp(ip,2), &
              excp(ip), um, bet)
          Else
            d1(3, ispin) = 0E0_dp
            d2(3, ispin) = 0E0_dp
            d2(5, ispin) = 0E0_dp
            vxcp(ip, 1) = 0E0_dp
            vxcp(ip, 2) = 0E0_dp
            excp(ip) = 0E0_dp
          End If
        End Do
      End Do

      Return
    End Subroutine

    Subroutine fpexcpbe(ro, rol, ro1, ro2, xr, s, c, v1, v2, exc, um, bet)
      Use mod_datatypes, Only: dp
!----------------------------------------------------------------------
!     driver routine for PBE GGA subroutines.
!     based on excpbe.f in Munich code (version on 20 Dec 2009)
!     coded by M. Ogura, Jun. 2011, Munich
!----------------------------------------------------------------------
      Implicit None

      Real (Kind=dp) :: c, exc, s, v1, v2, xr
      Real (Kind=dp) :: ro(2), rol(2), ro1(3, 2), ro2(5, 2)
      Real (Kind=dp) :: um, bet

      Real (Kind=dp) :: conf, conrs, d, drv1, drv2, drv2s, drv3, drv4, ec, ex, &
        fk, g, pi, rs, sk, ss, thrd, thrd2, tt, uu, vcdn, vcup, vv, vx, vxcdn, &
        vxcup, ww, x, xd, xu, y, z, zet
      Integer :: jsp, llda

      pi = 4E0_dp*atan(1E0_dp)
      thrd = 1E0_dp/3E0_dp
      thrd2 = 2E0_dp/3E0_dp
      conf = (3E0_dp*pi**2)**thrd
      conrs = (3E0_dp/(4E0_dp*pi))**thrd
      llda = 0
      exc = 0E0_dp
      vxcup = 0E0_dp
      vxcdn = 0E0_dp
      If (ro(1)>1E-12_dp .And. ro(2)>1E-12_dp) Then

!     ---begin the spin loop for exchange
        If (ro(1)<=1E-6_dp .Or. ro(2)<=1E-6_dp) llda = 1
        Do jsp = 1, 2
          d = 2E0_dp*ro(jsp)
          fk = conf*d**thrd
          x = ro1(1, jsp)
          y = ro1(2, jsp)
          z = ro1(3, jsp)
          drv1 = sqrt(x**2+y**2+z**2)*2E0_dp
          If (abs(drv1)<1E-8_dp) Then
            drv2 = 0E0_dp
          Else
            drv2s = 2E0_dp*y*z*ro2(4, jsp) + (z**2-y**2)*ro2(5, jsp) - &
              c/xr*y*(z**2+y**2)
            drv2 = x**2*ro2(1, jsp) + 2E0_dp*x*y*ro2(2, jsp) + &
              2E0_dp*x*z*ro2(3, jsp) - x*y**2/xr - x*z**2/xr - y**2*rol(jsp)
            If (abs(drv2s)>=1E-10_dp) drv2 = drv2 + drv2s/s
            drv2 = drv2/drv1*8E0_dp
          End If
          drv3 = ro2(1, jsp) + 2E0_dp/xr*x - rol(jsp)
          drv3 = drv3*2E0_dp
          ss = drv1/(d*2E0_dp*fk)
          uu = drv2/(d**2*(2E0_dp*fk)**3)
          vv = drv3/(d*(2E0_dp*fk)**2)
          Call excpbex(d, ss, uu, vv, ex, vx, llda, um)
          exc = exc + ex*(d/2E0_dp)/(ro(1)+ro(2))
          If (jsp==1) vxcup = vx
          If (jsp==2) vxcdn = vx
        End Do

!     ---correlation
        d = ro(1) + ro(2)
        zet = (ro(1)-ro(2))/d
        rs = conrs/d**thrd
        fk = 1.91915829E0_dp/rs
        sk = sqrt(4E0_dp*fk/pi)
        g = ((1E0_dp+zet)**thrd2+(1E0_dp-zet)**thrd2)/2E0_dp
        x = ro1(1, 1) + ro1(1, 2)
        y = ro1(2, 1) + ro1(2, 2)
        z = ro1(3, 1) + ro1(3, 2)
        drv1 = sqrt(x**2+y**2+z**2)
        If (drv1<1E-8_dp) Then
          drv2 = 0E0_dp
        Else
          drv2s = 2E0_dp*y*z*(ro2(4,1)+ro2(4,2)) + (z**2-y**2)*(ro2(5,1)+ro2(5 &
            ,2)) - c/xr*y*(z**2+y**2)
          drv2 = x**2*(ro2(1,1)+ro2(1,2)) + 2E0_dp*x*y*(ro2(2,1)+ro2(2,2)) + &
            2E0_dp*x*z*(ro2(3,1)+ro2(3,2)) - x*y**2/xr - x*z**2/xr - &
            y**2*(rol(1)+rol(2))
          If (abs(drv2s)>=1E-10_dp) drv2 = drv2 + drv2s/s
          drv2 = drv2/drv1
        End If
        drv3 = ro2(1, 1) + ro2(1, 2) + 2E0_dp/xr*x - rol(1) - rol(2)
        drv4 = x*(ro1(1,1)-ro1(1,2)-zet*x) + y*(ro1(2,1)-ro1(2,2)-zet*y) + &
          z*(ro1(3,1)-ro1(3,2)-zet*z)
        tt = drv1/(d*2E0_dp*sk*g)
        uu = drv2/(d**2*(2E0_dp*sk*g)**3)
        vv = drv3/(d*(2E0_dp*sk*g)**2)
        ww = drv4/(d**2*(2E0_dp*sk*g)**2)
        Call excpbec(rs, zet, tt, uu, vv, ww, ec, vcup, vcdn, llda, bet)
        exc = exc + ec
        vxcup = vxcup + vcup
        vxcdn = vxcdn + vcdn

      End If
!     ---convert from h to ry
      exc = 2E0_dp*exc
      xu = 2E0_dp*vxcup
      xd = 2E0_dp*vxcdn

      v1 = xu
      v2 = xd
    End Subroutine

    Subroutine excpbex(rho, s, u, v, ex, vx, llda, um)
      Use mod_datatypes, Only: dp
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burke's modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
!  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
!   (for U,V, see PW86(24))
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
!----------------------------------------------------------------------
! References:
! [a] J.P. Perdew, K. Burke, and M. Ernzerhof,
!     Phys. Rev. Lett. 77, 3865 (1996).
! [b] J.P. Perdew and Y. Wang,
!     Phys. Rev. B33, 8800 (1986); B40, 3399 (1989) (E).
!----------------------------------------------------------------------
! Formulas:
!    e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
! e_x[PBE]=e_x[unif]*FxPBE(s)
! FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13)
!----------------------------------------------------------------------
!
!  All input and output is in atomic units!
!
!  Modifications by: E. Engel
!  Last revision:    May 9, 2001
!engel
      Implicit None

! PARAMETER definitions
      Real (Kind=dp), Parameter :: thrd = 1.D0/3.D0
      Real (Kind=dp), Parameter :: thrd4 = 4.D0/3.D0
      Real (Kind=dp), Parameter :: ax = -0.738558766382022405884230032680836D0
      Real (Kind=dp), Parameter :: uk = 0.8040D0
      Real (Kind=dp) :: ul

! Dummy arguments
      Real (Kind=dp) :: ex, rho, s, u, v, vx, um
      Integer :: llda

! Local variables
      Real (Kind=dp) :: exunif, fs, fss, fxpbe, p0, s2

!----------------------------------------------------------------------
!    Define UL with via UM and UK
      ul = um/uk

!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif = ax*rho**thrd
      If (llda==1) Then
        ex = exunif
        vx = ex*thrd4
        Return
      End If
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! construct PBE enhancement factor
      s2 = s*s
      p0 = 1.E0_dp + ul*s2
      fxpbe = 1E0_dp + uk - uk/p0
      ex = exunif*fxpbe
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first and second derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
!  Fss=d Fs/ds
      fs = 2.E0_dp*uk*ul/(p0*p0)
      fss = -4.E0_dp*ul*s*fs/p0
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate potential from [b](24)
      vx = exunif*(thrd4*fxpbe-(u-thrd4*s2*s)*fss-v*fs)
    End Subroutine

    Subroutine excpbec(rs, zeta, t, uu, vv, ww, ec, vcup, vcdn, llda, bet)
      Use mod_datatypes, Only: dp
!engel
!  This subroutine evaluates the correlation energy per particle and
!  spin-up and spin-dn correlation potentials within the Perdew-Burke-
!  Ernzerhof GGA. It is a slightly modified version of K. Burke's
!  official PBE subroutine.
!
!  Input:  RS   = WIGNER-SEITZ RADIUS = ( 3 / (4*PI*(DUP+DDN)) )**(1/3)
!          ZETA = RELATIVE SPIN POLARIZATION = (DUP-DDN)/(DUP+DDN)
!          T    = ABS(GRAD D) / ( (2*SK*G) * D )
!          UU   = (GRAD D)*GRAD(ABS(GRAD D)) / ( (2*SK*G)**3 * D**2 )
!          VV   = (LAPLACIAN D) / ( (2*SK*G)**2 * D )
!          WW   = (GRAD D)*(GRAD ZETA) / ( (2*SK*G)**2 * D )
!  where:  FK   = LOCAL FERMI MOMENTUM = (3*PI**2*(DUP+DDN))**(1/3)
!          SK   = LOCAL SCREENING MOMENTUM = (4*FK/PI)**(1/2)

!  Output: EC   = correlation energy per particle
!          VCUP = spin-up correlation potential
!          VCDN = spin-dn correlation potential
!
!  All input and output is in atomic units!
!
! References:
! [a] J.P. Perdew, K. Burke, and M. Ernzerhof,
!     Phys. Rev. Lett. 77, 3865 (1996).
! [b] J. P. Perdew, K. Burke, and Y. Wang,
!     Phys. Rev. B54, 16533 (1996).
! [c] J. P. Perdew and Y. Wang,
!     Phys. Rev. B45, 13244 (1992).
!
!
!  Last revision:    May 9, 2001
!  Written by:       K. Burke, May 14, 1996.
!  Modifications by: E. Engel
!engel
      Implicit None

! PARAMETER definitions
      Real (Kind=dp), Parameter :: thrd = 1.E0/3.D0
      Real (Kind=dp), Parameter :: thrdm = -thrd
      Real (Kind=dp), Parameter :: thrd2 = 2.E0*thrd
      Real (Kind=dp), Parameter :: sixthm = thrdm/2.E0
      Real (Kind=dp), Parameter :: thrd4 = 4.E0*thrd
      Real (Kind=dp), Parameter :: gam = 0.5198420997897463295344212145565D0
      Real (Kind=dp), Parameter :: fzz = 8.D0/(9.D0*gam)
      Real (Kind=dp), Parameter :: gamma = &
        0.03109069086965489503494086371273D0
      Real (Kind=dp), Parameter :: eta = 1.D-12

! Dummy arguments
      Real (Kind=dp) :: ec, rs, t, uu, vcdn, vcup, vv, ww, zeta, bet
      Integer :: llda

! Local variables
      Real (Kind=dp) :: alfm, alfrsm, b, b2, bec, bg, comm, ecrs, eczeta, ep, &
        eprs, eu, eurs, f, fac, fact0, fact1, fact2, fact3, fact5, fz, g, g3, &
        g4, gz, h, hb, hbt, hrs, hrst, ht, htt, hz, hzt, pon, pref, q4, q5, &
        q8, q9, rsthrd, rtrs, t2, t4, t6, z4, delt
      External :: excgcor2

! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at
!          |zeta|=1.
!----------------------------------------------------------------------
!    Define DELT via BET and GAMMA
      delt = bet/gamma

!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      If (rs<3.E5_dp) Then
        rtrs = sqrt(rs)
        Call excgcor2(0.0310907E0_dp, 0.21370E0_dp, 7.5957E0_dp, 3.5876E0_dp, &
          1.6382E0_dp, 0.49294E0_dp, rtrs, eu, eurs)
        Call excgcor2(0.01554535E0_dp, 0.20548E0_dp, 14.1189E0_dp, &
          6.1977E0_dp, 3.3662E0_dp, 0.62517E0_dp, rtrs, ep, eprs)
        Call excgcor2(0.0168869E0_dp, 0.11125E0_dp, 10.357E0_dp, 3.6231E0_dp, &
          0.88026E0_dp, 0.49671E0_dp, rtrs, alfm, alfrsm)
        z4 = zeta**4
        f = ((1.E0_dp+zeta)**thrd4+(1.E0_dp-zeta)**thrd4-2.E0_dp)/gam
        ec = eu*(1.E0_dp-f*z4) + ep*f*z4 - alfm*f*(1.E0_dp-z4)/fzz
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZETA=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
        ecrs = eurs*(1.E0_dp-f*z4) + eprs*f*z4 - alfrsm*f*(1.E0_dp-z4)/fzz
        fz = thrd4*((1.E0_dp+zeta)**thrd-(1.E0_dp-zeta)**thrd)/gam
        eczeta = 4.E0_dp*(zeta**3)*f*(ep-eu+alfm/fzz) + &
          fz*(z4*ep-z4*eu-(1.E0_dp-z4)*alfm/fzz)
        comm = ec - rs*ecrs/3.E0_dp - zeta*eczeta
        vcup = comm + eczeta
        vcdn = comm - eczeta
        If (llda==1) Return
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
        g = ((1.E0_dp+zeta)**thrd2+(1.E0_dp-zeta)**thrd2)/2.E0_dp
        g3 = g**3
        pon = -ec/(g3*gamma)
        b = delt/(exp(pon)-1.E0_dp)
        b2 = b*b
        t2 = t*t
        t4 = t2*t2
        q4 = 1.E0_dp + b*t2
        q5 = 1.E0_dp + b*t2 + b2*t4
        h = g3*(bet/delt)*log(1.E0_dp+delt*q4*t2/q5)
        ec = ec + h
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
        g4 = g3*g
        t6 = t4*t2
        rsthrd = rs/3.E0_dp
        gz = (((1.E0_dp+zeta)**2+eta)**sixthm-((1.E0_dp-zeta)**2+eta)**sixthm) &
          /3.E0_dp
        fac = delt/b + 1.E0_dp
        bg = -3.E0_dp*b2*ec*fac/(bet*g4)
        bec = b2*fac/(bet*g3)
        q8 = q5*q5 + delt*q4*q5*t2
        q9 = 1.E0_dp + 2.E0_dp*b*t2
        hb = -bet*g3*b*t6*(2.E0_dp+b*t2)/q8
        hrs = -rsthrd*hb*bec*ecrs
        fact0 = 2.E0_dp*delt - 6.E0_dp*b
        fact1 = q5*q9 + q4*q9*q9
        hbt = 2.E0_dp*bet*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
        hrst = rsthrd*t2*hbt*bec*ecrs
        hz = 3.E0_dp*gz*h/g + hb*(bg*gz+bec*eczeta)
        ht = 2.E0_dp*bet*g3*q9/q8
        hzt = 3.E0_dp*gz*ht/g + hbt*(bg*gz+bec*eczeta)
        fact2 = q4*q5 + b*t2*(q4*q9+q5)
        fact3 = 2.E0_dp*b*q5*q9 + delt*fact2
        htt = 4.E0_dp*bet*g3*t*(2.E0_dp*b/q8-(q9*fact3/q8)/q8)
        comm = h + hrs + hrst + t2*ht/6.E0_dp + 7.E0_dp*t2*t*htt/6.E0_dp
        pref = hz - gz*t2*ht/g
        fact5 = gz*(2.E0_dp*ht+t*htt)/g
        comm = comm - pref*zeta - uu*htt - vv*ht - ww*(hzt-fact5)
        vcup = vcup + comm + pref
        vcdn = vcdn + comm - pref
      Else
        vcup = 0.E0_dp
        vcdn = 0.E0_dp
      End If
    End Subroutine
!*==excgcor2.f    processed by SPAG 6.55Rc at 08:17 on 20 Dec 2009

    Subroutine excgcor2(a, a1, b1, b2, b3, b4, rtrs, gg, ggrs)
      Use mod_datatypes, Only: dp
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
      Implicit None

! Dummy arguments
      Real (Kind=dp) :: a, a1, b1, b2, b3, b4, gg, ggrs, rtrs

! Local variables
      Real (Kind=dp) :: q0, q1, q2, q3

      q0 = -2.E0_dp*a*(1.E0_dp+a1*rtrs*rtrs)
      q1 = 2.E0_dp*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
      q2 = log(1.E0_dp+1.E0_dp/q1)
      gg = q0*q2
      q3 = a*(b1/rtrs+2.E0_dp*b2+rtrs*(3.E0_dp*b3+4.E0_dp*b4*rtrs))
      ggrs = -2.E0_dp*a*a1*q2 - q0*q3/(q1*(1.E0_dp+q1))
    End Subroutine
