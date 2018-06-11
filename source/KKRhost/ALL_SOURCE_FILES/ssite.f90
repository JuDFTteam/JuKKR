    Subroutine ssite(iwrregwf, iwrirrwf, nfilcbwf, calcint, getirrsol, soctl, &
      ctl, eryd, p, ihyper, iprint, ikm1lin, ikm2lin, nlq, nkmq, nlinq, nt, &
      nkm, iqat, tsst, msst, tsstlin, dzz, dzj, szz, szj, ozz, ozj, bzz, bzj, &
      qzz, qzj, tzz, tzj, vt, bt, at, z, nucleus, r, drdi, r2drdi, jws, imt, &
      ameopo, lopt, solver, cgc, ozzs, ozjs, nlmax, nqmax, linmax, nrmax, &
      nmmax, ntmax, nkmmax, nkmpmax, nlamax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   * ASSIGN QUANTUM NUMBERS AND CALL ROUTINE TO SOLVE                 *
!   * 8 COUPLED PARTIAL DIFFERENTIAL RADIAL DIRAC EQUATIONS. THE       *
!   * RESULTING WAVEFUNCTIONS ARE USED TO CALCULATE T-MATRICES IN      *
!   * THE KAPPA-MU REPRESENTATION                                      *
!   *                                                                  *
!   * + CALCULATION OF THE RADIAL INTEGRALS                            *
!   *   [ G1*G2 + F1*F2 ] R**2 DR                                      *
!   *                                                                  *
!   * FOR IHYPER <> 0 :                                                *
!   * CALCULATION OF THE HYPERFINE MATRIX ELEMENTS                     *
!   *                                                                  *
!   * RYD-UNITS USED THROUGHOUT                                        *
!   *                                                                  *
!   * NOTE: to save storage force  JG/JF  and  PR/QR  to share the     *
!   *       same storage by corresponding argument list in CALL ....   *
!   *                                                                  *
!   * 28/10/94  HE  tidy up,  P,Q used in <DIRAC> instead of g,f       *
!   * 05/10/96  HE  workspace for wavefunctions and matrices           *
!   *               is allocated dynamically !!!                       *
!   * 07/02/05  VP  few changes connected to the calculation of orbital*
!   *               polarisation                                       *
!   ********************************************************************


      Implicit None

!PARAMETER definitions
      Complex (Kind=dp) :: ci
      Parameter (ci=(0.0E0_dp,1.0E0_dp))

      Real (Kind=dp) :: f1, e0, a0, cautog

! conversion factor for hyperfine fields from A.U. to GAUSS
!                                 electron charge     in esu
!                                 Bohr-radius         in cm

      Parameter (f1=1.0E0_dp, e0=1.6021892E-19_dp*2.997930E+09_dp, &
        a0=0.52917706E-08_dp, cautog=e0/(a0*a0))

! Dummy arguments
      Logical :: calcint, getirrsol
      Complex (Kind=dp) :: eryd, p
      Integer :: ihyper, iprint, iwrirrwf, iwrregwf, linmax, nfilcbwf, nkm, &
        nkmmax, nlamax, nlmax, nmmax, nqmax, nrmax, nt, ntmax, nucleus
      Integer :: nkmpmax
      Character (Len=10) :: solver
      Real (Kind=dp) :: ameopo(nkmmax, nkmmax, nlamax, 3), &
        at(nrmax, nlamax, 3, ntmax), bt(nrmax, ntmax), ctl(ntmax, nlmax), &
        drdi(nrmax, nmmax), r(nrmax, nmmax), r2drdi(nrmax, nmmax), &
        soctl(ntmax, nlmax), vt(nrmax, ntmax)
      Real (Kind=dp) :: cgc(nkmpmax, 2)
      Complex (Kind=dp) :: ozjs(linmax, ntmax, 2), ozzs(linmax, ntmax, 2)
      Complex (Kind=dp) :: bzj(linmax, ntmax), bzz(linmax, ntmax), &
        dzj(linmax, ntmax), dzz(linmax, ntmax), msst(nkmmax, nkmmax, ntmax), &
        ozj(linmax, ntmax), ozz(linmax, ntmax), qzj(linmax, ntmax), &
        qzz(linmax, ntmax), szj(linmax, ntmax), szz(linmax, ntmax), &
        tsst(nkmmax, nkmmax, ntmax), tsstlin(linmax, ntmax), &
        tzj(linmax, ntmax), tzz(linmax, ntmax)
      Integer :: ikm1lin(linmax), ikm2lin(linmax), imt(ntmax), &
        iqat(nqmax, ntmax), jws(nmmax), lopt(ntmax), nkmq(nqmax), &
        nlinq(nqmax), nlq(nqmax), z(ntmax), muem05

! Local variables
      Complex (Kind=dp) :: a(2, 2), arg, b1, b2, cint(nrmax), crsq, csum, det, &
        dxp(2, 2), f11, f12, f21, f22, g11, g11p, g12, g12p, g21, g21p, g22, &
        g22p, gam(2, 2), gaminv(2, 2), hl, hlb1, hlb2, jf(nrmax, 2, 2), &
        jg(nrmax, 2, 2), jl, jlb1, jlb2, jlp, maux(nkmmax, nkmmax), &
        msst2(2, 2), nl, nlb1, nlb2, nlp, norm, pfac, pi(2, 2, nrmax), &
        pr(2, 2, nrmax), qi(2, 2, nrmax), qr(2, 2, nrmax), rmehf(2, 2), &
        rmehf1(2, 2), rmehf2(2, 2), sig(2, 2), tsst2(2, 2), xsst2(2, 2), &
        zf(nrmax, 2, 2), zfjf(2, 2), zfzf(2, 2), zg(nrmax, 2, 2), zgjg(2, 2), &
        zgzg(2, 2)
      Real (Kind=dp) :: ap(2, 2, nrmax), aq(2, 2, nrmax), c, cff(2, 2), &
        cfg(2, 2), cg1, cg2, cg4, cg5, cg8, cgf(2, 2), cgg(2, 2), ch(2, 2), &
        csqr, ctf(2, 2), ctg(2, 2), dovr(nrmax), drovrn(nrmax), mj, r1m(2, 2), &
        rkd(2, 2), rnuc, sk1, sk2, tdia1, tdia2, toff
      Real (Kind=dp) :: cog(2, 2, 2), cof(2, 2, 2)
      Complex (Kind=dp) :: cdjlzdz, cdnlzdz, cjlz, cnlz
      Real (Kind=dp) :: dble, dsqrt
      Integer :: i, i1, i2, i3, i5, ikm1, ikm2, il, im, in, info, &
        ipiv(nkmmax), iq, isk1, isk2, it, j, jlim, jtop, k, k1, k2, kap1, &
        kap2, kc, l, l1, lb1, lb2, lin, n, nsol, imkm1, imkm2, is, imj
      Integer :: ikapmue
      Integer :: isign, nint
      Real (Kind=dp) :: rnuctab
      Logical :: wronski

      Data r1m/1.0E0_dp, 0.0E0_dp, 0.0E0_dp, 1.0E0_dp/
      Data rkd/1.0E0_dp, 0.0E0_dp, 0.0E0_dp, -1.0E0_dp/
!     DATA RKD / 1.0D0, 0.0D0, 0.0D0, 1.0D0 /

      Call cinit(ntmax*linmax, dzz)
      Call cinit(ntmax*linmax, dzj)
      Call cinit(ntmax*linmax, szz)
      Call cinit(ntmax*linmax, szj)
      Call cinit(ntmax*linmax, ozz)
      Call cinit(ntmax*linmax, ozj)
      Call cinit(ntmax*linmax, bzz)
      Call cinit(ntmax*linmax, bzj)
      Call cinit(ntmax*linmax, qzz)
      Call cinit(ntmax*linmax, qzj)
      Call cinit(ntmax*linmax, tzz)
      Call cinit(ntmax*linmax, tzj)

      wronski = .True.
      wronski = .False.
!------------------------------------------------------------------------

      c = ctl(1, 1)
      csqr = c*c

!     calculate relativistic momentum

      p = sqrt(eryd*(1.0E0_dp+eryd/csqr))

      pfac = p/(1.0E0_dp+eryd/csqr)


! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      Do it = 1, nt

        If (iprint>0) Write (1337, '(A,I3,A,10F10.4)') ' SOLVER ', it, solver, &
          (soctl(it,il), il=1, nlmax)

        iq = iqat(1, it)
        im = imt(it)
        jtop = jws(im)
        If (nucleus/=0) Then
          rnuc = rnuctab(z(it))
          in = 1
          Do While (r(in,im)<rnuc)
            in = in + 1
          End Do
!        RLIM=R(IN,IM)
          jlim = in
          If (mod(jlim,2)==0) jlim = jlim - 1
          rnuc = r(jlim, im)
        End If
        Do i = 1, jtop
          dovr(i) = drdi(i, im)/r(i, im)
          If (nucleus/=0) drovrn(i) = (r(i,im)/rnuc)**3*drdi(i, im)
        End Do

        lin = 0

! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
        Do l = 0, (nlq(iq)-1)
          il = l + 1
          c = ctl(it, il)
          csqr = c*c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          kap1 = -l - 1
          kap2 = l
          If (l==0) kap2 = kap1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          isk1 = isign(1, kap1)
          isk2 = isign(1, kap2)
          sk1 = dble(isk1)
          sk2 = dble(isk2)
          l1 = l
          lb1 = l - isk1
          lb2 = l - isk2

          arg = p*r(jtop, im)
          jl = cjlz(l1, arg)
          jlb1 = cjlz(lb1, arg)
          jlb2 = cjlz(lb2, arg)
          nl = cnlz(l1, arg)
          nlb1 = cnlz(lb1, arg)
          nlb2 = cnlz(lb2, arg)
          hl = jl + ci*nl
          hlb1 = jlb1 + ci*nlb1
          hlb2 = jlb2 + ci*nlb2

          If (solver(1:7)=='ABM-SOC') Then
            jlp = cdjlzdz(l1, arg, 1)*p
            nlp = cdnlzdz(l1, arg, 1)*p
          End If

! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!DO MJ = -(DBLE(L)+0.5D0), + (DBLE(L)+0.5D0),1.0D0
          Do imj = 1, 2*l + 1
            mj = -(dble(l)+0.5E0_dp) + dble(imj-1)

!------------------------------------------------------------------------
! NO COUPLING FOR:  ABS(MUE)= J   +  J=L+1/2 == KAP=-L-1
            If (abs(mj)>=dble(l)) Then
              nsol = 1
            Else
              nsol = 2
            End If
!------------------------------------------------------------------------
            muem05 = nint(mj-0.5E0_dp)
            ikm1 = ikapmue(kap1, muem05)
            ikm2 = ikapmue(kap2, muem05)
            imkm1 = ikapmue(-kap1, muem05)
            imkm2 = ikapmue(-kap2, muem05)
!------------------------------------------------------------------------
            i5 = nrmax*2*2
            Call cinit(i5, zg)
            Call cinit(i5, zf)
            Call cinit(i5, jg)
            Call cinit(i5, jf)
            Call cinit(i5, pr)
            Call cinit(i5, qr)
            Call cinit(i5, pi)
            Call cinit(i5, qi)

            If (solver(1:2)=='BS') Then
              Call dirbs(getirrsol, ctl(it,il), eryd, l, mj, kap1, kap2, p, &
                cg1, cg2, cg4, cg5, cg8, vt(1,it), bt(1,it), z(it), nucleus, &
                r(1,im), drdi(1,im), dovr, jtop, pr, qr, pi, qi, zg, zf)
            Else If (solver=='ABM-BI') Then
              Stop &
                ' < DIRABMBI > : Not implemented. Set SOLVER=BS in inputcard'
!                   CALL DIRABMBI(GETIRRSOL,CTL(IT,IL),IT,ERYD,L,MJ,KAP1,
!      &                          KAP2,P,CG1,CG2,CG4,CG5,CG8,AMEOPC,
!      &                          AMEOPO,VT(1,IT),BT(1,IT),AT,Z(IT),
!      &                          NUCLEUS,R(1,IM),DRDI(1,IM),DOVR,JTOP,PR,
!      &                          QR,PI,QI,ZG,ZF,AP,AQ,NTMAX,NLAMAX,
!      &                          NKMMAX,NRMAX)
            Else If (solver(1:6)=='ABM-OP') Then
              Call dirabmop(getirrsol, ctl(it,il), it, eryd, l, mj, kap1, &
                kap2, p, cg1, cg2, cg4, cg5, cg8, ameopo, vt(1,it), bt(1,it), &
                at, z(it), nucleus, r(1,im), drdi(1,im), dovr, jtop, pr, qr, &
                pi, qi, zg, zf, ap, aq, lopt(it), ntmax, nlamax, nkmmax, &
                nrmax)
            Else If (solver=='ABM-SOC   ') Then
              Call dirabmsoc(getirrsol, ctl(it,il), soctl(it,il), it, eryd, l, &
                mj, kap1, kap2, p, cg1, cg2, cg4, cg5, cg8, vt(1,it), &
                bt(1,it), z(it), nucleus, r(1,im), drdi(1,im), dovr, jtop, &
                dxp, pr, qr, pi, qi, zg, zf, nrmax)
            Else If (solver=='ABM-SOC-II') Then
              Call dirabmsoc2(getirrsol, ctl(it,il), soctl(it,il), it, eryd, &
                l, mj, kap1, kap2, p, cg1, cg2, cg4, cg5, cg8, vt(1,it), &
                bt(1,it), z(it), nucleus, r(1,im), drdi(1,im), dovr, jtop, &
                dxp, pr, qr, pi, qi, zg, zf, nrmax)
            Else
              Write (6, *) 'No solver found for: ', solver
              Stop
            End If

!  wavefunctions at the muffin-tin-radius

            n = jtop

            g11 = pr(1, 1, n)/r(n, im)
            f11 = qr(1, 1, n)/(r(n,im)*c)
            g21 = pr(2, 1, n)/r(n, im)
            f21 = qr(2, 1, n)/(r(n,im)*c)
            g22 = pr(2, 2, n)/r(n, im)
            f22 = qr(2, 2, n)/(r(n,im)*c)
            g12 = pr(1, 2, n)/r(n, im)
            f12 = qr(1, 2, n)/(r(n,im)*c)

            If (solver(1:7)=='ABM-SOC') Then
              g11p = (dxp(1,1)/drdi(n,im)-g11)/r(n, im)
              g21p = (dxp(2,1)/drdi(n,im)-g21)/r(n, im)
              g12p = (dxp(1,2)/drdi(n,im)-g12)/r(n, im)
              g22p = (dxp(2,2)/drdi(n,im)-g22)/r(n, im)
            End If

! ------- the minor component for the soc-manipulated wf is meaningless

            If (solver(1:7)=='ABM-SOC') Then
              Call cinit(2*2*nrmax, qr)
              Call cinit(2*2*nrmax, qi)
            End If

!      COSD= NL * C * F11 - PFAC * SK1 * NLB1 * G11
!      SIND= JL * C * F11 - PFAC * SK1 * JLB1 * G11

! -------------------------------------------------------------------
!       T-SS  CONSTRUCTED USING EXPRESSIONS FROM H.E. + B.L.G. (1988)
! -------------------------------------------------------------------

            nl = (hl-jl)/ci
            nlb1 = (hlb1-jlb1)/ci
            nlb2 = (hlb2-jlb2)/ci


            If (solver(1:7)=='ABM-SOC') Then
              gam(1, 1) = jl*g11p - jlp*g11
              gam(2, 1) = jl*g21p - jlp*g21
              gam(1, 2) = jl*g12p - jlp*g12
              gam(2, 2) = jl*g22p - jlp*g22

              sig(1, 1) = nl*g11p - nlp*g11
              sig(2, 1) = nl*g21p - nlp*g21
              sig(1, 2) = nl*g12p - nlp*g12
              sig(2, 2) = nl*g22p - nlp*g22
            Else
              gam(1, 1) = jl*c*f11 - pfac*sk1*jlb1*g11
              gam(2, 1) = jl*c*f21 - pfac*sk2*jlb2*g21
              gam(1, 2) = jl*c*f12 - pfac*sk1*jlb1*g12
              gam(2, 2) = jl*c*f22 - pfac*sk2*jlb2*g22

              sig(1, 1) = nl*c*f11 - pfac*sk1*nlb1*g11
              sig(2, 1) = nl*c*f21 - pfac*sk2*nlb2*g21
              sig(1, 2) = nl*c*f12 - pfac*sk1*nlb1*g12
              sig(2, 2) = nl*c*f22 - pfac*sk2*nlb2*g22
            End If

            Call zcopy(nsol*nsol, gam, 1, gaminv, 1)
            Call zgetrf(nsol, nsol, gaminv, 2, ipiv, info)
            Call zgetri(nsol, gaminv, 2, ipiv, maux, 2*2, info)

            Do i2 = 1, nsol
              Do i1 = 1, nsol
                csum = 0.0E0_dp
                Do i3 = 1, nsol
                  csum = csum + sig(i1, i3)*gaminv(i3, i2)
                End Do
                xsst2(i1, i2) = p*csum
              End Do
            End Do

            Do i1 = 1, nsol
              Do i2 = 1, nsol
                msst2(i1, i2) = -xsst2(i1, i2)
              End Do
              msst2(i1, i1) = msst2(i1, i1) + ci*p
            End Do

            Call zcopy(nsol*nsol, msst2, 1, tsst2, 1)

            Call zgetrf(nsol, nsol, tsst2, 2, ipiv, info)
            Call zgetri(nsol, tsst2, 2, ipiv, maux, 2*2, info)

            If (iprint>=3) Write (1337, 100) it, l, mj, &
              ((tsst2(i1,i2),i2=1,nsol), i1=1, nsol)
!------------------------------------------------------------------------

!   COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION

            cgg(1, 1) = cg1
            cgg(1, 2) = cg2
            cgg(2, 1) = cg2
            cgg(2, 2) = cg4
            Call rinit(4, cgf)
            cgf(1, 1) = cg5
            cgf(2, 2) = cg8


!   COEFFICIENTS TO CALCULATE THE SPIN  DIPOLAR MOMENT TZ

            tdia1 = 2*mj/dble((2*l1+1)*(2*lb1+1))
            tdia2 = 2*mj/dble((2*l1+1)*(2*lb2+1))
            toff = -sqrt((l1+0.5E0_dp)**2-mj**2)/dble(2*l1+1)

            ctg(1, 1) = 0.5E0_dp*(cg1-3.0E0_dp*tdia1)
            ctg(1, 2) = 0.5E0_dp*(cg2-3.0E0_dp*toff)
            ctg(2, 1) = 0.5E0_dp*(cg2-3.0E0_dp*toff)
            ctg(2, 2) = 0.5E0_dp*(cg4-3.0E0_dp*tdia2)
            Call rinit(4, ctf)


!   COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION

            cfg(1, 1) = mj*(kap1+1.0E0_dp)/(kap1+0.5E0_dp)
            cfg(2, 2) = mj*(kap2+1.0E0_dp)/(kap2+0.5E0_dp)
            cfg(1, 2) = 0.5E0_dp*dsqrt(1.0E0_dp-(mj/(kap1+0.5E0_dp))**2)
            cfg(2, 1) = cfg(1, 2)
            Call rinit(4, cff)
            cff(1, 1) = mj*(-kap1+1.0E0_dp)/(-kap1+0.5E0_dp)
            cff(2, 2) = mj*(-kap2+1.0E0_dp)/(-kap2+0.5E0_dp)

!-----------------------------------------------------------------------
!   COEFFICIENTS TO CALCULATE THE ORBITAL POLARISATION

            Call rinit(4*2, cog)
            Call rinit(4*2, cof)
            Do is = 1, 2
              cog(1, 1, is) = cgc(ikm1, is)*cgc(ikm1, is)*dble(muem05-is+2)
              cof(1, 1, is) = cgc(imkm1, is)*cgc(imkm1, is)*dble(muem05-is+2)
            End Do

            If (nsol==2) Then
              Do is = 1, 2
                cog(2, 2, is) = cgc(ikm2, is)*cgc(ikm2, is)*dble(muem05-is+2)
                cof(2, 2, is) = cgc(imkm2, is)*cgc(imkm2, is)* &
                  dble(muem05-is+2)

                cog(1, 2, is) = cgc(ikm1, is)*cgc(ikm2, is)*dble(muem05-is+2)
                cog(2, 1, is) = cog(1, 2, is)
              End Do
            End If
!-----------------------------------------------------------------------
!   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
!        THE FACTOR  I  HAS BEEN OMITTED

            ch(1, 1) = 4.0E0_dp*kap1*mj/(4.0E0_dp*kap1*kap1-1.0E0_dp)
            ch(2, 2) = 4.0E0_dp*kap2*mj/(4.0E0_dp*kap2*kap2-1.0E0_dp)
            If (nsol==2) Then
              ch(1, 2) = dsqrt(0.25E0_dp-(mj/dble(kap1-kap2))**2)
              ch(2, 1) = ch(1, 2)
            End If
!-----------------------------------------------------------------------
!ALCULATE RADIAL INTEGRALS  UP TO   OR RWS(JTOP=JWS)


            If (nsol==1) Then
!====================================================================
! NO COUPLING TO OTHER SCATTERING CHANNELS
! REGULAR PART    Z*Z    Z = (GRA,FRA)

              norm = (p*nl-jl*xsst2(1,1))/g11

              Do i = 1, jtop
                zg(i, 1, 1) = (pr(1,1,i)/r(i,im))*norm
                zf(i, 1, 1) = (qr(1,1,i)/r(i,im)/c)*norm
                jg(i, 1, 1) = pi(1, 1, i)/r(i, im)
                jf(i, 1, 1) = qi(1, 1, i)/r(i, im)/c
              End Do


              If (iwrregwf/=0) Write (nfilcbwf, Rec=ikm1+(it-1)*nkm) it, l, &
                mj, nsol, 'REG', kap1, ikm1, (zg(i,1,1), zf(i,1,1), i=1, jtop)

              If (iwrirrwf/=0) Write (nfilcbwf, Rec=ikm1+(it-1+nt)*nkm) it, l, &
                mj, nsol, 'IRR', kap1, ikm1, (jg(i,1,1), jf(i,1,1), i=1, jtop)

!============================================== NO COUPLING = END ===
            Else
!====================================================================
! COUPLING OF TWO SCATTERING CHANNELS
!   Z(K1,K2):  INDEX 1: SPIN-ANGULAR CHARACTER
!              INDEX 2: BOUNDARY CONDITION


              det = g11*g22 - g12*g21

!OEFFICIENTS TO GET:   Z(K1,K1)  Z(K2,K1)
              b1 = p*nl - xsst2(1, 1)*jl
              b2 = -xsst2(2, 1)*jl
              a(1, 1) = (g22*b1-g12*b2)/det
              a(2, 1) = (g11*b2-g21*b1)/det

!OEFFICIENTS TO GET:   Z(K1,K2)  Z(K2,K2)
              b1 = -xsst2(1, 2)*jl
              b2 = p*nl - xsst2(2, 2)*jl
              a(1, 2) = (g22*b1-g12*b2)/det
              a(2, 2) = (g11*b2-g21*b1)/det

!ALCULATE FUNCTIONS: Z(K1,K1), Z(K2,K1), Z(K1,K2), Z(K2,K2)
              Do k = 1, nsol
                Do i = 1, jtop
                  zg(i, 1, k) = (pr(1,1,i)*a(1,k)+pr(1,2,i)*a(2,k))/r(i, im)
                  zf(i, 1, k) = (qr(1,1,i)*a(1,k)+qr(1,2,i)*a(2,k))/r(i, im)/c

                  zg(i, 2, k) = (pr(2,1,i)*a(1,k)+pr(2,2,i)*a(2,k))/r(i, im)
                  zf(i, 2, k) = (qr(2,1,i)*a(1,k)+qr(2,2,i)*a(2,k))/r(i, im)/c
                End Do
              End Do
              Do k = 1, nsol
                kc = 3 - k
                Do i = 1, jtop
                  jg(i, k, k) = pi(k, k, i)/r(i, im)
                  jf(i, k, k) = qi(k, k, i)/r(i, im)/c
                  jg(i, kc, k) = pi(kc, k, i)/r(i, im)
                  jf(i, kc, k) = qi(kc, k, i)/r(i, im)/c
                End Do
              End Do

!-----------------------------------------------------------------------
              If (iwrregwf/=0) Then
! solution 1
                Write (nfilcbwf, Rec=ikm1+(it-1)*nkm) it, l, mj, nsol, 'REG', &
                  kap1, ikm1, (zg(i,1,1), zf(i,1,1), i=1, jtop), kap2, ikm2, &
                  (zg(i,2,1), zf(i,2,1), i=1, jtop)

! solution 2
                Write (nfilcbwf, Rec=ikm2+(it-1)*nkm) it, l, mj, nsol, 'REG', &
                  kap2, ikm2, (zg(i,2,2), zf(i,2,2), i=1, jtop), kap1, ikm1, &
                  (zg(i,1,2), zf(i,1,2), i=1, jtop)
              End If

              If (iwrirrwf/=0) Then
! solution 1
                Write (nfilcbwf, Rec=ikm1+(it-1+nt)*nkm) it, l, mj, nsol, &
                  'IRR', kap1, ikm1, (jg(i,1,1), jf(i,1,1), i=1, jtop), kap2, &
                  ikm2, (jg(i,2,1), jf(i,2,1), i=1, jtop)

! solution 2
                Write (nfilcbwf, Rec=ikm2+(it-1+nt)*nkm) it, l, mj, nsol, &
                  'IRR', kap2, ikm2, (jg(i,2,2), jf(i,2,2), i=1, jtop), kap1, &
                  ikm1, (jg(i,1,2), jf(i,1,2), i=1, jtop)

              End If
!================================================= COUPLING = END ===
            End If



!ALCULATE SUM OF INTEGRALS TO BE MULTIPLIED TO   TAU(K1,K2)
            Do k1 = 1, nsol
              Do k2 = 1, nsol

                lin = lin + 1
                tsstlin(lin, it) = tsst2(k1, k2)
                If (calcint) Then
! REGULAR PART    Z*Z

                  Call cintabr(zg(1,1,k1), zg(1,1,k2), zgzg, zf(1,1,k1), &
                    zf(1,1,k2), zfzf, r2drdi(1,im), nsol, nsol, jtop, nrmax)

                  Call sumupint(dzz(lin,it), f1, zgzg, r1m, f1, zfzf, r1m, &
                    nsol)
                  Call sumupint(szz(lin,it), f1, zgzg, cgg, -f1, zfzf, cgf, &
                    nsol)
                  Call sumupint(ozz(lin,it), f1, zgzg, cfg, -f1, zfzf, cff, &
                    nsol)
! ----------------------------------------------------------------------
                  Do is = 1, 2
                    Call sumupint(ozzs(lin,it,is), f1, zgzg, cog(1,1,is), -f1, &
                      zfzf, cof(1,1,is), nsol)
                  End Do
! ----------------------------------------------------------------------
                  Call sumupint(qzz(lin,it), f1, zgzg, rkd, f1, zfzf, rkd, &
                    nsol)
                  Call sumupint(tzz(lin,it), f1, zgzg, ctg, -f1, zfzf, ctf, &
                    nsol)

!       write(66,'(3I3,2e16.7)') it,nsol,lin,DZZ(LIN,IT)
!       write(66,'(4e16.7)') ((ZGZG(ii,jj),ii=1,nsol),jj=1,nsol)
!       write(66,'(4e16.7)') ((ZFZF(ii,jj),ii=1,nsol),jj=1,nsol)


!-----------------------------------------------------------------------
                  If (ihyper==1) Then
                    Call cinthff(zg(1,1,k1), zf(1,1,k1), zg(1,1,k2), &
                      zf(1,1,k2), rmehf, nsol, nsol, jtop, cint, r(1,im), &
                      drdi(1,im), nrmax)

                    If (nucleus/=0) Then
! calculates integrals inside nucleus but up to now only
! approximately because jlim is not the nuclear radius
! the same arguments are valid for the irregular parts below
                      Call cinthff(zg(1,1,k1), zf(1,1,k1), zg(1,1,k2), &
                        zf(1,1,k2), rmehf1, nsol, nsol, jlim, cint, r(1,im), &
                        drdi(1,im), nrmax)
                      Call cinthff(zg(1,1,k1), zf(1,1,k1), zg(1,1,k2), &
                        zf(1,1,k2), rmehf2, nsol, nsol, jlim, cint, r(1,im), &
                        drovrn, nrmax)
                      Do i = 1, nsol
                        Do j = 1, nsol
                          rmehf(i, j) = rmehf(i, j) - rmehf1(i, j) + &
                            rmehf2(i, j)
                        End Do
                      End Do
                    End If
!                !end of nucleus.eq.0
                    Call sumupint(bzz(lin,it), cautog, rmehf, ch, 0.0E0_dp, &
                      rmehf, ch, nsol)
                  End If
!-----------------------------------------------------------------------

! IRREGULAR PART    Z*J
! THE  ENDING  A (B)  STANDS FOR THE DOMINATING (DIMINATED)
! SET OF SPIN-ANGULAR-CHAR:  I.E.  J==J(A,A)  FOR R>RMT

                  If (k1==k2) Then

                    Call cintabr(zg(1,1,k1), jg(1,1,k1), zgjg, zf(1,1,k1), &
                      jf(1,1,k1), zfjf, r2drdi(1,im), nsol, nsol, jtop, nrmax)

                    Call sumupint(dzj(lin,it), f1, zgjg, r1m, f1, zfjf, r1m, &
                      nsol)
                    Call sumupint(szj(lin,it), f1, zgjg, cgg, -f1, zfjf, cgf, &
                      nsol)
                    Call sumupint(ozj(lin,it), f1, zgjg, cfg, -f1, zfjf, cff, &
                      nsol)
! ----------------------------------------------------------------------
                    Do is = 1, 2
                      Call sumupint(ozjs(lin,it,is), f1, zgjg, cog(1,1,is), &
                        -f1, zfjf, cof(1,1,is), nsol)
                    End Do
! ----------------------------------------------------------------------
                    Call sumupint(qzj(lin,it), f1, zgjg, rkd, f1, zfjf, rkd, &
                      nsol)
                    Call sumupint(tzj(lin,it), f1, zgjg, ctg, -f1, zfjf, ctf, &
                      nsol)

!-----------------------------------------------------------------------
                    If (ihyper==1) Then
                      Call cinthff(zg(1,1,k1), zf(1,1,k1), jg(1,1,k1), &
                        jf(1,1,k1), rmehf, nsol, nsol, jtop, cint, r(1,im), &
                        drdi(1,im), nrmax)
                      If (nucleus/=0) Then
! calculates integrals inside nucleus but up to now only
! approximately because jlim is not the nuclear radius
! the same arguments are valid for the irregular parts below
                        Call cinthff(zg(1,1,k1), zf(1,1,k1), jg(1,1,k1), &
                          jf(1,1,k1), rmehf1, nsol, nsol, jlim, cint, r(1,im), &
                          drdi(1,im), nrmax)
                        Call cinthff(zg(1,1,k1), zf(1,1,k1), jg(1,1,k1), &
                          jf(1,1,k1), rmehf2, nsol, nsol, jlim, cint, r(1,im), &
                          drovrn, nrmax)
                        Do i = 1, nsol
                          Do j = 1, nsol
                            rmehf(i, j) = rmehf(i, j) - rmehf1(i, j) + &
                              rmehf2(i, j)
                          End Do
                        End Do
                      End If
!                !end of nucleus.eq.0

                      Call sumupint(bzj(lin,it), cautog, rmehf, ch, 0.0E0_dp, &
                        rmehf, ch, nsol)
                    End If
                  End If
!-----------------------------------------------------------------------
                End If
!           ! OF IF (.CALCINT.)
              End Do
            End Do

!heck WRONSKI-relationship

            If (wronski) Then
              Do i = 1, jtop, 40
                crsq = c*r(i, im)**2
                Write (1337, 110, Advance='no') it, l, nint(2*mj), i, &
                  r(i, im), 1.0E0_dp - (zf(i,1,1)*jg(i,1,1)-zg(i,1,1)*jf(i,1,1 &
                  )+zf(i,2,1)*jg(i,2,1)-zg(i,2,1)*jf(i,2,1))*crsq
                If (nsol==2) Then
                  Write (1337, 120) 1.0E0_dp - (zf(i,1,2)*jg(i,1,2)-zg(i,1,2)* &
                    jf(i,1,2)+zf(i,2,2)*jg(i,2,2)-zg(i,2,2)*jf(i,2,2))*crsq, &
                    1.0E0_dp - (zf(i,1,2)*jg(i,1,1)-zg(i,1,2)*jf(i,1,1)+zf(i,2 &
                    ,2)*jg(i,2,1)-zg(i,2,2)*jf(i,2,1))*crsq, &
                    1.0E0_dp - (zf(i,1,1)*jg(i,1,2)-zg(i,1,1)*jf(i,1,2)+ &
                    zf(i,2,1)*jg(i,2,2)-zg(i,2,1)*jf(i,2,2))*crsq
                Else
                  Write (1337, *)
                End If
              End Do
            End If

          End Do
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        End Do
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

        Call cinit(nkmmax*nkmmax, tsst(1,1,it))

        Do lin = 1, nlinq(iq)
          i1 = ikm1lin(lin)
          i2 = ikm2lin(lin)
          tsst(i1, i2, it) = tsstlin(lin, it)
        End Do

        Do j = 1, nkmq(iq)
          Call zcopy(nkmq(iq), tsst(1,j,it), 1, msst(1,j,it), 1)
        End Do

        Call zgetrf(nkmq(iq), nkmq(iq), msst(1,1,it), nkmmax, ipiv, info)
        Call zgetri(nkmq(iq), msst(1,1,it), nkmmax, ipiv, maux, nkmmax*nkmmax, &
          info)

        If (iprint>=4) Then
          Do lin = 1, nlinq(iq)
            i1 = ikm1lin(lin)
            i2 = ikm2lin(lin)
            Write (1337, 130) it, lin
            Write (1337, 130) it, i1, i2, ' DZZ ', dzz(lin, it), dzj(lin, it)
            Write (1337, 130) it, i1, i2, ' SZZ ', szz(lin, it), szj(lin, it)
            Write (1337, 130) it, i1, i2, ' OZZ ', ozz(lin, it), ozj(lin, it)
          End Do
        End If
      End Do
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

100   Format ('  t-ss(TLM)', 2I3, F5.1, 2E14.5, 2X, 2E14.5, :, /, 22X, 2E14.5, &
        2X, 2E14.5)
110   Format (' IT=', I2, ' L=', I2, ' MJ=', I2, '/2', I7, F10.6, 1(2X,2F9.6))
120   Format (3(2X,2F9.6))
130   Format (' IT=', I2, 2I3, A, 2X, 2E14.5, 2X, 2E14.5)
    End Subroutine

    Subroutine readwfun(nfil, it, l, mj, nsol, sreg, sirr, ikm1, kap1, ikm2, &
      kap2, nt, nkm, zg, zf, jg, jf, jtop, nrmax)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   *  reread the wave functions written by  <SSITE>  or  <CORE>       *
!   *                                                                  *
!   ********************************************************************
      Implicit None

! Dummy arguments
      Integer :: ikm1, ikm2, it, jtop, kap1, kap2, l, nfil, nkm, nrmax, nsol, &
        nt
      Real (Kind=dp) :: mj
      Character (Len=3) :: sirr, sreg
      Complex (Kind=dp) :: jf(nrmax, 2, 2), jg(nrmax, 2, 2), zf(nrmax, 2, 2), &
        zg(nrmax, 2, 2)

! Local variables
      Integer :: i, iflag, ikmin(2), itp, k, kapin(2), lp, nsolp
      Real (Kind=dp) :: mjp
      Character (Len=3) :: strp

      iflag = 0
!-----------------------------------------------------------------------
!                                    REGULAR wave function -- solution 1
      If (sreg=='REG' .Or. sreg=='COR') Then
        Read (nfil, Rec=ikm1+(it-1)*nkm) itp, lp, mjp, nsolp, strp, &
          (kapin(k), ikmin(k), (zg(i,k,1),zf(i,k,1),i=1,jtop), k=1, nsol)
        If (itp/=it .Or. lp/=l .Or. abs(mjp-mj)>0.001E0_dp .Or. &
          nsolp/=nsol .Or. strp/='REG') iflag = iflag + 1
        If (kap1/=kapin(1)) iflag = iflag + 1
        If (ikm1/=ikmin(1)) iflag = iflag + 1
        If (nsol>1) Then
          If (kap2/=kapin(2)) iflag = iflag + 1
          If (ikm2/=ikmin(2)) iflag = iflag + 1
        End If
      End If

!-----------------------------------------------------------------------
!                                  IRREGULAR wave function -- solution 1
      If (sirr=='IRR') Then
        Read (nfil, Rec=ikm1+(it-1+nt)*nkm) itp, lp, mjp, nsolp, strp, &
          (kapin(k), ikmin(k), (jg(i,k,1),jf(i,k,1),i=1,jtop), k=1, nsol)
        If (itp/=it .Or. lp/=l .Or. abs(mjp-mj)>0.001E0_dp .Or. &
          nsolp/=nsol .Or. strp/='IRR') iflag = iflag + 1
        If (kap1/=kapin(1)) iflag = iflag + 1
        If (ikm1/=ikmin(1)) iflag = iflag + 1
        If (nsol>1) Then
          If (kap2/=kapin(2)) iflag = iflag + 1
          If (ikm2/=ikmin(2)) iflag = iflag + 1
        End If
      End If

      If (nsol==2) Then
!-----------------------------------------------------------------------
!                                    REGULAR wave function -- solution 2
        If (sreg=='REG' .Or. sreg=='COR') Then
          Read (nfil, Rec=ikm2+(it-1)*nkm) itp, lp, mjp, nsolp, strp, &
            (kapin(k), ikmin(k), (zg(i,k,2),zf(i,k,2),i=1,jtop), k=2, 1, -1)
          If (itp/=it .Or. lp/=l .Or. abs(mjp-mj)>0.001E0_dp .Or. &
            nsolp/=nsol .Or. strp/='REG') iflag = iflag + 1
        End If

!-----------------------------------------------------------------------
!                                  IRREGULAR wave function -- solution 2
        If (sirr=='IRR') Then
          Read (nfil, Rec=ikm2+(it-1+nt)*nkm) itp, lp, mjp, nsolp, strp, &
            (kapin(k), ikmin(k), (jg(i,k,2),jf(i,k,2),i=1,jtop), k=2, 1, -1)
          If (itp/=it .Or. lp/=l .Or. abs(mjp-mj)>0.001E0_dp .Or. &
            nsolp/=nsol .Or. strp/='IRR') iflag = iflag + 1
        End If

      End If

!-----------------------------------------------------------------------

      If (iflag>0) Then
        Write (*, 100) iflag, it, l, mj, sreg, sirr
        Stop ' in <READWFUN>'
      End If
100   Format (/, /, 1X, 79('*'), /, 10X, 'error reading the wave functions', &
        ' IFLAG = 1', /, 10X, 'for  IT=', I2, ' L=', I2, ' MJ=', F4.1, &
        ' KEYS=', A, A)
    End Subroutine
