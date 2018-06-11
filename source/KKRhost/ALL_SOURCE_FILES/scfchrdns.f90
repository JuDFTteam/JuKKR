    Subroutine scfchrdns(nfilcbwf, r2drdi, jws, imt, shftef, totdos, muespn, &
      mueorb, irel, iprint, nt, nl, nkm, eryd, we, efermi, iecurr, netab, dos, &
      smt, omt, hff, dosi, smti, omti, hffi, dosm, dosl0, dosint, smtm, smtl0, &
      smtint, omtm, omtl0, omtint, hffm, hffl0, hffint, bcor, bcors, dzz, dzj, &
      szz, szj, ozz, ozj, bzz, bzj, ozzs, ozjs, omtls0, tautlin, nvaltot, &
      txtt, conc, nat, rhochr, rhospn, rhoorb, qel, gdia, gmdia, goff, ntmax, &
      nlmax, nmuemax, linmax, nrmax, nmmax, nkmmax, eband, ebandt)
      Use mod_datatypes, Only: dp
!   ********************************************************************
!   *                                                                  *
!   * SUBROUTINE TO CALCULATE THE  CHARGE, SPIN  AND  ORBITAL DENSITY  *
!   *                  WITHIN AN ATOMIC CELL                           *
!   *                                                                  *
!   * 12/03/96 HE                                                      *
!   ********************************************************************
      Implicit None

! PARAMETER definitions
      Integer :: ntmaxchk
      Parameter (ntmaxchk=10)
      Complex (Kind=dp) :: c0
      Parameter (c0=(0.0E0_dp,0.0E0_dp))
      Real (Kind=dp) :: pi
      Parameter (pi=3.141592653589793238462643E0_dp)

! Dummy arguments
      Complex (Kind=dp) :: eband, eryd, we
      Real (Kind=dp) :: efermi, mueorb, muespn, shftef, totdos, nvaltot
      Integer :: iecurr, iprint, irel, linmax, netab, nfilcbwf, nkm, nkmmax, &
        nl, nlmax, nmmax, nmuemax, nrmax, nt, ntmax
      Real (Kind=dp) :: bcor(ntmax), bcors(ntmax), conc(ntmax), dos(ntmax), &
        dosi(ntmax), gdia(nkmmax), gmdia(nkmmax), goff(nkmmax), hff(ntmax), &
        hffi(ntmax), omt(ntmax), omti(ntmax), qel(ntmax), &
        r2drdi(nrmax, nmmax), rhochr(nrmax, ntmax), rhoorb(nrmax, ntmax), &
        rhospn(nrmax, ntmax), smt(ntmax), smti(ntmax)
      Complex (Kind=dp) :: bzj(linmax, ntmax), bzz(linmax, ntmax), &
        dosint(nlmax, ntmax), dosl0(nlmax, ntmax), dosm(nmuemax), &
        dzj(linmax, ntmax), dzz(linmax, ntmax), ebandt(ntmax), &
        hffint(nlmax, ntmax), hffl0(nlmax, ntmax), hffm(nmuemax), &
        omtint(nlmax, ntmax), omtl0(nlmax, ntmax), omtm(nmuemax), &
        ozj(linmax, ntmax), ozz(linmax, ntmax), smtint(nlmax, ntmax), &
        smtl0(nlmax, ntmax), smtm(nmuemax), szj(linmax, ntmax), &
        szz(linmax, ntmax), tautlin(linmax, ntmax)
      Complex (Kind=dp) :: ozzs(linmax, ntmax, 2), ozjs(linmax, ntmax, 2), &
        omtls0(nlmax, ntmax, 2)
      Integer :: imt(ntmax), jws(nmmax), nat(ntmax)
      Character (Len=4) :: txtt(ntmax)

! Local variables
      Real (Kind=dp) :: aux, bdum(3), cff(2, 2), cfg(2, 2), cgf(2, 2), &
        cgg(2, 2), chko(ntmaxchk), chkq(ntmaxchk), chks(ntmaxchk), defermi, &
        dq, mj, mjmax, mjmin, r1m(2, 2), rint(nrmax), totnos
      Real (Kind=dp) :: dble, dsqrt
      Complex (Kind=dp) :: dosl, hffl, jf(nrmax, 2, 2), jg(nrmax, 2, 2), omtl, &
        smtl, wds, wof, wog, wsf, wsg, wt, zf(nrmax, 2, 2), zfjf, zfzf, &
        zg(nrmax, 2, 2), zgjg, zgzg
      Complex (Kind=dp) :: omtls(2), omtms(nmuemax, 2)
      Integer :: i, iflag, ikm1, ikm2, il, im, is, it, jj, jtop, k1, k2, ka, &
        kap1, kap2, kb, l, lin, lmax, mm, mue, nsol
      Integer :: ikapmue, imj
      Integer :: nint

      Save :: chko, chkq, chks

      Data r1m/1.0E0_dp, 0.0E0_dp, 0.0E0_dp, 1.0E0_dp/

      If (iecurr==1) Then

        Do it = 1, nt
          Do il = 1, nl
            dosint(il, it) = c0
            smtint(il, it) = c0
            omtint(il, it) = c0
            hffint(il, it) = c0
          End Do
          ebandt(it) = c0
        End Do

! ----------------------------- account for spin degeneracy for IREL <=1
        If (irel<=1) Then
          Do it = 1, nt
            im = imt(it)
            jtop = jws(im)
            Do i = 1, jtop
              rhochr(i, it) = rhochr(i, it)/2.0E0_dp
              rhospn(i, it) = 0.0E0_dp
              rhoorb(i, it) = 0.0E0_dp
            End Do
          End Do
        End If

        If (iprint>0) Then
          If (nt>ntmaxchk) Stop '<SCFCHRDNS> NT > NTMAXCHK'
          Do it = 1, nt
            im = imt(it)
            jtop = jws(im)
            Do i = 1, jtop
              rint(i) = rhochr(i, it)*r2drdi(i, im)
            End Do
            Call rintsimp(rint, jtop, chkq(it))
            Do i = 1, jtop
              rint(i) = rhospn(i, it)*r2drdi(i, im)
            End Do
            Call rintsimp(rint, jtop, chks(it))
            Do i = 1, jtop
              rint(i) = rhoorb(i, it)*r2drdi(i, im)
            End Do
            Call rintsimp(rint, jtop, chko(it))
          End Do
        End If
      End If

100   Continue
      totnos = 0.0E0_dp
      totdos = 0.0E0_dp
      muespn = 0.0E0_dp
      mueorb = 0.0E0_dp


! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      Do it = 1, nt
        im = imt(it)
        jtop = jws(im)

        lmax = nl - 1
        lin = 0

        dos(it) = 0.0E0_dp
        smt(it) = 0.0E0_dp
        omt(it) = 0.0E0_dp
        hff(it) = 0.0E0_dp
        dosi(it) = 0.0E0_dp
        smti(it) = 0.0E0_dp
        omti(it) = 0.0E0_dp
        hffi(it) = 0.0E0_dp

! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
        Do l = 0, lmax
          il = l + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          kap1 = -l - 1
          kap2 = l
          If (l==0) kap2 = kap1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          dosl = 0.0E0_dp
          smtl = 0.0E0_dp
          omtl = 0.0E0_dp
          hffl = 0.0E0_dp
          omtls(1) = omtl
          omtls(2) = omtl

          If (irel>1) Then
            mjmax = dble(l) + 0.5E0_dp
          Else
            mjmax = dble(l)
          End If
          mjmin = -mjmax
          mue = 0

! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
!DO MJ = MJMIN,MJMAX,1.0D0
          Do imj = nint(mjmin), nint(mjmax), 1
            mj = dble(imj)
            mue = mue + 1
            dosm(mue) = 0.0E0_dp
            smtm(mue) = 0.0E0_dp
            omtm(mue) = 0.0E0_dp
            hffm(mue) = 0.0E0_dp

            Do is = 1, 2
              omtms(mue, is) = 0.0E0_dp
            End Do

            If (irel<=1) Then
              nsol = 1
!           no coupling for:  abs(mue)= j   +  j=l+1/2 == kap=-l-1
            Else If (abs(mj)>dble(l)) Then
              nsol = 1
            Else
              nsol = 2
            End If
!------------------------------------------------------------------------
            ikm1 = ikapmue(kap1, nint(mj-0.5E0_dp))
            ikm2 = ikapmue(kap2, nint(mj-0.5E0_dp))
!------------------------------------------------------------------------
            If (irel<=1) Then
              ikm1 = il
              ikm2 = il
              If (nkm/=nl**2) Write (1337, 110) nkm
            End If


!   COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION

            cgg(1, 1) = gdia(ikm1)
            cgg(1, 2) = goff(ikm1)
            cgg(2, 1) = goff(ikm1)
            cgg(2, 2) = gdia(ikm2)
            Call rinit(4, cgf)
            cgf(1, 1) = gmdia(ikm1)
            cgf(2, 2) = gmdia(ikm2)

!   COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION

            cfg(1, 1) = mj*(kap1+1.0E0_dp)/(kap1+0.5E0_dp)
            cfg(2, 2) = mj*(kap2+1.0E0_dp)/(kap2+0.5E0_dp)
            cfg(1, 2) = 0.5E0_dp*dsqrt(1.0E0_dp-(mj/(kap1+0.5E0_dp))**2)
            cfg(2, 1) = cfg(1, 2)
            Call rinit(4, cff)
            cff(1, 1) = mj*(-kap1+1.0E0_dp)/(-kap1+0.5E0_dp)
            cff(2, 2) = mj*(-kap2+1.0E0_dp)/(-kap2+0.5E0_dp)

! -------------------------------------------------- read wave functions

            Call readwfun(nfilcbwf, it, l, mj, nsol, 'REG', 'IRR', ikm1, kap1, &
              ikm2, kap2, nt, nkm, zg, zf, jg, jf, jtop, nrmax)

            Do k1 = 1, nsol
              Do k2 = 1, nsol
                lin = lin + 1
                wt = -tautlin(lin, it)/pi
                dosm(mue) = dosm(mue) + wt*dzz(lin, it)
                smtm(mue) = smtm(mue) + wt*szz(lin, it)
                omtm(mue) = omtm(mue) + wt*ozz(lin, it)
                hffm(mue) = hffm(mue) + wt*bzz(lin, it)

                Do is = 1, 2
                  omtms(mue, is) = omtms(mue, is) + wt*ozzs(lin, it, is)
                End Do

                Do ka = 1, nsol
                  Do kb = 1, nsol
                    wds = we*wt*r1m(ka, kb)
                    wsg = we*wt*cgg(ka, kb)
                    wsf = we*wt*cgf(ka, kb)
                    wog = we*wt*cfg(ka, kb)
                    wof = we*wt*cff(ka, kb)
                    Do i = 1, jtop
                      zgzg = zg(i, ka, k1)*zg(i, kb, k2)
                      zfzf = zf(i, ka, k1)*zf(i, kb, k2)
                      rhochr(i, it) = rhochr(i, it) + aimag(wds*zgzg+wds*zfzf)
                      rhospn(i, it) = rhospn(i, it) + aimag(wsg*zgzg-wsf*zfzf)
                      rhoorb(i, it) = rhoorb(i, it) + aimag(wog*zgzg-wof*zfzf)
                    End Do
                  End Do
                End Do

!    NO IRREGULAR CONTRIBUTIONS TO THE BACKSCATTERING TERMS
                If (k1==k2) Then
                  dosm(mue) = dosm(mue) + dzj(lin, it)/pi
                  smtm(mue) = smtm(mue) + szj(lin, it)/pi
                  omtm(mue) = omtm(mue) + ozj(lin, it)/pi
                  hffm(mue) = hffm(mue) + bzj(lin, it)/pi

                  Do is = 1, 2
                    omtms(mue, is) = omtms(mue, is) + ozjs(lin, it, is)/pi
                  End Do

                  Do ka = 1, nsol
                    Do kb = 1, nsol
                      wds = we*r1m(ka, kb)/pi
                      wsg = we*cgg(ka, kb)/pi
                      wsf = we*cgf(ka, kb)/pi
                      wog = we*cfg(ka, kb)/pi
                      wof = we*cff(ka, kb)/pi
                      Do i = 1, jtop
                        zgjg = zg(i, ka, k1)*jg(i, kb, k2)
                        zfjf = zf(i, ka, k1)*jf(i, kb, k2)
                        rhochr(i, it) = rhochr(i, it) + &
                          aimag(wds*zgjg+wds*zfjf)
                        rhospn(i, it) = rhospn(i, it) + &
                          aimag(wsg*zgjg-wsf*zfjf)
                        rhoorb(i, it) = rhoorb(i, it) + &
                          aimag(wog*zgjg-wof*zfjf)
                      End Do
                    End Do
                  End Do
                End If
              End Do
            End Do


            dosl = dosl + dosm(mue)
            smtl = smtl + smtm(mue)
            omtl = omtl + omtm(mue)
            hffl = hffl + hffm(mue)

            Do is = 1, 2
              omtls(is) = omtls(is) + omtms(mue, is)
            End Do

          End Do
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

          ebandt(it) = ebandt(it) + we*dosl*eryd

          dosint(il, it) = dosint(il, it) + we*dosl
          smtint(il, it) = smtint(il, it) + we*smtl
          omtint(il, it) = omtint(il, it) + we*omtl
          hffint(il, it) = hffint(il, it) + we*hffl

          dosl0(il, it) = dosl
          smtl0(il, it) = smtl
          omtl0(il, it) = omtl
          hffl0(il, it) = hffl

          Do is = 1, 2
            omtls0(il, it, is) = omtls(is)
          End Do
! ----------------------------------------------------------------------
          If ((iprint>0) .And. (iecurr==netab)) Then

            If (irel>1) Then
              jj = 2*l + 2

              Write (1337, 150) iecurr, eryd, l, it, txtt(it), &
                'CRYSTAL TERMS       ', dosint(il, it), smtint(il, it), &
                omtint(il, it), (hffint(il,it)*1E-3_dp), dosl, smtl, omtl, &
                (hffl*1E-6_dp), ((-jj-1+2*mue), dosm(mue), smtm(mue), omtm(mue &
                ), (hffm(mue)*1E-6_dp), mue=1, jj)

            Else
              jj = 2*l + 1

              Write (1337, 240) iecurr, eryd, l, it, txtt(it), &
                'CRYSTAL TERMS       ', dosint(il, it), smtint(il, it), &
                omtint(il, it), (hffint(il,it)*1E-3_dp), dosl, smtl, omtl, &
                (hffl*1E-6_dp), ((-l-1+mm), dosm(mm), smtm(mm), omtm(mm), ( &
                hffm(mm)*1E-6_dp), mm=1, jj)

            End If
          End If
! ---------------------------------------------------------------------


          dos(it) = dos(it) + aimag(dosl)
          smt(it) = smt(it) + aimag(smtl)
          omt(it) = omt(it) + aimag(omtl)
          hff(it) = hff(it) + aimag(hffl)
          dosi(it) = dosi(it) + aimag(dosint(il,it))
          smti(it) = smti(it) + aimag(smtint(il,it))
          omti(it) = omti(it) + aimag(omtint(il,it))
          hffi(it) = hffi(it) + aimag(hffint(il,it))
        End Do
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL


        totnos = totnos + dosi(it)*conc(it)*nat(it)
        totdos = totdos + dos(it)*conc(it)*nat(it)
        muespn = muespn + smti(it)*conc(it)*nat(it)
        mueorb = mueorb + omti(it)*conc(it)*nat(it)

        If (iprint>0) Then

          Write (1337, 160) iecurr, eryd, it, txtt(it), dosi(it), smti(it), &
            omti(it), (hffi(it)*1E-3_dp), dos(it), smt(it), omt(it), &
            (hff(it)*1E-3_dp)

          If (it<nt) Then
            Write (1337, '(1X,79(''-''))')
          Else If ((iprint>0) .Or. (iecurr==netab)) Then
            Write (1337, 200) totdos, totnos, muespn, mueorb
          Else
            Write (1337, '('' '',79(''=''))')
          End If
        End If

      End Do
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

      If (iecurr==netab) Then

        eband = c0
        Do it = 1, nt
          eband = eband + ebandt(it)*conc(it)*nat(it)
        End Do

        If (irel>1) Then
          dq = nvaltot - totnos
        Else
          dq = nvaltot/2.0E0_dp - totnos
        End If

        defermi = dq/totdos

        If (abs(dq)>1E-06_dp) Then
          we = defermi
          efermi = efermi + defermi
          shftef = defermi

          Write (1337, '(/)')
          Write (1337, 220)(txtt(it), conc(it), it=1, nt)
          Write (1337, 230) dq, defermi, efermi

          Go To 100
        End If

        If (iprint>0) Then
          iflag = 0
          Do it = 1, nt
            im = imt(it)
            jtop = jws(im)
            Do i = 1, jtop
              rint(i) = rhochr(i, it)*r2drdi(i, im)
            End Do
            Call rintsimp(rint, jtop, aux)
            chkq(it) = aux - chkq(it)
            If (abs(chkq(it)-dosi(it))>1.0E-8_dp) Then
              iflag = 1
              Write (1337, 140) it, 'Q', dosi(it), chkq(it), dosi(it)/chkq(it)
            End If
            Do i = 1, jtop
              rint(i) = rhospn(i, it)*r2drdi(i, im)
            End Do
            Call rintsimp(rint, jtop, aux)
            chks(it) = aux - chks(it)
            If (abs(chks(it)-smti(it))>1.0E-8_dp) Then
              iflag = 1
              Write (1337, 140) it, 'S', smti(it), chks(it), smti(it)/chks(it)
            End If
            Do i = 1, jtop
              rint(i) = rhoorb(i, it)*r2drdi(i, im)
            End Do
            Call rintsimp(rint, jtop, aux)
            chko(it) = aux - chko(it)
            If (abs(chko(it)-omti(it))>1.0E-8_dp) Then
              iflag = 1
              Write (1337, 140) it, 'O', omti(it), chko(it), omti(it)/chko(it)
            End If
          End Do

          If (iflag==0) Then
            Write (1337, 120)
          Else
            Write (1337, 130)
          End If
        End If


        Do it = 1, nt

          Write (1337, 160)(iecurr+1), efermi, 0.0E0_dp, it, txtt(it)

          bdum(1) = bcors(it)*1E-3_dp
          bdum(2) = (bcor(it)-bcors(it))*1E-3_dp
          bdum(3) = bcor(it)*1E-3_dp

          Write (1337, 170)(aimag(dosl0(il,it)), aimag(dosint(il, &
            it)), aimag(smtl0(il,it)), aimag(smtint(il,it)), aimag(omtl0(il, &
            it)), aimag(omtint(il,it)), aimag(hffint(il, &
            it))*1E-3_dp, bdum(il), il=1, min(3,nl))
          If (nl>3) Write (1337, 180)(aimag(dosl0(il,it)), aimag(dosint(il, &
            it)), aimag(smtl0(il,it)), aimag(smtint(il,it)), aimag(omtl0(il, &
            it)), aimag(omtint(il,it)), aimag(hffint(il, &
            it))*1E-3_dp, il=4, nl)

          Write (1337, 190) dos(it), dosi(it), smt(it), smti(it), omt(it), &
            omti(it), (hffi(it)*1E-3_dp), ((hffi(it)+bcor(it))*1E-3_dp)

          If (it<nt) Then
            Write (1337, '(1X,79(''-''))')
          Else
            Write (1337, 210) totdos, totnos, muespn, mueorb, aimag(eband)
          End If

          im = imt(it)
          jtop = jws(im)

          Do i = 1, jtop
            rint(i) = rhochr(i, it)*r2drdi(i, im)
          End Do

          Call rintsimp(rint, jtop, qel(it))

! ----------------------------- account for spin degeneracy for IREL <=1
          If (irel<=1) Then
            qel(it) = qel(it)*2.0E0_dp
            Do i = 1, jtop
              rhochr(i, it) = rhochr(i, it)*2.0E0_dp
              rhospn(i, it) = 0.0E0_dp
              rhoorb(i, it) = 0.0E0_dp
            End Do
          End If

        End Do

      End If

110   Format ('warning in <SCFCHRDNS>:  IREL<=1 and  NL**2 <> NKM=', I5)
120   Format (/, 10X, 'integrals in <SCFCHRDNS> agree within 1D-8', /)
130   Format (/, 10X, '... integrals in <SCFCHRDNS>  NOT OK')
140   Format (' IT ', I3, 2X, A, 2X, F20.10, /, 12X, 4F20.10)
150   Format (/, I4, ' E=', 2F7.4, 3X, 'L=', I2, 3X, 'IT=', I2, 2X, A, 2X, &
        A20, /, 15X, 'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |', &
        '   B_tot   [kG]', /, ' INT(DE)  ', 2F8.3, 2X, 2F8.3, 2X, 2F8.3, &
        F10.1, F8.1, /, ' SUM(MJ)  ', 2F8.3, 2X, 2F8.3, 2X, 2F8.3, F10.1, &
        F8.1, 20(:,/,' MJ= ',I2,'/2 ',2F8.3,2X,2F8.3,2X,2F8.3,F10.1,F8.1))
160   Format (/, I4, ' E=', 2F7.4, 10X, 'IT=', I2, 2X, A, :, /, 15X, &
        'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |', &
        '   B_tot   [kG]', /, ' INT(DE) crystal  ', F8.3, 10X, F8.3, 10X, &
        F8.3, 10X, F8.1, /, ' TOTAL   crystal  ', F8.3, 10X, F8.3, 10X, F8.3, &
        10X, F8.1)
170   Format ('         DOS      NOS     P_spin   m_spin', &
        '    P_orb    m_orb    B_val      B_core', /, '  s ', 2F9.4, F10.4, &
        F9.4, F10.5, F9.5, F8.2, ' s  ', F8.2, :, /, '  p ', 2F9.4, F10.4, &
        F9.4, F10.5, F9.5, F8.2, ' ns ', F8.2, :, /, '  d ', 2F9.4, F10.4, &
        F9.4, F10.5, F9.5, F8.2, ' cor', F8.2)
180   Format ('  f ', 2F9.4, F10.4, F9.4, F10.5, F9.5, F8.2, :, /, '  g ', &
        2F9.4, F10.4, F9.4, F10.5, F9.5, F8.2)
190   Format (' sum', 2F9.4, F10.4, F9.4, F10.5, F9.5, F8.2, ' v+c', F8.2)
200   Format (' ', 79('-'), /, ' TDOS/NOS ', 2F8.3, ' MUE-SPIN:', F8.3, &
        '  MUE-ORB:', F8.3)
210   Format (' ', 79('-'), /, ' TOT', 2F9.4, 10X, F9.4, 10X, F9.5, /, &
        ' E_band', F17.6, ' [Ry]', /, ' ', 79('='))
220   Format ((' ',79('*'),/), /, ' KKR-run for: ', 15(A,F5.2))
230   Format (/, ' results extrapolated to corrected FERMI - ENERGY:', /, &
        ' CHARGE MISFIT     ', F9.5, /, ' E_F CORRECTION    ', F9.5, /, &
        ' NEW FERMI ENERGY  ', F9.5, /)
240   Format (/, I4, ' E=', 2F7.4, 3X, 'L=', I2, 3X, 'IT=', I2, 2X, A, 2X, &
        A20, /, 15X, 'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |', &
        '   B_tot   [kG]', /, ' INT(DE)  ', 2F8.3, 2X, 2F8.3, 2X, 2F8.3, &
        F10.1, F8.1, /, ' SUM(ML)  ', 2F8.3, 2X, 2F8.3, 2X, 2F8.3, F10.1, &
        F8.1, 20(:,/,' ML= ',I2,'   ',2F8.3,2X,2F8.3,2X,2F8.3,F10.1,F8.1))
    End Subroutine
