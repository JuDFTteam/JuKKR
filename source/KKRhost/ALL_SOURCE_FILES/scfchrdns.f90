module mod_scfchrdns

contains

subroutine scfchrdns(nfilcbwf, r2drdi, jws, imt, shftef, totdos, muespn, &
  mueorb, irel, iprint, nt, nl, nkm, eryd, we, efermi, iecurr, netab, dos, &
  smt, omt, hff, dosi, smti, omti, hffi, dosm, dosl0, dosint, smtm, smtl0, &
  smtint, omtm, omtl0, omtint, hffm, hffl0, hffint, bcor, bcors, dzz, dzj, &
  szz, szj, ozz, ozj, bzz, bzj, ozzs, ozjs, omtls0, tautlin, nvaltot, txtt, &
  conc, nat, rhochr, rhospn, rhoorb, qel, gdia, gmdia, goff, ntmax, nlmax, &
  nmuemax, linmax, nrmax, nmmax, nkmmax, eband, ebandt)
  ! ********************************************************************
  ! *                                                                  *
  ! * SUBROUTINE TO CALCULATE THE  CHARGE, SPIN  AND  ORBITAL DENSITY  *
  ! *                  WITHIN AN ATOMIC CELL                           *
  ! *                                                                  *
  ! * 12/03/96 HE                                                      *
  ! ********************************************************************
  use mod_datatypes, only: dp
  use mod_ssite, only: readwfun
  use mod_rintsimp
  use mod_ikapmue
  use mod_rinit
  implicit none

  ! PARAMETER definitions
  integer :: ntmaxchk
  parameter (ntmaxchk=10)
  complex (kind=dp) :: c0
  parameter (c0=(0.0e0_dp,0.0e0_dp))
  real (kind=dp) :: pi
  parameter (pi=3.141592653589793238462643e0_dp)

  ! Dummy arguments
  complex (kind=dp) :: eband, eryd, we
  real (kind=dp) :: efermi, mueorb, muespn, shftef, totdos, nvaltot
  integer :: iecurr, iprint, irel, linmax, netab, nfilcbwf, nkm, nkmmax, nl, &
    nlmax, nmmax, nmuemax, nrmax, nt, ntmax
  real (kind=dp) :: bcor(ntmax), bcors(ntmax), conc(ntmax), dos(ntmax), &
    dosi(ntmax), gdia(nkmmax), gmdia(nkmmax), goff(nkmmax), hff(ntmax), &
    hffi(ntmax), omt(ntmax), omti(ntmax), qel(ntmax), r2drdi(nrmax, nmmax), &
    rhochr(nrmax, ntmax), rhoorb(nrmax, ntmax), rhospn(nrmax, ntmax), &
    smt(ntmax), smti(ntmax)
  complex (kind=dp) :: bzj(linmax, ntmax), bzz(linmax, ntmax), &
    dosint(nlmax, ntmax), dosl0(nlmax, ntmax), dosm(nmuemax), &
    dzj(linmax, ntmax), dzz(linmax, ntmax), ebandt(ntmax), &
    hffint(nlmax, ntmax), hffl0(nlmax, ntmax), hffm(nmuemax), &
    omtint(nlmax, ntmax), omtl0(nlmax, ntmax), omtm(nmuemax), &
    ozj(linmax, ntmax), ozz(linmax, ntmax), smtint(nlmax, ntmax), &
    smtl0(nlmax, ntmax), smtm(nmuemax), szj(linmax, ntmax), &
    szz(linmax, ntmax), tautlin(linmax, ntmax)
  complex (kind=dp) :: ozzs(linmax, ntmax, 2), ozjs(linmax, ntmax, 2), &
    omtls0(nlmax, ntmax, 2)
  integer :: imt(ntmax), jws(nmmax), nat(ntmax)
  character (len=4) :: txtt(ntmax)

  ! Local variables
  real (kind=dp) :: aux, bdum(3), cff(2, 2), cfg(2, 2), cgf(2, 2), cgg(2, 2), &
    chko(ntmaxchk), chkq(ntmaxchk), chks(ntmaxchk), defermi, dq, mj, mjmax, &
    mjmin, r1m(2, 2), rint(nrmax), totnos
  complex (kind=dp) :: dosl, hffl, jf(nrmax, 2, 2), jg(nrmax, 2, 2), omtl, &
    smtl, wds, wof, wog, wsf, wsg, wt, zf(nrmax, 2, 2), zfjf, zfzf, &
    zg(nrmax, 2, 2), zgjg, zgzg
  complex (kind=dp) :: omtls(2), omtms(nmuemax, 2)
  integer :: i, iflag, ikm1, ikm2, il, im, is, it, jj, jtop, k1, k2, ka, kap1, &
    kap2, kb, l, lin, lmax, mm, mue, nsol
  integer :: nint

  save :: chko, chkq, chks

  data r1m/1.0e0_dp, 0.0e0_dp, 0.0e0_dp, 1.0e0_dp/

  if (iecurr==1) then

    do it = 1, nt
      do il = 1, nl
        dosint(il, it) = c0
        smtint(il, it) = c0
        omtint(il, it) = c0
        hffint(il, it) = c0
      end do
      ebandt(it) = c0
    end do

    ! ----------------------------- account for spin degeneracy for IREL <=1
    if (irel<=1) then
      do it = 1, nt
        im = imt(it)
        jtop = jws(im)
        do i = 1, jtop
          rhochr(i, it) = rhochr(i, it)/2.0e0_dp
          rhospn(i, it) = 0.0e0_dp
          rhoorb(i, it) = 0.0e0_dp
        end do
      end do
    end if

    if (iprint>0) then
      if (nt>ntmaxchk) stop '<SCFCHRDNS> NT > NTMAXCHK'
      do it = 1, nt
        im = imt(it)
        jtop = jws(im)
        do i = 1, jtop
          rint(i) = rhochr(i, it)*r2drdi(i, im)
        end do
        call rintsimp(rint, jtop, chkq(it))
        do i = 1, jtop
          rint(i) = rhospn(i, it)*r2drdi(i, im)
        end do
        call rintsimp(rint, jtop, chks(it))
        do i = 1, jtop
          rint(i) = rhoorb(i, it)*r2drdi(i, im)
        end do
        call rintsimp(rint, jtop, chko(it))
      end do
    end if
  end if

100 continue
  totnos = 0.0e0_dp
  totdos = 0.0e0_dp
  muespn = 0.0e0_dp
  mueorb = 0.0e0_dp


  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
  do it = 1, nt
    im = imt(it)
    jtop = jws(im)

    lmax = nl - 1
    lin = 0

    dos(it) = 0.0e0_dp
    smt(it) = 0.0e0_dp
    omt(it) = 0.0e0_dp
    hff(it) = 0.0e0_dp
    dosi(it) = 0.0e0_dp
    smti(it) = 0.0e0_dp
    omti(it) = 0.0e0_dp
    hffi(it) = 0.0e0_dp

    ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    do l = 0, lmax
      il = l + 1
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
      kap1 = -l - 1
      kap2 = l
      if (l==0) kap2 = kap1
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

      dosl = 0.0e0_dp
      smtl = 0.0e0_dp
      omtl = 0.0e0_dp
      hffl = 0.0e0_dp
      omtls(1) = omtl
      omtls(2) = omtl

      if (irel>1) then
        mjmax = real(l, kind=dp) + 0.5e0_dp
      else
        mjmax = real(l, kind=dp)
      end if
      mjmin = -mjmax
      mue = 0

      ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      do mj = mjmin,mjmax,1.0_dp
        mue = mue + 1
        dosm(mue) = 0.0e0_dp
        smtm(mue) = 0.0e0_dp
        omtm(mue) = 0.0e0_dp
        hffm(mue) = 0.0e0_dp

        do is = 1, 2
          omtms(mue, is) = 0.0e0_dp
        end do

        if (irel<=1) then
          nsol = 1
          ! no coupling for:  abs(mue)= j   +  j=l+1/2 == kap=-l-1
        else if (abs(mj)>real(l, kind=dp)) then
          nsol = 1
        else
          nsol = 2
        end if
        ! ------------------------------------------------------------------------
        ikm1 = ikapmue(kap1, nint(mj-0.5e0_dp))
        ikm2 = ikapmue(kap2, nint(mj-0.5e0_dp))
        ! ------------------------------------------------------------------------
        if (irel<=1) then
          ikm1 = il
          ikm2 = il
          if (nkm/=nl**2) write (1337, 110) nkm
        end if


        ! COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION

        cgg(1, 1) = gdia(ikm1)
        cgg(1, 2) = goff(ikm1)
        cgg(2, 1) = goff(ikm1)
        cgg(2, 2) = gdia(ikm2)
        call rinit(4, cgf)
        cgf(1, 1) = gmdia(ikm1)
        cgf(2, 2) = gmdia(ikm2)

        ! COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION

        cfg(1, 1) = mj*(kap1+1.0e0_dp)/(kap1+0.5e0_dp)
        cfg(2, 2) = mj*(kap2+1.0e0_dp)/(kap2+0.5e0_dp)
        cfg(1, 2) = 0.5e0_dp*sqrt(1.0e0_dp-(mj/(kap1+0.5e0_dp))**2)
        cfg(2, 1) = cfg(1, 2)
        call rinit(4, cff)
        cff(1, 1) = mj*(-kap1+1.0e0_dp)/(-kap1+0.5e0_dp)
        cff(2, 2) = mj*(-kap2+1.0e0_dp)/(-kap2+0.5e0_dp)

        ! -------------------------------------------------- read wave
        ! functions

        call readwfun(nfilcbwf, it, l, mj, nsol, 'REG', 'IRR', ikm1, kap1, &
          ikm2, kap2, nt, nkm, zg, zf, jg, jf, jtop, nrmax)

        do k1 = 1, nsol
          do k2 = 1, nsol
            lin = lin + 1
            wt = -tautlin(lin, it)/pi
            dosm(mue) = dosm(mue) + wt*dzz(lin, it)
            smtm(mue) = smtm(mue) + wt*szz(lin, it)
            omtm(mue) = omtm(mue) + wt*ozz(lin, it)
            hffm(mue) = hffm(mue) + wt*bzz(lin, it)

            do is = 1, 2
              omtms(mue, is) = omtms(mue, is) + wt*ozzs(lin, it, is)
            end do

            do ka = 1, nsol
              do kb = 1, nsol
                wds = we*wt*r1m(ka, kb)
                wsg = we*wt*cgg(ka, kb)
                wsf = we*wt*cgf(ka, kb)
                wog = we*wt*cfg(ka, kb)
                wof = we*wt*cff(ka, kb)
                do i = 1, jtop
                  zgzg = zg(i, ka, k1)*zg(i, kb, k2)
                  zfzf = zf(i, ka, k1)*zf(i, kb, k2)
                  rhochr(i, it) = rhochr(i, it) + aimag(wds*zgzg+wds*zfzf)
                  rhospn(i, it) = rhospn(i, it) + aimag(wsg*zgzg-wsf*zfzf)
                  rhoorb(i, it) = rhoorb(i, it) + aimag(wog*zgzg-wof*zfzf)
                end do
              end do
            end do

            ! NO IRREGULAR CONTRIBUTIONS TO THE BACKSCATTERING TERMS
            if (k1==k2) then
              dosm(mue) = dosm(mue) + dzj(lin, it)/pi
              smtm(mue) = smtm(mue) + szj(lin, it)/pi
              omtm(mue) = omtm(mue) + ozj(lin, it)/pi
              hffm(mue) = hffm(mue) + bzj(lin, it)/pi

              do is = 1, 2
                omtms(mue, is) = omtms(mue, is) + ozjs(lin, it, is)/pi
              end do

              do ka = 1, nsol
                do kb = 1, nsol
                  wds = we*r1m(ka, kb)/pi
                  wsg = we*cgg(ka, kb)/pi
                  wsf = we*cgf(ka, kb)/pi
                  wog = we*cfg(ka, kb)/pi
                  wof = we*cff(ka, kb)/pi
                  do i = 1, jtop
                    zgjg = zg(i, ka, k1)*jg(i, kb, k2)
                    zfjf = zf(i, ka, k1)*jf(i, kb, k2)
                    rhochr(i, it) = rhochr(i, it) + aimag(wds*zgjg+wds*zfjf)
                    rhospn(i, it) = rhospn(i, it) + aimag(wsg*zgjg-wsf*zfjf)
                    rhoorb(i, it) = rhoorb(i, it) + aimag(wog*zgjg-wof*zfjf)
                  end do
                end do
              end do
            end if
          end do
        end do


        dosl = dosl + dosm(mue)
        smtl = smtl + smtm(mue)
        omtl = omtl + omtm(mue)
        hffl = hffl + hffm(mue)

        do is = 1, 2
          omtls(is) = omtls(is) + omtms(mue, is)
        end do

      end do
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

      do is = 1, 2
        omtls0(il, it, is) = omtls(is)
      end do
      ! ----------------------------------------------------------------------
      if ((iprint>0) .and. (iecurr==netab)) then

        if (irel>1) then
          jj = 2*l + 2

          write (1337, 150) iecurr, eryd, l, it, txtt(it), &
            'CRYSTAL TERMS       ', dosint(il, it), smtint(il, it), &
            omtint(il, it), (hffint(il,it)*1e-3_dp), dosl, smtl, omtl, &
            (hffl*1e-6_dp), ((-jj-1+2*mue), dosm(mue), smtm(mue), omtm(mue), ( &
            hffm(mue)*1e-6_dp), mue=1, jj)

        else
          jj = 2*l + 1

          write (1337, 240) iecurr, eryd, l, it, txtt(it), &
            'CRYSTAL TERMS       ', dosint(il, it), smtint(il, it), &
            omtint(il, it), (hffint(il,it)*1e-3_dp), dosl, smtl, omtl, &
            (hffl*1e-6_dp), ((-l-1+mm), dosm(mm), smtm(mm), omtm(mm), (hffm(mm &
            )*1e-6_dp), mm=1, jj)

        end if
      end if
      ! ---------------------------------------------------------------------


      dos(it) = dos(it) + aimag(dosl)
      smt(it) = smt(it) + aimag(smtl)
      omt(it) = omt(it) + aimag(omtl)
      hff(it) = hff(it) + aimag(hffl)
      dosi(it) = dosi(it) + aimag(dosint(il,it))
      smti(it) = smti(it) + aimag(smtint(il,it))
      omti(it) = omti(it) + aimag(omtint(il,it))
      hffi(it) = hffi(it) + aimag(hffint(il,it))
    end do
    ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL


    totnos = totnos + dosi(it)*conc(it)*nat(it)
    totdos = totdos + dos(it)*conc(it)*nat(it)
    muespn = muespn + smti(it)*conc(it)*nat(it)
    mueorb = mueorb + omti(it)*conc(it)*nat(it)

    if (iprint>0) then

      write (1337, 160) iecurr, eryd, it, txtt(it), dosi(it), smti(it), &
        omti(it), (hffi(it)*1e-3_dp), dos(it), smt(it), omt(it), &
        (hff(it)*1e-3_dp)

      if (it<nt) then
        write (1337, '(1X,79(''-''))')
      else if ((iprint>0) .or. (iecurr==netab)) then
        write (1337, 200) totdos, totnos, muespn, mueorb
      else
        write (1337, '('' '',79(''=''))')
      end if
    end if

  end do
  ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

  if (iecurr==netab) then

    eband = c0
    do it = 1, nt
      eband = eband + ebandt(it)*conc(it)*nat(it)
    end do

    if (irel>1) then
      dq = nvaltot - totnos
    else
      dq = nvaltot/2.0e0_dp - totnos
    end if

    defermi = dq/totdos

    if (abs(dq)>1e-06_dp) then
      we = defermi
      efermi = efermi + defermi
      shftef = defermi

      write (1337, '(/)')
      write (1337, 220)(txtt(it), conc(it), it=1, nt)
      write (1337, 230) dq, defermi, efermi

      go to 100
    end if

    if (iprint>0) then
      iflag = 0
      do it = 1, nt
        im = imt(it)
        jtop = jws(im)
        do i = 1, jtop
          rint(i) = rhochr(i, it)*r2drdi(i, im)
        end do
        call rintsimp(rint, jtop, aux)
        chkq(it) = aux - chkq(it)
        if (abs(chkq(it)-dosi(it))>1.0e-8_dp) then
          iflag = 1
          write (1337, 140) it, 'Q', dosi(it), chkq(it), dosi(it)/chkq(it)
        end if
        do i = 1, jtop
          rint(i) = rhospn(i, it)*r2drdi(i, im)
        end do
        call rintsimp(rint, jtop, aux)
        chks(it) = aux - chks(it)
        if (abs(chks(it)-smti(it))>1.0e-8_dp) then
          iflag = 1
          write (1337, 140) it, 'S', smti(it), chks(it), smti(it)/chks(it)
        end if
        do i = 1, jtop
          rint(i) = rhoorb(i, it)*r2drdi(i, im)
        end do
        call rintsimp(rint, jtop, aux)
        chko(it) = aux - chko(it)
        if (abs(chko(it)-omti(it))>1.0e-8_dp) then
          iflag = 1
          write (1337, 140) it, 'O', omti(it), chko(it), omti(it)/chko(it)
        end if
      end do

      if (iflag==0) then
        write (1337, 120)
      else
        write (1337, 130)
      end if
    end if


    do it = 1, nt

      write (1337, 160)(iecurr+1), efermi, 0.0e0_dp, it, txtt(it)

      bdum(1) = bcors(it)*1e-3_dp
      bdum(2) = (bcor(it)-bcors(it))*1e-3_dp
      bdum(3) = bcor(it)*1e-3_dp

      write (1337, 170)(aimag(dosl0(il,it)), aimag(dosint(il, &
        it)), aimag(smtl0(il,it)), aimag(smtint(il,it)), aimag(omtl0(il, &
        it)), aimag(omtint(il,it)), aimag(hffint(il, &
        it))*1e-3_dp, bdum(il), il=1, min(3,nl))
      if (nl>3) write (1337, 180)(aimag(dosl0(il,it)), aimag(dosint(il, &
        it)), aimag(smtl0(il,it)), aimag(smtint(il,it)), aimag(omtl0(il, &
        it)), aimag(omtint(il,it)), aimag(hffint(il,it))*1e-3_dp, il=4, nl)

      write (1337, 190) dos(it), dosi(it), smt(it), smti(it), omt(it), &
        omti(it), (hffi(it)*1e-3_dp), ((hffi(it)+bcor(it))*1e-3_dp)

      if (it<nt) then
        write (1337, '(1X,79(''-''))')
      else
        write (1337, 210) totdos, totnos, muespn, mueorb, aimag(eband)
      end if

      im = imt(it)
      jtop = jws(im)

      do i = 1, jtop
        rint(i) = rhochr(i, it)*r2drdi(i, im)
      end do

      call rintsimp(rint, jtop, qel(it))

      ! ----------------------------- account for spin degeneracy for IREL <=1
      if (irel<=1) then
        qel(it) = qel(it)*2.0e0_dp
        do i = 1, jtop
          rhochr(i, it) = rhochr(i, it)*2.0e0_dp
          rhospn(i, it) = 0.0e0_dp
          rhoorb(i, it) = 0.0e0_dp
        end do
      end if

    end do

  end if

110 format ('warning in <SCFCHRDNS>:  IREL<=1 and  NL**2 <> NKM=', i5)
120 format (/, 10x, 'integrals in <SCFCHRDNS> agree within 1D-8', /)
130 format (/, 10x, '... integrals in <SCFCHRDNS>  NOT OK')
140 format (' IT ', i3, 2x, a, 2x, f20.10, /, 12x, 4f20.10)
150 format (/, i4, ' E=', 2f7.4, 3x, 'L=', i2, 3x, 'IT=', i2, 2x, a, 2x, a20, &
    /, 15x, 'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |', &
    '   B_tot   [kG]', /, ' INT(DE)  ', 2f8.3, 2x, 2f8.3, 2x, 2f8.3, f10.1, &
    f8.1, /, ' SUM(MJ)  ', 2f8.3, 2x, 2f8.3, 2x, 2f8.3, f10.1, f8.1, &
    20(:,/,' MJ= ',i2,'/2 ',2f8.3,2x,2f8.3,2x,2f8.3,f10.1,f8.1))
160 format (/, i4, ' E=', 2f7.4, 10x, 'IT=', i2, 2x, a, :, /, 15x, &
    'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |', '   B_tot   [kG]', &
    /, ' INT(DE) crystal  ', f8.3, 10x, f8.3, 10x, f8.3, 10x, f8.1, /, &
    ' TOTAL   crystal  ', f8.3, 10x, f8.3, 10x, f8.3, 10x, f8.1)
170 format ('         DOS      NOS     P_spin   m_spin', &
    '    P_orb    m_orb    B_val      B_core', /, '  s ', 2f9.4, f10.4, f9.4, &
    f10.5, f9.5, f8.2, ' s  ', f8.2, :, /, '  p ', 2f9.4, f10.4, f9.4, f10.5, &
    f9.5, f8.2, ' ns ', f8.2, :, /, '  d ', 2f9.4, f10.4, f9.4, f10.5, f9.5, &
    f8.2, ' cor', f8.2)
180 format ('  f ', 2f9.4, f10.4, f9.4, f10.5, f9.5, f8.2, :, /, '  g ', &
    2f9.4, f10.4, f9.4, f10.5, f9.5, f8.2)
190 format (' sum', 2f9.4, f10.4, f9.4, f10.5, f9.5, f8.2, ' v+c', f8.2)
200 format (' ', 79('-'), /, ' TDOS/NOS ', 2f8.3, ' MUE-SPIN:', f8.3, &
    '  MUE-ORB:', f8.3)
210 format (' ', 79('-'), /, ' TOT', 2f9.4, 10x, f9.4, 10x, f9.5, /, &
    ' E_band', f17.6, ' [Ry]', /, ' ', 79('='))
220 format ((' ',79('*'),/), /, ' KKR-run for: ', 15(a,f5.2))
230 format (/, ' results extrapolated to corrected FERMI - ENERGY:', /, &
    ' CHARGE MISFIT     ', f9.5, /, ' E_F CORRECTION    ', f9.5, /, &
    ' NEW FERMI ENERGY  ', f9.5, /)
240 format (/, i4, ' E=', 2f7.4, 3x, 'L=', i2, 3x, 'IT=', i2, 2x, a, 2x, a20, &
    /, 15x, 'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |', &
    '   B_tot   [kG]', /, ' INT(DE)  ', 2f8.3, 2x, 2f8.3, 2x, 2f8.3, f10.1, &
    f8.1, /, ' SUM(ML)  ', 2f8.3, 2x, 2f8.3, 2x, 2f8.3, f10.1, f8.1, &
    20(:,/,' ML= ',i2,'   ',2f8.3,2x,2f8.3,2x,2f8.3,f10.1,f8.1))
end subroutine scfchrdns

end module mod_scfchrdns
