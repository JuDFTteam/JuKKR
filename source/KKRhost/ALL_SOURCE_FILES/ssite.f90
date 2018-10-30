!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Assign quantum numbers and call routine to solve 8 coupled differentail radial Dirac eqautions.
!> Author: 
!> Assign quantum numbers and call routine to solve 8 coupled differentail radial 
!> Dirac eqautions. The resulting wavefunctions are used to calculate t-matrices in
!> the \(\kappa-\mu\) representation.
!> Calculation of the radial integrals:
!> \begin{equation}
!> \left[ G_1 G_2 + F_1 F_2\right]r^2 dr
!> \end{equation}
!> for `IHYPER <> 0`. Calculation of the hyperfine matrix elements.
!------------------------------------------------------------------------------------
!> @note
!> 
!> * Ryd-units are used throughout.
!> * To save storage force `JG/JF` and `PR/QR` to share the same storage by 
!> corresponding argument list in `CALL`.
!> * 28/10/94 HE tidy up, `P`,`Q` used in `<DIRAC>` instead of `g`,`f`
!> * 05/10/96 HE workspace for wavefunctions and matrices is allocated dynamically
!> * 07/02/05 VP few changes connected to the calculation of orbital polarisation
!> @endnote
!------------------------------------------------------------------------------------
module mod_ssite

contains

  !-------------------------------------------------------------------------------
  !> Summary: Assign quantum numbers and call routine to solve 8 coupled differentail radial Dirac eqautions.
  !> Author: 
  !> Category: dirac, single-site, KKRhost
  !> Deprecated: False 
  !> Assign quantum numbers and call routine to solve 8 coupled differentail radial 
  !> Dirac eqautions. The resulting wavefunctions are used to calculate t-matrices in
  !> the \(\kappa-\mu\) representation.
  !> Calculation of the radial integrals:
  !> \begin{equation}
  !> \left[ G_1 G_2 + F_1 F_2\right]r^2 dr
  !> \end{equation}
  !> for `IHYPER <> 0`. Calculation of the hyperfine matrix elements.
  !-------------------------------------------------------------------------------
  !> @note 
  !>
  !> * Ryd-units are used throughout.
  !> * To save storage force `JG/JF` and `PR/QR` to share the same storage by 
  !> corresponding argument list in `CALL`.
  !> * 28/10/94 HE tidy up, `P`,`Q` used in `<DIRAC>` instead of `g`,`f`
  !> * 05/10/96 HE workspace for wavefunctions and matrices is allocated dynamically
  !> * 07/02/05 VP few changes connected to the calculation of orbital polarisation
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine ssite(iwrregwf,iwrirrwf,nfilcbwf,calcint,getirrsol,soctl,ctl,eryd,p,   &
    ihyper,iprint,ikm1lin,ikm2lin,nlq,nkmq,nlinq,nt,nkm,iqat,tsst,msst,tsstlin,dzz, &
    dzj,szz,szj,ozz,ozj,bzz,bzj,qzz,qzj,tzz,tzj,vt,bt,at,zat,nucleus,r,drdi,r2drdi, &
    jws,imt,ameopo,lopt,solver,cgc,ozzs,ozjs,nlmax,nqmax,linmax,nrmax,nmmax,ntmax,  &
    nkmmax,nkmpmax,nlamax)

    use :: mod_datatypes, only: dp
    use :: mod_cinit
    use :: mod_cdjlzdz
    use :: mod_cdnlzdz
    use :: mod_dirac_op, only: dirabmop
    use :: mod_cinthff
    use :: mod_cintabr
    use :: mod_cjlz
    use :: mod_dirac_soc2, only: dirabmsoc2
    use :: mod_dirac_soc, only: dirabmsoc
    use :: mod_dirac_bs, only: dirbs
    use :: mod_cnlz
    use :: mod_ikapmue
    use :: mod_rinit
    use :: mod_sumupint
    use :: mod_rnuctab
    use :: mod_constants, only: ci

    implicit none

    real (kind=dp) :: f1, e0, a0, cautog

    ! conversion factor for hyperfine fields from A.U. to GAUSS
    ! electron charge     in esu
    ! Bohr-radius         in cm

    parameter (f1=1.0e0_dp, e0=1.6021892e-19_dp*2.997930e+09_dp, a0=0.52917706e-08_dp, cautog=e0/(a0*a0))

    ! Dummy arguments
    logical :: calcint, getirrsol
    complex (kind=dp) :: eryd, p
    integer :: ihyper, iprint, iwrirrwf, iwrregwf, linmax, nfilcbwf, nkm, nkmmax, nlamax, nlmax, nmmax, nqmax, nrmax, nt, ntmax, nucleus
    integer :: nkmpmax
    character (len=10) :: solver
    real (kind=dp) :: ameopo(nkmmax, nkmmax, nlamax, 3), at(nrmax, nlamax, 3, ntmax), bt(nrmax, ntmax), ctl(ntmax, nlmax), drdi(nrmax, nmmax), r(nrmax, nmmax), &
      r2drdi(nrmax, nmmax), soctl(ntmax, nlmax), vt(nrmax, ntmax)
    real (kind=dp) :: cgc(nkmpmax, 2)
    complex (kind=dp) :: ozjs(linmax, ntmax, 2), ozzs(linmax, ntmax, 2)
    complex (kind=dp) :: bzj(linmax, ntmax), bzz(linmax, ntmax), dzj(linmax, ntmax), dzz(linmax, ntmax), msst(nkmmax, nkmmax, ntmax), ozj(linmax, ntmax), ozz(linmax, ntmax), &
      qzj(linmax, ntmax), qzz(linmax, ntmax), szj(linmax, ntmax), szz(linmax, ntmax), tsst(nkmmax, nkmmax, ntmax), tsstlin(linmax, ntmax), tzj(linmax, ntmax), tzz(linmax, ntmax)
    integer :: ikm1lin(linmax), ikm2lin(linmax), imt(ntmax), iqat(nqmax, ntmax), jws(nmmax), lopt(ntmax), nkmq(nqmax), nlinq(nqmax), nlq(nqmax), zat(ntmax), muem05

    ! Local variables
    complex (kind=dp) :: a(2, 2), arg, b1, b2, cint(nrmax), crsq, csum, det, dxp(2, 2), f11, f12, f21, f22, g11, g11p, g12, g12p, g21, g21p, g22, g22p, gam(2, 2), gaminv(2, 2), hl, &
      hlb1, hlb2, jf(nrmax, 2, 2), jg(nrmax, 2, 2), jl, jlb1, jlb2, jlp, maux(nkmmax, nkmmax), msst2(2, 2), nl, nlb1, nlb2, nlp, norm, pfac, pi(2, 2, nrmax), pr(2, 2, nrmax), &
      qi(2, 2, nrmax), qr(2, 2, nrmax), rmehf(2, 2), rmehf1(2, 2), rmehf2(2, 2), sig(2, 2), tsst2(2, 2), xsst2(2, 2), zf(nrmax, 2, 2), zfjf(2, 2), zfzf(2, 2), zg(nrmax, 2, 2), &
      zgjg(2, 2), zgzg(2, 2)
    real (kind=dp) :: ap(2, 2, nrmax), aq(2, 2, nrmax), c, cff(2, 2), cfg(2, 2), cg1, cg2, cg4, cg5, cg8, cgf(2, 2), cgg(2, 2), ch(2, 2), csqr, ctf(2, 2), ctg(2, 2), dovr(nrmax), &
      drovrn(nrmax), mj, r1m(2, 2), rkd(2, 2), rnuc, sk1, sk2, tdia1, tdia2, toff
    real (kind=dp) :: cog(2, 2, 2), cof(2, 2, 2)
    integer :: i, i1, i2, i3, i5, ikm1, ikm2, il, im, in, info, ipiv(nkmmax), iq, isk1, isk2, it, j, jlim, jtop, k, k1, k2, kap1, kap2, kc, l, l1, lb1, lb2, lin, n, nsol, imkm1, &
      imkm2, is
    integer :: isign, nint
    logical :: wronski

    data r1m/1.0e0_dp, 0.0e0_dp, 0.0e0_dp, 1.0e0_dp/
    data rkd/1.0e0_dp, 0.0e0_dp, 0.0e0_dp, -1.0e0_dp/
    ! DATA RKD / 1.0D0, 0.0D0, 0.0D0, 1.0D0 /

    call cinit(ntmax*linmax, dzz)
    call cinit(ntmax*linmax, dzj)
    call cinit(ntmax*linmax, szz)
    call cinit(ntmax*linmax, szj)
    call cinit(ntmax*linmax, ozz)
    call cinit(ntmax*linmax, ozj)
    call cinit(ntmax*linmax, bzz)
    call cinit(ntmax*linmax, bzj)
    call cinit(ntmax*linmax, qzz)
    call cinit(ntmax*linmax, qzj)
    call cinit(ntmax*linmax, tzz)
    call cinit(ntmax*linmax, tzj)

    wronski = .true.
    wronski = .false.
    ! ------------------------------------------------------------------------

    c = ctl(1, 1)
    csqr = c*c

    ! calculate relativistic momentum

    p = sqrt(eryd*(1.0e0_dp+eryd/csqr))

    pfac = p/(1.0e0_dp+eryd/csqr)


    ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    do it = 1, nt

      if (iprint>0) write (1337, '(A,I3,A,10F10.4)') ' SOLVER ', it, solver, (soctl(it,il), il=1, nlmax)

      iq = iqat(1, it)
      im = imt(it)
      jtop = jws(im)
      if (nucleus/=0) then
        rnuc = rnuctab(zat(it))
        in = 1
        do while (r(in,im)<rnuc)
          in = in + 1
        end do
        ! RLIM=R(IN,IM)
        jlim = in
        if (mod(jlim,2)==0) jlim = jlim - 1
        rnuc = r(jlim, im)
      end if
      do i = 1, jtop
        dovr(i) = drdi(i, im)/r(i, im)
        if (nucleus/=0) drovrn(i) = (r(i,im)/rnuc)**3*drdi(i, im)
      end do

      lin = 0

      ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      do l = 0, (nlq(iq)-1)
        il = l + 1
        c = ctl(it, il)
        csqr = c*c

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        kap1 = -l - 1
        kap2 = l
        if (l==0) kap2 = kap1
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        isk1 = isign(1, kap1)
        isk2 = isign(1, kap2)
        sk1 = real(isk1, kind=dp)
        sk2 = real(isk2, kind=dp)
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

        if (solver(1:7)=='ABM-SOC') then
          jlp = cdjlzdz(l1, arg, 1)*p
          nlp = cdnlzdz(l1, arg, 1)*p
        end if

        ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
        do mj = -(real(l,kind=dp)+0.5_dp), +(real(l,kind=dp)+0.5_dp), 1.0_dp

          ! ------------------------------------------------------------------------
          ! NO COUPLING FOR:  ABS(MUE)= J   +  J=L+1/2 == KAP=-L-1
          if (abs(mj)>=real(l,kind=dp)) then
            nsol = 1
          else
            nsol = 2
          end if
          ! ------------------------------------------------------------------------
          muem05 = nint(mj-0.5e0_dp)
          ikm1 = ikapmue(kap1, muem05)
          ikm2 = ikapmue(kap2, muem05)
          imkm1 = ikapmue(-kap1, muem05)
          imkm2 = ikapmue(-kap2, muem05)
          ! ------------------------------------------------------------------------
          i5 = nrmax*2*2
          call cinit(i5, zg)
          call cinit(i5, zf)
          call cinit(i5, jg)
          call cinit(i5, jf)
          call cinit(i5, pr)
          call cinit(i5, qr)
          call cinit(i5, pi)
          call cinit(i5, qi)

          if (solver(1:2)=='BS') then
            call dirbs(getirrsol, ctl(it,il), eryd, l, mj, kap1, kap2, p, cg1, cg2, cg4, cg5, cg8, vt(1,it), bt(1,it), zat(it), nucleus, r(1,im), drdi(1,im), dovr, jtop, pr, qr, &
              pi, qi, zg, zf)
          else if (solver=='ABM-BI') then
            stop ' < DIRABMBI > : Not implemented. Set SOLVER=BS in inputcard'
            ! CALL DIRABMBI(GETIRRSOL,CTL(IT,IL),IT,ERYD,L,MJ,KAP1,
            ! &                          KAP2,P,CG1,CG2,CG4,CG5,CG8,AMEOPC,
            ! &                          AMEOPO,VT(1,IT),BT(1,IT),AT,Z(IT),
            ! &
            ! NUCLEUS,R(1,IM),DRDI(1,IM),DOVR,JTOP,PR,
            ! &                          QR,PI,QI,ZG,ZF,AP,AQ,NTMAX,NLAMAX,
            ! &                          NKMMAX,NRMAX)
          else if (solver(1:6)=='ABM-OP') then
            call dirabmop(getirrsol, ctl(it,il), it, eryd, l, mj, kap1, kap2, p, cg1, cg2, cg4, cg5, cg8, ameopo, vt(1,it), bt(1,it), at, zat(it), nucleus, r(1,im), drdi(1,im), &
              dovr, jtop, pr, qr, pi, qi, zg, zf, ap, aq, lopt(it), ntmax, nlamax, nkmmax, nrmax)
          else if (solver=='ABM-SOC   ') then
            call dirabmsoc(getirrsol, ctl(it,il), soctl(it,il), it, eryd, l, mj, kap1, kap2, p, cg1, cg2, cg4, cg5, cg8, vt(1,it), bt(1,it), zat(it), nucleus, r(1,im), drdi(1,im), &
              dovr, jtop, dxp, pr, qr, pi, qi, zg, zf, nrmax)
          else if (solver=='ABM-SOC-II') then
            call dirabmsoc2(getirrsol, ctl(it,il), soctl(it,il), it, eryd, l, mj, kap1, kap2, p, cg1, cg2, cg4, cg5, cg8, vt(1,it), bt(1,it), zat(it), nucleus, r(1,im), drdi(1,im), &
              dovr, jtop, dxp, pr, qr, pi, qi, zg, zf, nrmax)
          else
            write (6, *) 'No solver found for: ', solver
            stop
          end if

          ! wavefunctions at the muffin-tin-radius

          n = jtop

          g11 = pr(1, 1, n)/r(n, im)
          f11 = qr(1, 1, n)/(r(n,im)*c)
          g21 = pr(2, 1, n)/r(n, im)
          f21 = qr(2, 1, n)/(r(n,im)*c)
          g22 = pr(2, 2, n)/r(n, im)
          f22 = qr(2, 2, n)/(r(n,im)*c)
          g12 = pr(1, 2, n)/r(n, im)
          f12 = qr(1, 2, n)/(r(n,im)*c)

          if (solver(1:7)=='ABM-SOC') then
            g11p = (dxp(1,1)/drdi(n,im)-g11)/r(n, im)
            g21p = (dxp(2,1)/drdi(n,im)-g21)/r(n, im)
            g12p = (dxp(1,2)/drdi(n,im)-g12)/r(n, im)
            g22p = (dxp(2,2)/drdi(n,im)-g22)/r(n, im)
          end if

          ! ------- the minor component for the soc-manipulated wf is
          ! meaningless

          if (solver(1:7)=='ABM-SOC') then
            call cinit(2*2*nrmax, qr)
            call cinit(2*2*nrmax, qi)
          end if

          ! COSD= NL * C * F11 - PFAC * SK1 * NLB1 * G11
          ! SIND= JL * C * F11 - PFAC * SK1 * JLB1 * G11

          ! -------------------------------------------------------------------
          ! T-SS  CONSTRUCTED USING EXPRESSIONS FROM H.E. + B.L.G. (1988)
          ! -------------------------------------------------------------------

          nl = (hl-jl)/ci
          nlb1 = (hlb1-jlb1)/ci
          nlb2 = (hlb2-jlb2)/ci


          if (solver(1:7)=='ABM-SOC') then
            gam(1, 1) = jl*g11p - jlp*g11
            gam(2, 1) = jl*g21p - jlp*g21
            gam(1, 2) = jl*g12p - jlp*g12
            gam(2, 2) = jl*g22p - jlp*g22

            sig(1, 1) = nl*g11p - nlp*g11
            sig(2, 1) = nl*g21p - nlp*g21
            sig(1, 2) = nl*g12p - nlp*g12
            sig(2, 2) = nl*g22p - nlp*g22
          else
            gam(1, 1) = jl*c*f11 - pfac*sk1*jlb1*g11
            gam(2, 1) = jl*c*f21 - pfac*sk2*jlb2*g21
            gam(1, 2) = jl*c*f12 - pfac*sk1*jlb1*g12
            gam(2, 2) = jl*c*f22 - pfac*sk2*jlb2*g22

            sig(1, 1) = nl*c*f11 - pfac*sk1*nlb1*g11
            sig(2, 1) = nl*c*f21 - pfac*sk2*nlb2*g21
            sig(1, 2) = nl*c*f12 - pfac*sk1*nlb1*g12
            sig(2, 2) = nl*c*f22 - pfac*sk2*nlb2*g22
          end if

          call zcopy(nsol*nsol, gam, 1, gaminv, 1)
          call zgetrf(nsol, nsol, gaminv, 2, ipiv, info)
          call zgetri(nsol, gaminv, 2, ipiv, maux, 2*2, info)

          do i2 = 1, nsol
            do i1 = 1, nsol
              csum = 0.0e0_dp
              do i3 = 1, nsol
                csum = csum + sig(i1, i3)*gaminv(i3, i2)
              end do
              xsst2(i1, i2) = p*csum
            end do
          end do

          do i1 = 1, nsol
            do i2 = 1, nsol
              msst2(i1, i2) = -xsst2(i1, i2)
            end do
            msst2(i1, i1) = msst2(i1, i1) + ci*p
          end do

          call zcopy(nsol*nsol, msst2, 1, tsst2, 1)

          call zgetrf(nsol, nsol, tsst2, 2, ipiv, info)
          call zgetri(nsol, tsst2, 2, ipiv, maux, 2*2, info)

          if (iprint>=3) write (1337, 100) it, l, mj, ((tsst2(i1,i2),i2=1,nsol), i1=1, nsol)
          ! ------------------------------------------------------------------------

          ! COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION

          cgg(1, 1) = cg1
          cgg(1, 2) = cg2
          cgg(2, 1) = cg2
          cgg(2, 2) = cg4
          call rinit(4, cgf)
          cgf(1, 1) = cg5
          cgf(2, 2) = cg8


          ! COEFFICIENTS TO CALCULATE THE SPIN  DIPOLAR MOMENT TZ

          tdia1 = 2*mj/real((2*l1+1)*(2*lb1+1), kind=dp)
          tdia2 = 2*mj/real((2*l1+1)*(2*lb2+1), kind=dp)
          toff = -sqrt((l1+0.5e0_dp)**2-mj**2)/real(2*l1+1, kind=dp)

          ctg(1, 1) = 0.5e0_dp*(cg1-3.0e0_dp*tdia1)
          ctg(1, 2) = 0.5e0_dp*(cg2-3.0e0_dp*toff)
          ctg(2, 1) = 0.5e0_dp*(cg2-3.0e0_dp*toff)
          ctg(2, 2) = 0.5e0_dp*(cg4-3.0e0_dp*tdia2)
          call rinit(4, ctf)


          ! COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION

          cfg(1, 1) = mj*(kap1+1.0e0_dp)/(kap1+0.5e0_dp)
          cfg(2, 2) = mj*(kap2+1.0e0_dp)/(kap2+0.5e0_dp)
          cfg(1, 2) = 0.5e0_dp*sqrt(1.0e0_dp-(mj/(kap1+0.5e0_dp))**2)
          cfg(2, 1) = cfg(1, 2)
          call rinit(4, cff)
          cff(1, 1) = mj*(-kap1+1.0e0_dp)/(-kap1+0.5e0_dp)
          cff(2, 2) = mj*(-kap2+1.0e0_dp)/(-kap2+0.5e0_dp)

          ! -----------------------------------------------------------------------
          ! COEFFICIENTS TO CALCULATE THE ORBITAL POLARISATION

          call rinit(4*2, cog)
          call rinit(4*2, cof)
          do is = 1, 2
            cog(1, 1, is) = cgc(ikm1, is)*cgc(ikm1, is)*real(muem05-is+2, kind=dp)
            cof(1, 1, is) = cgc(imkm1, is)*cgc(imkm1, is)*real(muem05-is+2, kind=dp)
          end do

          if (nsol==2) then
            do is = 1, 2
              cog(2, 2, is) = cgc(ikm2, is)*cgc(ikm2, is)*real(muem05-is+2, kind=dp)
              cof(2, 2, is) = cgc(imkm2, is)*cgc(imkm2, is)*real(muem05-is+2, kind=dp)

              cog(1, 2, is) = cgc(ikm1, is)*cgc(ikm2, is)*real(muem05-is+2, kind=dp)
              cog(2, 1, is) = cog(1, 2, is)
            end do
          end if
          ! -----------------------------------------------------------------------
          ! ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
          ! THE FACTOR  I  HAS BEEN OMITTED

          ch(1, 1) = 4.0e0_dp*kap1*mj/(4.0e0_dp*kap1*kap1-1.0e0_dp)
          ch(2, 2) = 4.0e0_dp*kap2*mj/(4.0e0_dp*kap2*kap2-1.0e0_dp)
          if (nsol==2) then
            ch(1, 2) = sqrt(0.25e0_dp-(mj/real(kap1-kap2,kind=dp))**2)
            ch(2, 1) = ch(1, 2)
          end if
          ! -----------------------------------------------------------------------
          ! ALCULATE RADIAL INTEGRALS  UP TO   OR RWS(JTOP=JWS)


          if (nsol==1) then
            ! ====================================================================
            ! NO COUPLING TO OTHER SCATTERING CHANNELS
            ! REGULAR PART    Z*Z    Z = (GRA,FRA)

            norm = (p*nl-jl*xsst2(1,1))/g11

            do i = 1, jtop
              zg(i, 1, 1) = (pr(1,1,i)/r(i,im))*norm
              zf(i, 1, 1) = (qr(1,1,i)/r(i,im)/c)*norm
              jg(i, 1, 1) = pi(1, 1, i)/r(i, im)
              jf(i, 1, 1) = qi(1, 1, i)/r(i, im)/c
            end do


            if (iwrregwf/=0) write (nfilcbwf, rec=ikm1+(it-1)*nkm) it, l, mj, nsol, 'REG', kap1, ikm1, (zg(i,1,1), zf(i,1,1), i=1, jtop)

            if (iwrirrwf/=0) write (nfilcbwf, rec=ikm1+(it-1+nt)*nkm) it, l, mj, nsol, 'IRR', kap1, ikm1, (jg(i,1,1), jf(i,1,1), i=1, jtop)

            ! ============================================== NO COUPLING = END
            ! ===
          else
            ! ====================================================================
            ! COUPLING OF TWO SCATTERING CHANNELS
            ! Z(K1,K2):  INDEX 1: SPIN-ANGULAR CHARACTER
            ! INDEX 2: BOUNDARY CONDITION


            det = g11*g22 - g12*g21

            ! OEFFICIENTS TO GET:   Z(K1,K1)  Z(K2,K1)
            b1 = p*nl - xsst2(1, 1)*jl
            b2 = -xsst2(2, 1)*jl
            a(1, 1) = (g22*b1-g12*b2)/det
            a(2, 1) = (g11*b2-g21*b1)/det

            ! OEFFICIENTS TO GET:   Z(K1,K2)  Z(K2,K2)
            b1 = -xsst2(1, 2)*jl
            b2 = p*nl - xsst2(2, 2)*jl
            a(1, 2) = (g22*b1-g12*b2)/det
            a(2, 2) = (g11*b2-g21*b1)/det

            ! ALCULATE FUNCTIONS: Z(K1,K1), Z(K2,K1), Z(K1,K2), Z(K2,K2)
            do k = 1, nsol
              do i = 1, jtop
                zg(i, 1, k) = (pr(1,1,i)*a(1,k)+pr(1,2,i)*a(2,k))/r(i, im)
                zf(i, 1, k) = (qr(1,1,i)*a(1,k)+qr(1,2,i)*a(2,k))/r(i, im)/c

                zg(i, 2, k) = (pr(2,1,i)*a(1,k)+pr(2,2,i)*a(2,k))/r(i, im)
                zf(i, 2, k) = (qr(2,1,i)*a(1,k)+qr(2,2,i)*a(2,k))/r(i, im)/c
              end do
            end do
            do k = 1, nsol
              kc = 3 - k
              do i = 1, jtop
                jg(i, k, k) = pi(k, k, i)/r(i, im)
                jf(i, k, k) = qi(k, k, i)/r(i, im)/c
                jg(i, kc, k) = pi(kc, k, i)/r(i, im)
                jf(i, kc, k) = qi(kc, k, i)/r(i, im)/c
              end do
            end do

            ! -----------------------------------------------------------------------
            if (iwrregwf/=0) then
              ! solution 1
              write (nfilcbwf, rec=ikm1+(it-1)*nkm) it, l, mj, nsol, 'REG', kap1, ikm1, (zg(i,1,1), zf(i,1,1), i=1, jtop), kap2, ikm2, (zg(i,2,1), zf(i,2,1), i=1, jtop)

              ! solution 2
              write (nfilcbwf, rec=ikm2+(it-1)*nkm) it, l, mj, nsol, 'REG', kap2, ikm2, (zg(i,2,2), zf(i,2,2), i=1, jtop), kap1, ikm1, (zg(i,1,2), zf(i,1,2), i=1, jtop)
            end if

            if (iwrirrwf/=0) then
              ! solution 1
              write (nfilcbwf, rec=ikm1+(it-1+nt)*nkm) it, l, mj, nsol, 'IRR', kap1, ikm1, (jg(i,1,1), jf(i,1,1), i=1, jtop), kap2, ikm2, (jg(i,2,1), jf(i,2,1), i=1, jtop)

              ! solution 2
              write (nfilcbwf, rec=ikm2+(it-1+nt)*nkm) it, l, mj, nsol, 'IRR', kap2, ikm2, (jg(i,2,2), jf(i,2,2), i=1, jtop), kap1, ikm1, (jg(i,1,2), jf(i,1,2), i=1, jtop)

            end if
            ! ================================================= COUPLING = END
            ! ===
          end if



          ! ALCULATE SUM OF INTEGRALS TO BE MULTIPLIED TO   TAU(K1,K2)
          do k1 = 1, nsol
            do k2 = 1, nsol

              lin = lin + 1
              tsstlin(lin, it) = tsst2(k1, k2)
              if (calcint) then
                ! REGULAR PART    Z*Z

                call cintabr(zg(1,1,k1), zg(1,1,k2), zgzg, zf(1,1,k1), zf(1,1,k2), zfzf, r2drdi(1,im), nsol, nsol, jtop, nrmax)

                call sumupint(dzz(lin,it), f1, zgzg, r1m, f1, zfzf, r1m, nsol)
                call sumupint(szz(lin,it), f1, zgzg, cgg, -f1, zfzf, cgf, nsol)
                call sumupint(ozz(lin,it), f1, zgzg, cfg, -f1, zfzf, cff, nsol)
                ! ----------------------------------------------------------------------
                do is = 1, 2
                  call sumupint(ozzs(lin,it,is), f1, zgzg, cog(1,1,is), -f1, zfzf, cof(1,1,is), nsol)
                end do
                ! ----------------------------------------------------------------------
                call sumupint(qzz(lin,it), f1, zgzg, rkd, f1, zfzf, rkd, nsol)
                call sumupint(tzz(lin,it), f1, zgzg, ctg, -f1, zfzf, ctf, nsol)

                ! write(66,'(3I3,2e16.7)') it,nsol,lin,DZZ(LIN,IT)
                ! write(66,'(4e16.7)') ((ZGZG(ii,jj),ii=1,nsol),jj=1,nsol)
                ! write(66,'(4e16.7)') ((ZFZF(ii,jj),ii=1,nsol),jj=1,nsol)


                ! -----------------------------------------------------------------------
                if (ihyper==1) then
                  call cinthff(zg(1,1,k1), zf(1,1,k1), zg(1,1,k2), zf(1,1,k2), rmehf, nsol, nsol, jtop, cint, r(1,im), drdi(1,im), nrmax)

                  if (nucleus/=0) then
                    ! calculates integrals inside nucleus but up to now only
                    ! approximately because jlim is not the nuclear radius
                    ! the same arguments are valid for the irregular parts below
                    call cinthff(zg(1,1,k1), zf(1,1,k1), zg(1,1,k2), zf(1,1,k2), rmehf1, nsol, nsol, jlim, cint, r(1,im), drdi(1,im), nrmax)
                    call cinthff(zg(1,1,k1), zf(1,1,k1), zg(1,1,k2), zf(1,1,k2), rmehf2, nsol, nsol, jlim, cint, r(1,im), drovrn, nrmax)
                    do i = 1, nsol
                      do j = 1, nsol
                        rmehf(i, j) = rmehf(i, j) - rmehf1(i, j) + rmehf2(i, j)
                      end do
                    end do
                  end if
                  ! !end of nucleus.eq.0
                  call sumupint(bzz(lin,it), cautog, rmehf, ch, 0.0e0_dp, rmehf, ch, nsol)
                end if
                ! -----------------------------------------------------------------------

                ! IRREGULAR PART    Z*J
                ! THE  ENDING  A (B)  STANDS FOR THE DOMINATING (DIMINATED)
                ! SET OF SPIN-ANGULAR-CHAR:  I.E.  J==J(A,A)  FOR R>RMT

                if (k1==k2) then

                  call cintabr(zg(1,1,k1), jg(1,1,k1), zgjg, zf(1,1,k1), jf(1,1,k1), zfjf, r2drdi(1,im), nsol, nsol, jtop, nrmax)

                  call sumupint(dzj(lin,it), f1, zgjg, r1m, f1, zfjf, r1m, nsol)
                  call sumupint(szj(lin,it), f1, zgjg, cgg, -f1, zfjf, cgf, nsol)
                  call sumupint(ozj(lin,it), f1, zgjg, cfg, -f1, zfjf, cff, nsol)
                  ! ----------------------------------------------------------------------
                  do is = 1, 2
                    call sumupint(ozjs(lin,it,is), f1, zgjg, cog(1,1,is), -f1, zfjf, cof(1,1,is), nsol)
                  end do
                  ! ----------------------------------------------------------------------
                  call sumupint(qzj(lin,it), f1, zgjg, rkd, f1, zfjf, rkd, nsol)
                  call sumupint(tzj(lin,it), f1, zgjg, ctg, -f1, zfjf, ctf, nsol)

                  ! -----------------------------------------------------------------------
                  if (ihyper==1) then
                    call cinthff(zg(1,1,k1), zf(1,1,k1), jg(1,1,k1), jf(1,1,k1), rmehf, nsol, nsol, jtop, cint, r(1,im), drdi(1,im), nrmax)
                    if (nucleus/=0) then
                      ! calculates integrals inside nucleus but up to now only
                      ! approximately because jlim is not the nuclear radius
                      ! the same arguments are valid for the irregular parts
                      ! below
                      call cinthff(zg(1,1,k1), zf(1,1,k1), jg(1,1,k1), jf(1,1,k1), rmehf1, nsol, nsol, jlim, cint, r(1,im), drdi(1,im), nrmax)
                      call cinthff(zg(1,1,k1), zf(1,1,k1), jg(1,1,k1), jf(1,1,k1), rmehf2, nsol, nsol, jlim, cint, r(1,im), drovrn, nrmax)
                      do i = 1, nsol
                        do j = 1, nsol
                          rmehf(i, j) = rmehf(i, j) - rmehf1(i, j) + rmehf2(i, j)
                        end do
                      end do
                    end if
                    ! !end of nucleus.eq.0

                    call sumupint(bzj(lin,it), cautog, rmehf, ch, 0.0e0_dp, rmehf, ch, nsol)
                  end if
                end if
                ! -----------------------------------------------------------------------
              end if
              ! ! OF IF (.CALCINT.)
            end do
          end do

          ! heck WRONSKI-relationship

          if (wronski) then
            do i = 1, jtop, 40
              crsq = c*r(i, im)**2
              write (1337, 110, advance='no') it, l, nint(2*mj), i, r(i, im), 1.0e0_dp - (zf(i,1,1)*jg(i,1,1)-zg(i,1,1)*jf(i,1,1)+zf(i,2,1)*jg(i,2,1)-zg(i,2,1)*jf(i,2,1))*crsq
              if (nsol==2) then
                write (1337, 120) 1.0e0_dp - (zf(i,1,2)*jg(i,1,2)-zg(i,1,2)*jf(i,1,2)+zf(i,2,2)*jg(i,2,2)-zg(i,2,2)*jf(i,2,2))*crsq, &
                  1.0e0_dp - (zf(i,1,2)*jg(i,1,1)-zg(i,1,2)*jf(i,1,1)+zf(i,2,2)*jg(i,2,1)-zg(i,2,2)*jf(i,2,1))*crsq, 1.0e0_dp - (zf(i,1,1)*jg(i,1,2)-zg(i,1,1)*jf(i,1,2)+zf(i,2,1)* &
                  jg(i,2,2)-zg(i,2,1)*jf(i,2,2))*crsq
              else
                write (1337, *)
              end if
            end do
          end if

        end do
        ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      end do
      ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

      call cinit(nkmmax*nkmmax, tsst(1,1,it))

      do lin = 1, nlinq(iq)
        i1 = ikm1lin(lin)
        i2 = ikm2lin(lin)
        tsst(i1, i2, it) = tsstlin(lin, it)
      end do

      do j = 1, nkmq(iq)
        call zcopy(nkmq(iq), tsst(1,j,it), 1, msst(1,j,it), 1)
      end do

      call zgetrf(nkmq(iq), nkmq(iq), msst(1,1,it), nkmmax, ipiv, info)
      call zgetri(nkmq(iq), msst(1,1,it), nkmmax, ipiv, maux, nkmmax*nkmmax, info)

      if (iprint>=4) then
        do lin = 1, nlinq(iq)
          i1 = ikm1lin(lin)
          i2 = ikm2lin(lin)
          write (1337, 130) it, lin
          write (1337, 130) it, i1, i2, ' DZZ ', dzz(lin, it), dzj(lin, it)
          write (1337, 130) it, i1, i2, ' SZZ ', szz(lin, it), szj(lin, it)
          write (1337, 130) it, i1, i2, ' OZZ ', ozz(lin, it), ozj(lin, it)
        end do
      end if
    end do
    ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

100 format ('  t-ss(TLM)', 2i3, f5.1, 2e14.5, 2x, 2e14.5, :, /, 22x, 2e14.5, 2x, 2e14.5)
110 format (' IT=', i2, ' L=', i2, ' MJ=', i2, '/2', i7, f10.6, 1(2x,2f9.6))
120 format (3(2x,2f9.6))
130 format (' IT=', i2, 2i3, a, 2x, 2e14.5, 2x, 2e14.5)
  end subroutine ssite

  !-------------------------------------------------------------------------------
  !> Summary: Re-read the wave functions written by `<SSITE>` or `<CORE>`
  !> Author: Who wrote this subroutine
  !> Category: dirac, input-output, KKRhost
  !> Deprecated: False 
  !> Re-read the wave functions written by `<SSITE>` or `<CORE>`
  !-------------------------------------------------------------------------------
  subroutine readwfun(nfil,it,l,mj,nsol,sreg,sirr,ikm1,kap1,ikm2,kap2,nt,nkm,zg,zf, &
    jg,jf,jtop,nrmax)

    use :: mod_datatypes, only: dp

    implicit none

    ! Dummy arguments
    integer :: ikm1, ikm2, it, jtop, kap1, kap2, l, nfil, nkm, nrmax, nsol, nt
    real (kind=dp) :: mj
    character (len=3) :: sirr, sreg
    complex (kind=dp) :: jf(nrmax, 2, 2), jg(nrmax, 2, 2), zf(nrmax, 2, 2), zg(nrmax, 2, 2)

    ! Local variables
    integer :: i, iflag, ikmin(2), itp, k, kapin(2), lp, nsolp
    real (kind=dp) :: mjp
    character (len=3) :: strp

    iflag = 0
    ! -----------------------------------------------------------------------
    ! REGULAR wave function -- solution 1
    if (sreg=='REG' .or. sreg=='COR') then
      read (nfil, rec=ikm1+(it-1)*nkm) itp, lp, mjp, nsolp, strp, (kapin(k), ikmin(k), (zg(i,k,1),zf(i,k,1),i=1,jtop), k=1, nsol)
      if (itp/=it .or. lp/=l .or. abs(mjp-mj)>0.001e0_dp .or. nsolp/=nsol .or. strp/='REG') iflag = iflag + 1
      if (kap1/=kapin(1)) iflag = iflag + 1
      if (ikm1/=ikmin(1)) iflag = iflag + 1
      if (nsol>1) then
        if (kap2/=kapin(2)) iflag = iflag + 1
        if (ikm2/=ikmin(2)) iflag = iflag + 1
      end if
    end if

    ! -----------------------------------------------------------------------
    ! IRREGULAR wave function -- solution 1
    if (sirr=='IRR') then
      read (nfil, rec=ikm1+(it-1+nt)*nkm) itp, lp, mjp, nsolp, strp, (kapin(k), ikmin(k), (jg(i,k,1),jf(i,k,1),i=1,jtop), k=1, nsol)
      if (itp/=it .or. lp/=l .or. abs(mjp-mj)>0.001e0_dp .or. nsolp/=nsol .or. strp/='IRR') iflag = iflag + 1
      if (kap1/=kapin(1)) iflag = iflag + 1
      if (ikm1/=ikmin(1)) iflag = iflag + 1
      if (nsol>1) then
        if (kap2/=kapin(2)) iflag = iflag + 1
        if (ikm2/=ikmin(2)) iflag = iflag + 1
      end if
    end if

    if (nsol==2) then
      ! -----------------------------------------------------------------------
      ! REGULAR wave function -- solution 2
      if (sreg=='REG' .or. sreg=='COR') then
        read (nfil, rec=ikm2+(it-1)*nkm) itp, lp, mjp, nsolp, strp, (kapin(k), ikmin(k), (zg(i,k,2),zf(i,k,2),i=1,jtop), k=2, 1, -1)
        if (itp/=it .or. lp/=l .or. abs(mjp-mj)>0.001e0_dp .or. nsolp/=nsol .or. strp/='REG') iflag = iflag + 1
      end if

      ! -----------------------------------------------------------------------
      ! IRREGULAR wave function -- solution 2
      if (sirr=='IRR') then
        read (nfil, rec=ikm2+(it-1+nt)*nkm) itp, lp, mjp, nsolp, strp, (kapin(k), ikmin(k), (jg(i,k,2),jf(i,k,2),i=1,jtop), k=2, 1, -1)
        if (itp/=it .or. lp/=l .or. abs(mjp-mj)>0.001e0_dp .or. nsolp/=nsol .or. strp/='IRR') iflag = iflag + 1
      end if

    end if

    ! -----------------------------------------------------------------------

    if (iflag>0) then
      write (*, 100) iflag, it, l, mj, sreg, sirr
      stop ' in <READWFUN>'
    end if
100 format (/, /, 1x, 79('*'), /, 10x, 'error reading the wave functions', ' IFLAG = 1', /, 10x, 'for  IT=', i2, ' L=', i2, ' MJ=', f4.1, ' KEYS=', a, a)
  end subroutine readwfun

end module mod_ssite
